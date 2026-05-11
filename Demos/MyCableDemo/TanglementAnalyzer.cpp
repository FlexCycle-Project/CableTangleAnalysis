#include "TanglementAnalyzer.h"

#include "KymoKnotAnalyzer.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>

namespace
{
    static Real cross2D(const Vector2r& a, const Vector2r& b)
    {
        return a.x() * b.y() - a.y() * b.x();
    }

    static Real orient2D(const Vector2r& a, const Vector2r& b, const Vector2r& c)
    {
        return cross2D(b - a, c - a);
    }

    static bool onSegment2D(const Vector2r& a, const Vector2r& b, const Vector2r& p)
    {
        const Real eps = static_cast<Real>(1e-9);
        return p.x() >= std::min(a.x(), b.x()) - eps &&
               p.x() <= std::max(a.x(), b.x()) + eps &&
               p.y() >= std::min(a.y(), b.y()) - eps &&
               p.y() <= std::max(a.y(), b.y()) + eps;
    }

    static bool segmentsIntersect2D(
        const Vector2r& a0,
        const Vector2r& a1,
        const Vector2r& b0,
        const Vector2r& b1)
    {
        const Real eps = static_cast<Real>(1e-9);
        const Real o1 = orient2D(a0, a1, b0);
        const Real o2 = orient2D(a0, a1, b1);
        const Real o3 = orient2D(b0, b1, a0);
        const Real o4 = orient2D(b0, b1, a1);

        auto sign = [eps](Real value) -> int
        {
            if (value > eps)
                return 1;
            if (value < -eps)
                return -1;
            return 0;
        };

        const int s1 = sign(o1);
        const int s2 = sign(o2);
        const int s3 = sign(o3);
        const int s4 = sign(o4);

        if (s1 == 0 && s2 == 0 && s3 == 0 && s4 == 0)
            return false;

        if (s1 == 0 && onSegment2D(a0, a1, b0)) return true;
        if (s2 == 0 && onSegment2D(a0, a1, b1)) return true;
        if (s3 == 0 && onSegment2D(b0, b1, a0)) return true;
        if (s4 == 0 && onSegment2D(b0, b1, a1)) return true;

        return (s1 * s2 < 0) && (s3 * s4 < 0);
    }

    static bool pointInPolygon2D(
        const std::vector<Vector2r>& polygon,
        const Vector2r& point)
    {
        if (polygon.size() < 3)
            return false;

        bool inside = false;
        const Real eps = static_cast<Real>(1e-9);
        for (size_t i = 0, j = polygon.size() - 1; i < polygon.size(); j = i++)
        {
            const Vector2r& a = polygon[j];
            const Vector2r& b = polygon[i];
            const Real cross = std::abs(
                (b.x() - a.x()) * (point.y() - a.y()) -
                (b.y() - a.y()) * (point.x() - a.x()));
            if (cross <= eps &&
                point.x() >= std::min(a.x(), b.x()) - eps &&
                point.x() <= std::max(a.x(), b.x()) + eps &&
                point.y() >= std::min(a.y(), b.y()) - eps &&
                point.y() <= std::max(a.y(), b.y()) + eps)
            {
                return true;
            }

            const bool intersect =
                ((a.y() > point.y()) != (b.y() > point.y())) &&
                (point.x() < (b.x() - a.x()) * (point.y() - a.y()) /
                    ((b.y() - a.y()) + static_cast<Real>(1e-12)) + a.x());
            if (intersect)
                inside = !inside;
        }

        return inside;
    }

    static std::vector<TanglementAnalyzer::CrossingEvent> findCrossings(
        const std::vector<std::vector<Vector3r>>& ropePieces)
    {
        std::vector<TanglementAnalyzer::CrossingEvent> crossingEvents;
        if (ropePieces.empty())
            return crossingEvents;

        for (size_t pieceIndex = 0; pieceIndex < ropePieces.size(); ++pieceIndex)
        {
            const std::vector<Vector3r>& piece = ropePieces[pieceIndex];
            if (piece.size() < 2)
                continue;

            const int lastEdge = static_cast<int>(piece.size()) - 2;
            for (int edgeA = 0; edgeA <= lastEdge; ++edgeA)
            {
                for (int edgeB = edgeA + 2; edgeB <= lastEdge; ++edgeB)
                {
                    const Vector2r a0(
                        piece[static_cast<size_t>(edgeA)].x(),
                        piece[static_cast<size_t>(edgeA)].z());
                    const Vector2r a1(
                        piece[static_cast<size_t>(edgeA + 1)].x(),
                        piece[static_cast<size_t>(edgeA + 1)].z());
                    const Vector2r b0(
                        piece[static_cast<size_t>(edgeB)].x(),
                        piece[static_cast<size_t>(edgeB)].z());
                    const Vector2r b1(
                        piece[static_cast<size_t>(edgeB + 1)].x(),
                        piece[static_cast<size_t>(edgeB + 1)].z());

                    if (!segmentsIntersect2D(a0, a1, b0, b1))
                        continue;

                    const Vector2r ad = a1 - a0;
                    const Vector2r bd = b1 - b0;
                    const Vector2r diff = b0 - a0;
                    const Real denom = ad.x() * bd.y() - ad.y() * bd.x();
                    if (std::abs(denom) <= static_cast<Real>(1e-10))
                        continue;

                    const Real ua = (diff.x() * bd.y() - diff.y() * bd.x()) / denom;
                    const Real ub = (diff.x() * ad.y() - diff.y() * ad.x()) / denom;
                    if (ua <= static_cast<Real>(1e-8) || ua >= static_cast<Real>(1.0) - static_cast<Real>(1e-8) ||
                        ub <= static_cast<Real>(1e-8) || ub >= static_cast<Real>(1.0) - static_cast<Real>(1e-8))
                    {
                        continue;
                    }

                    const Vector3r firstPoint3D =
                        piece[static_cast<size_t>(edgeA)] +
                        ua * (piece[static_cast<size_t>(edgeA + 1)] - piece[static_cast<size_t>(edgeA)]);
                    const Vector3r secondPoint3D =
                        piece[static_cast<size_t>(edgeB)] +
                        ub * (piece[static_cast<size_t>(edgeB + 1)] - piece[static_cast<size_t>(edgeB)]);

                    TanglementAnalyzer::CrossingEvent event;
                    event.order = static_cast<Real>(0.5) *
                        (static_cast<Real>(edgeA) + static_cast<Real>(edgeB));
                    event.pieceIndex = pieceIndex;
                    event.firstEdge = edgeA;
                    event.secondEdge = edgeB;
                    event.heightDelta = firstPoint3D.y() - secondPoint3D.y();
                    event.firstOverSecond = event.heightDelta > static_cast<Real>(0.0);
                    crossingEvents.push_back(event);
                }
            }
        }

        return crossingEvents;
    }

    static std::vector<TanglementAnalyzer::CrossingEvent> filterCrossingsForRegion(
        const std::vector<TanglementAnalyzer::CrossingEvent>& allCrossingEvents,
        size_t pieceIndex,
        const TanglementAnalyzer::LocalKnotRegion& region)
    {
        std::vector<TanglementAnalyzer::CrossingEvent> crossingEvents;
        const int regionStart = region.start;
        const int regionEnd = region.end;
        for (const TanglementAnalyzer::CrossingEvent& event : allCrossingEvents)
        {
            if (event.pieceIndex != pieceIndex)
                continue;
            const int firstEdge = std::min(event.firstEdge, event.secondEdge);
            const int secondEdge = std::max(event.firstEdge, event.secondEdge);
            if ((firstEdge < regionStart || firstEdge > regionEnd) &&
                (secondEdge < regionStart || secondEdge > regionEnd))
                continue;
            crossingEvents.push_back(event);
        }

        return crossingEvents;
    }

    static bool areMutuallyIntersectingHoops(
        const TanglementAnalyzer::HoopCandidate& a,
        const TanglementAnalyzer::HoopCandidate& b)
    {
        int aStart = std::min(a.startEdge, a.endEdge);
        int aEnd = std::max(a.startEdge, a.endEdge);
        int bStart = std::min(b.startEdge, b.endEdge);
        int bEnd = std::max(b.startEdge, b.endEdge);

        const bool overlap = aStart <= bEnd && bStart <= aEnd;
        const bool aContainsB = aStart <= bStart && aEnd >= bEnd;
        const bool bContainsA = bStart <= aStart && bEnd >= aEnd;
        return overlap && !aContainsB && !bContainsA;
    }

    static std::pair<int, int> scoreHoopPair(
        const std::vector<TanglementAnalyzer::CrossingEvent>& regionCrossings,
        const std::vector<Vector3r>& piece,
        const TanglementAnalyzer::HoopCandidate& hoopA,
        const TanglementAnalyzer::HoopCandidate& hoopB)
    {
        const auto edgeInHoop = [](int edgeIndex, const TanglementAnalyzer::HoopCandidate& hoop) -> bool
        {
            const int hoopStart = std::min(hoop.startEdge, hoop.endEdge);
            const int hoopEnd = std::max(hoop.startEdge, hoop.endEdge);
            return edgeIndex >= hoopStart && edgeIndex <= hoopEnd;
        };

        const auto edgeInExclusiveTargetPart = [](int edgeIndex,
                                                  const TanglementAnalyzer::HoopCandidate& sourceHoop,
                                                  const TanglementAnalyzer::HoopCandidate& targetHoop) -> bool
        {
            const int s0 = std::min(sourceHoop.startEdge, sourceHoop.endEdge);
            const int s1 = std::max(sourceHoop.startEdge, sourceHoop.endEdge);
            const int t0 = std::min(targetHoop.startEdge, targetHoop.endEdge);
            const int t1 = std::max(targetHoop.startEdge, targetHoop.endEdge);

            return edgeIndex >= t0 && edgeIndex <= t1 &&
                   !(edgeIndex > s0 && edgeIndex < s1);
        };

        const auto scoreDirectional = [&](const TanglementAnalyzer::HoopCandidate& sourceHoop,
                                          const TanglementAnalyzer::HoopCandidate& targetHoop) -> int
        {
            std::vector<Vector2r> sourceHoopPolygon;
            const int sourceStart = std::min(sourceHoop.startEdge, sourceHoop.endEdge);
            const int sourceEnd = std::max(sourceHoop.startEdge, sourceHoop.endEdge);
            const int sourceLastVertex = std::min(sourceEnd + 1, static_cast<int>(piece.size()) - 1);
            sourceHoopPolygon.reserve(static_cast<size_t>(sourceLastVertex - sourceStart + 1));
            for (int vertexIndex = sourceStart; vertexIndex <= sourceLastVertex; ++vertexIndex)
            {
                const Vector3r& point3D = piece[static_cast<size_t>(vertexIndex)];
                sourceHoopPolygon.emplace_back(point3D.x(), point3D.z());
            }

            std::vector<std::pair<Real, bool>> states;
            states.reserve(regionCrossings.size());

            for (const TanglementAnalyzer::CrossingEvent& event : regionCrossings)
            {
                const bool firstInSource = edgeInHoop(event.firstEdge, sourceHoop);
                const bool secondInSource = edgeInHoop(event.secondEdge, sourceHoop);
                const bool firstInTargetOnly =
                    edgeInExclusiveTargetPart(event.firstEdge, sourceHoop, targetHoop);
                const bool secondInTargetOnly =
                    edgeInExclusiveTargetPart(event.secondEdge, sourceHoop, targetHoop);

                if (firstInSource && secondInTargetOnly)
                {
                    states.emplace_back(static_cast<Real>(event.firstEdge), event.firstOverSecond);
                }
                else if (secondInSource && firstInTargetOnly)
                {
                    states.emplace_back(static_cast<Real>(event.secondEdge), !event.firstOverSecond);
                }
            }

            std::sort(
                states.begin(),
                states.end(),
                [](const std::pair<Real, bool>& lhs, const std::pair<Real, bool>& rhs)
                {
                    return lhs.first < rhs.first;
                });

            int score = 0;
            for (size_t i = 0; i + 1 < states.size(); ++i)
            {
                const bool firstState = states[i].second;
                const bool secondState = states[i + 1].second;
                if (secondState != firstState)
                {
                    const Real midpointEdge =
                        static_cast<Real>(0.5) * (states[i].first + states[i + 1].first);
                    const Real clampedEdgeIndex = std::max(
                        static_cast<Real>(0.0),
                        std::min(midpointEdge, static_cast<Real>(piece.size() - 1)));
                    const size_t midpointIndex = static_cast<size_t>(
                        std::round(clampedEdgeIndex));
                    const Vector3r& point3D = piece[midpointIndex];
                    const Vector2r midpoint2D(point3D.x(), point3D.z());
                    if (pointInPolygon2D(sourceHoopPolygon, midpoint2D))
                        ++score;
                }
            }

            return score;
        };

        const int scoreA = scoreDirectional(hoopA, hoopB);
        const int scoreB = scoreDirectional(hoopB, hoopA);
        return std::make_pair(scoreA, scoreB);
    }

}

namespace
{
    static std::vector<int> buildProblematicHoopMiddleHighlights(
        const std::vector<Vector3r>& piece,
        const std::vector<TanglementAnalyzer::CrossingEvent>& allCrossingEvents,
        size_t pieceIndex,
        const TanglementAnalyzer::LocalKnotRegion& knotRegion,
        int pieceEdgeCount,
        std::vector<TanglementAnalyzer::HoopCandidate>* outHoops,
        std::vector<int>* outTargetHoopMiddleEdges)
    {
        std::vector<int> middleEdges;
        if (outHoops != nullptr)
            outHoops->clear();
        if (outTargetHoopMiddleEdges != nullptr)
            outTargetHoopMiddleEdges->clear();

        if (pieceEdgeCount < 3 || allCrossingEvents.empty())
            return middleEdges;

        const std::vector<TanglementAnalyzer::CrossingEvent> crossingEvents =
            filterCrossingsForRegion(allCrossingEvents, pieceIndex, knotRegion);
        if (crossingEvents.empty())
            return middleEdges;
        std::vector<TanglementAnalyzer::HoopCandidate> hoops;
        hoops.reserve(crossingEvents.size());
        for (const TanglementAnalyzer::CrossingEvent& event : crossingEvents)
        {
            TanglementAnalyzer::HoopCandidate hoop;
            hoop.startEdge = std::min(event.firstEdge, event.secondEdge);
            hoop.endEdge = std::max(event.firstEdge, event.secondEdge);
            hoop.order = event.order;
            hoops.push_back(hoop);
        }

        const auto hoopContains = [](const TanglementAnalyzer::HoopCandidate& outer,
                                      const TanglementAnalyzer::HoopCandidate& inner) -> bool
        {
            const int outerStart = std::min(outer.startEdge, outer.endEdge);
            const int outerEnd = std::max(outer.startEdge, outer.endEdge);
            const int innerStart = std::min(inner.startEdge, inner.endEdge);
            const int innerEnd = std::max(inner.startEdge, inner.endEdge);
            return outerStart <= innerStart && outerEnd >= innerEnd;
        };

        std::vector<int> pairSumScores(hoops.size(), 0);
        std::vector<int> blackPointScores(hoops.size(), 0);
        // Nested hoops stay in the data. For each source hoop, pick the best
        // left-side and right-side representative that still intersects it,
        // then score those representative sides against the source once.
        for (size_t sourceIndex = 0; sourceIndex < hoops.size(); ++sourceIndex)
        {
            const TanglementAnalyzer::HoopCandidate& sourceHoop = hoops[sourceIndex];
            const int s0 = std::min(sourceHoop.startEdge, sourceHoop.endEdge);
            const int s1 = std::max(sourceHoop.startEdge, sourceHoop.endEdge);

            size_t bestLeftIndex = std::numeric_limits<size_t>::max();
            int bestLeftCoverage = std::numeric_limits<int>::min();
            size_t bestRightIndex = std::numeric_limits<size_t>::max();
            int bestRightCoverage = std::numeric_limits<int>::min();

            for (size_t targetIndex = 0; targetIndex < hoops.size(); ++targetIndex)
            {
                if (targetIndex == sourceIndex)
                    continue;

                const TanglementAnalyzer::HoopCandidate& targetHoop = hoops[targetIndex];
                if (!areMutuallyIntersectingHoops(sourceHoop, targetHoop))
                    continue;

                const int t0 = std::min(targetHoop.startEdge, targetHoop.endEdge);
                const int t1 = std::max(targetHoop.startEdge, targetHoop.endEdge);

                if (t0 < s0 && t1 >= s0)
                {
                    const int leftCoverage = s0 - t0;
                    if (leftCoverage > bestLeftCoverage)
                    {
                        bestLeftCoverage = leftCoverage;
                        bestLeftIndex = targetIndex;
                    }
                }

                if (t1 > s1 && t0 <= s1)
                {
                    const int rightCoverage = t1 - s1;
                    if (rightCoverage > bestRightCoverage)
                    {
                        bestRightCoverage = rightCoverage;
                        bestRightIndex = targetIndex;
                    }
                }
            }

            if (bestLeftIndex != std::numeric_limits<size_t>::max() &&
                bestRightIndex != std::numeric_limits<size_t>::max() &&
                bestLeftIndex == bestRightIndex)
            {
                const std::pair<int, int> pairScores =
                    scoreHoopPair(crossingEvents, piece, sourceHoop, hoops[bestLeftIndex]);
                if (pairScores.first != 0)
                    pairSumScores[sourceIndex] += pairScores.first;
                continue;
            }

            if (bestLeftIndex != std::numeric_limits<size_t>::max())
            {
                const std::pair<int, int> pairScores =
                    scoreHoopPair(crossingEvents, piece, sourceHoop, hoops[bestLeftIndex]);
                if (pairScores.first != 0)
                    pairSumScores[sourceIndex] += pairScores.first;
                if (pairScores.first != 0 && pairScores.second != 0)
                    blackPointScores[sourceIndex] += pairScores.first;
            }

            if (bestRightIndex != std::numeric_limits<size_t>::max())
            {
                const std::pair<int, int> pairScores =
                    scoreHoopPair(crossingEvents, piece, sourceHoop, hoops[bestRightIndex]);
                if (pairScores.first != 0)
                    pairSumScores[sourceIndex] += pairScores.first;
                if (pairScores.first != 0 && pairScores.second != 0)
                    blackPointScores[sourceIndex] += pairScores.first;
            }
        }

        for (size_t i = 0; i < hoops.size(); ++i)
        {
            hoops[i].pairSumScore = pairSumScores[i];
        }

        if (outHoops != nullptr)
            *outHoops = hoops;

        for (size_t i = 0; i < hoops.size(); ++i)
        {
            if (blackPointScores[i] == 0)
                continue;

            const int hoopStart = std::max(0, std::min(hoops[i].startEdge, hoops[i].endEdge));
            const int hoopEnd = std::min(
                std::max(0, pieceEdgeCount - 1),
                std::max(hoops[i].startEdge, hoops[i].endEdge));
            for (int edgeIndex = hoopStart; edgeIndex <= hoopEnd; ++edgeIndex)
                middleEdges.push_back(edgeIndex);
        }

        std::sort(middleEdges.begin(), middleEdges.end());
        middleEdges.erase(std::unique(middleEdges.begin(), middleEdges.end()), middleEdges.end());

        if (outTargetHoopMiddleEdges != nullptr)
        {
            for (size_t i = 0; i < hoops.size(); ++i)
            {
                if (blackPointScores[i] == 0)
                    continue;

                int targetMiddleEdge = static_cast<int>(
                    std::round(
                        static_cast<Real>(0.5) *
                        (static_cast<Real>(hoops[i].startEdge) +
                         static_cast<Real>(hoops[i].endEdge))));
                targetMiddleEdge = std::max(0, targetMiddleEdge);
                targetMiddleEdge = std::min(std::max(0, pieceEdgeCount - 1), targetMiddleEdge);
                outTargetHoopMiddleEdges->push_back(targetMiddleEdge);
            }

            std::sort(outTargetHoopMiddleEdges->begin(), outTargetHoopMiddleEdges->end());
            outTargetHoopMiddleEdges->erase(
                std::unique(
                    outTargetHoopMiddleEdges->begin(),
                    outTargetHoopMiddleEdges->end()),
                outTargetHoopMiddleEdges->end());
        }

        return middleEdges;
    }

}

namespace
{
    static TanglementAnalyzer::TanglementMetricsSnapshot buildTanglementMetricsSnapshot(
        const std::vector<TanglementAnalyzer::CrossingEvent>& allCrossingEvents,
        const std::vector<std::vector<Vector3r>>& ropePieces,
        const std::vector<SelfCollisionContact>& contacts,
        const std::vector<PBD::KymoKnotAnalysisResult>& kymoResults)
    {
        TanglementAnalyzer::TanglementMetricsSnapshot snapshot;
        snapshot.pieceCount = ropePieces.size();
        snapshot.selfCollisionCount = contacts.size();

        Real longestPieceSpanFraction = static_cast<Real>(0.0);
        size_t piecesWithKnot = 0;

        for (const std::vector<Vector3r>& piece : ropePieces)
            snapshot.totalPointCount += piece.size();

        std::vector<unsigned int> touchedEdges;
        touchedEdges.reserve(contacts.size() * 2);
        for (const SelfCollisionContact& contact : contacts)
        {
            touchedEdges.push_back(contact.edgeA);
            touchedEdges.push_back(contact.edgeB);
        }
        std::sort(touchedEdges.begin(), touchedEdges.end());
        touchedEdges.erase(
            std::unique(touchedEdges.begin(), touchedEdges.end()),
            touchedEdges.end());
        snapshot.contactParticipationCount = touchedEdges.size();

        snapshot.gravityCrossings = static_cast<int>(allCrossingEvents.size());

        for (size_t pieceIndex = 0; pieceIndex < ropePieces.size(); ++pieceIndex)
        {
            const std::vector<Vector3r>& piece = ropePieces[pieceIndex];
            const size_t pieceEdges = piece.size() > 1 ? piece.size() - 1 : 0;
            if (pieceIndex >= kymoResults.size() || pieceEdges == 0)
                continue;

            const PBD::KymoKnotAnalysisResult& result = kymoResults[pieceIndex];
            if (!result.available || !result.success || result.segments.empty())
                continue;

            ++piecesWithKnot;
            size_t pieceLongestSegment = 0;
            for (const PBD::KymoKnotSegment& segment : result.segments)
            {
                const size_t segmentLength =
                    static_cast<size_t>(std::max(0, segment.length));
                pieceLongestSegment = std::max(pieceLongestSegment, segmentLength);
            }

            const Real pieceSpanFraction = static_cast<Real>(pieceLongestSegment) /
                static_cast<Real>(std::max<size_t>(1, pieceEdges));
            longestPieceSpanFraction += pieceSpanFraction;
        }

        snapshot.knotSpanFraction = piecesWithKnot > 0
            ? longestPieceSpanFraction / static_cast<Real>(piecesWithKnot)
            : static_cast<Real>(0.0);

        return snapshot;
    }
}

TanglementAnalyzer::TanglementAnalysisResult TanglementAnalyzer::analyzeTanglement(
    const std::vector<std::vector<Vector3r>>& ropePieces,
    const std::vector<SelfCollisionContact>& contacts,
    PBD::KymoKnotSearchMode searchMode)
{
    TanglementAnalysisResult result;
    result.knotRegions.resize(ropePieces.size());
    result.crossingEdges.resize(ropePieces.size());
    result.crossingCenterEdges.resize(ropePieces.size());
    result.crossingMiddleEdges.resize(ropePieces.size());
    result.scoredHoops.resize(ropePieces.size());
    result.problematicHoopMiddleEdges.resize(ropePieces.size());
    result.targetHoopMiddleEdges.resize(ropePieces.size());
    result.kymoResults.reserve(ropePieces.size());

    const std::vector<CrossingEvent> allCrossingEvents = findCrossings(ropePieces);

    for (size_t pieceIndex = 0; pieceIndex < ropePieces.size(); ++pieceIndex)
    {
        const std::vector<Vector3r>& cablePoints = ropePieces[pieceIndex];
        const PBD::KymoKnotAnalysisResult knotResult =
            PBD::KymoKnotAnalyzer::analyze(
                PBD::CableGeometry{ cablePoints },
                searchMode);

        result.kymoResults.push_back(knotResult);
        if (!knotResult.available || !knotResult.success || knotResult.segments.empty())
            continue;

        for (const PBD::KymoKnotSegment& segment : knotResult.segments)
        {
            const int localStart = std::max(0, segment.start);
            const int localEnd = std::min(
                static_cast<int>(cablePoints.size()) - 1,
                segment.end);
            if (localEnd < localStart)
                continue;
            result.knotRegions[pieceIndex].push_back(
                LocalKnotRegion{
                    localStart,
                    localEnd,
                    static_cast<size_t>(std::max(1, segment.length))});
            break;
        }
    }

    for (size_t pieceIndex = 0; pieceIndex < ropePieces.size(); ++pieceIndex)
    {
        const std::vector<Vector3r>& piece = ropePieces[pieceIndex];
        if (result.knotRegions[pieceIndex].empty())
            continue;

        const LocalKnotRegion& knotRegion = result.knotRegions[pieceIndex].front();
        const std::vector<CrossingEvent> regionCrossings =
            filterCrossingsForRegion(
                allCrossingEvents,
                pieceIndex,
                knotRegion);
        if (!regionCrossings.empty())
        {
            std::vector<int> edges;
            std::vector<int> centerEdges;
            std::vector<int> middleEdges;

            for (const CrossingEvent& event : regionCrossings)
            {
                edges.push_back(event.firstEdge);
                edges.push_back(event.secondEdge);
                centerEdges.push_back(static_cast<int>(std::round(event.order)));
            }

            std::sort(edges.begin(), edges.end());
            edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
            std::sort(centerEdges.begin(), centerEdges.end());
            centerEdges.erase(std::unique(centerEdges.begin(), centerEdges.end()), centerEdges.end());

            for (size_t i = 0; i + 1 < edges.size(); ++i)
            {
                int centerEdge = static_cast<int>(
                    std::round(
                        static_cast<Real>(0.5) *
                        (static_cast<Real>(edges[i]) +
                         static_cast<Real>(edges[i + 1]))));
                centerEdge = std::max(0, centerEdge);
                centerEdge = std::min(static_cast<int>(piece.size()) - 2, centerEdge);
                middleEdges.push_back(centerEdge);
            }

            std::sort(middleEdges.begin(), middleEdges.end());
            middleEdges.erase(std::unique(middleEdges.begin(), middleEdges.end()), middleEdges.end());

            result.crossingEdges[pieceIndex] = std::move(edges);
            result.crossingCenterEdges[pieceIndex] = std::move(centerEdges);
            result.crossingMiddleEdges[pieceIndex] = std::move(middleEdges);
        }

        result.problematicHoopMiddleEdges[pieceIndex] =
            buildProblematicHoopMiddleHighlights(
                piece,
                allCrossingEvents,
                pieceIndex,
                knotRegion,
                static_cast<int>(piece.size() > 1 ? piece.size() - 1 : 0),
                &result.scoredHoops[pieceIndex],
                &result.targetHoopMiddleEdges[pieceIndex]);
    }

    result.metrics = buildTanglementMetricsSnapshot(
        allCrossingEvents,
        ropePieces,
        contacts,
        result.kymoResults);
    return result;
}

namespace
{
    static void printHoopScoresForPiece(
        const std::vector<TanglementAnalyzer::HoopCandidate>& hoops)
    {
        bool printedAnyScore = false;
        for (size_t i = 0; i < hoops.size(); ++i)
        {
            if (hoops[i].pairSumScore == 0)
                continue;

            if (!printedAnyScore)
            {
                std::cout << "    hoop scores:" << std::endl;
                printedAnyScore = true;
            }
            std::cout << "      (" << i << ") pairSumScore=" << hoops[i].pairSumScore
                      << std::endl;
        }
    }

}

void TanglementAnalyzer::printTanglementReport(
    const std::vector<std::vector<Vector3r>>& ropePieces,
    const TanglementAnalysisResult& result)
{
    const std::vector<CrossingEvent> allCrossingEvents = findCrossings(ropePieces);

    std::cout << "Tanglement metrics:" << std::endl;
    std::cout << "  global: pieces=" << result.metrics.pieceCount
              << " points=" << result.metrics.totalPointCount
              << " selfCollisions=" << result.metrics.selfCollisionCount
              << " (" << result.metrics.contactParticipationCount << ")"
              << " crossings=" << result.metrics.gravityCrossings
              << std::endl;
    std::cout << "  global: knotSpanFraction=" << result.metrics.knotSpanFraction
              << std::endl;

    for (size_t pieceIndex = 0; pieceIndex < ropePieces.size(); ++pieceIndex)
    {
        std::cout << "  piece " << pieceIndex << " crossings:";
        bool printedAnyCrossings = false;
        for (const CrossingEvent& event : allCrossingEvents)
        {
            if (event.pieceIndex != pieceIndex)
                continue;

            const int lo = std::min(event.firstEdge, event.secondEdge);
            const int hi = std::max(event.firstEdge, event.secondEdge);
            const char* overUnder =
                std::abs(event.heightDelta) <= static_cast<Real>(1e-8)
                    ? "unclear"
                    : (event.firstOverSecond ? "first-over" : "second-over");
            std::cout << " [enc=" << lo << "-" << hi
                      << ", edges=" << event.firstEdge
                      << "," << event.secondEdge
                      << ", " << overUnder
                      << "]";
            printedAnyCrossings = true;
        }
        if (!printedAnyCrossings)
            std::cout << " none";
        std::cout << std::endl;
    }

    for (size_t pieceIndex = 0; pieceIndex < ropePieces.size(); ++pieceIndex)
    {
        if (result.knotRegions[pieceIndex].empty())
            continue;

        std::cout << "  piece " << pieceIndex << std::endl;
        printHoopScoresForPiece(result.scoredHoops[pieceIndex]);
    }

    for (size_t pieceIndex = 0; pieceIndex < ropePieces.size(); ++pieceIndex)
    {
        const std::vector<Vector3r>& cablePoints = ropePieces[pieceIndex];
        if (cablePoints.size() < 2)
        {
            std::cout << "  piece " << pieceIndex << ": too short" << std::endl;
            continue;
        }

        std::cout << "  piece " << pieceIndex << ": ";
        const PBD::KymoKnotAnalysisResult& knotResult = result.kymoResults[pieceIndex];
        std::cout << PBD::KymoKnotAnalyzer::describeAnalysisResult(knotResult, &cablePoints)
                  << std::endl;

        if (pieceIndex < result.scoredHoops.size() && !result.scoredHoops[pieceIndex].empty())
        {
            std::cout << "    hoop scores:";
            for (const TanglementAnalyzer::HoopCandidate& hoop : result.scoredHoops[pieceIndex])
            {
                if (hoop.pairSumScore == 0)
                    continue;

                std::cout << " ["
                          << std::min(hoop.startEdge, hoop.endEdge)
                          << "-"
                          << std::max(hoop.startEdge, hoop.endEdge)
                          << ": pairSumScore=" << hoop.pairSumScore
                          << "]";
            }
            std::cout << std::endl;
        }
    }
}
