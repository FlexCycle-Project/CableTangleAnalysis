#include "TanglementAnalyzer.h"

#include "KymoKnotAnalyzer.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>
#include <unordered_map>
#include <unordered_set>

namespace
{
    struct HitSample
    {
        // edge is a continuous edge-coordinate along the cable piece; state is
        // the over/under relation at that crossing. Hit lists use state changes
        // to infer whether a cable strand threads through a hoop.
        Real edge = static_cast<Real>(0.0);
        bool state = false;
    };

    struct HitListBuckets
    {
        std::vector<HitSample> decreasingSideHits;
        std::vector<HitSample> increasingSideHits;
        std::vector<HitSample> cableHits;
    };

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
        // Projection crossings are detected in the XZ plane; the Y coordinate
        // is used later only to decide over/under ordering.
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
        // Ray-casting containment test for projected hoop polygons. Points on
        // the boundary count as inside so borderline threadings remain visible.
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
        // Build the projected crossing diagram for all pieces. Same-piece
        // adjacent edges are skipped, but cross-piece edges are all considered.
        std::vector<TanglementAnalyzer::CrossingEvent> crossingEvents;
        if (ropePieces.empty())
            return crossingEvents;

        for (size_t firstPieceIndex = 0; firstPieceIndex < ropePieces.size(); ++firstPieceIndex)
        {
            const std::vector<Vector3r>& firstPiece = ropePieces[firstPieceIndex];
            if (firstPiece.size() < 2)
                continue;

            const int firstLastEdge = static_cast<int>(firstPiece.size()) - 2;
            for (size_t secondPieceIndex = firstPieceIndex; secondPieceIndex < ropePieces.size(); ++secondPieceIndex)
            {
                const std::vector<Vector3r>& secondPiece = ropePieces[secondPieceIndex];
                if (secondPiece.size() < 2)
                    continue;

                const int secondLastEdge = static_cast<int>(secondPiece.size()) - 2;
                for (int edgeA = 0; edgeA <= firstLastEdge; ++edgeA)
                {
                    const int secondEdgeStart =
                        firstPieceIndex == secondPieceIndex ? edgeA + 2 : 0;
                    for (int edgeB = secondEdgeStart; edgeB <= secondLastEdge; ++edgeB)
                    {
                        if (firstPieceIndex == secondPieceIndex && edgeB <= edgeA + 1)
                            continue;

                        const Vector2r a0(
                            firstPiece[static_cast<size_t>(edgeA)].x(),
                            firstPiece[static_cast<size_t>(edgeA)].z());
                        const Vector2r a1(
                            firstPiece[static_cast<size_t>(edgeA + 1)].x(),
                            firstPiece[static_cast<size_t>(edgeA + 1)].z());
                        const Vector2r b0(
                            secondPiece[static_cast<size_t>(edgeB)].x(),
                            secondPiece[static_cast<size_t>(edgeB)].z());
                        const Vector2r b1(
                            secondPiece[static_cast<size_t>(edgeB + 1)].x(),
                            secondPiece[static_cast<size_t>(edgeB + 1)].z());

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
                            firstPiece[static_cast<size_t>(edgeA)] +
                            ua * (firstPiece[static_cast<size_t>(edgeA + 1)] - firstPiece[static_cast<size_t>(edgeA)]);
                        const Vector3r secondPoint3D =
                            secondPiece[static_cast<size_t>(edgeB)] +
                            ub * (secondPiece[static_cast<size_t>(edgeB + 1)] - secondPiece[static_cast<size_t>(edgeB)]);

                        // Store both the projected crossing location in edge
                        // coordinates and the 3D height ordering.
                        TanglementAnalyzer::CrossingEvent event;
                        event.order = static_cast<Real>(0.5) *
                            (static_cast<Real>(edgeA) + static_cast<Real>(edgeB));
                        event.pieceIndex = firstPieceIndex;
                        event.otherPieceIndex = secondPieceIndex;
                        event.firstEdge = edgeA;
                        event.secondEdge = edgeB;
                        event.heightDelta = firstPoint3D.y() - secondPoint3D.y();
                        event.firstOverSecond = event.heightDelta > static_cast<Real>(0.0);
                        crossingEvents.push_back(event);
                    }
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
            if (event.pieceIndex != pieceIndex || event.otherPieceIndex != pieceIndex)
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

    static int scoreHoopAgainstCable(
        const std::vector<TanglementAnalyzer::CrossingEvent>& allCrossingEvents,
        size_t sourcePieceIndex,
        const std::vector<Vector3r>& sourcePiece,
        const TanglementAnalyzer::HoopCandidate& sourceHoop,
        size_t targetPieceIndex,
        const std::vector<Vector3r>& targetPiece,
        HitListBuckets* outHitLists = nullptr)
    {
        // A hoop is the subchain between two self-crossing edges. The score
        // counts how many strands enter/exit this projected loop with changing
        // over-under state, which is the heuristic signal for threading.
        const auto edgeInHoop = [](int edgeIndex, const TanglementAnalyzer::HoopCandidate& hoop) -> bool
        {
            const int hoopStart = std::min(hoop.startEdge, hoop.endEdge);
            const int hoopEnd = std::max(hoop.startEdge, hoop.endEdge);
            return edgeIndex >= hoopStart && edgeIndex <= hoopEnd;
        };

        std::vector<Vector2r> sourceHoopPolygon;
        const int sourceStart = std::min(sourceHoop.startEdge, sourceHoop.endEdge);
        const int sourceEnd = std::max(sourceHoop.startEdge, sourceHoop.endEdge);
        const int sourceCableStart = 0;
        const int sourceCableEnd = std::max(0, static_cast<int>(sourcePiece.size()) - 1);
        const int sourceLastVertex = std::min(sourceEnd + 1, static_cast<int>(sourcePiece.size()) - 1);
        sourceHoopPolygon.reserve(static_cast<size_t>(sourceLastVertex - sourceStart + 1));
        for (int vertexIndex = sourceStart; vertexIndex <= sourceLastVertex; ++vertexIndex)
        {
            const Vector3r& point3D = sourcePiece[static_cast<size_t>(vertexIndex)];
            sourceHoopPolygon.emplace_back(point3D.x(), point3D.z());
        }

        // Same-piece hits are split by which side of the hoop they approach
        // from. Cross-piece hits use one shared list because the other cable
        // has no meaningful "before/after this hoop" edge interval.
        const auto samplePointOnPiece = [](const std::vector<Vector3r>& piece, Real edgePos) -> Vector3r
        {
            if (piece.empty())
                return Vector3r::Zero();

            const Real clampedEdgePos = std::max(
                static_cast<Real>(0.0),
                std::min(edgePos, static_cast<Real>(piece.size() - 1)));
            const int lowerIndex = static_cast<int>(std::floor(clampedEdgePos));
            const int upperIndex = std::min(lowerIndex + 1, static_cast<int>(piece.size()) - 1);
            const Real alpha = clampedEdgePos - static_cast<Real>(lowerIndex);
            return (static_cast<Real>(1.0) - alpha) * piece[static_cast<size_t>(lowerIndex)] +
                alpha * piece[static_cast<size_t>(upperIndex)];
        };

        std::vector<HitSample> increasingSideHits;
        std::vector<HitSample> decreasingSideHits;
        std::vector<HitSample> cableHits;
        std::vector<HitSample> wholeCableHits;
        increasingSideHits.reserve(allCrossingEvents.size());
        decreasingSideHits.reserve(allCrossingEvents.size());
        cableHits.reserve(allCrossingEvents.size());
        wholeCableHits.reserve(allCrossingEvents.size());

        for (const TanglementAnalyzer::CrossingEvent& event : allCrossingEvents)
        {
            int sourceEdge = -1;
            int targetEdge = -1;
            bool sourceOverTarget = false;

            if (event.pieceIndex == sourcePieceIndex && event.otherPieceIndex == targetPieceIndex)
            {
                sourceEdge = event.firstEdge;
                targetEdge = event.secondEdge;
                sourceOverTarget = event.firstOverSecond;
            }
            else if (event.pieceIndex == targetPieceIndex && event.otherPieceIndex == sourcePieceIndex)
            {
                sourceEdge = event.secondEdge;
                targetEdge = event.firstEdge;
                sourceOverTarget = !event.firstOverSecond;
            }
            else
            {
                continue;
            }

            const bool sourceInHoop = edgeInHoop(sourceEdge, sourceHoop);
            const bool targetInHoop = edgeInHoop(targetEdge, sourceHoop);
            if (sourceInHoop == targetInHoop)
                continue;

            // Only crossings between the hoop interior and exterior can be
            // threadings. Crossings wholly inside or wholly outside the hoop
            // do not indicate a strand passing through the loop boundary.
            const int outsideEdge = sourceInHoop ? targetEdge : sourceEdge;
            const bool outsideEdgeOverInside = sourceInHoop ? sourceOverTarget : !sourceOverTarget;
            const bool outsideEdgeIsInsideSourceHoop =
                outsideEdge >= sourceStart && outsideEdge <= sourceEnd;
            if (outsideEdgeIsInsideSourceHoop)
                continue;

            const HitSample sample{
                static_cast<Real>(outsideEdge),
                outsideEdgeOverInside
            };

            if (targetPieceIndex == sourcePieceIndex)
            {
                if (outsideEdge >= sourceCableStart && outsideEdge < sourceStart)
                    decreasingSideHits.push_back(sample);
                else if (outsideEdge > sourceEnd && outsideEdge <= sourceCableEnd)
                    increasingSideHits.push_back(sample);
            }
            else
            {
                cableHits.push_back(sample);
            }
        }

        const auto scoreHitList = [&](std::vector<HitSample>& hitList, bool descending) -> int
        {
            std::sort(
                hitList.begin(),
                hitList.end(),
                [descending](const HitSample& lhs, const HitSample& rhs)
                {
                    return descending ? lhs.edge > rhs.edge : lhs.edge < rhs.edge;
                });

            int score = 0;
            for (size_t i = 0; i + 1 < hitList.size(); ++i)
            {
                if (hitList[i + 1].state == hitList[i].state)
                    continue;

                // A change in over/under state between consecutive hits is a
                // candidate threading. The midpoint test rejects changes whose
                // intervening strand does not pass through the hoop interior.
                const Real midpointEdge =
                    static_cast<Real>(0.5) * (hitList[i].edge + hitList[i + 1].edge);
                const Vector3r point3D = samplePointOnPiece(targetPiece, midpointEdge);
                const Vector2r midpoint2D(point3D.x(), point3D.z());
                if (pointInPolygon2D(sourceHoopPolygon, midpoint2D))
                    ++score;
            }
            return score;
        };

        int score = 0;
        const int decreasingScore = scoreHitList(decreasingSideHits, true);
        const int increasingScore = scoreHitList(increasingSideHits, false);
        const int targetCableScore = scoreHitList(cableHits, false);
        score += decreasingScore;
        score += increasingScore;
        score += targetCableScore;
        if (outHitLists != nullptr)
        {
            outHitLists->decreasingSideHits = std::move(decreasingSideHits);
            outHitLists->increasingSideHits = std::move(increasingSideHits);
            outHitLists->cableHits = std::move(cableHits);
        }
        return score;
    }

}

    void TanglementAnalyzer::printHoopScoreHitLists(
        const std::vector<std::vector<Vector3r>>& ropePieces,
        size_t sourcePieceIndex,
        const HoopCandidate& sourceHoop)
    {
    const std::vector<CrossingEvent> allCrossingEvents = findCrossings(ropePieces);
    if (sourcePieceIndex >= ropePieces.size())
        return;

    const std::vector<Vector3r>& sourcePiece = ropePieces[sourcePieceIndex];
    const auto cableLabel = [](size_t pieceIndex) -> std::string
    {
        return "cable " + std::to_string(pieceIndex + 1);
    };

    std::cout << "    hoop hitlists:" << std::endl;
    std::cout << "      source " << cableLabel(sourcePieceIndex)
              << " hoop [" << sourceHoop.startEdge << "," << sourceHoop.endEdge << "]"
              << std::endl;

    int totalScore = 0;
    for (size_t targetPieceIndex = 0; targetPieceIndex < ropePieces.size(); ++targetPieceIndex)
    {
        const std::vector<Vector3r>& targetPiece = ropePieces[targetPieceIndex];
        HitListBuckets hitLists;
        const int score = scoreHoopAgainstCable(
            allCrossingEvents,
            sourcePieceIndex,
            sourcePiece,
            sourceHoop,
            targetPieceIndex,
            targetPiece,
            &hitLists);
        totalScore += score;

        const auto printHits = [](const std::string& label, const std::vector<HitSample>& hits)
        {
            std::cout << "        " << label << ":";
            if (hits.empty())
            {
                std::cout << " empty";
            }
            else
            {
                for (const HitSample& hit : hits)
                    std::cout << " [" << hit.edge << "," << (hit.state ? 1 : 0) << "]";
            }
            std::cout << std::endl;
        };

        if (targetPieceIndex == sourcePieceIndex)
        {
            std::cout << "      source " << cableLabel(targetPieceIndex) << std::endl;
            printHits("ascendingHits", hitLists.increasingSideHits);
            printHits("descendingHits", hitLists.decreasingSideHits);
        }
        else
        {
            std::cout << "      target " << cableLabel(targetPieceIndex) << std::endl;
            printHits("hits", hitLists.cableHits);
        }
    }

    std::cout << "      total score=" << totalScore << std::endl;
}

namespace
{
    static std::vector<int> buildProblematicHoopMiddleHighlights(
        const std::vector<Vector3r>& piece,
        const std::vector<std::vector<Vector3r>>& ropePieces,
        const std::vector<TanglementAnalyzer::CrossingEvent>& allCrossingEvents,
        size_t pieceIndex,
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

        std::vector<TanglementAnalyzer::CrossingEvent> crossingEvents;
        crossingEvents.reserve(allCrossingEvents.size());
        for (const TanglementAnalyzer::CrossingEvent& event : allCrossingEvents)
        {
            if (event.pieceIndex == pieceIndex &&
                event.otherPieceIndex == pieceIndex)
            {
                crossingEvents.push_back(event);
            }
        }
        if (crossingEvents.empty())
            return middleEdges;
        std::vector<TanglementAnalyzer::HoopCandidate> hoops;
        hoops.reserve(crossingEvents.size());
        for (const TanglementAnalyzer::CrossingEvent& event : crossingEvents)
        {
            // Each same-piece crossing defines a candidate loop spanning the
            // two crossing edges. Nested hoops are kept; their scores can differ
            // because other strands may thread one loop but not the other.
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
        for (size_t sourceIndex = 0; sourceIndex < hoops.size(); ++sourceIndex)
        {
            const TanglementAnalyzer::HoopCandidate& sourceHoop = hoops[sourceIndex];
            for (size_t targetPieceIndex = 0; targetPieceIndex < ropePieces.size(); ++targetPieceIndex)
            {
                // Score the hoop against every piece, including itself. This
                // lets the same metric handle self-threading and multi-piece
                // cable threadings after cuts.
                const std::vector<Vector3r>& targetPiece = ropePieces[targetPieceIndex];
                if (targetPiece.size() < 2)
                    continue;
                pairSumScores[sourceIndex] += scoreHoopAgainstCable(
                    allCrossingEvents,
                    pieceIndex,
                    piece,
                    sourceHoop,
                    targetPieceIndex,
                    targetPiece);
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
            if (hoops[i].pairSumScore == 0)
                continue;

            const int midpointEdge = std::max(
                0,
                std::min(
                    static_cast<int>(
                        std::round(
                            static_cast<Real>(0.5) *
                            (static_cast<Real>(hoops[i].startEdge) +
                             static_cast<Real>(hoops[i].endEdge)))),
                    std::max(0, pieceEdgeCount - 1)));
            middleEdges.push_back(midpointEdge);
        }

        std::sort(middleEdges.begin(), middleEdges.end());
        middleEdges.erase(std::unique(middleEdges.begin(), middleEdges.end()), middleEdges.end());

        if (outTargetHoopMiddleEdges != nullptr)
        {
            for (size_t i = 0; i < hoops.size(); ++i)
            {
                if (hoops[i].pairSumScore == 0)
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

    static std::vector<Real> buildSegmentScoresForPiece(
        const std::vector<TanglementAnalyzer::HoopCandidate>& hoops,
        int pieceEdgeCount)
    {
        if (pieceEdgeCount <= 0)
            return {};

        std::vector<Real> segmentScores(static_cast<size_t>(pieceEdgeCount), static_cast<Real>(0.0));
        bool hasAnyScore = false;

        for (const TanglementAnalyzer::HoopCandidate& hoop : hoops)
        {
            if (hoop.pairSumScore == 0)
                continue;

            // Propagate each hoop score to its covered edge interval. This
            // converts loop-level threading evidence into an edge-level score
            // suitable for visualization or finer cut selection.
            const int hoopStart = std::max(0, std::min(hoop.startEdge, hoop.endEdge));
            const int hoopEnd = std::min(pieceEdgeCount - 1, std::max(hoop.startEdge, hoop.endEdge));
            for (int edgeIndex = hoopStart; edgeIndex <= hoopEnd; ++edgeIndex)
            {
                segmentScores[static_cast<size_t>(edgeIndex)] += static_cast<Real>(hoop.pairSumScore);
                hasAnyScore = true;
            }
        }

        return hasAnyScore ? segmentScores : std::vector<Real>();
    }

    struct SegmentRef
    {
        size_t pieceIndex = 0;
        size_t edgeIndex = 0;
        Vector3r p0 = Vector3r::Zero();
        Vector3r p1 = Vector3r::Zero();
    };

    struct Int3
    {
        int x, y, z;

        bool operator==(const Int3& other) const
        {
            return x == other.x && y == other.y && z == other.z;
        }
    };

    struct Int3Hash
    {
        std::size_t operator()(const Int3& c) const
        {
            const std::size_t h1 = std::hash<int>()(c.x);
            const std::size_t h2 = std::hash<int>()(c.y);
            const std::size_t h3 = std::hash<int>()(c.z);
            return h1 ^ (h2 * 73856093u) ^ (h3 * 19349663u);
        }
    };

    class DisjointSet
    {
    public:
        explicit DisjointSet(size_t n)
            : m_parent(n), m_rank(n, 0)
        {
            for (size_t i = 0; i < n; ++i)
                m_parent[i] = i;
        }

        size_t find(size_t x)
        {
            if (m_parent[x] == x)
                return x;
            m_parent[x] = find(m_parent[x]);
            return m_parent[x];
        }

        void unite(size_t a, size_t b)
        {
            a = find(a);
            b = find(b);
            if (a == b)
                return;

            if (m_rank[a] < m_rank[b])
                std::swap(a, b);

            m_parent[b] = a;
            if (m_rank[a] == m_rank[b])
                ++m_rank[a];
        }

    private:
        std::vector<size_t> m_parent;
        std::vector<unsigned char> m_rank;
    };

    static uint64_t makeOrderedSegmentPairKey(size_t a, size_t b)
    {
        const uint64_t lo = static_cast<uint64_t>(std::min(a, b));
        const uint64_t hi = static_cast<uint64_t>(std::max(a, b));
        return (lo << 32) | (hi & 0xffffffffu);
    }

    static std::vector<SegmentRef> collectCableSegments(
        const std::vector<std::vector<Vector3r>>& ropePieces)
    {
        std::vector<SegmentRef> segments;
        for (size_t pieceIndex = 0; pieceIndex < ropePieces.size(); ++pieceIndex)
        {
            const std::vector<Vector3r>& piece = ropePieces[pieceIndex];
            if (piece.size() < 2)
                continue;

            for (size_t edgeIndex = 0; edgeIndex + 1 < piece.size(); ++edgeIndex)
            {
                segments.push_back(
                    SegmentRef{
                        pieceIndex,
                        edgeIndex,
                        piece[edgeIndex],
                        piece[edgeIndex + 1] });
            }
        }

        return segments;
    }

    static std::vector<std::pair<size_t, size_t>> collectNearbySegmentPairsSpatialHash(
        const std::vector<SegmentRef>& segments,
        Real inflate,
        Real cellSize,
        const RopeConfig& ropeCfg,
        const SelfCollisionConfig& collisionCfg)
    {
        std::vector<std::pair<size_t, size_t>> pairs;
        if (segments.size() < 2)
            return pairs;

        // Metric broad phase, separate from physical collision. It uses a
        // larger influence radius to count near misses as geometric clutter.
        const Real eps = static_cast<Real>(1e-8);
        const Real safeCellSize = std::max(cellSize, std::max(inflate, eps));
        std::unordered_map<Int3, std::vector<size_t>, Int3Hash> grid;
        grid.reserve(segments.size() * 2);

        auto cellOf = [safeCellSize](const Vector3r& p) -> Int3
        {
            return Int3{
                static_cast<int>(std::floor(p.x() / safeCellSize)),
                static_cast<int>(std::floor(p.y() / safeCellSize)),
                static_cast<int>(std::floor(p.z() / safeCellSize))
            };
        };

        for (size_t segmentIndex = 0; segmentIndex < segments.size(); ++segmentIndex)
        {
            const SegmentRef& segment = segments[segmentIndex];
            Vector3r bmin = segment.p0.cwiseMin(segment.p1);
            Vector3r bmax = segment.p0.cwiseMax(segment.p1);
            bmin.array() -= inflate;
            bmax.array() += inflate;

            const Int3 cmin = cellOf(bmin);
            const Int3 cmax = cellOf(bmax);
            for (int x = cmin.x; x <= cmax.x; ++x)
            {
                for (int y = cmin.y; y <= cmax.y; ++y)
                {
                    for (int z = cmin.z; z <= cmax.z; ++z)
                        grid[Int3{x, y, z}].push_back(segmentIndex);
                }
            }
        }

        std::unordered_set<uint64_t> uniquePairs;
        uniquePairs.reserve(segments.size() * 8);

        for (const auto& entry : grid)
        {
            const std::vector<size_t>& cellSegments = entry.second;
            for (size_t a = 0; a < cellSegments.size(); ++a)
            {
                for (size_t b = a + 1; b < cellSegments.size(); ++b)
                {
                    const size_t first = cellSegments[a];
                    const size_t second = cellSegments[b];
                    if (first == second)
                        continue;

                    const SegmentRef& segA = segments[first];
                    const SegmentRef& segB = segments[second];
                    if (segA.pieceIndex == segB.pieceIndex)
                    {
                        const Real arcDist = static_cast<Real>(
                            std::abs(static_cast<int>(segA.edgeIndex) - static_cast<int>(segB.edgeIndex))) *
                            ropeCfg.dx;
                        if (arcDist < collisionCfg.excludeArcLength)
                            continue;
                    }

                    uniquePairs.insert(makeOrderedSegmentPairKey(first, second));
                }
            }
        }

        pairs.reserve(uniquePairs.size());
        for (uint64_t key : uniquePairs)
        {
            pairs.emplace_back(
                static_cast<size_t>(key >> 32),
                static_cast<size_t>(key & 0xffffffffu));
        }

        return pairs;
    }

}

namespace
{
    static TanglementAnalyzer::TanglementMetricsSnapshot buildTanglementMetricsSnapshot(
        const std::vector<TanglementAnalyzer::CrossingEvent>& allCrossingEvents,
        const std::vector<std::vector<Vector3r>>& ropePieces,
        const std::vector<SelfCollisionContact>& contacts,
        const std::vector<PBD::KymoKnotAnalysisResult>& kymoResults,
        const RopeConfig& ropeCfg,
        const SelfCollisionConfig& collisionCfg)
    {
        // This snapshot combines topology-adjacent information (crossings,
        // KymoKnot span) with geometric density information (close segment
        // pairs, contact participation, connected close-pair clusters).
        TanglementAnalyzer::TanglementMetricsSnapshot snapshot;
        snapshot.pieceCount = ropePieces.size();
        snapshot.selfCollisionCount = contacts.size();

        const std::vector<SegmentRef> segments = collectCableSegments(ropePieces);
        snapshot.totalSegmentCount = segments.size();

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

        // Proximity metrics use a soft penalty inside the influence distance.
        // Squaring the normalized deficit emphasizes very close segment pairs.
        const Real influenceDistance =
            std::max(static_cast<Real>(2.0) * collisionCfg.ropeRadius,
                     static_cast<Real>(4.0) * collisionCfg.ropeRadius);
        const std::vector<std::pair<size_t, size_t>> candidatePairs =
            collectNearbySegmentPairsSpatialHash(
                segments,
                influenceDistance,
                collisionCfg.hashCellSize,
                ropeCfg,
                collisionCfg);

        DisjointSet componentSets(segments.size());
        std::vector<char> segmentTouched(segments.size(), 0);
        Real samePieceProximitySum = static_cast<Real>(0.0);
        Real crossPieceProximitySum = static_cast<Real>(0.0);
        Real minimumPairDistance = std::numeric_limits<Real>::infinity();
        size_t samePieceClosePairCount = 0;
        size_t crossPieceClosePairCount = 0;

        for (const std::pair<size_t, size_t>& pair : candidatePairs)
        {
            const size_t first = pair.first;
            const size_t second = pair.second;
            const SegmentRef& segA = segments[first];
            const SegmentRef& segB = segments[second];

            Real s = static_cast<Real>(0.0);
            Real t = static_cast<Real>(0.0);
            Vector3r c1 = Vector3r::Zero();
            Vector3r c2 = Vector3r::Zero();
            closestPtSegmentSegment(segA.p0, segA.p1, segB.p0, segB.p1, s, t, c1, c2);

            const Real dist = (c2 - c1).norm();
            minimumPairDistance = std::min(minimumPairDistance, dist);
            if (dist >= influenceDistance)
                continue;

            const Real normalizedDeficit =
                (influenceDistance - dist) /
                std::max(influenceDistance, static_cast<Real>(1e-12));
            const Real pairPenalty = normalizedDeficit * normalizedDeficit;

            if (segA.pieceIndex == segB.pieceIndex)
            {
                ++samePieceClosePairCount;
                samePieceProximitySum += pairPenalty;
            }
            else
            {
                ++crossPieceClosePairCount;
                crossPieceProximitySum += pairPenalty;
            }

            segmentTouched[first] = 1;
            segmentTouched[second] = 1;
            componentSets.unite(first, second);
        }

        snapshot.minimumPairDistance = minimumPairDistance;
        snapshot.closePairCount = samePieceClosePairCount + crossPieceClosePairCount;
        snapshot.samePieceClosePairCount = samePieceClosePairCount;
        snapshot.crossPieceClosePairCount = crossPieceClosePairCount;
        snapshot.closeSegmentCount = static_cast<size_t>(
            std::count(segmentTouched.begin(), segmentTouched.end(), static_cast<char>(1)));
        snapshot.samePieceProximityScore = segments.empty()
            ? static_cast<Real>(0.0)
            : samePieceProximitySum / static_cast<Real>(segments.size());
        snapshot.crossPieceProximityScore = segments.empty()
            ? static_cast<Real>(0.0)
            : crossPieceProximitySum / static_cast<Real>(segments.size());

        struct ComponentStats
        {
            size_t segmentCount = 0;
            std::unordered_set<size_t> pieceIndices;
        };

        std::unordered_map<size_t, ComponentStats> components;
        for (size_t segmentIndex = 0; segmentIndex < segments.size(); ++segmentIndex)
        {
            if (!segmentTouched[segmentIndex])
                continue;

            const size_t root = componentSets.find(segmentIndex);
            components[root].segmentCount += 1;
            components[root].pieceIndices.insert(segments[segmentIndex].pieceIndex);
        }

        Real clusterPenaltySum = static_cast<Real>(0.0);
        for (const auto& entry : components)
        {
            // Only multi-piece close-pair components contribute to the cluster
            // penalty; single-piece folding is counted separately as same-piece
            // proximity and is less important for inter-cable separation.
            const size_t segmentCount = entry.second.segmentCount;
            const size_t pieceCount = entry.second.pieceIndices.size();
            if (segmentCount < 2 || pieceCount < 2)
                continue;

            ++snapshot.entangledComponentCount;
            ++snapshot.multiPieceComponentCount;
            clusterPenaltySum += static_cast<Real>(segmentCount - 1);
        }

        snapshot.clusterPenalty = segments.empty()
            ? static_cast<Real>(0.0)
            : clusterPenaltySum / static_cast<Real>(segments.size());
        const Real samePieceWeight = static_cast<Real>(0.1);
        const Real componentWeight = static_cast<Real>(0.25);
        snapshot.entanglementScore =
            snapshot.crossPieceProximityScore +
            samePieceWeight * snapshot.samePieceProximityScore +
            componentWeight * snapshot.clusterPenalty;
        snapshot.proximityScore = snapshot.entanglementScore;

        for (size_t pieceIndex = 0; pieceIndex < ropePieces.size(); ++pieceIndex)
        {
            // KymoKnot reports localized knot spans per piece. The catalog uses
            // the average longest span fraction over pieces that contain knots.
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

namespace
{
    static TanglementAnalyzer::TanglementAnalysisResult analyzeTanglementInternal(
        const std::vector<std::vector<Vector3r>>& ropePieces,
        const std::vector<SelfCollisionContact>& contacts,
        const RopeConfig& ropeCfg,
        const SelfCollisionConfig& collisionCfg,
        PBD::KymoKnotSearchMode searchMode,
        bool runKymoKnot)
{
    using AnalysisResult = TanglementAnalyzer::TanglementAnalysisResult;
    using CrossingEvent = TanglementAnalyzer::CrossingEvent;

    AnalysisResult result;
    result.knotRegions.resize(ropePieces.size());
    result.crossingEdges.resize(ropePieces.size());
    result.crossingCenterEdges.resize(ropePieces.size());
    result.crossingMiddleEdges.resize(ropePieces.size());
    result.scoredHoops.resize(ropePieces.size());
    result.problematicHoopMiddleEdges.resize(ropePieces.size());
    result.targetHoopMiddleEdges.resize(ropePieces.size());
    result.segmentScores.resize(ropePieces.size());
    result.kymoResults.reserve(ropePieces.size());

    const std::vector<CrossingEvent> allCrossingEvents = findCrossings(ropePieces);

    if (runKymoKnot)
    {
        // Full analysis path: run KymoKnot to localize actual knot cores. The
        // geometry-only path skips this expensive step for evaluation loops
        // where only projected threading counts are needed.
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
                    TanglementAnalyzer::LocalKnotRegion{
                        localStart,
                        localEnd,
                        static_cast<size_t>(std::max(1, segment.length))});
                break;
            }
        }
    }

    for (size_t pieceIndex = 0; pieceIndex < ropePieces.size(); ++pieceIndex)
    {
        // Crossing/hoop features are computed from the projected geometry for
        // both full and geometry-only analysis.
        const std::vector<Vector3r>& piece = ropePieces[pieceIndex];
        std::vector<CrossingEvent> pieceCrossings;
        pieceCrossings.reserve(allCrossingEvents.size());
        for (const CrossingEvent& event : allCrossingEvents)
        {
            if (event.pieceIndex == pieceIndex &&
                event.otherPieceIndex == pieceIndex)
            {
                pieceCrossings.push_back(event);
            }
        }
        if (!pieceCrossings.empty())
        {
            std::vector<int> edges;
            std::vector<int> centerEdges;
            std::vector<int> middleEdges;

            for (const CrossingEvent& event : pieceCrossings)
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
                ropePieces,
                allCrossingEvents,
                pieceIndex,
                static_cast<int>(piece.size() > 1 ? piece.size() - 1 : 0),
                &result.scoredHoops[pieceIndex],
                &result.targetHoopMiddleEdges[pieceIndex]);

        result.segmentScores[pieceIndex] = buildSegmentScoresForPiece(
            result.scoredHoops[pieceIndex],
            static_cast<int>(piece.size() > 1 ? piece.size() - 1 : 0));
    }

    result.metrics = buildTanglementMetricsSnapshot(
        allCrossingEvents,
        ropePieces,
        contacts,
        result.kymoResults,
        ropeCfg,
        collisionCfg);
    result.allCrossingEvents = allCrossingEvents;
    return result;
}
}

TanglementAnalyzer::TanglementAnalysisResult TanglementAnalyzer::analyzeTanglement(
    const std::vector<std::vector<Vector3r>>& ropePieces,
    const std::vector<SelfCollisionContact>& contacts,
    const RopeConfig& ropeCfg,
    const SelfCollisionConfig& collisionCfg,
    PBD::KymoKnotSearchMode searchMode)
{
    return analyzeTanglementInternal(
        ropePieces,
        contacts,
        ropeCfg,
        collisionCfg,
        searchMode,
        true);
}

TanglementAnalyzer::TanglementAnalysisResult TanglementAnalyzer::analyzeTanglementGeometryOnly(
    const std::vector<std::vector<Vector3r>>& ropePieces,
    const std::vector<SelfCollisionContact>& contacts,
    const RopeConfig& ropeCfg,
    const SelfCollisionConfig& collisionCfg)
{
    return analyzeTanglementInternal(
        ropePieces,
        contacts,
        ropeCfg,
        collisionCfg,
        PBD::KymoKnotSearchMode::TopDown,
        false);
}

TanglementAnalyzer::TanglementMetricsSnapshot TanglementAnalyzer::evaluateContinuousEntanglementMetrics(
    const std::vector<std::vector<Vector3r>>& ropePieces,
    const std::vector<SelfCollisionContact>& contacts,
    const RopeConfig& ropeCfg,
    const SelfCollisionConfig& collisionCfg)
{
    const std::vector<CrossingEvent> noCrossings;
    const std::vector<PBD::KymoKnotAnalysisResult> noKymoResults;
    return buildTanglementMetricsSnapshot(
        noCrossings,
        ropePieces,
        contacts,
        noKymoResults,
        ropeCfg,
        collisionCfg);
}

TanglementAnalyzer::ThreadingMetricsSnapshot TanglementAnalyzer::computeThreadingMetrics(
    const TanglementAnalysisResult& result)
{
    // Threading count is the sum of positive hoop scores. It is intentionally
    // not a topological invariant; it measures how many scored loop-threading
    // relationships remain in the current projected geometry.
    ThreadingMetricsSnapshot snapshot;
    Real totalPositiveScore = static_cast<Real>(0.0);

    for (const std::vector<HoopCandidate>& pieceHoops : result.scoredHoops)
    {
        for (const HoopCandidate& hoop : pieceHoops)
        {
            if (hoop.pairSumScore <= 0)
                continue;

            ++snapshot.threadedHoopCount;
            snapshot.totalThreadingCount += hoop.pairSumScore;
            snapshot.maxHoopScore = std::max(snapshot.maxHoopScore, hoop.pairSumScore);
            totalPositiveScore += static_cast<Real>(hoop.pairSumScore);
        }
    }

    if (snapshot.threadedHoopCount > 0)
        snapshot.averageHoopScore = totalPositiveScore / static_cast<Real>(snapshot.threadedHoopCount);

    return snapshot;
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
    const ThreadingMetricsSnapshot threading = computeThreadingMetrics(result);
    std::cout << "  global: pieces=" << result.metrics.pieceCount
              << " points=" << result.metrics.totalPointCount
              << " segments=" << result.metrics.totalSegmentCount
              << " selfCollisions=" << result.metrics.selfCollisionCount
              << " (" << result.metrics.contactParticipationCount << ")"
              << " crossings=" << result.metrics.gravityCrossings
              << std::endl;
    std::cout << "  global: entanglementScore=" << result.metrics.entanglementScore
              << " proximity=" << result.metrics.proximityScore
              << " cross=" << result.metrics.crossPieceProximityScore
              << " same=" << result.metrics.samePieceProximityScore
              << " clusterPenalty=" << result.metrics.clusterPenalty
              << " closePairs=" << result.metrics.closePairCount
              << " [same=" << result.metrics.samePieceClosePairCount
              << ", cross=" << result.metrics.crossPieceClosePairCount << "]"
              << " closeSegments=" << result.metrics.closeSegmentCount
              << " components=" << result.metrics.entangledComponentCount
              << " multiPieceComponents=" << result.metrics.multiPieceComponentCount
              << std::endl;
    std::cout << "  global: minimumPairDistance=";
    if (std::isfinite(result.metrics.minimumPairDistance))
        std::cout << result.metrics.minimumPairDistance;
    else
        std::cout << "n/a";
    std::cout << std::endl;
    std::cout << "  threading: total=" << threading.totalThreadingCount
              << " threadedHoops=" << threading.threadedHoopCount
              << " max=" << threading.maxHoopScore
              << " average=" << threading.averageHoopScore
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
