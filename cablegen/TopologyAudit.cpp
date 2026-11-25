#include "TopologyAudit.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace
{
    struct EdgeRef
    {
        size_t pieceIndex = 0;
        int edgeIndex = -1;
        Vector3r p0 = Vector3r::Zero();
        Vector3r p1 = Vector3r::Zero();
        Real minX = static_cast<Real>(0.0);
        Real maxX = static_cast<Real>(0.0);
        Real minZ = static_cast<Real>(0.0);
        Real maxZ = static_cast<Real>(0.0);
    };

    using CrossingKey = std::tuple<size_t, int, size_t, int>;

    static CrossingKey makeCrossingKey(const TopologyAudit::Crossing& crossing)
    {
        return std::make_tuple(
            crossing.pieceA,
            crossing.edgeA,
            crossing.pieceB,
            crossing.edgeB);
    }

    static std::uint64_t packCell(int cx, int cz)
    {
        return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(cx)) << 32) |
               static_cast<std::uint64_t>(static_cast<std::uint32_t>(cz));
    }

    static std::uint64_t packPair(size_t a, size_t b)
    {
        const size_t lo = std::min(a, b);
        const size_t hi = std::max(a, b);
        return (static_cast<std::uint64_t>(lo) << 32) |
               static_cast<std::uint64_t>(hi & 0xffffffffu);
    }

    // Intersection of two XZ-projected segments. Endpoint intersections are
    // rejected to avoid double-counting a crossing at a shared particle.
    static bool intersectProjectedSegments(
        const Vector2r& a0,
        const Vector2r& a1,
        const Vector2r& b0,
        const Vector2r& b1,
        Real& ua,
        Real& ub,
        Real parallelEpsilon,
        Real endpointEpsilon)
    {
        const Real ax = a1.x() - a0.x();
        const Real az = a1.y() - a0.y();
        const Real bx = b1.x() - b0.x();
        const Real bz = b1.y() - b0.y();
        const Real denom = ax * bz - az * bx;
        if (std::abs(denom) <= parallelEpsilon)
            return false;

        const Real dx = b0.x() - a0.x();
        const Real dz = b0.y() - a0.y();
        ua = (dx * bz - dz * bx) / denom;
        ub = (dx * az - dz * ax) / denom;

        return ua > endpointEpsilon && ua < static_cast<Real>(1.0) - endpointEpsilon &&
               ub > endpointEpsilon && ub < static_cast<Real>(1.0) - endpointEpsilon;
    }

    static std::vector<EdgeRef> collectEdges(
        const std::vector<std::vector<Vector3r>>& pieces)
    {
        std::vector<EdgeRef> edges;
        for (size_t pieceIndex = 0; pieceIndex < pieces.size(); ++pieceIndex)
        {
            const std::vector<Vector3r>& piece = pieces[pieceIndex];
            if (piece.size() < 2)
                continue;

            for (size_t edgeIndex = 0; edgeIndex + 1 < piece.size(); ++edgeIndex)
            {
                EdgeRef edge;
                edge.pieceIndex = pieceIndex;
                edge.edgeIndex = static_cast<int>(edgeIndex);
                edge.p0 = piece[edgeIndex];
                edge.p1 = piece[edgeIndex + 1];
                edge.minX = std::min(edge.p0.x(), edge.p1.x());
                edge.maxX = std::max(edge.p0.x(), edge.p1.x());
                edge.minZ = std::min(edge.p0.z(), edge.p1.z());
                edge.maxZ = std::max(edge.p0.z(), edge.p1.z());
                edges.push_back(edge);
            }
        }
        return edges;
    }

    static std::vector<std::pair<size_t, size_t>> broadPhaseProjectedPairs(
        const std::vector<EdgeRef>& edges,
        Real cellSize)
    {
        std::vector<std::pair<size_t, size_t>> pairs;
        if (edges.size() < 2)
            return pairs;

        const Real safeCell = std::max(cellSize, static_cast<Real>(1e-6));
        std::unordered_map<std::uint64_t, std::vector<size_t>> grid;
        grid.reserve(edges.size() * 2);

        for (size_t edgeIndex = 0; edgeIndex < edges.size(); ++edgeIndex)
        {
            const EdgeRef& edge = edges[edgeIndex];
            const int minCx = static_cast<int>(std::floor(edge.minX / safeCell));
            const int maxCx = static_cast<int>(std::floor(edge.maxX / safeCell));
            const int minCz = static_cast<int>(std::floor(edge.minZ / safeCell));
            const int maxCz = static_cast<int>(std::floor(edge.maxZ / safeCell));

            for (int cx = minCx; cx <= maxCx; ++cx)
            {
                for (int cz = minCz; cz <= maxCz; ++cz)
                    grid[packCell(cx, cz)].push_back(edgeIndex);
            }
        }

        std::unordered_set<std::uint64_t> uniquePairs;
        uniquePairs.reserve(edges.size() * 6);
        for (const auto& entry : grid)
        {
            const std::vector<size_t>& bucket = entry.second;
            for (size_t a = 0; a < bucket.size(); ++a)
            {
                for (size_t b = a + 1; b < bucket.size(); ++b)
                {
                    const size_t first = bucket[a];
                    const size_t second = bucket[b];
                    if (first == second)
                        continue;

                    const EdgeRef& edgeA = edges[first];
                    const EdgeRef& edgeB = edges[second];
                    if (edgeA.pieceIndex == edgeB.pieceIndex &&
                        std::abs(edgeA.edgeIndex - edgeB.edgeIndex) <= 1)
                    {
                        continue;
                    }

                    uniquePairs.insert(packPair(first, second));
                }
            }
        }

        pairs.reserve(uniquePairs.size());
        for (std::uint64_t key : uniquePairs)
        {
            pairs.emplace_back(
                static_cast<size_t>(key >> 32),
                static_cast<size_t>(key & 0xffffffffu));
        }
        return pairs;
    }
}

namespace TopologyAudit
{
    std::vector<Crossing> computeCrossings(
        const std::vector<std::vector<Vector3r>>& pieces,
        const Settings& settings)
    {
        std::vector<Crossing> crossings;
        const std::vector<EdgeRef> edges = collectEdges(pieces);
        if (edges.size() < 2)
            return crossings;

        const std::vector<std::pair<size_t, size_t>> pairs =
            broadPhaseProjectedPairs(edges, settings.cellSize);

        crossings.reserve(pairs.size() / 2);
        for (const std::pair<size_t, size_t>& pair : pairs)
        {
            const EdgeRef& edgeA = edges[pair.first];
            const EdgeRef& edgeB = edges[pair.second];

            const Vector2r a0(edgeA.p0.x(), edgeA.p0.z());
            const Vector2r a1(edgeA.p1.x(), edgeA.p1.z());
            const Vector2r b0(edgeB.p0.x(), edgeB.p0.z());
            const Vector2r b1(edgeB.p1.x(), edgeB.p1.z());

            Real ua = static_cast<Real>(0.0);
            Real ub = static_cast<Real>(0.0);
            if (!intersectProjectedSegments(
                    a0,
                    a1,
                    b0,
                    b1,
                    ua,
                    ub,
                    settings.parallelEpsilon,
                    settings.endpointEpsilon))
            {
                continue;
            }

            const Vector3r firstPoint =
                edgeA.p0 + ua * (edgeA.p1 - edgeA.p0);
            const Vector3r secondPoint =
                edgeB.p0 + ub * (edgeB.p1 - edgeB.p0);

            Crossing crossing;
            crossing.pieceA = edgeA.pieceIndex;
            crossing.edgeA = edgeA.edgeIndex;
            crossing.pieceB = edgeB.pieceIndex;
            crossing.edgeB = edgeB.edgeIndex;
            crossing.x = static_cast<Real>(0.5) * (firstPoint.x() + secondPoint.x());
            crossing.z = static_cast<Real>(0.5) * (firstPoint.z() + secondPoint.z());
            crossing.heightDelta = firstPoint.y() - secondPoint.y();
            crossing.firstOverSecond = crossing.heightDelta > static_cast<Real>(0.0);
            crossings.push_back(crossing);
        }

        std::sort(
            crossings.begin(),
            crossings.end(),
            [](const Crossing& a, const Crossing& b)
            {
                return makeCrossingKey(a) < makeCrossingKey(b);
            });

        return crossings;
    }

    void Auditor::reset()
    {
        previous_.clear();
        hasPrevious_ = false;
    }

    void Auditor::begin(
        const std::vector<std::vector<Vector3r>>& pieces,
        const Settings& settings)
    {
        settings_ = settings;
        previous_ = computeCrossings(pieces, settings_);
        hasPrevious_ = true;
    }

    StepResult Auditor::audit(
        unsigned int step,
        const std::vector<std::vector<Vector3r>>& pieces)
    {
        StepResult result;
        const std::vector<Crossing> current = computeCrossings(pieces, settings_);
        result.crossingCount = current.size();

        if (!hasPrevious_)
        {
            previous_ = current;
            hasPrevious_ = true;

            for (const Crossing& crossing : current)
            {
                const Real absHeight = std::abs(crossing.heightDelta);
                result.minAbsHeight = std::min(result.minAbsHeight, absHeight);
                if (absHeight < settings_.contactEpsilon)
                    ++result.centerlineContactCount;
                if (absHeight < static_cast<Real>(2.0) * settings_.surfaceRadius)
                    ++result.surfaceOverlapCount;
            }
            return result;
        }

        std::map<CrossingKey, Crossing> previousByKey;
        std::map<CrossingKey, Crossing> currentByKey;
        for (const Crossing& crossing : previous_)
            previousByKey[makeCrossingKey(crossing)] = crossing;
        for (const Crossing& crossing : current)
            currentByKey[makeCrossingKey(crossing)] = crossing;

        for (const auto& entry : currentByKey)
        {
            const Crossing& crossing = entry.second;
            const Real absHeight = std::abs(crossing.heightDelta);
            result.minAbsHeight = std::min(result.minAbsHeight, absHeight);

            if (absHeight < settings_.contactEpsilon)
                ++result.centerlineContactCount;
            if (absHeight < static_cast<Real>(2.0) * settings_.surfaceRadius)
                ++result.surfaceOverlapCount;

            const auto previousIt = previousByKey.find(entry.first);
            if (previousIt == previousByKey.end())
            {
                // A crossing born with (almost) no vertical separation is a
                // topology touch/pass-through, not a Reidemeister crossing
                // created by two strands moving apart.
                if (absHeight < settings_.contactEpsilon)
                {
                    Violation violation;
                    violation.type = "BIRTH_AT_CONTACT";
                    violation.step = step;
                    violation.current = crossing;
                    violation.hasCurrent = true;
                    violation.currentHeight = crossing.heightDelta;
                    result.violations.push_back(violation);
                }
                continue;
            }

            const Real previousHeight = previousIt->second.heightDelta;
            if ((previousHeight > static_cast<Real>(0.0) &&
                 crossing.heightDelta < static_cast<Real>(0.0)) ||
                (previousHeight < static_cast<Real>(0.0) &&
                 crossing.heightDelta > static_cast<Real>(0.0)))
            {
                Violation violation;
                violation.type = "CROSSING_FLIP";
                violation.step = step;
                violation.previous = previousIt->second;
                violation.current = crossing;
                violation.hasPrevious = true;
                violation.hasCurrent = true;
                violation.previousHeight = previousHeight;
                violation.currentHeight = crossing.heightDelta;
                result.violations.push_back(violation);
            }
            else if (absHeight < settings_.contactEpsilon)
            {
                Violation violation;
                violation.type = "CENTERLINE_CONTACT";
                violation.step = step;
                violation.previous = previousIt->second;
                violation.current = crossing;
                violation.hasPrevious = true;
                violation.hasCurrent = true;
                violation.previousHeight = previousHeight;
                violation.currentHeight = crossing.heightDelta;
                result.violations.push_back(violation);
            }
        }

        for (const auto& entry : previousByKey)
        {
            if (currentByKey.find(entry.first) != currentByKey.end())
                continue;

            const Real absHeight = std::abs(entry.second.heightDelta);
            if (absHeight < settings_.contactEpsilon)
            {
                Violation violation;
                violation.type = "DEATH_AT_CONTACT";
                violation.step = step;
                violation.previous = entry.second;
                violation.hasPrevious = true;
                violation.previousHeight = entry.second.heightDelta;
                result.violations.push_back(violation);
            }
        }

        previous_ = current;
        return result;
    }
}
