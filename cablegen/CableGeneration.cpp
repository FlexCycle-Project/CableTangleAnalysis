#include "CableGeneration.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>

namespace CableGeneration
{
std::vector<Vector3r> generateRandomCableMessPolyline(
    const CableLayoutConfig& layoutCfg,
    const RopeConfig& ropeCfg,
    const RandomCableConfig& randomCfg,
    Real groundY)
{
    const Real minHeight = groundY + randomCfg.minHeightAboveGround;
    const Vector3r startPoint(
        -layoutCfg.endpointX,
        minHeight,
        static_cast<Real>(0.0));
    const Vector3r endPoint(
        layoutCfg.endpointX,
        minHeight,
        static_cast<Real>(0.0));

    std::random_device rd;
    const unsigned int seed = (randomCfg.seed != 0) ? randomCfg.seed : rd();
    std::mt19937 rng(seed);
    std::uniform_real_distribution<Real> unit01(
        static_cast<Real>(0.0),
        static_cast<Real>(1.0));

    auto chooseTurn = [&]() -> int
    {
        const Real totalWeight =
            randomCfg.straightChance +
            randomCfg.leftChance +
            randomCfg.rightChance;
        if (totalWeight <= static_cast<Real>(1e-12))
            return 0;

        const Real roll = unit01(rng) * totalWeight;
        if (roll < randomCfg.straightChance)
            return 0;
        if (roll < randomCfg.straightChance + randomCfg.leftChance)
            return 1;
        return -1;
    };

    auto headingToDirection = [](int heading) -> Vector3r
    {
        switch (heading & 3)
        {
        case 0: return Vector3r( static_cast<Real>(1.0), static_cast<Real>(0.0), static_cast<Real>(0.0));
        case 1: return Vector3r( static_cast<Real>(0.0), static_cast<Real>(0.0), static_cast<Real>(1.0));
        case 2: return Vector3r(-static_cast<Real>(1.0), static_cast<Real>(0.0), static_cast<Real>(0.0));
        default: return Vector3r( static_cast<Real>(0.0), static_cast<Real>(0.0),-static_cast<Real>(1.0));
        }
    };

    const int controlCount = std::max(2, randomCfg.controlPointCount);
    std::vector<Vector3r> middlePoints;
    middlePoints.reserve(static_cast<size_t>(controlCount));

    const Vector3r middleBase =
        static_cast<Real>(0.5) * (startPoint + endPoint);
    int heading = 0;
    Vector3r current = middleBase;
    const Real stepLength = randomCfg.controlStepLength;

    for (int i = 0; i < controlCount; ++i)
    {
        const int turn = chooseTurn();
        if (turn > 0)
            heading = (heading + 1) & 3;
        else if (turn < 0)
            heading = (heading + 3) & 3;

        current += stepLength * headingToDirection(heading);
        current.z() -= randomCfg.centerPull * current.z();
        current.z() = std::max(
            -randomCfg.maxLateralOffset,
            std::min(randomCfg.maxLateralOffset, current.z()));
        current.y() = minHeight;
        middlePoints.push_back(current);
    }

    Vector3r centroid(0.0, 0.0, 0.0);
    for (const Vector3r& p : middlePoints)
        centroid += p;
    centroid /= static_cast<Real>(middlePoints.size());

    const Vector3r targetMiddle = middleBase;
    for (Vector3r& p : middlePoints)
    {
        p += targetMiddle - centroid;
        p.y() = minHeight;
    }

    std::vector<Vector3r> controlPoints;
    const int bridgeSamples = std::max(2, controlCount / 4);
    controlPoints.reserve(
        static_cast<size_t>(middlePoints.size()) + 2u * static_cast<size_t>(bridgeSamples) + 2u);

    auto appendBridge = [&](const Vector3r& from, const Vector3r& to, bool includeFirstPoint)
    {
        for (int i = includeFirstPoint ? 0 : 1; i <= bridgeSamples; ++i)
        {
            const Real t = static_cast<Real>(i) / static_cast<Real>(bridgeSamples);
            controlPoints.push_back((static_cast<Real>(1.0) - t) * from + t * to);
        }
    };

    appendBridge(startPoint, middlePoints.front(), true);
    for (size_t i = 0; i < middlePoints.size(); ++i)
    {
        if (i == 0)
            continue;
        controlPoints.push_back(middlePoints[i]);
    }
    appendBridge(middlePoints.back(), endPoint, false);

    smoothPolyline(
        controlPoints,
        randomCfg.smoothingIterations,
        randomCfg.smoothingAlpha);

    controlPoints.front() = startPoint;
    controlPoints.back() = endPoint;

    Real minY = std::numeric_limits<Real>::max();
    for (const Vector3r& p : controlPoints)
        minY = std::min(minY, p.y());
    if (minY < minHeight)
    {
        const Real lift = minHeight - minY;
        for (Vector3r& p : controlPoints)
            p.y() += lift;
        controlPoints.front() = startPoint;
        controlPoints.back() = endPoint;
    }

    const Real pathLength = [&controlPoints]() -> Real
    {
        Real total = static_cast<Real>(0.0);
        for (size_t i = 1; i < controlPoints.size(); ++i)
            total += (controlPoints[i] - controlPoints[i - 1]).norm();
        return total;
    }();

    const int particleCount = std::max(
        2,
        static_cast<int>(std::ceil(pathLength / ropeCfg.dx)) + 1);

    return resamplePolylineByLength(controlPoints, particleCount);
}

void smoothPolyline(
    std::vector<Vector3r>& points,
    int iterations,
    Real alpha)
{
    if (points.size() < 3 || iterations <= 0)
        return;

    alpha = clamp01(alpha);

    for (int iter = 0; iter < iterations; ++iter)
    {
        const std::vector<Vector3r> oldPoints = points;
        for (size_t i = 1; i + 1 < points.size(); ++i)
        {
            const Vector3r neighborAverage =
                static_cast<Real>(0.5) * (oldPoints[i - 1] + oldPoints[i + 1]);
            points[i] = oldPoints[i] + alpha * (neighborAverage - oldPoints[i]);
        }

        points.front() = oldPoints.front();
        points.back() = oldPoints.back();
    }
}

std::vector<Vector3r> resamplePolylineByLength(
    const std::vector<Vector3r>& polyline,
    int wantedCount)
{
    std::vector<Vector3r> result;
    if (polyline.size() < 2 || wantedCount < 2)
        return result;

    std::vector<Real> cumulativeLength(polyline.size(), static_cast<Real>(0.0));
    for (size_t i = 1; i < polyline.size(); ++i)
    {
        cumulativeLength[i] =
            cumulativeLength[i - 1] +
            (polyline[i] - polyline[i - 1]).norm();
    }

    const Real totalLength = cumulativeLength.back();
    if (totalLength < static_cast<Real>(1e-8))
    {
        result.assign(static_cast<size_t>(wantedCount), polyline.front());
        return result;
    }

    result.reserve(static_cast<size_t>(wantedCount));
    const Real spacing = totalLength / static_cast<Real>(wantedCount - 1);
    size_t segment = 0;

    for (int i = 0; i < wantedCount; ++i)
    {
        const Real targetLength =
            std::min(totalLength, spacing * static_cast<Real>(i));

        while (segment + 1 < cumulativeLength.size() &&
               cumulativeLength[segment + 1] < targetLength)
        {
            ++segment;
        }

        if (segment + 1 >= polyline.size())
        {
            result.push_back(polyline.back());
            continue;
        }

        const Real segmentLength =
            cumulativeLength[segment + 1] - cumulativeLength[segment];
        if (segmentLength < static_cast<Real>(1e-12))
        {
            result.push_back(polyline[segment]);
            continue;
        }

        const Real localT =
            (targetLength - cumulativeLength[segment]) / segmentLength;
        result.push_back(
            (static_cast<Real>(1.0) - localT) * polyline[segment] +
            localT * polyline[segment + 1]);
    }

    return result;
}

}
