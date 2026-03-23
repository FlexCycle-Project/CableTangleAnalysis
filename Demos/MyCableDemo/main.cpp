#include "Common/Common.h"
#include "Demos/Visualization/MiniGL.h"
#include "Simulation/TimeManager.h"
#include "Simulation/SimulationModel.h"
#include "Demos/Common/DemoBase.h"
#include "Simulation/Simulation.h"

#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <cstdint>
#include <limits>

#if defined(_DEBUG) && !defined(EIGEN_ALIGN)
#define new DEBUG_NEW
#endif

using namespace PBD;
using namespace Eigen;
using namespace std;

// Rod material, discretization, and rendering settings.
struct RopeConfig
{
    Real stretchStiffness = static_cast<Real>(1.0);
    Real shearStiffnessX = static_cast<Real>(1.0);
    Real shearStiffnessY = static_cast<Real>(1.0);
    Real bendStiffnessX = static_cast<Real>(0.005);
    Real bendStiffnessY = static_cast<Real>(0.005);
    Real twistStiffness = static_cast<Real>(0.02);
    Real dx = static_cast<Real>(0.05);
    Real renderRadius = static_cast<Real>(0.06);
};

// Per-step simulation controls and endpoint pulling.
struct PhysicsConfig
{
    Real timeStep = static_cast<Real>(0.005);
    Real pullForce = static_cast<Real>(500.0);
    Real pullStep = static_cast<Real>(10.0);
    Real velocityRecomputeBlend = static_cast<Real>(0.05);
    bool usePullA = true;
    bool usePullB = true;
};

// Direct controls for the random sphere-mess mode.
struct RandomShapeConfig
{
    unsigned int seed = 4;
    unsigned int controlPointCount = 50;
    int smoothIterations = 1;
    Real smoothAlpha = static_cast<Real>(0.30);
    Real endpointX = static_cast<Real>(10.0);
    Real height = static_cast<Real>(6.0);
    Real sphereRadius = static_cast<Real>(2.0);
    Real generatorStepSize = static_cast<Real>(1);
};

// Tuning for the custom rope self-collision pass.
struct SelfCollisionConfig
{
    bool enabled = true;
    Real ropeRadius = static_cast<Real>(0.06);
    Real compliance = static_cast<Real>(1e-8);
    unsigned int iterations = 4;
    Real excludeArcLength = static_cast<Real>(0.9);
    Real normalDamping = static_cast<Real>(0.2);
    Real tangentialDamping = static_cast<Real>(0.08);
    Real maxCorrection = static_cast<Real>(0.015);
    Real hashCellSize = static_cast<Real>(0.25);
};

// Single source of truth for the demo.
struct AppConfig
{
    RopeConfig rope;
    PhysicsConfig physics;
    RandomShapeConfig random;
    Real groundY = static_cast<Real>(0.0);
    SelfCollisionConfig selfCollision;
};

// Integer grid coordinate for the self-collision spatial hash.
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

// Velocity-side damping data saved during self-collision resolution.
struct VelocityPairCorrection
{
    unsigned int i0, i1, j0, j1;
    Real wP0, wP1, wQ0, wQ1;
    Vector3r nrm;
};

// Derived layout for random mode before we build the actual polyline.
struct RandomCableLayout
{
    Vector3r startPoint;
    Vector3r endPoint;
    Vector3r sphereEntry;
    Vector3r sphereExit;
};

static void buildModel();
static void reset();
static void timeStep();
static void render();

static std::vector<Vector3r> generateRandomCablePolyline(
    const RandomShapeConfig& randomCfg,
    const RopeConfig& ropeCfg);

static RandomCableLayout makeRandomCableLayout(const RandomShapeConfig& randomCfg);

static void addCableModel(
    SimulationModel* model,
    std::vector<Vector3r>& ropePoints,
    const RopeConfig& ropeCfg);

static std::vector<Vector3r> generateRandomPolylineInSphere(
    const Vector3r& start,
    const Vector3r& end,
    int controlCount,
    Real sphereRadius,
    Real generatorStepSize,
    unsigned int seed);

static void smoothPolyline(
    std::vector<Vector3r>& points,
    int iterations,
    Real alpha);

static std::vector<Vector3r> resamplePolylineByLength(
    const std::vector<Vector3r>& polyline,
    int wantedCount);

static void applyGroundConstraint(Real groundY);
static void applyPullForces(SimulationModel* model, Real dt, const PhysicsConfig& physicsCfg);
static void applySelfCollisionXPBD(
    SimulationModel* model,
    Real dt,
    const RopeConfig& ropeCfg,
    const SelfCollisionConfig& collisionCfg);

static Real clamp01(Real x);
static Real getInverseMass(const ParticleData& pd, unsigned int index);
static Vector3r clampVec(const Vector3r& v, Real maxNorm);
static uint64_t makeOrderedEdgePairKey(unsigned int a, unsigned int b);

static void closestPtSegmentSegment(
    const Vector3r& p1, const Vector3r& q1,
    const Vector3r& p2, const Vector3r& q2,
    Real& s, Real& t,
    Vector3r& c1, Vector3r& c2);

static void recomputeVelocitiesFromPositions(
    SimulationModel* model,
    const std::vector<Vector3r>& oldPositions,
    Real dt,
    Real blend);

static std::vector<std::pair<unsigned int, unsigned int>>
collectSelfCollisionCandidatesSpatialHash(
    const ParticleData& particles,
    Real inflate,
    Real cellSize);

static DemoBase* base = nullptr;
static AppConfig g_cfg;
static unsigned int g_firstParticle = 0;
static unsigned int g_lastParticle = 0;

int main(int argc, char** argv)
{
    std::cout << "Random sphere mess demo" << std::endl;

    base = new DemoBase();

    base->init(argc, argv, "Random Sphere Mess Demo");

    SimulationModel* model = new SimulationModel();
    model->init();
    Simulation::getCurrent()->setModel(model);

    buildModel();

    MiniGL::setClientIdleFunc(timeStep);
    MiniGL::setClientSceneFunc(render);
    MiniGL::addKeyFunc('r', reset);
    MiniGL::addKeyFunc('1', []() { g_cfg.physics.usePullA = !g_cfg.physics.usePullA; });
    MiniGL::addKeyFunc('2', []() { g_cfg.physics.usePullB = !g_cfg.physics.usePullB; });
    MiniGL::addKeyFunc('e', []() { g_cfg.physics.pullForce += g_cfg.physics.pullStep; });
    MiniGL::addKeyFunc('d', []() {
        g_cfg.physics.pullForce =
            std::max(static_cast<Real>(0.0), g_cfg.physics.pullForce - g_cfg.physics.pullStep);
    });

    MiniGL::setViewport(
        40.0f, 0.1f, 500.0f,
        Vector3r(0.0, 50.0, 0.0),
        Vector3r(0.0, 0.0, 0.01));

    MiniGL::mainLoop();

    delete Simulation::getCurrent();
    delete base;
    delete model;
    return 0;
}

static void reset()
{
    Simulation::getCurrent()->reset();
    base->getSelectedParticles().clear();

    Simulation::getCurrent()->getModel()->cleanup();
    buildModel();
}

static void buildModel()
{
    TimeManager::getCurrent()->setTimeStepSize(g_cfg.physics.timeStep);

    SimulationModel* model = Simulation::getCurrent()->getModel();
    std::vector<Vector3r> ropePoints =
        generateRandomCablePolyline(g_cfg.random, g_cfg.rope);

    addCableModel(model, ropePoints, g_cfg.rope);

    std::cout << "Built cable with "
              << model->getParticles().getNumberOfParticles()
              << " particles." << std::endl;
}

static std::vector<Vector3r> generateRandomCablePolyline(
    const RandomShapeConfig& randomCfg,
    const RopeConfig& ropeCfg)
{
    const RandomCableLayout layout = makeRandomCableLayout(randomCfg);

    // Build a coarse random path inside the sphere, smooth it, then resample
    // the full path to the final rope particle count.
    std::vector<Vector3r> middlePoints = generateRandomPolylineInSphere(
        layout.sphereEntry,
        layout.sphereExit,
        static_cast<int>(randomCfg.controlPointCount),
        randomCfg.sphereRadius,
        randomCfg.generatorStepSize,
        randomCfg.seed);

    smoothPolyline(
        middlePoints,
        randomCfg.smoothIterations,
        randomCfg.smoothAlpha);

    if (!middlePoints.empty())
    {
        middlePoints.front() = layout.sphereEntry;
        middlePoints.back() = layout.sphereExit;
    }

    std::vector<Vector3r> controlPoints;
    controlPoints.reserve(middlePoints.size() + 2);
    controlPoints.push_back(layout.startPoint);
    for (const auto& p : middlePoints)
        controlPoints.push_back(p);
    controlPoints.push_back(layout.endPoint);

    Real totalLength = static_cast<Real>(0.0);
    for (size_t i = 1; i < controlPoints.size(); ++i)
        totalLength += (controlPoints[i] - controlPoints[i - 1]).norm();

    const int particleCount =
        std::max(2, static_cast<int>(std::ceil(totalLength / ropeCfg.dx)) + 1);

    std::vector<Vector3r> ropePoints =
        resamplePolylineByLength(controlPoints, particleCount);
    if (!ropePoints.empty())
    {
        ropePoints.front() = layout.startPoint;
        ropePoints.back() = layout.endPoint;
    }

    return ropePoints;
}

static RandomCableLayout makeRandomCableLayout(const RandomShapeConfig& randomCfg)
{
    const Vector3r axis(1.0, 0.0, 0.0);
    const Vector3r startPoint(
        -randomCfg.endpointX,
        randomCfg.height,
        static_cast<Real>(0.0));
    const Vector3r endPoint(randomCfg.endpointX, randomCfg.height, static_cast<Real>(0.0));
    const Vector3r sphereCenter(static_cast<Real>(0.0), randomCfg.height, static_cast<Real>(0.0));

    return RandomCableLayout{
        startPoint,
        endPoint,
        sphereCenter - randomCfg.sphereRadius * axis,
        sphereCenter + randomCfg.sphereRadius * axis
    };
}

static void addCableModel(
    SimulationModel* model,
    std::vector<Vector3r>& ropePoints,
    const RopeConfig& ropeCfg)
{
    // Convert the sampled centerline into orientation frames and then attach
    // the rod constraints expected by the PBD line model.
    const int particleCount = static_cast<int>(ropePoints.size());
    const int quaternionCount = particleCount - 1;

    std::vector<Quaternionr> quaternions(quaternionCount);
    Vector3r previousDirection(0, 0, 1);
    for (int i = 0; i < quaternionCount; ++i)
    {
        const Vector3r edge = ropePoints[i + 1] - ropePoints[i];
        Vector3r currentDirection;

        if (edge.squaredNorm() > static_cast<Real>(1e-12))
            currentDirection = edge.normalized();
        else
            currentDirection = previousDirection;

        const Quaternionr rotation =
            Quaternionr::FromTwoVectors(previousDirection, currentDirection);
        quaternions[i] = (i == 0) ? rotation : rotation * quaternions[i - 1];
        previousDirection = currentDirection;
    }

    std::vector<unsigned int> pointIndices(2 * particleCount - 2);
    std::vector<unsigned int> quatIndices(quaternionCount);
    for (int i = 0; i < particleCount - 1; ++i)
    {
        pointIndices[2 * i] = static_cast<unsigned int>(i);
        pointIndices[2 * i + 1] = static_cast<unsigned int>(i + 1);
    }
    for (int i = 0; i < quaternionCount; ++i)
        quatIndices[i] = static_cast<unsigned int>(i);

    model->addLineModel(
        static_cast<unsigned int>(particleCount),
        static_cast<unsigned int>(quaternionCount),
        ropePoints.data(),
        quaternions.data(),
        pointIndices.data(),
        quatIndices.data());

    ParticleData& particles = model->getParticles();
    const unsigned int firstParticle =
        particles.getNumberOfParticles() - static_cast<unsigned int>(particleCount);
    const unsigned int lastParticle = particles.getNumberOfParticles() - 1;

    for (unsigned int i = firstParticle; i <= lastParticle; ++i)
        particles.setMass(i, 1.0);

    g_firstParticle = firstParticle;
    g_lastParticle = lastParticle;

    OrientationData& orientations = model->getOrientations();
    const unsigned int firstQuat =
        orientations.getNumberOfQuaternions() - static_cast<unsigned int>(quaternionCount);
    const unsigned int lastQuat = orientations.getNumberOfQuaternions() - 1;

    for (unsigned int i = firstQuat; i <= lastQuat; ++i)
        orientations.setMass(i, 1.0);

    LineModel* lineModel = model->getLineModels().back();
    const unsigned int pointOffset = lineModel->getIndexOffset();
    const unsigned int quatOffset = lineModel->getIndexOffsetQuaternions();
    const LineModel::Edges& edges = lineModel->getEdges();

    for (unsigned int i = 0; i < edges.size(); ++i)
    {
        const unsigned int p0 = edges[i].m_vert[0] + pointOffset;
        const unsigned int p1 = edges[i].m_vert[1] + pointOffset;
        const unsigned int q = edges[i].m_quat + quatOffset;

        model->addStretchShearConstraint(
            p0,
            p1,
            q,
            ropeCfg.stretchStiffness,
            ropeCfg.shearStiffnessX,
            ropeCfg.shearStiffnessY);
    }

    for (unsigned int i = 0; i + 1 < edges.size(); ++i)
    {
        const unsigned int q0 = edges[i].m_quat + quatOffset;
        const unsigned int q1 = edges[i + 1].m_quat + quatOffset;

        model->addBendTwistConstraint(
            q0,
            q1,
            ropeCfg.twistStiffness,
            ropeCfg.bendStiffnessX,
            ropeCfg.bendStiffnessY);
    }
}

static std::vector<Vector3r> generateRandomPolylineInSphere(
    const Vector3r& start,
    const Vector3r& end,
    int controlCount,
    Real sphereRadius,
    Real generatorStepSize,
    unsigned int seed)
{
    std::mt19937 rng(seed);
    std::normal_distribution<Real> normal01(0.0, 1.0);
    const Vector3r sphereCenter = static_cast<Real>(0.5) * (start + end);

    auto randomUnitVector = [&]() -> Vector3r
    {
        Vector3r v;
        do
        {
            v = Vector3r(normal01(rng), normal01(rng), normal01(rng));
        }
        while (v.squaredNorm() < static_cast<Real>(1e-10));

        v.normalize();
        return v;
    };

    auto projectInsideSphere = [&](Vector3r p) -> Vector3r
    {
        const Vector3r d = p - sphereCenter;
        const Real n = d.norm();
        if (n > sphereRadius && n > static_cast<Real>(1e-12))
            p = sphereCenter + (sphereRadius / n) * d;
        return p;
    };

    auto validPoint = [&](const Vector3r& p) -> bool
    {
        return (p - sphereCenter).squaredNorm() <= sphereRadius * sphereRadius;
    };

    std::vector<Vector3r> points;
    points.reserve(controlCount);
    points.push_back(start);

    Vector3r prevDir = end - start;
    if (prevDir.squaredNorm() < static_cast<Real>(1e-12))
        prevDir = Vector3r(1, 0, 0);
    else
        prevDir.normalize();

    for (int i = 1; i < controlCount - 1; ++i)
    {
        const Vector3r current = points.back();
        const int remainingSegments = controlCount - 1 - i;
        Vector3r bestCandidate = current;
        Real bestScore = -std::numeric_limits<Real>::max();
        bool found = false;

        for (int trial = 0; trial < 60; ++trial)
        {
            Vector3r dirToEnd = end - current;
            if (dirToEnd.squaredNorm() > static_cast<Real>(1e-12))
                dirToEnd.normalize();
            else
                dirToEnd = prevDir;

            Vector3r radialDir = current - sphereCenter;
            if (radialDir.squaredNorm() > static_cast<Real>(1e-12))
                radialDir.normalize();
            else
                radialDir = randomUnitVector();

            Vector3r tangentDir = prevDir.cross(radialDir);
            if (tangentDir.squaredNorm() < static_cast<Real>(1e-12))
                tangentDir = randomUnitVector();
            else
                tangentDir.normalize();

            Vector3r randomDir = randomUnitVector();

            Vector3r dir =
                static_cast<Real>(0.18) * dirToEnd +
                static_cast<Real>(0.28) * prevDir +
                static_cast<Real>(0.34) * tangentDir +
                static_cast<Real>(0.20) * randomDir;

            if (dir.squaredNorm() < static_cast<Real>(1e-12))
                dir = dirToEnd;

            dir.normalize();

            Vector3r candidate = projectInsideSphere(current + generatorStepSize * dir);

            if (!validPoint(candidate))
                continue;

            const Real distToEnd = (end - candidate).norm();
            const Real maxReach =
                static_cast<Real>(remainingSegments) * generatorStepSize * static_cast<Real>(1.55);
            const Real minReach =
                static_cast<Real>(remainingSegments) * generatorStepSize * static_cast<Real>(0.10);

            if (distToEnd > maxReach)
                continue;

            Real minOldDist = std::numeric_limits<Real>::max();
            for (size_t k = 0; k < points.size(); ++k)
                minOldDist = std::min(minOldDist, (candidate - points[k]).norm());

            const Real radial =
                (candidate - sphereCenter).norm() /
                std::max(static_cast<Real>(1e-8), sphereRadius);

            Real score = static_cast<Real>(0.0);
            score += static_cast<Real>(0.35) * radial;
            score += static_cast<Real>(0.20) * std::max(static_cast<Real>(0.0), dir.dot(prevDir));
            score -= static_cast<Real>(0.10) * distToEnd;

            if (minOldDist < static_cast<Real>(0.45) * generatorStepSize)
                score -= static_cast<Real>(4.0);
            else if (minOldDist < static_cast<Real>(1.10) * generatorStepSize)
                score += static_cast<Real>(0.35);

            if (distToEnd < minReach)
                score -= static_cast<Real>(1.5);

            if (score > bestScore)
            {
                bestScore = score;
                bestCandidate = candidate;
                found = true;
            }
        }

        if (!found)
        {
            Vector3r fallbackDir = end - current;
            if (fallbackDir.squaredNorm() < static_cast<Real>(1e-12))
                fallbackDir = prevDir;
            else
                fallbackDir.normalize();

            bestCandidate = projectInsideSphere(current + generatorStepSize * fallbackDir);
        }

        Vector3r newDir = bestCandidate - current;
        if (newDir.squaredNorm() > static_cast<Real>(1e-12))
            prevDir = newDir.normalized();

        points.push_back(bestCandidate);
    }

    points.push_back(end);
    return points;
}

static void smoothPolyline(
    std::vector<Vector3r>& points,
    int iterations,
    Real alpha)
{
    if (points.size() < 3)
        return;

    // Lightweight Laplacian smoothing with pinned endpoints.
    const Vector3r firstLocked = points.front();
    const Vector3r lastLocked = points.back();

    for (int iter = 0; iter < iterations; ++iter)
    {
        const std::vector<Vector3r> oldPoints = points;

        for (size_t i = 1; i + 1 < points.size(); ++i)
        {
            const Vector3r neighborAverage =
                static_cast<Real>(0.5) * (oldPoints[i - 1] + oldPoints[i + 1]);

            points[i] = oldPoints[i] + alpha * (neighborAverage - oldPoints[i]);
        }

        points.front() = firstLocked;
        points.back() = lastLocked;
    }
}

static std::vector<Vector3r> resamplePolylineByLength(
    const std::vector<Vector3r>& polyline,
    int wantedCount)
{
    std::vector<Vector3r> result;

    if (polyline.size() < 2 || wantedCount < 2)
        return result;

    std::vector<Real> cumulativeLength(polyline.size(), 0.0);
    for (size_t i = 1; i < polyline.size(); ++i)
    {
        cumulativeLength[i] =
            cumulativeLength[i - 1] +
            (polyline[i] - polyline[i - 1]).norm();
    }

    const Real totalLength = cumulativeLength.back();
    if (totalLength < static_cast<Real>(1e-8))
    {
        result.assign(wantedCount, polyline.front());
        return result;
    }

    // Uniform arc-length spacing keeps rod particles distributed evenly even
    // when the source polyline is very non-uniform.
    const Real spacing = totalLength / static_cast<Real>(wantedCount - 1);
    result.reserve(wantedCount);

    size_t currentSegment = 0;
    for (int k = 0; k < wantedCount; ++k)
    {
        const Real targetDistance = spacing * static_cast<Real>(k);
        while (currentSegment + 1 < cumulativeLength.size() &&
               cumulativeLength[currentSegment + 1] < targetDistance)
        {
            currentSegment++;
        }

        if (currentSegment + 1 >= polyline.size())
        {
            result.push_back(polyline.back());
            continue;
        }

        const Real s0 = cumulativeLength[currentSegment];
        const Real s1 = cumulativeLength[currentSegment + 1];
        const Real segmentLength = std::max(static_cast<Real>(1e-8), s1 - s0);
        const Real t = (targetDistance - s0) / segmentLength;

        const Vector3r p =
            (static_cast<Real>(1.0) - t) * polyline[currentSegment] +
            t * polyline[currentSegment + 1];
        result.push_back(p);
    }

    return result;
}

static void timeStep()
{
    if (base->getValue<bool>(DemoBase::PAUSE))
        return;

    SimulationModel* model = Simulation::getCurrent()->getModel();
    ParticleData& particles = model->getParticles();

    const unsigned int stepsPerFrame =
        base->getValue<unsigned int>(DemoBase::NUM_STEPS_PER_RENDER);
    const Real dt = TimeManager::getCurrent()->getTimeStepSize();

    for (unsigned int step = 0; step < stepsPerFrame; ++step)
    {
        if (g_cfg.physics.usePullA || g_cfg.physics.usePullB)
            applyPullForces(model, dt, g_cfg.physics);

        Simulation::getCurrent()->getTimeStep()->step(*model);

        if (g_cfg.selfCollision.enabled)
        {
            std::vector<Vector3r> oldPositions(particles.getNumberOfParticles());
            for (unsigned int i = 0; i < particles.getNumberOfParticles(); ++i)
                oldPositions[i] = particles.getPosition(i);

            for (unsigned int k = 0; k < g_cfg.selfCollision.iterations; ++k)
                applySelfCollisionXPBD(model, dt, g_cfg.rope, g_cfg.selfCollision);

            applyGroundConstraint(g_cfg.groundY);

            recomputeVelocitiesFromPositions(
                model,
                oldPositions,
                dt,
                g_cfg.physics.velocityRecomputeBlend);
        }
        else
        {
            applyGroundConstraint(g_cfg.groundY);
        }

        base->step();
    }
}

static void applyPullForces(SimulationModel* model, Real dt, const PhysicsConfig& physicsCfg)
{
    ParticleData& particles = model->getParticles();

    // Endpoint pulls are applied directly as velocity changes each step.
    auto applyForce = [&](unsigned int particleIndex, const Vector3r& forceDir)
    {
        if (particleIndex >= particles.getNumberOfParticles())
            return;
        if (particles.getMass(particleIndex) <= 0.0)
            return;

        Vector3r& velocity = particles.getVelocity(particleIndex);
        velocity += dt * physicsCfg.pullForce * forceDir;
    };

    if (physicsCfg.usePullA)
        applyForce(g_lastParticle, Vector3r(1.0, 0.0, 0.0));

    if (physicsCfg.usePullB)
        applyForce(g_firstParticle, Vector3r(-1.0, 0.0, 0.0));
}

static Real clamp01(Real x)
{
    if (x < static_cast<Real>(0.0))
        return static_cast<Real>(0.0);
    if (x > static_cast<Real>(1.0))
        return static_cast<Real>(1.0);
    return x;
}

static Real getInverseMass(const ParticleData& pd, unsigned int index)
{
    return (pd.getMass(index) > 0.0)
        ? static_cast<Real>(1.0) / pd.getMass(index)
        : static_cast<Real>(0.0);
}

static void closestPtSegmentSegment(
    const Vector3r& p1, const Vector3r& q1,
    const Vector3r& p2, const Vector3r& q2,
    Real& s, Real& t,
    Vector3r& c1, Vector3r& c2)
{
    const Real eps = static_cast<Real>(1e-8);

    const Vector3r d1 = q1 - p1;
    const Vector3r d2 = q2 - p2;
    const Vector3r r = p1 - p2;

    const Real a = d1.dot(d1);
    const Real e = d2.dot(d2);
    const Real f = d2.dot(r);

    if (a <= eps && e <= eps)
    {
        s = static_cast<Real>(0.0);
        t = static_cast<Real>(0.0);
        c1 = p1;
        c2 = p2;
        return;
    }

    const Real c = d1.dot(r);

    if (a <= eps)
    {
        s = static_cast<Real>(0.0);
        t = clamp01(f / e);
    }
    else if (e <= eps)
    {
        t = static_cast<Real>(0.0);
        s = clamp01(-c / a);
    }
    else
    {
        const Real b = d1.dot(d2);
        const Real denom = a * e - b * b;

        if (std::abs(denom) > eps)
            s = clamp01((b * f - c * e) / denom);
        else
            s = static_cast<Real>(0.0);

        t = (b * s + f) / e;

        if (t < static_cast<Real>(0.0))
        {
            t = static_cast<Real>(0.0);
            s = clamp01(-c / a);
        }
        else if (t > static_cast<Real>(1.0))
        {
            t = static_cast<Real>(1.0);
            s = clamp01((b - c) / a);
        }
    }

    c1 = p1 + s * d1;
    c2 = p2 + t * d2;
}

static void recomputeVelocitiesFromPositions(
    SimulationModel* model,
    const std::vector<Vector3r>& oldPositions,
    Real dt,
    Real blend)
{
    ParticleData& pd = model->getParticles();
    const unsigned int n = pd.getNumberOfParticles();
    if (oldPositions.size() != n || dt <= static_cast<Real>(1e-12))
        return;

    blend = std::max(static_cast<Real>(0.0), std::min(static_cast<Real>(1.0), blend));

    // Blend physical velocities toward the actually observed post-correction motion.
    for (unsigned int i = 0; i < n; ++i)
    {
        if (pd.getMass(i) <= 0.0)
            continue;

        const Vector3r vPos = (pd.getPosition(i) - oldPositions[i]) / dt;
        Vector3r& v = pd.getVelocity(i);
        v = (static_cast<Real>(1.0) - blend) * v + blend * vPos;
    }
}

static std::vector<std::pair<unsigned int, unsigned int>>
collectSelfCollisionCandidatesSpatialHash(
    const ParticleData& pd,
    Real inflate,
    Real cellSize)
{
    const unsigned int n = pd.getNumberOfParticles();
    if (n < 4)
        return {};

    const unsigned int nEdges = n - 1;
    std::unordered_map<Int3, std::vector<unsigned int>, Int3Hash> grid;
    grid.reserve(nEdges * 2);

    // Hash each edge AABB into a 3D grid, then only test pairs that share cells.
    auto cellOf = [cellSize](const Vector3r& p) -> Int3
    {
        return Int3{
            static_cast<int>(std::floor(p.x() / cellSize)),
            static_cast<int>(std::floor(p.y() / cellSize)),
            static_cast<int>(std::floor(p.z() / cellSize))
        };
    };

    for (unsigned int e = 0; e < nEdges; ++e)
    {
        const Vector3r p0 = pd.getPosition(e);
        const Vector3r p1 = pd.getPosition(e + 1);

        Vector3r bmin = p0.cwiseMin(p1);
        Vector3r bmax = p0.cwiseMax(p1);
        bmin.array() -= inflate;
        bmax.array() += inflate;

        const Int3 cmin = cellOf(bmin);
        const Int3 cmax = cellOf(bmax);

        for (int x = cmin.x; x <= cmax.x; ++x)
        {
            for (int y = cmin.y; y <= cmax.y; ++y)
            {
                for (int z = cmin.z; z <= cmax.z; ++z)
                    grid[Int3{x, y, z}].push_back(e);
            }
        }
    }

    std::unordered_set<uint64_t> uniquePairs;
    uniquePairs.reserve(nEdges * 8);

    for (const auto& entry : grid)
    {
        const std::vector<unsigned int>& edges = entry.second;
        for (size_t a = 0; a < edges.size(); ++a)
        {
            for (size_t b = a + 1; b < edges.size(); ++b)
            {
                const unsigned int i = edges[a];
                const unsigned int j = edges[b];

                if (i == j)
                    continue;
                if (j <= i + 1 && j >= i)
                    continue;

                uniquePairs.insert(makeOrderedEdgePairKey(i, j));
            }
        }
    }

    std::vector<std::pair<unsigned int, unsigned int>> pairs;
    pairs.reserve(uniquePairs.size());
    for (uint64_t key : uniquePairs)
    {
        pairs.emplace_back(
            static_cast<unsigned int>(key >> 32),
            static_cast<unsigned int>(key & 0xffffffffu));
    }

    return pairs;
}

static uint64_t makeOrderedEdgePairKey(unsigned int a, unsigned int b)
{
    const uint64_t lo = std::min<uint64_t>(a, b);
    const uint64_t hi = std::max<uint64_t>(a, b);
    return (lo << 32) | hi;
}

static Vector3r clampVec(const Vector3r& v, Real maxNorm)
{
    const Real len = v.norm();
    if (len > maxNorm && len > static_cast<Real>(1e-12))
        return v * (maxNorm / len);
    return v;
}

static void applySelfCollisionXPBD(
    SimulationModel* model,
    Real dt,
    const RopeConfig& ropeCfg,
    const SelfCollisionConfig& collisionCfg)
{
    ParticleData& pd = model->getParticles();
    const unsigned int n = pd.getNumberOfParticles();
    if (n < 4)
        return;

    const Real minDist = static_cast<Real>(2.0) * collisionCfg.ropeRadius;
    const Real eps = static_cast<Real>(1e-8);
    const Real alpha = collisionCfg.compliance / (dt * dt);

    static std::vector<Vector3r> posDelta;
    static std::vector<unsigned int> contactCount;
    static std::vector<VelocityPairCorrection> velocityCorrections;

    if (posDelta.size() < n)
        posDelta.resize(n);
    if (contactCount.size() < n)
        contactCount.resize(n);

    for (unsigned int i = 0; i < n; ++i)
    {
        posDelta[i].setZero();
        contactCount[i] = 0;
    }

    velocityCorrections.clear();
    velocityCorrections.reserve(n * 2);

    const std::vector<std::pair<unsigned int, unsigned int>> candidates =
        collectSelfCollisionCandidatesSpatialHash(
            pd,
            minDist,
            collisionCfg.hashCellSize);

    // First solve position-level separation, then damp the corresponding
    // relative contact velocities.
    for (const auto& pair : candidates)
    {
        const unsigned int i = pair.first;
        const unsigned int j = pair.second;

        const Real arcDist = static_cast<Real>(std::abs((int)j - (int)i)) * ropeCfg.dx;
        if (arcDist < collisionCfg.excludeArcLength)
            continue;

        const unsigned int i0 = i;
        const unsigned int i1 = i + 1;
        const unsigned int j0 = j;
        const unsigned int j1 = j + 1;

        const Vector3r p0 = pd.getPosition(i0);
        const Vector3r p1 = pd.getPosition(i1);
        const Vector3r q0 = pd.getPosition(j0);
        const Vector3r q1 = pd.getPosition(j1);

        Real s, t;
        Vector3r c1, c2;
        closestPtSegmentSegment(p0, p1, q0, q1, s, t, c1, c2);

        Vector3r delta = c2 - c1;
        Real dist = delta.norm();

        Vector3r nrm;
        if (dist > eps)
        {
            nrm = delta / dist;
        }
        else
        {
            nrm = (p1 - p0).cross(q1 - q0);
            if (nrm.squaredNorm() < eps)
                nrm = Vector3r(1, 0, 0);
            else
                nrm.normalize();

            dist = static_cast<Real>(0.0);
        }

        const Real C = dist - minDist;
        if (C >= static_cast<Real>(0.0))
            continue;

        const Real wP0 = static_cast<Real>(1.0) - s;
        const Real wP1 = s;
        const Real wQ0 = static_cast<Real>(1.0) - t;
        const Real wQ1 = t;

        const Real invMassP0 = getInverseMass(pd, i0);
        const Real invMassP1 = getInverseMass(pd, i1);
        const Real invMassQ0 = getInverseMass(pd, j0);
        const Real invMassQ1 = getInverseMass(pd, j1);

        if (invMassP0 + invMassP1 + invMassQ0 + invMassQ1 <= eps)
            continue;

        const Vector3r gradP0 = -wP0 * nrm;
        const Vector3r gradP1 = -wP1 * nrm;
        const Vector3r gradQ0 = wQ0 * nrm;
        const Vector3r gradQ1 = wQ1 * nrm;

        const Real denom =
            invMassP0 * gradP0.squaredNorm() +
            invMassP1 * gradP1.squaredNorm() +
            invMassQ0 * gradQ0.squaredNorm() +
            invMassQ1 * gradQ1.squaredNorm() +
            alpha;
        if (denom <= eps)
            continue;

        const Real deltaLambda = -C / denom;

        Vector3r dp0 = clampVec(invMassP0 * deltaLambda * gradP0, collisionCfg.maxCorrection);
        Vector3r dp1 = clampVec(invMassP1 * deltaLambda * gradP1, collisionCfg.maxCorrection);
        Vector3r dq0 = clampVec(invMassQ0 * deltaLambda * gradQ0, collisionCfg.maxCorrection);
        Vector3r dq1 = clampVec(invMassQ1 * deltaLambda * gradQ1, collisionCfg.maxCorrection);

        if (invMassP0 > 0.0)
        {
            posDelta[i0] += dp0;
            contactCount[i0] += 1;
        }
        if (invMassP1 > 0.0)
        {
            posDelta[i1] += dp1;
            contactCount[i1] += 1;
        }
        if (invMassQ0 > 0.0)
        {
            posDelta[j0] += dq0;
            contactCount[j0] += 1;
        }
        if (invMassQ1 > 0.0)
        {
            posDelta[j1] += dq1;
            contactCount[j1] += 1;
        }

        velocityCorrections.push_back({ i0, i1, j0, j1, wP0, wP1, wQ0, wQ1, nrm });
    }

    for (unsigned int i = 0; i < n; ++i)
    {
        if (pd.getMass(i) <= 0.0 || contactCount[i] == 0)
            continue;

        const Real averageScale =
            static_cast<Real>(0.5) / static_cast<Real>(contactCount[i]);
        Vector3r corr = clampVec(averageScale * posDelta[i], collisionCfg.maxCorrection);
        pd.setPosition(i, pd.getPosition(i) + corr);
    }

    for (const auto& vc : velocityCorrections)
    {
        const Real invMassP0 = getInverseMass(pd, vc.i0);
        const Real invMassP1 = getInverseMass(pd, vc.i1);
        const Real invMassQ0 = getInverseMass(pd, vc.j0);
        const Real invMassQ1 = getInverseMass(pd, vc.j1);

        if (invMassP0 + invMassP1 + invMassQ0 + invMassQ1 <= eps)
            continue;

        const Vector3r vP =
            vc.wP0 * pd.getVelocity(vc.i0) +
            vc.wP1 * pd.getVelocity(vc.i1);
        const Vector3r vQ =
            vc.wQ0 * pd.getVelocity(vc.j0) +
            vc.wQ1 * pd.getVelocity(vc.j1);

        const Vector3r relV = vQ - vP;
        const Real vn = relV.dot(vc.nrm);

        auto applyImpulse = [&](const Vector3r& impulse)
        {
            if (invMassP0 > 0.0) pd.getVelocity(vc.i0) -= (vc.wP0 * invMassP0) * impulse;
            if (invMassP1 > 0.0) pd.getVelocity(vc.i1) -= (vc.wP1 * invMassP1) * impulse;
            if (invMassQ0 > 0.0) pd.getVelocity(vc.j0) += (vc.wQ0 * invMassQ0) * impulse;
            if (invMassQ1 > 0.0) pd.getVelocity(vc.j1) += (vc.wQ1 * invMassQ1) * impulse;
        };

        if (vn < static_cast<Real>(0.0) && collisionCfg.normalDamping > static_cast<Real>(0.0))
            applyImpulse(-collisionCfg.normalDamping * vn * vc.nrm);

        const Vector3r vt = relV - relV.dot(vc.nrm) * vc.nrm;
        if (vt.squaredNorm() > eps && collisionCfg.tangentialDamping > static_cast<Real>(0.0))
            applyImpulse(-collisionCfg.tangentialDamping * vt);
    }
}

static void render()
{
    base->render();

    SimulationModel* model = Simulation::getCurrent()->getModel();
    const ParticleData& particles = model->getParticles();

    float blue[4] = {0.0f, 0.0f, 0.8f, 1.0f};

    // Render the rod centerline as spheres plus cylinders.
    for (unsigned int m = 0; m < model->getLineModels().size(); ++m)
    {
        LineModel* lineModel = model->getLineModels()[m];
        const unsigned int pointOffset = lineModel->getIndexOffset();

        for (unsigned int e = 0; e < lineModel->getEdges().size(); ++e)
        {
            const unsigned int i0 = lineModel->getEdges()[e].m_vert[0] + pointOffset;
            const unsigned int i1 = lineModel->getEdges()[e].m_vert[1] + pointOffset;

            const Vector3r& p0 = particles.getPosition(i0);
            const Vector3r& p1 = particles.getPosition(i1);

            MiniGL::drawSphere(p0, 0.07f, blue);
            if (e == lineModel->getEdges().size() - 1)
                MiniGL::drawSphere(p1, 0.07f, blue);

            MiniGL::drawCylinder(p0, p1, blue, static_cast<float>(g_cfg.rope.renderRadius));
        }
    }
}

static void applyGroundConstraint(Real groundY)
{
    SimulationModel* model = Simulation::getCurrent()->getModel();
    ParticleData& particles = model->getParticles();

    // Clamp particles to the ground plane and remove downward velocity.
    for (unsigned int i = 0; i < particles.getNumberOfParticles(); ++i)
    {
        Vector3r p = particles.getPosition(i);
        if (p.y() < groundY)
        {
            p.y() = groundY;
            particles.setPosition(i, p);

            Vector3r& v = particles.getVelocity(i);
            if (v.y() < 0.0)
                v.y() = 0.0;
        }
    }
}
