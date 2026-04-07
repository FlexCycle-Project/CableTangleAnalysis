#include "Common/Common.h"
#include "Demos/Visualization/MiniGL.h"
#include "Simulation/TimeManager.h"
#include "Simulation/SimulationModel.h"
#include "Simulation/Constraints.h"
#include "Demos/Common/DemoBase.h"
#include "Simulation/Simulation.h"
#include "Utils/FileSystem.h"
#include "KymoKnotAnalyzer.h"

#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <cstdint>
#include <limits>
#include <sstream>
#include <utility>

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
    Real dx = static_cast<Real>(0.1);
    Real renderRadius = static_cast<Real>(0.06);
};

// Per-step simulation controls and endpoint pulling.
struct PhysicsConfig
{
    Real timeStep = static_cast<Real>(0.005);
    Real pullSpring = static_cast<Real>(1000.0);
    Real pullDamping = static_cast<Real>(2.0);
    Real pullMaxAccel = static_cast<Real>(1000.0);
    Real targetStep = static_cast<Real>(0.05);
    Real velocityRecomputeBlend = static_cast<Real>(0.05);
};

// Slow bend/twist rest-state creep makes the rod behave more like a cable.
struct RestShapeCreepConfig
{
    bool enabled = true;
    unsigned int stepInterval = 20;
    Real creep = static_cast<Real>(0.01);
    Real activationThreshold = static_cast<Real>(1e-4);
};

// Endpoints for the base cable layout.
struct CableLayoutConfig
{
    Real endpointX = static_cast<Real>(10.0);
    Real height = static_cast<Real>(6.0);
};

// Randomized control-walk settings for the mess generator.
struct RandomCableConfig
{
    bool enabled = true;
    unsigned int seed = 0;
    int controlPointCount = 30;
    Real controlStepLength = static_cast<Real>(1.5);
    Real straightChance = static_cast<Real>(0);
    Real leftChance = static_cast<Real>(0.5);
    Real rightChance = static_cast<Real>(0.5);
    Real centerPull = static_cast<Real>(0.15);
    Real maxLateralOffset = static_cast<Real>(4.5);
    int smoothingIterations = 1;
    Real smoothingAlpha = static_cast<Real>(0.08);
    Real minHeightAboveGround = static_cast<Real>(5.0);
};

// Tuning for the custom rope self-collision pass.
struct SelfCollisionConfig
{
    bool enabled = true;
    Real ropeRadius = static_cast<Real>(0.08);
    Real compliance = static_cast<Real>(1e-10);
    unsigned int iterations = 6;
    Real excludeArcLength = static_cast<Real>(0.9);
    Real normalDamping = static_cast<Real>(0.2);
    Real tangentialDamping = static_cast<Real>(0.08);
    Real maxCorrection = static_cast<Real>(0.05);
    Real hashCellSize = static_cast<Real>(0.2);
};

// Single source of truth for the demo.
struct AppConfig
{
    RopeConfig rope;
    PhysicsConfig physics;
    RestShapeCreepConfig restCreep;
    CableLayoutConfig layout;
    RandomCableConfig randomCable;
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

struct SelfCollisionContact
{
    unsigned int edgeA;
    unsigned int edgeB;
};

static void buildModel();
static void reset();
static void regenerateRandomCable();
static void timeStep();
static void render();

static std::vector<Vector3r> generateStraightCablePolyline(
    const CableLayoutConfig& layoutCfg,
    const RopeConfig& ropeCfg);
static std::vector<Vector3r> generateRandomCableMessPolyline(
    const CableLayoutConfig& layoutCfg,
    const RopeConfig& ropeCfg,
    const RandomCableConfig& randomCfg,
    Real groundY);
static void smoothPolyline(
    std::vector<Vector3r>& points,
    int iterations,
    Real alpha);
static std::vector<Vector3r> resamplePolylineByLength(
    const std::vector<Vector3r>& polyline,
    int wantedCount);

static void addCableModel(
    SimulationModel* model,
    std::vector<Vector3r>& ropePoints,
    const RopeConfig& ropeCfg);

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

static std::vector<SelfCollisionContact> detectSelfCollisionContacts(
    SimulationModel* model,
    const RopeConfig& ropeCfg,
    const SelfCollisionConfig& collisionCfg);
static unsigned int countSelfCollisionContactClusters(
    const std::vector<SelfCollisionContact>& contacts);
static void printSelfCollisionContacts(
    const std::vector<SelfCollisionContact>& contacts);
static void printTanglementMetrics(
    SimulationModel* model,
    const RopeConfig& ropeCfg,
    const SelfCollisionConfig& collisionCfg);
static void toggleKymoKnotSearchMode();
static const char* getKymoKnotSearchModeName(PBD::KymoKnotSearchMode mode);
static std::vector<Vector3r> getCurrentCablePoints(SimulationModel* model);
static bool isLocalIndexInKymoKnotSegment(
    int localIndex,
    const PBD::KymoKnotSegment& segment);
static bool isLocalIndexHighlighted(int localIndex);
static bool isLocalEdgeHighlighted(int localEdgeIndex);
static void cacheKymoKnotResult(const PBD::KymoKnotAnalysisResult& result);
static void getHeightColor(Real y, float* color);

static std::string getCableSnapshotDirectory();
static std::string getCableSnapshotListPath();
static std::string getCableSnapshotPath();
static std::vector<std::string> getAvailableSnapshotPaths();
static bool loadCableSnapshot(
    const std::string& snapshotPath,
    std::vector<Vector3r>& ropePoints);
static void applyBendTwistRestStateCreep(
    SimulationModel* model,
    const RestShapeCreepConfig& creepCfg);
static std::vector<Vector3r> generateMathematicalKnotPolyline(
    const std::string& knotName);
static bool writeCableSnapshotFile(
    const std::string& snapshotPath,
    const std::vector<Vector3r>& ropePoints);
static void ensureCableSnapshotLibrary();
static void printSnapshotRecognitionReport(
    const std::string& snapshotPath,
    const std::vector<Vector3r>& ropePoints);
static void auditSnapshotLibrary();
static std::string getSnapshotDisplayName(const std::string& snapshotPath);

static DemoBase* base = nullptr;
static AppConfig g_cfg;
static unsigned int g_firstParticle = 0;
static unsigned int g_lastParticle = 0;
static Vector3r g_targetA(0.0, 0.0, 0.0);
static Vector3r g_targetB(0.0, 0.0, 0.0);
static bool g_useSavedCableSnapshot = false;
static unsigned int g_restCreepCounter = 0;
static int g_snapshotCycleIndex = -1;
static PBD::KymoKnotAnalysisResult g_lastKymoKnotResult;
static bool g_hasKymoKnotResult = false;
static std::string g_currentSnapshotPath;
static PBD::KymoKnotSearchMode g_kymoKnotSearchMode = PBD::KymoKnotSearchMode::TopDown;

int main(int argc, char** argv)
{
    base = new DemoBase();

    base->init(argc, argv, "Random Cable Mess Demo");

    SimulationModel* model = new SimulationModel();
    model->init();
    Simulation::getCurrent()->setModel(model);

    buildModel();

    // Keyboard input for this demo lives here.
    MiniGL::setClientIdleFunc(timeStep);
    MiniGL::setClientSceneFunc(render);
    MiniGL::addKeyFunc('r', reset);
    MiniGL::addKeyFunc('t', []() {
        printTanglementMetrics(
            Simulation::getCurrent()->getModel(),
            g_cfg.rope,
            g_cfg.selfCollision);
    });
    MiniGL::addKeyFunc('b', []() {
        toggleKymoKnotSearchMode();
    });
    MiniGL::addKeyFunc('g', []() {
        const std::vector<SelfCollisionContact> contacts =
            detectSelfCollisionContacts(
                Simulation::getCurrent()->getModel(),
                g_cfg.rope,
                g_cfg.selfCollision);
        printSelfCollisionContacts(contacts);
    });
    MiniGL::addKeyFunc('u', []() {
        const std::vector<std::string> snapshotPaths = getAvailableSnapshotPaths();
        if (snapshotPaths.empty())
            return;

        g_snapshotCycleIndex = (g_snapshotCycleIndex + 1) %
            static_cast<int>(snapshotPaths.size());
        g_useSavedCableSnapshot = true;
        reset();
    });
    MiniGL::addKeyFunc('p', []() {
        g_useSavedCableSnapshot = false;
        g_snapshotCycleIndex = -1;
        reset();
    });
    MiniGL::addKeyFunc('o', []() {
        regenerateRandomCable();
    });
    MiniGL::addKeyFunc('z', []() { MiniGL::move(0, 0, -0.5); });
    MiniGL::addKeyFunc('x', []() { MiniGL::move(0, 0, 0.5); });

    // E/D/F/S/C/V now move the magenta target.
    MiniGL::addKeyFunc('e', []() { g_targetB.z() += g_cfg.physics.targetStep; });
    MiniGL::addKeyFunc('d', []() { g_targetB.z() -= g_cfg.physics.targetStep; });
    MiniGL::addKeyFunc('f', []() { g_targetB.x() -= g_cfg.physics.targetStep; });
    MiniGL::addKeyFunc('s', []() { g_targetB.x() += g_cfg.physics.targetStep; });
    MiniGL::addKeyFunc('c', []() { g_targetB.y() -= g_cfg.physics.targetStep; });
    MiniGL::addKeyFunc('v', []() { g_targetB.y() += g_cfg.physics.targetStep; });

    // I/K/L/J/N/M now move the yellow target.
    MiniGL::addKeyFunc('i', []() { g_targetA.z() += g_cfg.physics.targetStep; });
    MiniGL::addKeyFunc('k', []() { g_targetA.z() -= g_cfg.physics.targetStep; });
    MiniGL::addKeyFunc('l', []() { g_targetA.x() -= g_cfg.physics.targetStep; });
    MiniGL::addKeyFunc('j', []() { g_targetA.x() += g_cfg.physics.targetStep; });
    MiniGL::addKeyFunc('n', []() { g_targetA.y() -= g_cfg.physics.targetStep; });
    MiniGL::addKeyFunc('m', []() { g_targetA.y() += g_cfg.physics.targetStep; });

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

static void regenerateRandomCable()
{
    g_useSavedCableSnapshot = false;
    g_snapshotCycleIndex = -1;
    g_cfg.randomCable.enabled = true;

    std::random_device rd;
    unsigned int newSeed = rd();
    if (newSeed == 0)
        newSeed = 1;
    g_cfg.randomCable.seed = newSeed;

    std::cout << "Regenerating random cable with seed " << g_cfg.randomCable.seed << std::endl;
    reset();
}

static void buildModel()
{
    TimeManager::getCurrent()->setTimeStepSize(g_cfg.physics.timeStep);

    SimulationModel* model = Simulation::getCurrent()->getModel();
    std::vector<Vector3r> ropePoints;
    g_currentSnapshotPath.clear();
    if (!g_useSavedCableSnapshot || !loadCableSnapshot(getCableSnapshotPath(), ropePoints))
    {
        if (g_cfg.randomCable.enabled)
        {
            ropePoints = generateRandomCableMessPolyline(
                g_cfg.layout,
                g_cfg.rope,
                g_cfg.randomCable,
                g_cfg.groundY);
        }
        else
        {
            ropePoints = generateStraightCablePolyline(g_cfg.layout, g_cfg.rope);
        }
        if (g_useSavedCableSnapshot)
            g_useSavedCableSnapshot = false;
    }

    g_hasKymoKnotResult = false;
    g_lastKymoKnotResult = PBD::KymoKnotAnalysisResult();

    addCableModel(model, ropePoints, g_cfg.rope);

    if (g_firstParticle < model->getParticles().getNumberOfParticles() &&
        g_lastParticle < model->getParticles().getNumberOfParticles())
    {
        g_targetA = model->getParticles().getPosition(g_firstParticle);
        g_targetB = model->getParticles().getPosition(g_lastParticle);
    }
}

static std::vector<Vector3r> generateStraightCablePolyline(
    const CableLayoutConfig& layoutCfg,
    const RopeConfig& ropeCfg)
{
    const Vector3r startPoint(
        -layoutCfg.endpointX,
        layoutCfg.height,
        static_cast<Real>(0.0));
    const Vector3r endPoint(
        layoutCfg.endpointX,
        layoutCfg.height,
        static_cast<Real>(0.0));
    const Real totalLength = (endPoint - startPoint).norm();
    const int particleCount =
        std::max(2, static_cast<int>(std::ceil(totalLength / ropeCfg.dx)) + 1);

    std::vector<Vector3r> ropePoints;
    ropePoints.reserve(particleCount);

    for (int i = 0; i < particleCount; ++i)
    {
        const Real t =
            static_cast<Real>(i) / static_cast<Real>(std::max(1, particleCount - 1));
        ropePoints.push_back(
            (static_cast<Real>(1.0) - t) * startPoint + t * endPoint);
    }

    return ropePoints;
}

static std::vector<Vector3r> generateRandomCableMessPolyline(
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

static void smoothPolyline(
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

static std::vector<Vector3r> resamplePolylineByLength(
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

static std::string getCableSnapshotDirectory()
{
    return Utilities::FileSystem::normalizePath(
        base->getOutputPath() + "/saved_cable");
}

static std::string getCableSnapshotListPath()
{
    return Utilities::FileSystem::normalizePath(
        getCableSnapshotDirectory() + "/snapshotlist.txt");
}

static std::string getCableSnapshotPath()
{
    const std::vector<std::string> snapshotPaths = getAvailableSnapshotPaths();
    if (snapshotPaths.empty())
        return Utilities::FileSystem::normalizePath(getCableSnapshotDirectory() + "/snapshot.txt");

    int index = g_snapshotCycleIndex;
    if (index < 0 || index >= static_cast<int>(snapshotPaths.size()))
        index = 0;
    return snapshotPaths[static_cast<size_t>(index)];
}

static std::vector<std::string> getAvailableSnapshotPaths()
{
    ensureCableSnapshotLibrary();

    std::vector<std::string> snapshotPaths;
    const std::string listPath = getCableSnapshotListPath();
    if (Utilities::FileSystem::isFile(listPath))
    {
        std::ifstream in(listPath.c_str());
        std::string line;
        while (std::getline(in, line))
        {
            const std::size_t commentPos = line.find('#');
            if (commentPos != std::string::npos)
                line = line.substr(0, commentPos);

            auto trim = [](std::string& value)
            {
                const std::size_t begin = value.find_first_not_of(" \t\r\n");
                if (begin == std::string::npos)
                {
                    value.clear();
                    return;
                }

                const std::size_t end = value.find_last_not_of(" \t\r\n");
                value = value.substr(begin, end - begin + 1);
            };

            trim(line);
            if (line.empty())
                continue;

            std::string candidate = Utilities::FileSystem::normalizePath(line);
            if (Utilities::FileSystem::isRelativePath(candidate))
                candidate = Utilities::FileSystem::normalizePath(getCableSnapshotDirectory() + "/" + candidate);

            if (Utilities::FileSystem::isFile(candidate))
                snapshotPaths.push_back(candidate);
        }

        if (!snapshotPaths.empty())
            return snapshotPaths;
    }

    std::vector<std::string> fileNames;
    const std::string snapshotDirectory = getCableSnapshotDirectory();

    if (!Utilities::FileSystem::getFilesInDirectory(snapshotDirectory, fileNames))
        return snapshotPaths;

    std::sort(fileNames.begin(), fileNames.end());
    for (const std::string& fileName : fileNames)
    {
        if (fileName == "." || fileName == ".." || fileName == "snapshotlist.txt")
            continue;

        if (Utilities::FileSystem::getFileExt(fileName) != "txt")
            continue;

        const std::string fullPath =
            Utilities::FileSystem::normalizePath(snapshotDirectory + "/" + fileName);
        if (Utilities::FileSystem::isFile(fullPath))
            snapshotPaths.push_back(fullPath);
    }

    return snapshotPaths;
}

static bool loadCableSnapshot(
    const std::string& snapshotPath,
    std::vector<Vector3r>& ropePoints)
{
    ropePoints.clear();
    g_currentSnapshotPath.clear();

    if (base == nullptr)
        return false;

    std::ifstream in(snapshotPath);
    if (!in.is_open())
        return false;

    size_t count = 0;
    in >> count;
    if (!in.good() || count < 2)
        return false;

    ropePoints.reserve(count);
    for (size_t i = 0; i < count; ++i)
    {
        Real x, y, z;
        in >> x >> y >> z;
        if (!in.good())
        {
            ropePoints.clear();
            return false;
        }
        ropePoints.emplace_back(x, y, z);
    }
    g_currentSnapshotPath = snapshotPath;
    return true;
}

static std::string getSnapshotDisplayName(const std::string& snapshotPath)
{
    if (snapshotPath.empty())
        return "<none>";
    return Utilities::FileSystem::getFileNameWithExt(snapshotPath);
}

static bool writeCableSnapshotFile(
    const std::string& snapshotPath,
    const std::vector<Vector3r>& ropePoints)
{
    std::ofstream out(snapshotPath.c_str(), std::ios::out | std::ios::trunc);
    if (!out.is_open())
        return false;

    out << ropePoints.size() << "\n";
    out << std::fixed << std::setprecision(9);
    for (const Vector3r& p : ropePoints)
        out << p.x() << " " << p.y() << " " << p.z() << "\n";
    return true;
}

static std::vector<Vector3r> generateMathematicalKnotPolyline(
    const std::string& knotName)
{
    auto sampleTrefoil = [](Real t) -> Vector3r
    {
        return Vector3r(
            std::sin(t) + static_cast<Real>(2.0) * std::sin(static_cast<Real>(2.0) * t),
            std::cos(t) - static_cast<Real>(2.0) * std::cos(static_cast<Real>(2.0) * t),
            -std::sin(static_cast<Real>(3.0) * t));
    };

    auto sampleFigureEight = [](Real t) -> Vector3r
    {
        return Vector3r(
            (static_cast<Real>(2.0) + std::cos(static_cast<Real>(2.0) * t)) * std::cos(static_cast<Real>(3.0) * t),
            (static_cast<Real>(2.0) + std::cos(static_cast<Real>(2.0) * t)) * std::sin(static_cast<Real>(3.0) * t),
            std::sin(static_cast<Real>(4.0) * t));
    };

    auto sampleUnknot = [](Real t) -> Vector3r
    {
        return Vector3r(
            std::cos(t),
            std::sin(t),
            static_cast<Real>(0.15) * std::sin(static_cast<Real>(2.0) * t));
    };

    auto sampleKnot = [&](const std::string& type, Real t) -> Vector3r
    {
        if (type == "figure8" || type == "figure-eight" || type == "4_1")
            return sampleFigureEight(t);
        if (type == "unknot")
            return sampleUnknot(t);
        return sampleTrefoil(t);
    };

    auto buildCore = [&](const std::string& type, int samples) -> std::vector<Vector3r>
    {
        std::vector<Vector3r> points;
        points.reserve(samples);
        for (int i = 0; i < samples; ++i)
        {
            const Real t =
                static_cast<Real>(2.0) * static_cast<Real>(std::acos(static_cast<Real>(-1.0))) *
                static_cast<Real>(i) / static_cast<Real>(samples);
            points.push_back(sampleKnot(type, t));
        }

        Vector3r center(0.0, 0.0, 0.0);
        for (const Vector3r& p : points)
            center += p;
        center /= static_cast<Real>(points.size());

        Real maxRadius = static_cast<Real>(0.0);
        for (Vector3r& p : points)
        {
            p -= center;
            maxRadius = std::max(maxRadius, p.norm());
        }

        if (maxRadius > static_cast<Real>(1e-12))
        {
            const Real scale = static_cast<Real>(3.0) / maxRadius;
            for (Vector3r& p : points)
                p *= scale;
        }

        return points;
    };

    auto appendOffset = [](std::vector<Vector3r>& dst, const std::vector<Vector3r>& src, const Vector3r& offset)
    {
        for (const Vector3r& p : src)
            dst.push_back(p + offset);
    };

    std::vector<Vector3r> ropePoints;
    const Real knotLift = static_cast<Real>(5.0);

    if (knotName == "multi" || knotName == "double" || knotName == "mixed")
    {
        const std::vector<Vector3r> trefoilCore = buildCore("trefoil", 300);
        const std::vector<Vector3r> figureEightCore = buildCore("figure8", 300);

        const Vector3r trefoilOffset(static_cast<Real>(-16.0), knotLift, static_cast<Real>(0.0));
        const Vector3r figureEightOffset(static_cast<Real>(16.0), knotLift, static_cast<Real>(0.0));

        const int leadCount = 24;
        const Real leadLength = static_cast<Real>(5.0);
        const int spacerCount = 80;

        std::vector<Vector3r> trefoilPlaced;
        appendOffset(trefoilPlaced, trefoilCore, trefoilOffset);
        std::vector<Vector3r> figureEightPlaced;
        appendOffset(figureEightPlaced, figureEightCore, figureEightOffset);

        const Vector3r leadStartDirection = (trefoilPlaced[1] - trefoilPlaced[0]).normalized();
        const Vector3r bridgeStart = trefoilPlaced.back();
        const Vector3r bridgeEnd = figureEightPlaced.front();
        const Vector3r endDirection = (figureEightPlaced.back() - figureEightPlaced[figureEightPlaced.size() - 2]).normalized();

        ropePoints.reserve(
            trefoilPlaced.size() +
            figureEightPlaced.size() +
            static_cast<size_t>(leadCount) * 2 +
            static_cast<size_t>(spacerCount));

        for (int i = leadCount; i >= 1; --i)
        {
            const Real alpha = static_cast<Real>(i) / static_cast<Real>(leadCount);
            ropePoints.push_back(trefoilPlaced.front() - leadStartDirection * (leadLength * alpha));
        }

        ropePoints.insert(ropePoints.end(), trefoilPlaced.begin(), trefoilPlaced.end());

        for (int i = 1; i <= spacerCount; ++i)
        {
            const Real alpha = static_cast<Real>(i) / static_cast<Real>(spacerCount + 1);
            ropePoints.push_back((static_cast<Real>(1.0) - alpha) * bridgeStart + alpha * bridgeEnd);
        }

        ropePoints.insert(ropePoints.end(), figureEightPlaced.begin(), figureEightPlaced.end());

        for (int i = 1; i <= leadCount; ++i)
        {
            const Real alpha = static_cast<Real>(i) / static_cast<Real>(leadCount);
            ropePoints.push_back(figureEightPlaced.back() + endDirection * (leadLength * alpha));
        }

        return ropePoints;
    }

    std::vector<Vector3r> corePoints = buildCore(knotName, 300);
    for (Vector3r& p : corePoints)
        p.y() += knotLift;

    const int leadCount = 18;
    const Real leadLength = static_cast<Real>(3.0);
    const Vector3r startDirection =
        (corePoints[1] - corePoints[0]).normalized();
    const Vector3r endDirection =
        (corePoints[corePoints.size() - 1] - corePoints[corePoints.size() - 2]).normalized();

    ropePoints.reserve(corePoints.size() + 2 * static_cast<size_t>(leadCount));

    for (int i = leadCount; i >= 1; --i)
    {
        const Real alpha = static_cast<Real>(i) / static_cast<Real>(leadCount);
        ropePoints.push_back(corePoints.front() - startDirection * (leadLength * alpha));
    }

    ropePoints.insert(ropePoints.end(), corePoints.begin(), corePoints.end());

    for (int i = 1; i <= leadCount; ++i)
    {
        const Real alpha = static_cast<Real>(i) / static_cast<Real>(leadCount);
        ropePoints.push_back(corePoints.back() + endDirection * (leadLength * alpha));
    }

    return ropePoints;
}

static void ensureCableSnapshotLibrary()
{
    static bool initialized = false;
    if (initialized || base == nullptr)
        return;
    initialized = true;

    const std::string snapshotDirectory = getCableSnapshotDirectory();
    Utilities::FileSystem::makeDirs(snapshotDirectory);

    const std::vector<std::pair<std::string, std::vector<Vector3r>>> sampleSnapshots =
    {
        { "mathematical_trefoil.txt", generateMathematicalKnotPolyline("trefoil") },
        { "mathematical_figure8.txt", generateMathematicalKnotPolyline("figure8") },
        { "mathematical_multi.txt", generateMathematicalKnotPolyline("multi") }
    };

    for (const auto& sample : sampleSnapshots)
    {
        const std::string snapshotPath =
            Utilities::FileSystem::normalizePath(snapshotDirectory + "/" + sample.first);
        writeCableSnapshotFile(snapshotPath, sample.second);
    }

    std::vector<std::string> fileNames;
    if (!Utilities::FileSystem::getFilesInDirectory(snapshotDirectory, fileNames))
        return;

    std::sort(fileNames.begin(), fileNames.end());

    std::vector<std::string> orderedNames;
    orderedNames.reserve(fileNames.size());

    const std::vector<std::string> preferredOrder =
    {
        "mathematical_trefoil.txt",
        "mathematical_figure8.txt",
        "mathematical_multi.txt"
    };

    for (const std::string& preferred : preferredOrder)
    {
        if (std::find(fileNames.begin(), fileNames.end(), preferred) != fileNames.end())
            orderedNames.push_back(preferred);
    }

    for (const std::string& fileName : fileNames)
    {
        if (fileName == "." || fileName == ".." || fileName == "snapshotlist.txt")
            continue;
        if (Utilities::FileSystem::getFileExt(fileName) != "txt")
            continue;
        if (std::find(orderedNames.begin(), orderedNames.end(), fileName) != orderedNames.end())
            continue;
        orderedNames.push_back(fileName);
    }

    std::ofstream listOut(getCableSnapshotListPath().c_str(), std::ios::out | std::ios::trunc);
    if (!listOut.is_open())
        return;

    listOut << "# Snapshot list used by the MyCableDemo snapshot cycler.\n";
    listOut << "# The first entries are synthetic mathematical knots.\n";
    for (const std::string& fileName : orderedNames)
        listOut << fileName << "\n";
}

static void printSnapshotRecognitionReport(
    const std::string& snapshotPath,
    const std::vector<Vector3r>& ropePoints)
{
    const std::string bridgeRoot =
        (base != nullptr)
            ? Utilities::FileSystem::normalizePath(base->getOutputPath() + "/kymoknot_bridge")
            : "kymoknot_bridge";
    const PBD::KymoKnotAnalysisResult knotResult =
        PBD::KymoKnotAnalyzer::analyzeLinear(
            ropePoints,
            bridgeRoot,
            g_kymoKnotSearchMode);

    std::cout << "Snapshot " << getSnapshotDisplayName(snapshotPath) << ": ";
    if (!knotResult.available)
    {
        g_hasKymoKnotResult = false;
        std::cout << knotResult.message << std::endl;
        return;
    }

    cacheKymoKnotResult(knotResult);
    std::cout << (knotResult.success ? "recognized" : "not recognized") << std::endl;
    if (!knotResult.success)
    {
        std::cout << "  KymoKnot: " << knotResult.message << std::endl;
        return;
    }

    for (size_t i = 0; i < knotResult.segments.size(); ++i)
    {
        const PBD::KymoKnotSegment& segment = knotResult.segments[i];
        std::cout << "  [" << i << "] "
                  << segment.knotType
                  << " start=" << segment.start
                  << " end=" << segment.end
                  << " length=" << segment.length
                  << std::endl;
    }
}

static void auditSnapshotLibrary()
{
    const std::vector<std::string> snapshotPaths = getAvailableSnapshotPaths();
    if (snapshotPaths.empty())
    {
        std::cout << "Snapshot audit: no snapshot files found." << std::endl;
        return;
    }

    std::cout << "Snapshot audit:" << std::endl;
    for (const std::string& snapshotPath : snapshotPaths)
    {
        std::vector<Vector3r> ropePoints;
        if (!loadCableSnapshot(snapshotPath, ropePoints))
        {
            std::cout << "  " << getSnapshotDisplayName(snapshotPath)
                      << ": failed to load" << std::endl;
            continue;
        }

        printSnapshotRecognitionReport(snapshotPath, ropePoints);
    }
}

static void applyBendTwistRestStateCreep(
    SimulationModel* model,
    const RestShapeCreepConfig& creepCfg)
{
    if (!creepCfg.enabled || model == nullptr)
        return;

    g_restCreepCounter++;
    if (creepCfg.stepInterval == 0 || (g_restCreepCounter % creepCfg.stepInterval) != 0)
        return;

    OrientationData& orientations = model->getOrientations();
    SimulationModel::ConstraintVector& constraints = model->getConstraints();

    for (Constraint* constraint : constraints)
    {
        if (constraint == nullptr || constraint->getTypeId() != BendTwistConstraint::TYPE_ID)
            continue;

        BendTwistConstraint* bendTwist = static_cast<BendTwistConstraint*>(constraint);
        const unsigned int q0Index = bendTwist->m_bodies[0];
        const unsigned int q1Index = bendTwist->m_bodies[1];
        if (q0Index >= orientations.getNumberOfQuaternions() ||
            q1Index >= orientations.getNumberOfQuaternions())
        {
            continue;
        }

        Quaternionr currentRelative =
            orientations.getQuaternion(q0Index).conjugate() *
            orientations.getQuaternion(q1Index);
        currentRelative.normalize();

        Quaternionr& restRelative = bendTwist->m_restDarbouxVector;
        if (restRelative.dot(currentRelative) < static_cast<Real>(0.0))
            currentRelative.coeffs() *= static_cast<Real>(-1.0);

        const Real mismatch =
            (restRelative.coeffs() - currentRelative.coeffs()).squaredNorm();
        if (mismatch < creepCfg.activationThreshold)
            continue;

        Eigen::Matrix<Real, 4, 1> blendedCoeffs =
            (static_cast<Real>(1.0) - creepCfg.creep) * restRelative.coeffs() +
            creepCfg.creep * currentRelative.coeffs();

        const Real norm = blendedCoeffs.norm();
        if (norm > static_cast<Real>(1e-12))
            restRelative.coeffs() = blendedCoeffs / norm;
    }
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

        applyBendTwistRestStateCreep(model, g_cfg.restCreep);

        base->step();
    }
}

static void applyPullForces(SimulationModel* model, Real dt, const PhysicsConfig& physicsCfg)
{
    ParticleData& particles = model->getParticles();

    // Targets act like damped springs at the velocity level, so they steer
    // the ends without hard-position pinning them.
    auto pullOneParticle = [&](unsigned int particleIndex, const Vector3r& target)
    {
        if (particleIndex >= particles.getNumberOfParticles())
            return;
        if (particles.getMass(particleIndex) <= 0.0)
            return;

        const Vector3r position = particles.getPosition(particleIndex);
        Vector3r& velocity = particles.getVelocity(particleIndex);
        Vector3r acceleration =
            physicsCfg.pullSpring * (target - position) -
            physicsCfg.pullDamping * velocity;

        const Real accelNorm = acceleration.norm();
        if (accelNorm > physicsCfg.pullMaxAccel)
            acceleration *= (physicsCfg.pullMaxAccel / accelNorm);

        velocity += dt * acceleration;
    };

    pullOneParticle(g_firstParticle, g_targetA);
    pullOneParticle(g_lastParticle, g_targetB);
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

static std::vector<SelfCollisionContact> detectSelfCollisionContacts(
    SimulationModel* model,
    const RopeConfig& ropeCfg,
    const SelfCollisionConfig& collisionCfg)
{
    std::vector<SelfCollisionContact> contacts;
    if (model == nullptr)
        return contacts;

    ParticleData& pd = model->getParticles();
    const unsigned int n = pd.getNumberOfParticles();
    if (n < 4)
        return contacts;

    const Real minDist = static_cast<Real>(2.0) * collisionCfg.ropeRadius;
    const Real eps = static_cast<Real>(1e-8);
    const std::vector<std::pair<unsigned int, unsigned int>> candidates =
        collectSelfCollisionCandidatesSpatialHash(
            pd,
            minDist,
            collisionCfg.hashCellSize);

    contacts.reserve(candidates.size());
    for (const auto& pair : candidates)
    {
        const unsigned int i = pair.first;
        const unsigned int j = pair.second;

        const Real arcDist = static_cast<Real>(std::abs((int)j - (int)i)) * ropeCfg.dx;
        if (arcDist < collisionCfg.excludeArcLength)
            continue;

        const Vector3r p0 = pd.getPosition(i);
        const Vector3r p1 = pd.getPosition(i + 1);
        const Vector3r q0 = pd.getPosition(j);
        const Vector3r q1 = pd.getPosition(j + 1);

        Real s, t;
        Vector3r c1, c2;
        closestPtSegmentSegment(p0, p1, q0, q1, s, t, c1, c2);

        const Vector3r delta = c2 - c1;
        const Real dist = delta.norm();
        if (dist - minDist < static_cast<Real>(0.0) || dist <= eps)
            contacts.push_back({ i, j });
    }

    return contacts;
}

static unsigned int countSelfCollisionContactClusters(
    const std::vector<SelfCollisionContact>& contacts)
{
    const size_t n = contacts.size();
    if (n == 0)
        return 0;

    std::vector<unsigned char> visited(n, 0);
    unsigned int clusterCount = 0;

    auto contactsShareEdge = [](const SelfCollisionContact& a, const SelfCollisionContact& b)
    {
        return a.edgeA == b.edgeA ||
               a.edgeA == b.edgeB ||
               a.edgeB == b.edgeA ||
               a.edgeB == b.edgeB;
    };

    for (size_t start = 0; start < n; ++start)
    {
        if (visited[start] != 0)
            continue;

        clusterCount++;
        std::vector<size_t> stack(1, start);
        visited[start] = 1;

        while (!stack.empty())
        {
            const size_t current = stack.back();
            stack.pop_back();

            for (size_t next = 0; next < n; ++next)
            {
                if (visited[next] != 0)
                    continue;
                if (!contactsShareEdge(contacts[current], contacts[next]))
                    continue;

                visited[next] = 1;
                stack.push_back(next);
            }
        }
    }

    return clusterCount;
}

static void printSelfCollisionContacts(
    const std::vector<SelfCollisionContact>& contacts)
{
    std::cout << "Self-collision contacts (" << contacts.size() << "):";
    if (contacts.empty())
    {
        std::cout << " none" << std::endl;
        return;
    }

    std::cout << std::endl;
    for (const SelfCollisionContact& contact : contacts)
    {
        std::cout << "  (" << contact.edgeA << ", " << contact.edgeB << ")"
                  << std::endl;
    }
}

static std::vector<Vector3r> getCurrentCablePoints(SimulationModel* model)
{
    std::vector<Vector3r> points;
    if (model == nullptr || g_lastParticle < g_firstParticle)
        return points;

    ParticleData& particles = model->getParticles();
    const unsigned int count = g_lastParticle - g_firstParticle + 1;
    points.reserve(count);
    for (unsigned int i = 0; i < count; ++i)
        points.push_back(particles.getPosition(g_firstParticle + i));

    return points;
}

static bool isLocalIndexInKymoKnotSegment(
    int localIndex,
    const PBD::KymoKnotSegment& segment)
{
    if (localIndex < 0 || segment.start < 0 || segment.end < 0)
        return false;

    if (segment.start <= segment.end)
        return localIndex >= segment.start && localIndex <= segment.end;

    return localIndex >= segment.start || localIndex <= segment.end;
}

static bool isLocalIndexHighlighted(int localIndex)
{
    if (!g_hasKymoKnotResult)
        return false;

    for (const PBD::KymoKnotSegment& segment : g_lastKymoKnotResult.segments)
    {
        if (isLocalIndexInKymoKnotSegment(localIndex, segment))
            return true;
    }
    return false;
}

static bool isLocalEdgeHighlighted(int localEdgeIndex)
{
    return isLocalIndexHighlighted(localEdgeIndex) ||
           isLocalIndexHighlighted(localEdgeIndex + 1);
}

static void cacheKymoKnotResult(const PBD::KymoKnotAnalysisResult& result)
{
    g_lastKymoKnotResult = result;
    g_hasKymoKnotResult = result.available && !result.segments.empty();
}

static const char* getKymoKnotSearchModeName(PBD::KymoKnotSearchMode mode)
{
    return (mode == PBD::KymoKnotSearchMode::BottomUp) ? "bottom-up" : "top-down";
}

static void toggleKymoKnotSearchMode()
{
    g_kymoKnotSearchMode =
        (g_kymoKnotSearchMode == PBD::KymoKnotSearchMode::TopDown)
            ? PBD::KymoKnotSearchMode::BottomUp
            : PBD::KymoKnotSearchMode::TopDown;

    std::cout << "KymoKnot search mode: "
              << getKymoKnotSearchModeName(g_kymoKnotSearchMode)
              << std::endl;
}

static void printTanglementMetrics(
    SimulationModel* model,
    const RopeConfig& ropeCfg,
    const SelfCollisionConfig& collisionCfg)
{
    const std::vector<SelfCollisionContact> contacts =
        detectSelfCollisionContacts(model, ropeCfg, collisionCfg);
    const unsigned int clusterCount = countSelfCollisionContactClusters(contacts);

    std::cout << "Tanglement metrics:"
              << " self-collisions=" << contacts.size()
              << ", collision-blobs=" << clusterCount
              << std::endl;

    const std::vector<Vector3r> cablePoints = getCurrentCablePoints(model);
    const std::string bridgeRoot =
        (base != nullptr)
            ? Utilities::FileSystem::normalizePath(base->getOutputPath() + "/kymoknot_bridge")
            : "kymoknot_bridge";
    const PBD::KymoKnotAnalysisResult knotResult =
        PBD::KymoKnotAnalyzer::analyzeLinear(
            cablePoints,
            bridgeRoot,
            g_kymoKnotSearchMode);

    if (!knotResult.available)
    {
        g_hasKymoKnotResult = false;
        std::cout << "KymoKnot: " << knotResult.message << std::endl;
        return;
    }

    cacheKymoKnotResult(knotResult);
    std::cout << "KymoKnot segments: " << knotResult.segments.size() << std::endl;
    for (size_t i = 0; i < knotResult.segments.size(); ++i)
    {
        const PBD::KymoKnotSegment& segment = knotResult.segments[i];
        std::cout << "  [" << i << "] "
                  << segment.knotType
                  << " start=" << segment.start
                  << " end=" << segment.end
                  << " length=" << segment.length
                  << std::endl;
    }
}

static void getHeightColor(Real y, float* color)
{
    const Real minY = g_cfg.groundY - static_cast<Real>(0.25);
    const Real maxY = std::max(
        g_cfg.layout.height + static_cast<Real>(2.0),
        g_cfg.groundY + static_cast<Real>(1.0));
    const Real t = clamp01((y - minY) / std::max(maxY - minY, static_cast<Real>(1e-6)));

    const Vector3r low(0.08, 0.35, 0.85);
    const Vector3r mid(0.15, 0.75, 0.55);
    const Vector3r high(0.95, 0.72, 0.12);

    Vector3r rgb;
    if (t < static_cast<Real>(0.5))
    {
        const Real localT = static_cast<Real>(2.0) * t;
        rgb = (static_cast<Real>(1.0) - localT) * low + localT * mid;
    }
    else
    {
        const Real localT = static_cast<Real>(2.0) * (t - static_cast<Real>(0.5));
        rgb = (static_cast<Real>(1.0) - localT) * mid + localT * high;
    }

    color[0] = static_cast<float>(rgb.x());
    color[1] = static_cast<float>(rgb.y());
    color[2] = static_cast<float>(rgb.z());
    color[3] = 1.0f;
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
    // `base->render()` is the built-in demo GUI + scene helper rendering.
    base->render();

    SimulationModel* model = Simulation::getCurrent()->getModel();
    const ParticleData& particles = model->getParticles();

    float yellow[4] = {1.0f, 0.9f, 0.1f, 1.0f};
    float magenta[4] = {0.9f, 0.1f, 0.9f, 1.0f};
    float segmentColor[4] = {0.0f, 0.0f, 0.8f, 1.0f};
    float knotColor[4] = {0.98f, 0.25f, 0.18f, 1.0f};
    float knotSphereColor[4] = {1.0f, 0.45f, 0.20f, 1.0f};

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
            const Real avgY = static_cast<Real>(0.5) * (p0.y() + p1.y());
            getHeightColor(avgY, segmentColor);

            const int localEdgeIndex =
                static_cast<int>(i0) - static_cast<int>(g_firstParticle);
            const bool highlightKnot = isLocalEdgeHighlighted(localEdgeIndex);
            const float* activeColor = highlightKnot ? knotColor : segmentColor;

            MiniGL::drawSphere(p0, 0.03f, highlightKnot ? knotSphereColor : segmentColor);
            if (e == lineModel->getEdges().size() - 1)
                MiniGL::drawSphere(p1, 0.03f, isLocalIndexHighlighted(static_cast<int>(i1) - static_cast<int>(g_firstParticle)) ? knotSphereColor : segmentColor);

            MiniGL::drawCylinder(
                p0,
                p1,
                activeColor,
                static_cast<float>(g_cfg.rope.renderRadius));
        }
    }

    MiniGL::drawSphere(g_targetA, 0.08f, yellow);
    MiniGL::drawSphere(g_targetB, 0.08f, magenta);
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
