#include "Common/Common.h"
#include "Demos/Visualization/MiniGL.h"
#include "Simulation/TimeManager.h"
#include "Simulation/SimulationModel.h"
#include "Simulation/Constraints.h"
#include "Demos/Common/DemoBase.h"
#include "Simulation/Simulation.h"
#include "Utils/FileSystem.h"

#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
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

// Endpoints for the straight cable layout.
struct CableLayoutConfig
{
    Real endpointX = static_cast<Real>(10.0);
    Real height = static_cast<Real>(6.0);
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
    RestShapeCreepConfig restCreep;
    CableLayoutConfig layout;
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

static void buildModel();
static void reset();
static void timeStep();
static void render();

static std::vector<Vector3r> generateStraightCablePolyline(
    const CableLayoutConfig& layoutCfg,
    const RopeConfig& ropeCfg);

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

static std::string getCableSnapshotPath();
static bool saveCurrentCableSnapshot();
static bool loadCableSnapshot(std::vector<Vector3r>& ropePoints);
static void applyBendTwistRestStateCreep(
    SimulationModel* model,
    const RestShapeCreepConfig& creepCfg);

static DemoBase* base = nullptr;
static AppConfig g_cfg;
static unsigned int g_firstParticle = 0;
static unsigned int g_lastParticle = 0;
static Vector3r g_targetA(0.0, 0.0, 0.0);
static Vector3r g_targetB(0.0, 0.0, 0.0);
static bool g_useSavedCableSnapshot = false;
static unsigned int g_restCreepCounter = 0;

int main(int argc, char** argv)
{
    std::cout << "Straight cable demo" << std::endl;

    base = new DemoBase();

    base->init(argc, argv, "Straight Cable Demo");

    SimulationModel* model = new SimulationModel();
    model->init();
    Simulation::getCurrent()->setModel(model);

    buildModel();

    // Keyboard input for this demo lives here.
    MiniGL::setClientIdleFunc(timeStep);
    MiniGL::setClientSceneFunc(render);
    MiniGL::addKeyFunc('r', reset);
    MiniGL::addKeyFunc('o', []() {
        if (saveCurrentCableSnapshot())
        {
            g_useSavedCableSnapshot = true;
            std::cout << "Saved cable snapshot. Reset now rebuilds from saved data." << std::endl;
        }
    });
    MiniGL::addKeyFunc('p', []() {
        g_useSavedCableSnapshot = false;
        reset();
        std::cout << "Snapshot replay disabled. Rebuilt straight cable." << std::endl;
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

static void buildModel()
{
    TimeManager::getCurrent()->setTimeStepSize(g_cfg.physics.timeStep);

    SimulationModel* model = Simulation::getCurrent()->getModel();
    std::vector<Vector3r> ropePoints;
    if (!g_useSavedCableSnapshot || !loadCableSnapshot(ropePoints))
    {
        ropePoints = generateStraightCablePolyline(g_cfg.layout, g_cfg.rope);
        if (g_useSavedCableSnapshot)
        {
            std::cout << "Saved cable snapshot missing or invalid. Falling back to straight cable."
                      << std::endl;
            g_useSavedCableSnapshot = false;
        }
    }

    addCableModel(model, ropePoints, g_cfg.rope);

    if (g_firstParticle < model->getParticles().getNumberOfParticles() &&
        g_lastParticle < model->getParticles().getNumberOfParticles())
    {
        g_targetA = model->getParticles().getPosition(g_firstParticle);
        g_targetB = model->getParticles().getPosition(g_lastParticle);
    }

    std::cout << "Built cable with "
              << model->getParticles().getNumberOfParticles()
              << " particles." << std::endl;
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

static std::string getCableSnapshotPath()
{
    return Utilities::FileSystem::normalizePath(
        base->getOutputPath() + "/saved_cable/snapshot.txt");
}

static bool saveCurrentCableSnapshot()
{
    if (base == nullptr)
        return false;

    SimulationModel* model = Simulation::getCurrent()->getModel();
    if (model == nullptr)
        return false;

    const ParticleData& particles = model->getParticles();
    if (g_firstParticle >= particles.getNumberOfParticles() ||
        g_lastParticle >= particles.getNumberOfParticles() ||
        g_firstParticle > g_lastParticle)
    {
        return false;
    }

    const std::string snapshotPath = getCableSnapshotPath();
    Utilities::FileSystem::makeDirs(Utilities::FileSystem::getFilePath(snapshotPath));

    std::ofstream out(snapshotPath);
    if (!out.is_open())
    {
        std::cout << "Failed to open snapshot file: " << snapshotPath << std::endl;
        return false;
    }

    out << (g_lastParticle - g_firstParticle + 1) << "\n";
    out << std::setprecision(17);
    for (unsigned int i = g_firstParticle; i <= g_lastParticle; ++i)
    {
        const Vector3r& p = particles.getPosition(i);
        out << p.x() << " " << p.y() << " " << p.z() << "\n";
    }

    std::cout << "Cable snapshot saved to: " << snapshotPath << std::endl;
    return true;
}

static bool loadCableSnapshot(std::vector<Vector3r>& ropePoints)
{
    ropePoints.clear();

    if (base == nullptr)
        return false;

    std::ifstream in(getCableSnapshotPath());
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

    return true;
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

    float blue[4] = {0.0f, 0.0f, 0.8f, 1.0f};
    float yellow[4] = {1.0f, 0.9f, 0.1f, 1.0f};
    float magenta[4] = {0.9f, 0.1f, 0.9f, 1.0f};

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

            MiniGL::drawSphere(p0, 0.03f, blue);
            if (e == lineModel->getEdges().size() - 1)
                MiniGL::drawSphere(p1, 0.03f, blue);

            MiniGL::drawCylinder(p0, p1, blue, static_cast<float>(g_cfg.rope.renderRadius));
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
