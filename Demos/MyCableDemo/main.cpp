#include "Common/Common.h"
#include "Demos/Visualization/MiniGL.h"
#include "Simulation/TimeManager.h"
#include "Simulation/SimulationModel.h"
#include "Simulation/Constraints.h"
#include "Demos/Common/DemoBase.h"
#include "Simulation/Simulation.h"
#include "Utils/FileSystem.h"
#include "CablePhysics.h"
#include "KymoKnotAnalyzer.h"

#include <GLFW/glfw3.h>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <random>
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

extern "C" void glfwGetCursorPos(GLFWwindow* window, double* xpos, double* ypos);

static void buildModel();
static void reset();
static void regenerateRandomCable();
static void togglePullForces();
static void adjustPullStrength(Real scale);
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
static void addCableModel(
    SimulationModel* model,
    std::vector<Vector3r>& ropePoints,
    const RopeConfig& ropeCfg,
    std::pair<unsigned int, unsigned int>* particleRangeOut);
static void initializeCablePieces();
static void rebuildCableFromPieces(
    SimulationModel* model,
    const std::vector<std::vector<Vector3r>>& ropePieces);
static std::vector<std::vector<Vector3r>> collectCurrentCablePieces(
    SimulationModel* model);
static bool splitCablePieceAtEdge(
    std::vector<std::vector<Vector3r>>& ropePieces,
    size_t pieceIndex,
    size_t edgeIndex);
static bool buildPickRay(
    const ::Vector2i& screenPos,
    Vector3r& rayStart,
    Vector3r& rayEnd);
static bool selectCableEdgeAtScreenPosition(
    const ::Vector2i& screenPos,
    size_t& outPieceIndex,
    size_t& outEdgeIndex,
    Real& outRayDepthT,
    Real& outEdgeT,
    Vector3r& outGrabbedPoint);
static bool beginMouseDragAtScreenPosition(const ::Vector2i& screenPos);
static void updateMouseDragTargetAtScreenPosition(const ::Vector2i& screenPos);
static void applyMouseDragForces(SimulationModel* model, Real dt);
static bool cutSelectedCableSegment();
static void endMouseDrag();
static bool handleMousePress(int button, int action, int mods);
static bool handleMouseMove(int x, int y);

static void printTanglementMetrics(
    SimulationModel* model,
    const RopeConfig& ropeCfg,
    const SelfCollisionConfig& collisionCfg);
static void toggleKymoKnotSearchMode();
static const char* getKymoKnotSearchModeName(PBD::KymoKnotSearchMode mode);
static bool isLocalIndexInKymoKnotSegment(
    int localIndex,
    const PBD::KymoKnotSegment& segment);
static bool isLocalIndexHighlighted(size_t pieceIndex, int localIndex);
static bool isLocalEdgeHighlighted(size_t pieceIndex, int localEdgeIndex);
static void cacheKymoKnotResult(const PBD::KymoKnotAnalysisResult& result);
static void cacheCableKymoKnotResults(
    const std::vector<PBD::KymoKnotAnalysisResult>& results);
static void getHeightColor(Real y, float* color);

static std::string getCableSnapshotDirectory();
static std::string getCableSnapshotListPath();
static std::string getCableSnapshotPath();
static std::vector<std::string> getAvailableSnapshotPaths();
static bool loadCableSnapshot(
    const std::string& snapshotPath,
    std::vector<Vector3r>& ropePoints);
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
static bool g_useSavedCableSnapshot = false;
unsigned int g_restCreepCounter = 0;
static int g_snapshotCycleIndex = -1;
static PBD::KymoKnotAnalysisResult g_lastKymoKnotResult;
static bool g_hasKymoKnotResult = false;
static std::vector<PBD::KymoKnotAnalysisResult> g_cachedCableKymoKnotResults;
static bool g_hasCableKymoKnotResults = false;
static std::string g_currentSnapshotPath;
static PBD::KymoKnotSearchMode g_kymoKnotSearchMode = PBD::KymoKnotSearchMode::TopDown;
static std::vector<std::vector<Vector3r>> g_cablePieces;
std::vector<std::pair<unsigned int, unsigned int>> g_cablePieceRanges;
bool g_pullEnabled = false;

struct MouseDragState
{
    bool selected = false;
    bool active = false;
    size_t pieceIndex = std::numeric_limits<size_t>::max();
    size_t edgeIndex = std::numeric_limits<size_t>::max();
    Real edgeT = static_cast<Real>(0.0);
    Real rayDepthT = static_cast<Real>(0.0);
    Vector3r grabbedPoint = Vector3r::Zero();
    Vector3r targetPoint = Vector3r::Zero();
};

static MouseDragState g_mouseDrag;
static const Real g_mouseDragSpring = static_cast<Real>(1800.0);
static const Real g_mouseDragDamping = static_cast<Real>(24.0);
static const Real g_mouseDragMaxAccel = static_cast<Real>(3000.0);

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
    MiniGL::addMousePressFunc(handleMousePress);
    MiniGL::addMouseMoveFunc(handleMouseMove);
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
    MiniGL::addKeyFunc('c', []() {
        if (!cutSelectedCableSegment())
            std::cout << "No selected cable segment to cut" << std::endl;
    });
    MiniGL::addKeyFunc('u', []() {
        const std::vector<std::string> snapshotPaths = getAvailableSnapshotPaths();
        if (snapshotPaths.empty())
            return;

        g_snapshotCycleIndex = (g_snapshotCycleIndex + 1) %
            static_cast<int>(snapshotPaths.size());
        g_useSavedCableSnapshot = true;
        g_cablePieces.clear();
        reset();
    });
    MiniGL::addKeyFunc('p', []() {
        g_useSavedCableSnapshot = false;
        g_snapshotCycleIndex = -1;
        g_cablePieces.clear();
        reset();
    });
    MiniGL::addKeyFunc('o', []() {
        g_cablePieces.clear();
        regenerateRandomCable();
    });
    MiniGL::addKeyFunc('*', []() {
        togglePullForces();
    });
    MiniGL::addKeyFunc('-', []() {
        adjustPullStrength(static_cast<Real>(0.85));
    });
    MiniGL::addKeyFunc('+', []() {
        adjustPullStrength(static_cast<Real>(1.15));
    });
    MiniGL::addKeyFunc('z', []() { MiniGL::move(0, 0, -1); });
    MiniGL::addKeyFunc('x', []() { MiniGL::move(0, 0, 1); });

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
    g_mouseDrag = MouseDragState();
    g_cachedCableKymoKnotResults.clear();
    g_hasCableKymoKnotResults = false;
    g_cablePieces.clear();
    Simulation::getCurrent()->reset();
    base->getSelectedParticles().clear();

    Simulation::getCurrent()->getModel()->cleanup();
    buildModel();
}

static void regenerateRandomCable()
{
    g_mouseDrag = MouseDragState();
    g_useSavedCableSnapshot = false;
    g_snapshotCycleIndex = -1;
    g_cfg.randomCable.enabled = true;
    g_cablePieces.clear();

    std::random_device rd;
    unsigned int newSeed = rd();
    if (newSeed == 0)
        newSeed = 1;
    g_cfg.randomCable.seed = newSeed;

    std::cout << "Regenerating random cable with seed " << g_cfg.randomCable.seed << std::endl;
    reset();
}

static void togglePullForces()
{
    g_pullEnabled = !g_pullEnabled;
    std::cout << "End pull force: " << (g_pullEnabled ? "on" : "off") << std::endl;
}

static void adjustPullStrength(Real scale)
{
    g_cfg.physics.pullSpring = std::max(
        static_cast<Real>(0.0),
        g_cfg.physics.pullSpring * scale);
    std::cout << "End pull strength: " << g_cfg.physics.pullSpring << std::endl;
}

static void initializeCablePieces()
{
    g_cablePieces.clear();

    std::vector<Vector3r> ropePoints;
    if (g_useSavedCableSnapshot && loadCableSnapshot(getCableSnapshotPath(), ropePoints))
    {
        g_cablePieces.push_back(std::move(ropePoints));
        return;
    }

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

    g_cablePieces.push_back(std::move(ropePoints));
}

static void rebuildCableFromPieces(
    SimulationModel* model,
    const std::vector<std::vector<Vector3r>>& ropePieces)
{
    if (model == nullptr)
        return;

    g_cablePieceRanges.clear();
    g_firstParticle = 0;
    g_lastParticle = 0;

    for (const std::vector<Vector3r>& piece : ropePieces)
    {
        if (piece.size() < 2)
            continue;

        std::vector<Vector3r> copy = piece;
        std::pair<unsigned int, unsigned int> range(0, 0);
        addCableModel(model, copy, g_cfg.rope, &range);
        g_cablePieceRanges.push_back(range);
    }

    if (!g_cablePieceRanges.empty())
    {
        g_firstParticle = g_cablePieceRanges.front().first;
        g_lastParticle = g_cablePieceRanges.back().second;
    }
}

static std::vector<std::vector<Vector3r>> collectCurrentCablePieces(
    SimulationModel* model)
{
    std::vector<std::vector<Vector3r>> ropePieces;
    if (model == nullptr)
        return ropePieces;

    const ParticleData& particles = model->getParticles();
    for (LineModel* lineModel : model->getLineModels())
    {
        if (lineModel == nullptr)
            continue;

        const LineModel::Edges& edges = lineModel->getEdges();
        if (edges.empty())
            continue;

        const unsigned int pointOffset = lineModel->getIndexOffset();
        const unsigned int pointCount = static_cast<unsigned int>(edges.size()) + 1u;
        if (pointOffset + pointCount > particles.getNumberOfParticles())
            continue;

        std::vector<Vector3r> piece;
        piece.reserve(pointCount);
        for (unsigned int i = 0; i < pointCount; ++i)
            piece.push_back(particles.getPosition(pointOffset + i));
        ropePieces.push_back(std::move(piece));
    }

    return ropePieces;
}

static bool splitCablePieceAtEdge(
    std::vector<std::vector<Vector3r>>& ropePieces,
    size_t pieceIndex,
    size_t edgeIndex)
{
    if (pieceIndex >= ropePieces.size())
        return false;

    std::vector<Vector3r>& piece = ropePieces[pieceIndex];
    if (piece.size() < 3 || edgeIndex + 1 >= piece.size())
        return false;

    std::vector<Vector3r> left(
        piece.begin(),
        piece.begin() + static_cast<std::ptrdiff_t>(edgeIndex + 1));
    std::vector<Vector3r> right(
        piece.begin() + static_cast<std::ptrdiff_t>(edgeIndex + 1),
        piece.end());
    if (left.size() < 2 || right.size() < 2)
        return false;

    piece = std::move(left);
    ropePieces.insert(
        ropePieces.begin() + static_cast<std::ptrdiff_t>(pieceIndex + 1),
        std::move(right));
    return true;
}

static bool buildPickRay(
    const ::Vector2i& screenPos,
    Vector3r& rayStart,
    Vector3r& rayEnd)
{
    GLint viewport[4];
    GLdouble mv[16];
    GLdouble pm[16];
    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_MODELVIEW_MATRIX, mv);
    glGetDoublev(GL_PROJECTION_MATRIX, pm);

    Matrix4r modelView;
    Matrix4r projection;
    for (int c = 0; c < 4; ++c)
    {
        for (int r = 0; r < 4; ++r)
        {
            modelView(r, c) = static_cast<Real>(mv[c * 4 + r]);
            projection(r, c) = static_cast<Real>(pm[c * 4 + r]);
        }
    }

    const Matrix4r invMVP = (projection * modelView).inverse();
    const Real x = (static_cast<Real>(screenPos[0] - viewport[0]) /
        static_cast<Real>(std::max(1, viewport[2])) * static_cast<Real>(2.0)) - static_cast<Real>(1.0);
    const Real y = (static_cast<Real>(viewport[3] - screenPos[1] - viewport[1]) /
        static_cast<Real>(std::max(1, viewport[3])) * static_cast<Real>(2.0)) - static_cast<Real>(1.0);

    const Vector4r nearClip(x, y, static_cast<Real>(-1.0), static_cast<Real>(1.0));
    const Vector4r farClip(x, y, static_cast<Real>(1.0), static_cast<Real>(1.0));
    Vector4r nearWorld = invMVP * nearClip;
    Vector4r farWorld = invMVP * farClip;

    if (std::abs(nearWorld.w()) < static_cast<Real>(1e-8) ||
        std::abs(farWorld.w()) < static_cast<Real>(1e-8))
    {
        return false;
    }

    nearWorld /= nearWorld.w();
    farWorld /= farWorld.w();
    rayStart = nearWorld.head<3>();
    rayEnd = farWorld.head<3>();
    return true;
}

static bool selectCableEdgeAtScreenPosition(
    const ::Vector2i& screenPos,
    size_t& outPieceIndex,
    size_t& outEdgeIndex,
    Real& outRayDepthT,
    Real& outEdgeT,
    Vector3r& outGrabbedPoint)
{
    SimulationModel* model = Simulation::getCurrent()->getModel();
    if (model == nullptr)
        return false;

    Vector3r rayStart(0.0, 0.0, 0.0);
    Vector3r rayEnd(0.0, 0.0, 0.0);
    if (!buildPickRay(screenPos, rayStart, rayEnd))
        return false;

    const std::vector<std::vector<Vector3r>> ropePieces = collectCurrentCablePieces(model);
    if (ropePieces.empty())
        return false;

    size_t bestPiece = std::numeric_limits<size_t>::max();
    size_t bestEdge = std::numeric_limits<size_t>::max();
    Real bestDistance = std::numeric_limits<Real>::max();
    Real bestRayT = static_cast<Real>(0.0);
    Real bestEdgeT = static_cast<Real>(0.0);
    Vector3r bestGrabbedPoint = Vector3r::Zero();

    for (size_t pieceIndex = 0; pieceIndex < ropePieces.size(); ++pieceIndex)
    {
        const std::vector<Vector3r>& piece = ropePieces[pieceIndex];
        for (size_t edgeIndex = 0; edgeIndex + 1 < piece.size(); ++edgeIndex)
        {
            Real s = static_cast<Real>(0.0);
            Real t = static_cast<Real>(0.0);
            Vector3r c1, c2;
            closestPtSegmentSegment(
                rayStart,
                rayEnd,
                piece[edgeIndex],
                piece[edgeIndex + 1],
                s,
                t,
                c1,
                c2);

            const Real distance = (c1 - c2).norm();
            if (distance < bestDistance)
            {
                bestDistance = distance;
                bestPiece = pieceIndex;
                bestEdge = edgeIndex;
                bestRayT = s;
                bestEdgeT = t;
                bestGrabbedPoint = c2;
            }
        }
    }

    const Real dragThreshold = std::max(
        static_cast<Real>(0.15),
        static_cast<Real>(2.0) * g_cfg.rope.renderRadius);
    if (bestPiece == std::numeric_limits<size_t>::max() || bestDistance > dragThreshold)
        return false;

    outPieceIndex = bestPiece;
    outEdgeIndex = bestEdge;
    outRayDepthT = bestRayT;
    outEdgeT = bestEdgeT;
    outGrabbedPoint = bestGrabbedPoint;
    return true;
}

static bool beginMouseDragAtScreenPosition(const ::Vector2i& screenPos)
{
    size_t pieceIndex = std::numeric_limits<size_t>::max();
    size_t edgeIndex = std::numeric_limits<size_t>::max();
    Real rayDepthT = static_cast<Real>(0.0);
    Real edgeT = static_cast<Real>(0.0);
    Vector3r grabbedPoint = Vector3r::Zero();
    if (!selectCableEdgeAtScreenPosition(
            screenPos,
            pieceIndex,
            edgeIndex,
            rayDepthT,
            edgeT,
            grabbedPoint))
    {
        return false;
    }

    g_mouseDrag.active = true;
    g_mouseDrag.selected = true;
    g_mouseDrag.pieceIndex = pieceIndex;
    g_mouseDrag.edgeIndex = edgeIndex;
    g_mouseDrag.edgeT = edgeT;
    g_mouseDrag.rayDepthT = rayDepthT;
    g_mouseDrag.grabbedPoint = grabbedPoint;
    g_mouseDrag.targetPoint = grabbedPoint;
    return true;
}

static void updateMouseDragTargetAtScreenPosition(const ::Vector2i& screenPos)
{
    if (!g_mouseDrag.active)
        return;

    Vector3r rayStart(0.0, 0.0, 0.0);
    Vector3r rayEnd(0.0, 0.0, 0.0);
    if (!buildPickRay(screenPos, rayStart, rayEnd))
        return;

    g_mouseDrag.targetPoint = rayStart + g_mouseDrag.rayDepthT * (rayEnd - rayStart);
}

static void applyMouseDragForces(SimulationModel* model, Real dt)
{
    if (!g_mouseDrag.active || model == nullptr)
        return;
    if (g_mouseDrag.pieceIndex >= g_cablePieceRanges.size())
        return;

    const std::pair<unsigned int, unsigned int> range =
        g_cablePieceRanges[g_mouseDrag.pieceIndex];

    if (range.second <= range.first)
        return;

    const unsigned int localEdgeIndex = static_cast<unsigned int>(g_mouseDrag.edgeIndex);
    const unsigned int i0 = range.first + localEdgeIndex;
    const unsigned int i1 = i0 + 1;

    ParticleData& particles = model->getParticles();
    if (i1 >= particles.getNumberOfParticles())
        return;

    const Real w0 = static_cast<Real>(1.0) - g_mouseDrag.edgeT;
    const Real w1 = g_mouseDrag.edgeT;

    const Vector3r p0 = particles.getPosition(i0);
    const Vector3r p1 = particles.getPosition(i1);
    const Vector3r v0 = particles.getVelocity(i0);
    const Vector3r v1 = particles.getVelocity(i1);

    const Vector3r currentPoint = w0 * p0 + w1 * p1;
    const Vector3r currentVelocity = w0 * v0 + w1 * v1;

    Vector3r acceleration =
        g_mouseDragSpring * (g_mouseDrag.targetPoint - currentPoint) -
        g_mouseDragDamping * currentVelocity;

    const Real accelNorm = acceleration.norm();
    if (accelNorm > g_mouseDragMaxAccel)
        acceleration *= (g_mouseDragMaxAccel / accelNorm);

    if (particles.getMass(i0) > 0.0)
        particles.getVelocity(i0) += dt * w0 * acceleration;
    if (particles.getMass(i1) > 0.0)
        particles.getVelocity(i1) += dt * w1 * acceleration;
}

static void endMouseDrag()
{
    g_mouseDrag.active = false;
}

static bool cutSelectedCableSegment()
{
    if (!g_mouseDrag.selected ||
        g_mouseDrag.pieceIndex == std::numeric_limits<size_t>::max() ||
        g_mouseDrag.edgeIndex == std::numeric_limits<size_t>::max())
    {
        return false;
    }

    SimulationModel* model = Simulation::getCurrent()->getModel();
    if (model == nullptr)
        return false;

    const size_t pieceIndex = g_mouseDrag.pieceIndex;
    const size_t edgeIndex = g_mouseDrag.edgeIndex;
    std::vector<std::vector<Vector3r>> splitPieces = collectCurrentCablePieces(model);
    if (!splitCablePieceAtEdge(splitPieces, pieceIndex, edgeIndex))
        return false;

    g_cablePieces = std::move(splitPieces);
    g_useSavedCableSnapshot = false;
    g_snapshotCycleIndex = -1;
    g_mouseDrag = MouseDragState();
    g_cachedCableKymoKnotResults.clear();
    g_hasCableKymoKnotResults = false;

    std::cout << "Cut cable piece " << pieceIndex
              << " at edge " << edgeIndex << std::endl;
    Simulation::getCurrent()->getModel()->cleanup();
    buildModel();
    return true;
}

static bool handleMousePress(int button, int action, int mods)
{
    (void)mods;
    if (button != 0)
        return false;

    if (action == 1)
    {
        double x = 0.0;
        double y = 0.0;
        glfwGetCursorPos(MiniGL::getWindow(), &x, &y);
        return beginMouseDragAtScreenPosition(::Vector2i(static_cast<int>(x), static_cast<int>(y)));
    }

    if (action == 0)
    {
        endMouseDrag();
        return true;
    }

    return false;
}

static bool handleMouseMove(int x, int y)
{
    if (!g_mouseDrag.active)
        return false;

    if (glfwGetMouseButton(MiniGL::getWindow(), GLFW_MOUSE_BUTTON_1) != GLFW_PRESS)
    {
        endMouseDrag();
        return false;
    }

    updateMouseDragTargetAtScreenPosition(::Vector2i(x, y));
    return true;
}

static void buildModel()
{
    TimeManager::getCurrent()->setTimeStepSize(g_cfg.physics.timeStep);

    SimulationModel* model = Simulation::getCurrent()->getModel();
    g_currentSnapshotPath.clear();
    if (g_cablePieces.empty())
        initializeCablePieces();

    g_hasKymoKnotResult = false;
    g_lastKymoKnotResult = PBD::KymoKnotAnalysisResult();

    rebuildCableFromPieces(model, g_cablePieces);
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

static void addCableModel(
    SimulationModel* model,
    std::vector<Vector3r>& ropePoints,
    const RopeConfig& ropeCfg)
{
    addCableModel(model, ropePoints, ropeCfg, nullptr);
}

static void addCableModel(
    SimulationModel* model,
    std::vector<Vector3r>& ropePoints,
    const RopeConfig& ropeCfg,
    std::pair<unsigned int, unsigned int>* particleRangeOut)
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
    if (particleRangeOut != nullptr)
        *particleRangeOut = std::make_pair(firstParticle, lastParticle);

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
        applyMouseDragForces(model, dt);
        if (g_pullEnabled)
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

static bool isLocalIndexHighlighted(size_t pieceIndex, int localIndex)
{
    if (!g_hasCableKymoKnotResults ||
        pieceIndex >= g_cachedCableKymoKnotResults.size())
        return false;

    const PBD::KymoKnotAnalysisResult& knotResult =
        g_cachedCableKymoKnotResults[pieceIndex];
    if (!knotResult.available)
        return false;

    for (const PBD::KymoKnotSegment& segment : knotResult.segments)
    {
        if (isLocalIndexInKymoKnotSegment(localIndex, segment))
            return true;
    }
    return false;
}

static bool isLocalEdgeHighlighted(size_t pieceIndex, int localEdgeIndex)
{
    return isLocalIndexHighlighted(pieceIndex, localEdgeIndex) ||
           isLocalIndexHighlighted(pieceIndex, localEdgeIndex + 1);
}

static void cacheKymoKnotResult(const PBD::KymoKnotAnalysisResult& result)
{
    g_lastKymoKnotResult = result;
    g_hasKymoKnotResult = result.available && !result.segments.empty();
}

static void cacheCableKymoKnotResults(
    const std::vector<PBD::KymoKnotAnalysisResult>& results)
{
    g_cachedCableKymoKnotResults = results;
    g_hasCableKymoKnotResults = false;

    for (const PBD::KymoKnotAnalysisResult& result : g_cachedCableKymoKnotResults)
    {
        if (result.available && !result.segments.empty())
        {
            g_hasCableKymoKnotResults = true;
            break;
        }
    }
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

    std::cout << "Tanglement metrics:"
              << " self-collisions=" << contacts.size()
              << std::endl;

    const std::string bridgeRoot =
        (base != nullptr)
            ? Utilities::FileSystem::normalizePath(base->getOutputPath() + "/kymoknot_bridge")
            : "kymoknot_bridge";
    const std::vector<std::vector<Vector3r>> ropePieces =
        collectCurrentCablePieces(model);

    if (ropePieces.empty())
    {
        g_cachedCableKymoKnotResults.clear();
        g_hasCableKymoKnotResults = false;
        g_hasKymoKnotResult = false;
        std::cout << "KymoKnot: no cable pieces" << std::endl;
        return;
    }

    std::vector<PBD::KymoKnotAnalysisResult> cableResults;
    cableResults.reserve(ropePieces.size());
    for (size_t pieceIndex = 0; pieceIndex < ropePieces.size(); ++pieceIndex)
    {
        const std::vector<Vector3r>& cablePoints = ropePieces[pieceIndex];
        if (cablePoints.size() < 2)
        {
            cableResults.push_back(PBD::KymoKnotAnalysisResult());
            continue;
        }

        const PBD::KymoKnotAnalysisResult knotResult =
            PBD::KymoKnotAnalyzer::analyzeLinear(
                cablePoints,
                bridgeRoot,
                g_kymoKnotSearchMode);
        cableResults.push_back(knotResult);

        std::cout << "Cable piece " << pieceIndex << ": ";
        if (!knotResult.available)
        {
            std::cout << knotResult.message << std::endl;
            continue;
        }

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

    cacheCableKymoKnotResults(cableResults);
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

static void render()
{
    // `base->render()` is the built-in demo GUI + scene helper rendering.
    base->render();

    SimulationModel* model = Simulation::getCurrent()->getModel();
    const ParticleData& particles = model->getParticles();

    float segmentColor[4] = {0.0f, 0.0f, 0.8f, 1.0f};
    float knotColor[4] = {0.98f, 0.25f, 0.18f, 1.0f};
    float knotSphereColor[4] = {1.0f, 0.45f, 0.20f, 1.0f};
    float selectedColor[4] = {0.98f, 0.82f, 0.10f, 1.0f};
    float selectedSphereColor[4] = {1.0f, 0.92f, 0.30f, 1.0f};

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

            const int localEdgeIndex = static_cast<int>(e);
            const bool highlightKnot = isLocalEdgeHighlighted(m, localEdgeIndex);
            const bool highlightSelected =
                g_mouseDrag.selected &&
                m == g_mouseDrag.pieceIndex &&
                e == g_mouseDrag.edgeIndex;
            const float* activeColor =
                highlightSelected ? selectedColor :
                (highlightKnot ? knotColor : segmentColor);

            MiniGL::drawSphere(
                p0,
                0.03f,
                highlightSelected ? selectedSphereColor :
                (highlightKnot ? knotSphereColor : segmentColor));
            if (e == lineModel->getEdges().size() - 1)
                MiniGL::drawSphere(
                    p1,
                    0.03f,
                    isLocalIndexHighlighted(
                        m,
                        localEdgeIndex + 1)
                        ? knotSphereColor
                        : segmentColor);

            MiniGL::drawCylinder(
                p0,
                p1,
                activeColor,
                static_cast<float>(g_cfg.rope.renderRadius));
        }
    }
}
