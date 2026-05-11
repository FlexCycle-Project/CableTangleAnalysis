#include "Common/Common.h"
#include "Demos/Visualization/MiniGL.h"
#include "Simulation/TimeManager.h"
#include "Simulation/SimulationModel.h"
#include "Simulation/Constraints.h"
#include "Demos/Common/DemoBase.h"
#include "Simulation/Simulation.h"
#include "CablePhysics.h"
#include "CableGeneration.h"
#include "CableVisualization.h"
#include "TanglementAnalyzer.h"
#include "KymoKnotAnalyzer.h"

#include <GLFW/glfw3.h>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <cstdint>
#include <limits>
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
static void regenerateMathematicalKnot();
static void toggleMathematicalKnotPreset();
static void togglePullForces();
static void adjustPullStrength(Real scale);
static void timeStep();
static void render();

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
static bool beginMouseDragAtScreenPosition(const ::Vector2i& screenPos);
static void updateMouseDragTargetAtScreenPosition(const ::Vector2i& screenPos);
static bool cutSelectedCableSegment();
static void endMouseDrag();
static bool handleMousePress(int button, int action, int mods);
static bool handleMouseMove(int x, int y);

static DemoBase* base = nullptr;
static AppConfig g_cfg;
static unsigned int g_firstParticle = 0;
static unsigned int g_lastParticle = 0;
static unsigned int g_restCreepCounter = 0;
static std::vector<PBD::KymoKnotAnalysisResult> g_cachedCableKymoKnotResults;
static bool g_hasCableKymoKnotResults = false;
static CableVisualization::VisualizationState g_visualizationState;
static PBD::KymoKnotSearchMode g_kymoKnotSearchMode = PBD::KymoKnotSearchMode::TopDown;
static std::vector<std::vector<Vector3r>> g_cablePieces;
static std::vector<std::pair<unsigned int, unsigned int>> g_cablePieceRanges;
static bool g_useMathematicalKnot = false;
static std::string g_mathematicalKnotName = "5_2";
bool g_pullEnabled = false;
static bool g_gravityEnabled = true;
static Vector3r g_gravityReference = Vector3r(0.0, -9.81, 0.0);
static void captureGravityReference();
static void applyGravityState();
static void toggleGravity();

static MouseDragState g_mouseDrag;
static CableVisualization::HoopMarkerSelection g_hoopMarkerSelection;
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
    captureGravityReference();
    applyGravityState();

    // Keyboard input for this demo lives here.
    MiniGL::setClientIdleFunc(timeStep);
    MiniGL::setClientSceneFunc(render);
    MiniGL::addMousePressFunc(handleMousePress);
    MiniGL::addMouseMoveFunc(handleMouseMove);
    MiniGL::addKeyFunc('r', reset);
    MiniGL::addKeyFunc('t', []() {
        SimulationModel* model = Simulation::getCurrent()->getModel();
        if (model == nullptr)
        {
            std::cout << "No simulation model available." << std::endl;
            return;
        }

        const std::vector<SelfCollisionContact> contacts =
            detectSelfCollisionContacts(
                model,
                g_cablePieceRanges,
                g_cfg.rope,
                g_cfg.selfCollision);
        const std::vector<std::vector<Vector3r>> ropePieces =
            collectCurrentCablePieces(model);
        const TanglementAnalyzer::TanglementAnalysisResult tanglement =
            TanglementAnalyzer::analyzeTanglement(
                ropePieces,
                contacts,
                g_kymoKnotSearchMode);

        g_cachedCableKymoKnotResults = tanglement.kymoResults;
        g_hasCableKymoKnotResults = false;
        for (const PBD::KymoKnotAnalysisResult& knotResult : g_cachedCableKymoKnotResults)
        {
            if (knotResult.available && knotResult.success && !knotResult.segments.empty())
            {
                g_hasCableKymoKnotResults = true;
                break;
            }
        }
        g_visualizationState = CableVisualization::buildVisualizationState(tanglement);

        TanglementAnalyzer::printTanglementReport(ropePieces, tanglement);
    });
    MiniGL::addKeyFunc('c', []() {
        if (!cutSelectedCableSegment())
            std::cout << "No selected cable segment to cut" << std::endl;
    });
    MiniGL::addKeyFunc('v', []() {
        if (g_mouseDrag.active && g_mouseDrag.selected)
        {
            g_mouseDrag.upwardPull = true;
            std::cout << "Upward drag enabled for selected cable segment" << std::endl;
        }
    });
    MiniGL::addKeyFunc('k', []() {
        toggleMathematicalKnotPreset();
    });
    MiniGL::addKeyFunc('o', []() {
        g_cablePieces.clear();
        regenerateRandomCable();
    });
    MiniGL::addKeyFunc('*', []() {
        togglePullForces();
    });
    MiniGL::addKeyFunc('g', []() {
        toggleGravity();
    });
    MiniGL::addKeyFunc('-', []() {
        adjustPullStrength(static_cast<Real>(0.85));
    });
    MiniGL::addKeyFunc('+', []() {
        adjustPullStrength(static_cast<Real>(1.15));
    });
    MiniGL::addKeyFunc('z', []() { MiniGL::move(0, 0, -1); });
    MiniGL::addKeyFunc('x', []() { MiniGL::move(0, 0, 2); });

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
    g_hoopMarkerSelection = CableVisualization::HoopMarkerSelection();
    g_cachedCableKymoKnotResults.clear();
    g_hasCableKymoKnotResults = false;
    g_visualizationState.clear();
    g_cablePieces.clear();
    g_restCreepCounter = 0;
    Simulation::getCurrent()->reset();
    base->getSelectedParticles().clear();

    Simulation::getCurrent()->getModel()->cleanup();
    buildModel();
    applyGravityState();
}

static void regenerateRandomCable()
{
    g_mouseDrag = MouseDragState();
    g_hoopMarkerSelection = CableVisualization::HoopMarkerSelection();
    g_useMathematicalKnot = false;
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

static void regenerateMathematicalKnot()
{
    g_mouseDrag = MouseDragState();
    g_hoopMarkerSelection = CableVisualization::HoopMarkerSelection();
    g_useMathematicalKnot = true;
    g_cfg.randomCable.enabled = false;
    g_cablePieces.clear();

    std::cout << "Regenerating mathematical knot: " << g_mathematicalKnotName << std::endl;
    reset();
}

static void toggleMathematicalKnotPreset()
{
    if (!g_useMathematicalKnot)
    {
        g_mathematicalKnotName = "5_2";
        regenerateMathematicalKnot();
        return;
    }

    g_useMathematicalKnot = false;
    g_mouseDrag = MouseDragState();
    g_hoopMarkerSelection = CableVisualization::HoopMarkerSelection();
    g_cfg.randomCable.enabled = false;
    g_cablePieces.clear();
    std::cout << "Regenerating straight cable" << std::endl;
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

static void captureGravityReference()
{
    if (!Simulation::hasCurrent())
        return;

    Real* grav = Simulation::getCurrent()->getVecValue<Real>(Simulation::GRAVITATION);
    if (grav == nullptr)
        return;

    g_gravityReference = Vector3r(grav[0], grav[1], grav[2]);
}

static void applyGravityState()
{
    if (!Simulation::hasCurrent())
        return;

    Vector3r gravity = g_gravityEnabled
        ? g_gravityReference
        : Vector3r::Zero();
    Simulation::getCurrent()->setVecValue<Real>(Simulation::GRAVITATION, gravity.data());
}

static void toggleGravity()
{
    g_gravityEnabled = !g_gravityEnabled;
    applyGravityState();
    std::cout << "Gravity: " << (g_gravityEnabled ? "on" : "off") << std::endl;
}

static void initializeCablePieces()
{
    g_cablePieces.clear();
    std::vector<Vector3r> ropePoints;
    if (g_useMathematicalKnot)
    {
        ropePoints = CableGeneration::generateMathematicalKnotPolyline(g_mathematicalKnotName);
    }
    else if (g_cfg.randomCable.enabled)
    {
        ropePoints = CableGeneration::generateRandomCableMessPolyline(
            g_cfg.layout,
            g_cfg.rope,
            g_cfg.randomCable,
            g_cfg.groundY);
    }
    else
    {
        ropePoints = CableGeneration::generateStraightCablePolyline(g_cfg.layout, g_cfg.rope);
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

static bool beginMouseDragAtScreenPosition(const ::Vector2i& screenPos)
{
    const std::vector<std::vector<Vector3r>> ropePieces =
        collectCurrentCablePieces(Simulation::getCurrent()->getModel());
    CableVisualization::HoopMarkerSelection hoopPick;
        if (CableVisualization::selectHoopMarkerAtScreenPosition(
            screenPos,
            ropePieces,
            g_cfg.rope.renderRadius,
            g_visualizationState,
            hoopPick))
    {
        g_hoopMarkerSelection = hoopPick;
        g_mouseDrag = MouseDragState();
        std::cout << "Selected hoop marker on cable piece " << hoopPick.pieceIndex;
        if (hoopPick.pieceIndex < g_visualizationState.scoredHoops.size() &&
            hoopPick.pieceIndex < ropePieces.size())
        {
            if (hoopPick.hoopIndex < g_visualizationState.scoredHoops[hoopPick.pieceIndex].size())
            {
                const TanglementAnalyzer::HoopCandidate& hoop =
                    g_visualizationState.scoredHoops[hoopPick.pieceIndex][hoopPick.hoopIndex];
                if (hoop.pairSumScore != 0)
                {
                    std::cout << " hoop " << hoopPick.hoopIndex
                              << " score=" << hoop.pairSumScore;
                }
            }
        }
        std::cout << std::endl;
        return true;
    }

    CableVisualization::CablePickResult pickResult;
    if (!CableVisualization::selectCableEdgeAtScreenPosition(
            screenPos,
            ropePieces,
            g_cfg.rope.renderRadius,
            pickResult))
    {
        return false;
    }

    g_hoopMarkerSelection = CableVisualization::HoopMarkerSelection();
    g_mouseDrag.active = true;
    g_mouseDrag.selected = true;
    g_mouseDrag.pieceIndex = pickResult.pieceIndex;
    g_mouseDrag.edgeIndex = pickResult.edgeIndex;
    g_mouseDrag.edgeT = pickResult.edgeT;
    g_mouseDrag.rayDepthT = pickResult.rayDepthT;
    g_mouseDrag.grabbedPoint = pickResult.grabbedPoint;
    g_mouseDrag.targetPoint = pickResult.grabbedPoint;
    g_mouseDrag.upwardPull = false;

    std::cout << "Selected cable piece " << pickResult.pieceIndex
              << " segment " << pickResult.edgeIndex << std::endl;
    return true;
}

static void updateMouseDragTargetAtScreenPosition(const ::Vector2i& screenPos)
{
    if (!g_mouseDrag.active)
        return;

    Vector3r rayStart(0.0, 0.0, 0.0);
    Vector3r rayEnd(0.0, 0.0, 0.0);
    if (!CableVisualization::buildPickRay(screenPos, rayStart, rayEnd))
        return;

    g_mouseDrag.targetPoint = rayStart + g_mouseDrag.rayDepthT * (rayEnd - rayStart);
}

static void endMouseDrag()
{
    g_mouseDrag.active = false;
    g_mouseDrag.upwardPull = false;
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
    g_mouseDrag = MouseDragState();
    g_hoopMarkerSelection = CableVisualization::HoopMarkerSelection();
    g_cachedCableKymoKnotResults.clear();
    g_hasCableKymoKnotResults = false;
    g_visualizationState.clear();

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
    if (g_cablePieces.empty())
        initializeCablePieces();

    rebuildCableFromPieces(model, g_cablePieces);
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
        applyMouseDragForces(
            model,
            dt,
            g_mouseDrag,
            g_cablePieceRanges,
            g_cfg.rope,
            g_mouseDragSpring,
            g_mouseDragDamping,
            g_mouseDragMaxAccel);
        if (g_pullEnabled)
            applyPullForces(model, dt, g_cablePieceRanges, g_cfg.physics);

        Simulation::getCurrent()->getTimeStep()->step(*model);

        if (g_cfg.selfCollision.enabled)
        {
            std::vector<Vector3r> oldPositions(particles.getNumberOfParticles());
            for (unsigned int i = 0; i < particles.getNumberOfParticles(); ++i)
                oldPositions[i] = particles.getPosition(i);

            for (unsigned int k = 0; k < g_cfg.selfCollision.iterations; ++k)
                applySelfCollision(
                    model,
                    dt,
                    g_cablePieceRanges,
                    g_cfg.rope,
                    g_cfg.selfCollision);

            applyGroundConstraintToModel(model, g_cfg.groundY);

            recomputeVelocitiesFromPositions(
                model,
                oldPositions,
                dt,
                g_cfg.physics.velocityRecomputeBlend);
        }
        else
        {
            applyGroundConstraintToModel(model, g_cfg.groundY);
        }

        applyBendTwistRestStateCreep(model, g_cfg.restCreep, g_restCreepCounter);

        base->step();
    }
}

static void render()
{
    base->render();
    const CableVisualization::HighlightSelection selection
    {
        g_mouseDrag.selected,
        g_mouseDrag.pieceIndex,
        g_mouseDrag.edgeIndex
    };

    CableVisualization::renderCable(
        Simulation::getCurrent()->getModel(),
        g_cfg.rope,
        selection,
        g_hoopMarkerSelection,
        g_visualizationState);
}

