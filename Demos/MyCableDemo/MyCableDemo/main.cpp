#include "Common/Common.h"
#include "Demos/Visualization/MiniGL.h"
#include "Demos/Visualization/Selection.h"
#include "Simulation/TimeManager.h"
#include <Eigen/Dense>
#include "Simulation/SimulationModel.h"
#include "Simulation/TimeStepController.h"
#include <iostream>
#include "Demos/Visualization/Visualization.h"
#include "Utils/Logger.h"
#include "Utils/Timing.h"
#include "Utils/FileSystem.h"
#include "Demos/Common/DemoBase.h"
#include "Simulation/Simulation.h"
#include "../Common/imguiParameters.h"

#include <vector>


// Enable memory leak detection
#if defined(_DEBUG) && !defined(EIGEN_ALIGN)
#define new DEBUG_NEW
#endif

using namespace PBD;
using namespace Eigen;
using namespace std;
using namespace Utilities;

static void timeStep();
static void buildModel();
static void createStraightCable(const Vector3r& position, const Matrix3r& orientation, Real length, int nPoints);
static void render();
static void reset();

// Lightweight ground contact
static void applyGroundConstraint(Real groundY);

// ----------------------------------------------------------------------------
// Globals / settings
// ----------------------------------------------------------------------------
static DemoBase* base = nullptr;

// cable discretization
static const unsigned int nParticles = 60;
static Real g_length = static_cast<Real>(6.0);

// visualization
static bool drawFrames = false;
static bool topDownCamera = true;
static Real g_renderRadius = static_cast<Real>(0.2);

// ground
static Real g_groundY = static_cast<Real>(0.0);
static bool g_enableGround = true;

// interactive end (kinematic target)
static unsigned int g_endParticle = 0;
static Vector3r g_endTarget(0, 0, 0);
static bool g_enableKinematicEnd = true;

// ----------------------------------------------------------------------------
// main
// ----------------------------------------------------------------------------
int main(int argc, char** argv)
{
    std::cout << "This demo shows a simple Cosserat rod (elastic cable) implementation using position-based dynamics." << std::endl;

        base = new DemoBase();
    base->init(argc, argv, "MyCableDemo (Cosserat rod)");

    SimulationModel* model = new SimulationModel();
    model->init();
    Simulation::getCurrent()->setModel(model);

    buildModel();

    base->createParameterGUI();

    // add additional demo parameters
    {
        // Draw frames
        auto* p = new imguiParameters::imguiBoolParameter();
        p->description = "Draw frames";
        p->label = "Draw frames";
        p->readOnly = false;
        p->getFct = []() -> bool { return drawFrames; };
        p->setFct = [](bool b) -> void { drawFrames = b; };
        imguiParameters::addParam("Visualization", "Elastic rods", p);
    }
    {
        // Top-down camera
        auto* p = new imguiParameters::imguiBoolParameter();
        p->description = "Lock camera to top-down view";
        p->label = "Top-down camera";
        p->readOnly = false;
        p->getFct = []() -> bool { return topDownCamera; };
        p->setFct = [](bool b) -> void { topDownCamera = b; };
        imguiParameters::addParam("Visualization", "Elastic rods", p);
    }
    {
        // Ground
        auto* p = new imguiParameters::imguiBoolParameter();
        p->description = "Enable simple ground constraint (y >= groundY)";
        p->label = "Enable ground";
        p->readOnly = false;
        p->getFct = []() -> bool { return g_enableGround; };
        p->setFct = [](bool b) -> void { g_enableGround = b; };
        imguiParameters::addParam("Simulation", "Elastic rods", p);
    }
    {
        // Kinematic end
        auto* p = new imguiParameters::imguiBoolParameter();
        p->description = "Make cable end kinematic (controlled by target)";
        p->label = "Kinematic end";
        p->readOnly = false;
        p->getFct = []() -> bool { return g_enableKinematicEnd; };
        p->setFct = [](bool b) -> void { g_enableKinematicEnd = b; };
        imguiParameters::addParam("Simulation", "Elastic rods", p);
    }

    // OpenGL callbacks
    MiniGL::setClientIdleFunc(timeStep);
    MiniGL::addKeyFunc('r', reset);

    // Fail-safe end control: IJKL + U/O (x,z,y)
    MiniGL::addKeyFunc('i', []() { g_endTarget.z() += 0.25; });
    MiniGL::addKeyFunc('k', []() { g_endTarget.z() -= 0.25; });
    MiniGL::addKeyFunc('j', []() { g_endTarget.x() -= 0.25; });
    MiniGL::addKeyFunc('l', []() { g_endTarget.x() += 0.25; });
    MiniGL::addKeyFunc('u', []() { g_endTarget.y() += 0.25; });
    MiniGL::addKeyFunc('o', []() { g_endTarget.y() -= 0.25; });

    MiniGL::setClientSceneFunc(render);

    // initial camera
    if (topDownCamera)
    {
        // Camera above +Y, looking down.
        MiniGL::setViewport(40.0f, 0.1f, 500.0f,
            Vector3r(0, 50, 0),
            Vector3r(0.0, 0, 0.01));
    }
    else
    {
        MiniGL::setViewport(40.0f, 0.1f, 500.0f,
            Vector3r(5.0, 10.0, 30.0),
            Vector3r(5.0, 0.0, 0.0));
    }

    // NOTE:
    // Most PBD demos already support mouse picking/dragging via DemoBase/Selection.
    // Try: left click a particle and drag (mouse spring). If it doesn't work in your build,
    // tell me and I'll give you a tiny patch that wires the mouse callbacks explicitly.

    MiniGL::mainLoop();

    Utilities::Timing::printAverageTimes();
    Utilities::Timing::printTimeSums();

    delete Simulation::getCurrent();
    delete base;
    delete model;

    return 0;
}

// ----------------------------------------------------------------------------
// Reset / timestep
// ----------------------------------------------------------------------------
static void reset()
{
    Utilities::Timing::printAverageTimes();
    Utilities::Timing::reset();

    Simulation::getCurrent()->reset();
    base->getSelectedParticles().clear();

    Simulation::getCurrent()->getModel()->cleanup();
    buildModel();
}

static void timeStep()
{
    const Real pauseAt = base->getValue<Real>(DemoBase::PAUSE_AT);
    if ((pauseAt > 0.0) && (pauseAt < TimeManager::getCurrent()->getTime()))
        base->setValue(DemoBase::PAUSE, true);

    if (base->getValue<bool>(DemoBase::PAUSE))
        return;

    SimulationModel* model = Simulation::getCurrent()->getModel();
    const unsigned int numSteps = base->getValue<unsigned int>(DemoBase::NUM_STEPS_PER_RENDER);

    for (unsigned int i = 0; i < numSteps; i++)
    {
        // apply kinematic end target before each substep
        if (g_enableKinematicEnd && g_endParticle < model->getParticles().getNumberOfParticles())
        {
            ParticleData& pd = model->getParticles();
            // setPosition exists in this framework; if your build complains, switch to pd.getPosition(...) = ...
            pd.setPosition(g_endParticle, g_endTarget);
        }

        START_TIMING("SimStep");
        Simulation::getCurrent()->getTimeStep()->step(*model);
        STOP_TIMING_AVG;

        if (g_enableGround)
            applyGroundConstraint(g_groundY);

        base->step();
    }
}

// ----------------------------------------------------------------------------
// Scene / model
// ----------------------------------------------------------------------------
static void buildModel()
{
    TimeManager::getCurrent()->setTimeStepSize(static_cast<Real>(0.005));

    // start above the ground
    createStraightCable(Vector3r(0, 6, 0), Matrix3r::Identity(), g_length, (int)nParticles);

    // place end target at initial end position
    SimulationModel* model = Simulation::getCurrent()->getModel();
    ParticleData& pd = model->getParticles();
    if (pd.getNumberOfParticles() > 0 && g_endParticle < pd.getNumberOfParticles())
        g_endTarget = pd.getPosition(g_endParticle);


    SimulationModel* m = Simulation::getCurrent()->getModel();
    std::cout << "[buildModel] lineModels=" << m->getLineModels().size()
        << " particles=" << m->getParticles().getNumberOfParticles()
        << " quats=" << m->getOrientations().getNumberOfQuaternions()
        << std::endl;
}

static void createStraightCable(const Vector3r& position, const Matrix3r& orientation, Real length, int nPoints)
{
    const int nQuaternions = nPoints - 1;
    vector<Vector3r> points(nPoints);
    vector<Quaternionr> quaternions(nQuaternions);

    const Real dx = length / static_cast<Real>(nPoints - 1);

    // init particles: straight line along local +z
    for (int i = 0; i < nPoints; i++)
    {
        Vector3r p(0.0, 0.0, dx * static_cast<Real>(i));
        points[i] = orientation * p + position;
    }

    // init quaternions: align frames to segment tangents
    Vector3r from(0, 0, 1);
    for (int i = 0; i < nQuaternions; i++)
    {
        Vector3r to = (points[i + 1] - points[i]).normalized();
        Quaternionr dq = Quaternionr::FromTwoVectors(from, to);
        if (i == 0) quaternions[i] = dq;
        else        quaternions[i] = dq * quaternions[i - 1];
        from = to;
    }

    vector<unsigned int> indices(2 * nPoints - 1);
    vector<unsigned int> indicesQuaternions(nQuaternions);

    for (int i = 0; i < nPoints - 1; i++)
    {
        indices[2 * i] = (unsigned int)i;
        indices[2 * i + 1] = (unsigned int)(i + 1);
    }
    for (int i = 0; i < nQuaternions; i++)
        indicesQuaternions[i] = (unsigned int)i;

    SimulationModel* model = Simulation::getCurrent()->getModel();
    model->addLineModel(nPoints, nQuaternions, &points[0], &quaternions[0], &indices[0], &indicesQuaternions[0]);

    // mass setup
    ParticleData& pd = model->getParticles();
    const int nPointsTotal = pd.getNumberOfParticles();

    // last inserted rod points are [nPointsTotal-nPoints, ..., nPointsTotal-1]
    const unsigned int first = (unsigned int)(nPointsTotal - nPoints);
    const unsigned int last = (unsigned int)(nPointsTotal - 1);

    for (unsigned int i = last; i > first; i--)
        pd.setMass(i, 1.0);

    // pin first point
    pd.setMass(first, 0.0);

    // make last point kinematic (draggable target)
    g_endParticle = last;
    pd.setMass(g_endParticle, 0.0);

    // orientation masses
    OrientationData& od = model->getOrientations();
    const unsigned int nQuaternionsTotal = od.getNumberOfQuaternions();
    const unsigned int firstQ = (unsigned int)(nQuaternionsTotal - nQuaternions);
    const unsigned int lastQ = (unsigned int)(nQuaternionsTotal - 1);

    for (unsigned int i = lastQ; i > firstQ; i--)
        od.setMass(i, 1.0);

    // pin first quaternion for stability
    od.setMass(firstQ, 0.0);

    // init constraints (Cosserat rod)
    const size_t rodNumber = model->getLineModels().size() - 1;
    const unsigned int offset = model->getLineModels()[rodNumber]->getIndexOffset();
    const unsigned int offsetQuaternions = model->getLineModels()[rodNumber]->getIndexOffsetQuaternions();
    const size_t nEdges = model->getLineModels()[rodNumber]->getEdges().size();
    const LineModel::Edges& edges = model->getLineModels()[rodNumber]->getEdges();

    // stretch/shear constraints
    for (unsigned int i = 0; i < nEdges; i++)
    {
        const unsigned int v1 = edges[i].m_vert[0] + offset;
        const unsigned int v2 = edges[i].m_vert[1] + offset;
        const unsigned int q1 = edges[i].m_quat + offsetQuaternions;
        model->addStretchShearConstraint(
            v1, v2, q1,
            model->getRodStretchingStiffness(),
            model->getRodShearingStiffnessX(),
            model->getRodShearingStiffnessY()
        );
    }

    // bend/twist constraints
    for (unsigned int i = 0; i < nEdges - 1; i++)
    {
        const unsigned int q1 = edges[i].m_quat + offsetQuaternions;
        const unsigned int q2 = edges[i + 1].m_quat + offsetQuaternions;
        model->addBendTwistConstraint(
            q1, q2,
            model->getRodTwistingStiffness(),
            model->getRodBendingStiffnessX(),
            model->getRodBendingStiffnessY()
        );
    }
}

// ----------------------------------------------------------------------------
// Rendering
// ----------------------------------------------------------------------------
static void renderLineModels()
{
    SimulationModel* model = Simulation::getCurrent()->getModel();
    const ParticleData& pd = model->getParticles();
    const OrientationData& od = model->getOrientations();

    float red[4] = { 0.8f, 0.0f, 0.0f, 1 };
    float green[4] = { 0.0f, 0.8f, 0.0f, 1 };
    float blue[4] = { 0.0f, 0.0f, 0.8f, 1 };

    for (unsigned int i = 0; i < model->getLineModels().size(); i++)
    {
        LineModel* lineModel = model->getLineModels()[i];

        for (unsigned int e = 0; e < lineModel->getEdges().size(); e++)
        {
            const unsigned int indexOffset = lineModel->getIndexOffset();
            const unsigned int indexOffsetQuaternions = lineModel->getIndexOffsetQuaternions();
            const unsigned int i1 = lineModel->getEdges()[e].m_vert[0] + indexOffset;
            const unsigned int i2 = lineModel->getEdges()[e].m_vert[1] + indexOffset;
            const unsigned int iq = lineModel->getEdges()[e].m_quat + indexOffsetQuaternions;

            const Vector3r& v1 = pd.getPosition(i1);
            const Vector3r& v2 = pd.getPosition(i2);
            const Quaternionr& q = od.getQuaternion(iq);

            // lightweight: tiny spheres + thin cylinder
            MiniGL::drawSphere(v1, 0.07f, blue);
            if (e == lineModel->getEdges().size() - 1)
                MiniGL::drawSphere(v2, 0.07f, blue);

            MiniGL::drawCylinder(v1, v2, blue, (float)g_renderRadius);

            if (drawFrames)
            {
                Vector3r vm = 0.5 * (v1 + v2);
                Real scale = static_cast<Real>(0.15);
                Vector3r d1 = q._transformVector(Vector3r(1, 0, 0)) * scale;
                Vector3r d2 = q._transformVector(Vector3r(0, 1, 0)) * scale;
                Vector3r d3 = q._transformVector(Vector3r(0, 0, 1)) * scale;
                MiniGL::drawCylinder(vm, vm + d1, red, 0.01f);
                MiniGL::drawCylinder(vm, vm + d2, green, 0.01f);
                MiniGL::drawCylinder(vm, vm + d3, blue, 0.01f);
            }
        }
    }

	float yellow[4] = { 0.9f, 0.9f, 0.0f, 1 };
	MiniGL::drawSphere(g_endTarget, 0.15f, yellow);
}

static void render()
{
    //std::cout << "drawing sphere" << endl;
	float white[4] = { 0.8f, 0.8f, 0.8f, 1 };
	MiniGL::drawSphere(Vector3r(0, 0, 0), 0.1f, white);

    base->render();
    renderLineModels();
}

// ----------------------------------------------------------------------------
// Helpers
// ----------------------------------------------------------------------------
static void applyGroundConstraint(Real groundY)
{
    SimulationModel* model = Simulation::getCurrent()->getModel();
    ParticleData& pd = model->getParticles();

    for (unsigned int i = 0; i < pd.getNumberOfParticles(); i++)
    {
        Vector3r p = pd.getPosition(i);
        if (p.y() < groundY)
        {
            p.y() = groundY;
            pd.setPosition(i, p);
        }
    }
}
