#ifndef MY_CABLE_DEMO_CABLE_PHYSICS_H
#define MY_CABLE_DEMO_CABLE_PHYSICS_H

#include "Common/Common.h"

#include <string>
#include <utility>
#include <vector>

namespace PBD
{
    class SimulationModel;
}

struct RopeConfig
{
    Real stretchStiffness = static_cast<Real>(1.0);
    Real shearStiffnessX = static_cast<Real>(1.0);
    Real shearStiffnessY = static_cast<Real>(1.0);
    Real bendStiffnessX = static_cast<Real>(0.003);
    Real bendStiffnessY = static_cast<Real>(0.003);
    Real twistStiffness = static_cast<Real>(0.012);
    Real dx = static_cast<Real>(0.1);
    Real renderRadius = static_cast<Real>(0.06);
};

struct PhysicsConfig
{
    Real timeStep = static_cast<Real>(0.005);
    Real pullSpring = static_cast<Real>(1000.0);
    Real pullDamping = static_cast<Real>(2.0);
    Real pullMaxAccel = static_cast<Real>(1000.0);
    Real velocityRecomputeBlend = static_cast<Real>(0.03);
};

struct RestShapeCreepConfig
{
    bool enabled = true;
    unsigned int stepInterval = 20;
    Real creep = static_cast<Real>(0.02);
    Real activationThreshold = static_cast<Real>(1e-4);
};

struct CableLayoutConfig
{
    Real endpointX = static_cast<Real>(10.0);
    Real height = static_cast<Real>(6.0);
};

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
    int smoothingIterations = 3;
    Real smoothingAlpha = static_cast<Real>(0.12);
    Real minHeightAboveGround = static_cast<Real>(5.0);
};

struct SelfCollisionConfig
{
    bool enabled = true;
    Real ropeRadius = static_cast<Real>(0.08);
    Real compliance = static_cast<Real>(1e-10);
    unsigned int iterations = 8;
    Real excludeArcLength = static_cast<Real>(0.9);
    Real normalDamping = static_cast<Real>(0.15);
    Real tangentialDamping = static_cast<Real>(0.05);
    Real maxCorrection = static_cast<Real>(0.08);
    Real hashCellSize = static_cast<Real>(0.2);
};

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

struct SelfCollisionContact
{
    unsigned int edgeA;
    unsigned int edgeB;
};

extern unsigned int g_restCreepCounter;
extern bool g_pullEnabled;
extern std::vector<std::pair<unsigned int, unsigned int>> g_cablePieceRanges;

Real clamp01(Real x);
void closestPtSegmentSegment(
    const Vector3r& p1, const Vector3r& q1,
    const Vector3r& p2, const Vector3r& q2,
    Real& s, Real& t,
    Vector3r& c1, Vector3r& c2);
void recomputeVelocitiesFromPositions(
    PBD::SimulationModel* model,
    const std::vector<Vector3r>& oldPositions,
    Real dt,
    Real blend);
void applyGroundConstraint(Real groundY);
void applyPullForces(PBD::SimulationModel* model, Real dt, const PhysicsConfig& physicsCfg);
void applySelfCollisionXPBD(
    PBD::SimulationModel* model,
    Real dt,
    const RopeConfig& ropeCfg,
    const SelfCollisionConfig& collisionCfg);
std::vector<SelfCollisionContact> detectSelfCollisionContacts(
    PBD::SimulationModel* model,
    const RopeConfig& ropeCfg,
    const SelfCollisionConfig& collisionCfg);
void applyBendTwistRestStateCreep(
    PBD::SimulationModel* model,
    const RestShapeCreepConfig& creepCfg);

#endif
