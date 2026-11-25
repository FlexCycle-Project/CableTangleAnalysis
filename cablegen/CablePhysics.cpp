#include "CablePhysics.h"

#include "Simulation/Constraints.h"
#include "Simulation/Simulation.h"
#include "Simulation/SimulationModel.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <unordered_set>

using namespace PBD;

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

struct VelocityPairCorrection
{
    // Contact damping is applied at segment closest points. These fields store
    // the four endpoint indices and barycentric weights needed to distribute
    // the relative velocity correction back to particles.
    unsigned int i0, i1, j0, j1;
    Real wP0, wP1, wQ0, wQ1;
    Vector3r nrm;
};

struct CableEdgeRef
{
    size_t pieceIndex;
    unsigned int localEdgeIndex;
    unsigned int startParticle;
};

static std::uint64_t makeEdgePairKey(unsigned int edgeA, unsigned int edgeB)
{
    const unsigned int a = std::min(edgeA, edgeB);
    const unsigned int b = std::max(edgeA, edgeB);
    return (static_cast<std::uint64_t>(a) << 32) | static_cast<std::uint64_t>(b);
}

Real clamp01(Real x)
{
    if (x < static_cast<Real>(0.0))
        return static_cast<Real>(0.0);
    if (x > static_cast<Real>(1.0))
        return static_cast<Real>(1.0);
    return x;
}

static Vector3r clampVec(const Vector3r& v, Real maxNorm)
{
    const Real len = v.norm();
    if (len > maxNorm && len > static_cast<Real>(1e-12))
        return v * (maxNorm / len);
    return v;
}

static Vector3r computeSpringDampedAcceleration(
    const Vector3r& springTerm,
    const Vector3r& currentVelocity,
    Real damping,
    Real maxAccel)
{
    Vector3r acceleration = springTerm - damping * currentVelocity;
    const Real accelNorm = acceleration.norm();
    if (accelNorm > maxAccel)
        acceleration *= (maxAccel / accelNorm);
    return acceleration;
}

static uint64_t makeOrderedEdgePairKey(unsigned int a, unsigned int b)
{
    const uint64_t lo = std::min<uint64_t>(a, b);
    const uint64_t hi = std::max<uint64_t>(a, b);
    return (lo << 32) | hi;
}

void closestPtSegmentSegment(
    const Vector3r& p1, const Vector3r& q1,
    const Vector3r& p2, const Vector3r& q2,
    Real& s, Real& t,
    Vector3r& c1, Vector3r& c2)
{
    // Standard closest-points query between two finite line segments. The
    // parameters s and t are barycentric coordinates on segment p1-q1 and
    // p2-q2, and are later reused to distribute contact corrections.
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

void recomputeVelocitiesFromPositions(
    SimulationModel* model,
    const std::vector<Vector3r>& oldPositions,
    Real dt,
    Real blend)
{
    // PBD modifies positions directly, so the final velocity must be rebuilt
    // from the actual corrected displacement. Blending keeps this correction
    // from injecting sharp velocity changes after contact projection.
    ParticleData& pd = model->getParticles();
    const unsigned int n = pd.getNumberOfParticles();
    if (oldPositions.size() != n || dt <= static_cast<Real>(1e-12))
        return;

    blend = std::max(static_cast<Real>(0.0), std::min(static_cast<Real>(1.0), blend));

    for (unsigned int i = 0; i < n; ++i)
    {
        const Vector3r vPos = (pd.getPosition(i) - oldPositions[i]) / dt;
        Vector3r& v = pd.getVelocity(i);
        v = (static_cast<Real>(1.0) - blend) * v + blend * vPos;
    }
}

void addCableModel(
    SimulationModel* model,
    std::vector<Vector3r>& ropePoints,
    const RopeConfig& ropeCfg)
{
    addCableModel(model, ropePoints, ropeCfg, nullptr);
}

void addCableModel(
    SimulationModel* model,
    std::vector<Vector3r>& ropePoints,
    const RopeConfig& ropeCfg,
    std::pair<unsigned int, unsigned int>* particleRangeOut)
{
    if (model == nullptr)
        return;

    const int particleCount = static_cast<int>(ropePoints.size());
    const int quaternionCount = particleCount - 1;

    std::vector<Quaternionr> quaternions(quaternionCount);
    Vector3r previousDirection(0, 0, 1);
    for (int i = 0; i < quaternionCount; ++i)
    {
        // Segment orientations store the local Cosserat frame. The cumulative
        // quaternion chain follows the cable tangent while avoiding Euler-angle
        // singularities.
        const Vector3r edge = ropePoints[i + 1] - ropePoints[i];
        Vector3r currentDirection;

        if (edge.squaredNorm() > static_cast<Real>(1e-12))
            currentDirection = edge.normalized();
        else
            currentDirection = previousDirection;

        const Quaternionr rotation =
            Quaternionr::FromTwoVectors(previousDirection, currentDirection);

        if (i == 0)
        {
            quaternions[i] = rotation;
        }
        else
        {
            quaternions[i] = rotation * quaternions[i - 1];
        }
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

    const size_t firstConstraintIndex = model->getConstraints().size();

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

    if (particleRangeOut != nullptr)
        *particleRangeOut = std::make_pair(firstParticle, lastParticle);

    LineModel* lineModel = model->getLineModels().back();
    const unsigned int pointOffset = lineModel->getIndexOffset();
    const unsigned int quatOffset = lineModel->getIndexOffsetQuaternions();
    const LineModel::Edges& edges = lineModel->getEdges();

    for (unsigned int i = 0; i < edges.size(); ++i)
    {
        // Stretch-shear constraints couple each segment length/direction with
        // its material frame, giving axial and shear resistance.
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
        // Bend-twist constraints compare neighboring material frames, which is
        // the discrete analogue of curvature and twist in a Cosserat rod.
        const unsigned int q0 = edges[i].m_quat + quatOffset;
        const unsigned int q1 = edges[i + 1].m_quat + quatOffset;

        model->addBendTwistConstraint(
            q0,
            q1,
            ropeCfg.twistStiffness,
            ropeCfg.bendStiffnessX,
            ropeCfg.bendStiffnessY);
    }

    setNewCableBendTwistRestStateToStraight(model, firstConstraintIndex);
}

void setNewCableBendTwistRestStateToStraight(
    SimulationModel* model,
    size_t firstConstraintIndex)
{
    if (model == nullptr)
        return;

    const SimulationModel::ConstraintVector& constraints = model->getConstraints();
    for (size_t i = firstConstraintIndex; i < constraints.size(); ++i)
    {
        // New randomly generated cables should not treat their initial random
        // shape as stress-free. Resetting bend-twist rest states to identity
        // makes the cable relax toward a straight natural shape.
        Constraint* constraint = constraints[i];
        if (constraint == nullptr ||
            constraint->getTypeId() != BendTwistConstraint::TYPE_ID)
        {
            continue;
        }

        BendTwistConstraint* bendTwist =
            static_cast<BendTwistConstraint*>(constraint);
        bendTwist->m_restDarbouxVector =
            Quaternionr(static_cast<Real>(1.0), static_cast<Real>(0.0), static_cast<Real>(0.0), static_cast<Real>(0.0));
    }
}

void applyMouseDragForces(
    SimulationModel* model,
    Real dt,
    const MouseDragState& dragState,
    const std::vector<std::pair<unsigned int, unsigned int>>& pieceRanges,
    const RopeConfig& ropeCfg,
    Real spring,
    Real damping,
    Real maxAccel)
{
    if (!dragState.active || model == nullptr)
        return;
    if (dragState.pieceIndex >= pieceRanges.size())
        return;

    const std::pair<unsigned int, unsigned int> range =
        pieceRanges[dragState.pieceIndex];

    if (range.second <= range.first)
        return;

    const unsigned int localEdgeIndex = static_cast<unsigned int>(dragState.edgeIndex);
    const unsigned int i0 = range.first + localEdgeIndex;
    const unsigned int i1 = i0 + 1;

    ParticleData& particles = model->getParticles();
    if (i1 >= particles.getNumberOfParticles())
        return;

    const Real w0 = static_cast<Real>(1.0) - dragState.edgeT;
    const Real w1 = dragState.edgeT;

    // Dragging acts at the selected point on an edge, not at a single particle.
    // The same barycentric weights distribute the spring-damper force to the
    // two endpoints, keeping the interaction continuous along the cable.
    const Vector3r p0 = particles.getPosition(i0);
    const Vector3r p1 = particles.getPosition(i1);
    const Vector3r v0 = particles.getVelocity(i0);
    const Vector3r v1 = particles.getVelocity(i1);

    const Vector3r currentPoint = w0 * p0 + w1 * p1;
    const Vector3r currentVelocity = w0 * v0 + w1 * v1;
    Vector3r targetPoint = dragState.targetPoint;
    if (dragState.upwardPull)
        targetPoint.y() = std::max(targetPoint.y(), currentPoint.y() + static_cast<Real>(0.1));

    const Vector3r acceleration = computeSpringDampedAcceleration(
        spring * (targetPoint - currentPoint),
        currentVelocity,
        damping,
        maxAccel);

    (void)ropeCfg;

    particles.getVelocity(i0) += dt * w0 * acceleration;
    particles.getVelocity(i1) += dt * w1 * acceleration;
}

Real applyGroundConstraintToModel(PBD::SimulationModel* model, Real groundY)
{
    if (model == nullptr)
        return static_cast<Real>(0.0);

    Real maxPenetration = static_cast<Real>(0.0);
    PBD::ParticleData& particles = model->getParticles();
    for (unsigned int i = 0; i < particles.getNumberOfParticles(); ++i)
    {
        Vector3r p = particles.getPosition(i);
        const Real penetration = groundY - p.y();
        if (penetration > maxPenetration)
            maxPenetration = penetration;

        if (p.y() < groundY)
        {
            p.y() = groundY;
            particles.setPosition(i, p);

            Vector3r& v = particles.getVelocity(i);
            if (v.y() < 0.0)
                v.y() = 0.0;
        }
    }
    return maxPenetration;
}

static std::vector<CableEdgeRef> collectCableEdges(
    const ParticleData& pd,
    const std::vector<std::pair<unsigned int, unsigned int>>& pieceRanges)
{
    // Keep both local edge indices and global particle indices. Local indices
    // are used for arc-length exclusion inside one cable piece; global indices
    // are needed to read/write particles in the simulation model.
    std::vector<CableEdgeRef> edges;
    if (pieceRanges.empty())
        return edges;

    for (size_t pieceIndex = 0; pieceIndex < pieceRanges.size(); ++pieceIndex)
    {
        const std::pair<unsigned int, unsigned int>& range = pieceRanges[pieceIndex];
        if (range.second <= range.first)
            continue;
        if (range.second >= pd.getNumberOfParticles())
            continue;

        const unsigned int edgeCount = range.second - range.first;
        for (unsigned int localEdgeIndex = 0; localEdgeIndex < edgeCount; ++localEdgeIndex)
        {
            edges.push_back(
                CableEdgeRef{
                    pieceIndex,
                    localEdgeIndex,
                    range.first + localEdgeIndex });
        }
    }

    return edges;
}

static std::vector<std::pair<unsigned int, unsigned int>>
collectSelfCollisionCandidatesSpatialHash(
    const ParticleData& pd,
    const std::vector<CableEdgeRef>& edgeRefs,
    Real inflate,
    Real cellSize)
{
    if (edgeRefs.size() < 2)
        return {};

    // Broad phase: insert each edge's inflated AABB into a uniform grid. Edges
    // sharing a grid cell become narrow-phase candidates, avoiding an O(n^2)
    // all-pairs segment test in dense cable configurations.
    std::unordered_map<Int3, std::vector<unsigned int>, Int3Hash> grid;
    grid.reserve(edgeRefs.size() * 2);

    auto cellOf = [cellSize](const Vector3r& p) -> Int3
    {
        return Int3{
            static_cast<int>(std::floor(p.x() / cellSize)),
            static_cast<int>(std::floor(p.y() / cellSize)),
            static_cast<int>(std::floor(p.z() / cellSize))
        };
    };

    for (unsigned int e = 0; e < edgeRefs.size(); ++e)
    {
        const unsigned int p0Index = edgeRefs[e].startParticle;
        const unsigned int p1Index = p0Index + 1;
        const Vector3r p0 = pd.getPosition(p0Index);
        const Vector3r p1 = pd.getPosition(p1Index);

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
    uniquePairs.reserve(edgeRefs.size() * 8);

    // The same edge pair can appear in multiple cells; uniquePairs collapses
    // duplicates before narrow-phase closest-point checks.
    for (const auto& entry : grid)
    {
        const std::vector<unsigned int>& edgeIndices = entry.second;
        for (size_t a = 0; a < edgeIndices.size(); ++a)
        {
            for (size_t b = a + 1; b < edgeIndices.size(); ++b)
            {
                const unsigned int i = edgeIndices[a];
                const unsigned int j = edgeIndices[b];

                if (i == j)
                    continue;
                if (edgeRefs[i].pieceIndex == edgeRefs[j].pieceIndex &&
                    edgeRefs[i].localEdgeIndex + 1 >= edgeRefs[j].localEdgeIndex &&
                    edgeRefs[j].localEdgeIndex + 1 >= edgeRefs[i].localEdgeIndex)
                {
                    continue;
                }

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

void applyBendTwistRestStateCreep(
    SimulationModel* model,
    const RestShapeCreepConfig& creepCfg,
    unsigned int& creepCounter)
{
    if (!creepCfg.enabled || model == nullptr)
        return;

    creepCounter++;
    if (creepCounter < creepCfg.activationDelaySteps)
        return;
    if (creepCounter == creepCfg.activationDelaySteps)
    {
        // After the initial settling phase, capture the current relative
        // segment orientations as the new bend-twist rest state. This gives
        // the cable deformation memory without changing stretch rest lengths.
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
            if (bendTwist->m_restDarbouxVector.dot(currentRelative) < static_cast<Real>(0.0))
                currentRelative.coeffs() *= static_cast<Real>(-1.0);

            bendTwist->m_restDarbouxVector = currentRelative;
        }

        return;
    }
    if (creepCfg.stepInterval == 0 || (creepCounter % creepCfg.stepInterval) != 0)
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

        // Slow blend toward sustained deformation. The quaternion sign check
        // above keeps both rotations on the same hemisphere before blending.
        Eigen::Matrix<Real, 4, 1> blendedCoeffs =
            (static_cast<Real>(1.0) - creepCfg.creep) * restRelative.coeffs() +
            creepCfg.creep * currentRelative.coeffs();

        const Real norm = blendedCoeffs.norm();
        if (norm > static_cast<Real>(1e-12))
            restRelative.coeffs() = blendedCoeffs / norm;

    }
}

void applyPullForces(
    SimulationModel* model,
    Real dt,
    const std::vector<std::pair<unsigned int, unsigned int>>& pieceRanges,
    const PhysicsConfig& physicsCfg)
{
    ParticleData& particles = model->getParticles();

    auto pullOneParticle = [&](unsigned int particleIndex, const Vector3r& direction)
    {
        if (particleIndex >= particles.getNumberOfParticles())
            return;
        Vector3r& velocity = particles.getVelocity(particleIndex);
        const Vector3r acceleration = computeSpringDampedAcceleration(
            physicsCfg.pullSpring * direction,
            velocity,
            physicsCfg.pullDamping,
            physicsCfg.pullMaxAccel);

        velocity += dt * acceleration;
    };

    for (const std::pair<unsigned int, unsigned int>& range : pieceRanges)
    {
        const unsigned int firstParticle = range.first;
        const unsigned int lastParticle = range.second;
        if (firstParticle >= particles.getNumberOfParticles() ||
            lastParticle >= particles.getNumberOfParticles() ||
            firstParticle >= lastParticle)
        {
            continue;
        }

        if (physicsCfg.pullLeftScale > static_cast<Real>(0.0))
        {
            pullOneParticle(
                firstParticle,
                Vector3r(
                    -physicsCfg.pullLeftScale,
                    static_cast<Real>(0.0),
                    static_cast<Real>(0.0)));
        }
        if (physicsCfg.pullRightScale > static_cast<Real>(0.0))
        {
            pullOneParticle(
                lastParticle,
                Vector3r(
                    physicsCfg.pullRightScale,
                    static_cast<Real>(0.0),
                    static_cast<Real>(0.0)));
        }
    }
}

static std::vector<SelfCollisionContact> detectSelfCollisionContactsFromParticles(
    const ParticleData& pd,
    const std::vector<std::pair<unsigned int, unsigned int>>& pieceRanges,
    const RopeConfig& ropeCfg,
    const SelfCollisionConfig& collisionCfg)
{
    std::vector<SelfCollisionContact> contacts;
    const unsigned int n = pd.getNumberOfParticles();
    if (n < 4)
        return contacts;

    const Real minDist = static_cast<Real>(2.0) * collisionCfg.ropeRadius;
    const Real eps = static_cast<Real>(1e-8);
    const std::vector<CableEdgeRef> edgeRefs = collectCableEdges(pd, pieceRanges);
    const std::vector<std::pair<unsigned int, unsigned int>> candidates =
        collectSelfCollisionCandidatesSpatialHash(
            pd,
            edgeRefs,
            minDist,
            collisionCfg.hashCellSize);

    // Narrow phase: confirm spatial-hash candidates by segment-segment
    // closest distance. Adjacent edges on the same piece are ignored because
    // they are connected by construction and are not self-collisions.
    contacts.reserve(candidates.size());
    for (const auto& pair : candidates)
    {
        const CableEdgeRef& edgeA = edgeRefs[pair.first];
        const CableEdgeRef& edgeB = edgeRefs[pair.second];
        const unsigned int i = edgeA.startParticle;
        const unsigned int j = edgeB.startParticle;

        if (edgeA.pieceIndex == edgeB.pieceIndex)
        {
            const Real arcDist =
                static_cast<Real>(
                    std::abs(
                        static_cast<int>(edgeB.localEdgeIndex) -
                        static_cast<int>(edgeA.localEdgeIndex))) *
                ropeCfg.dx;
            if (arcDist < collisionCfg.excludeArcLength)
                continue;
        }

        if (i + 1 >= n || j + 1 >= n)
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

std::vector<SelfCollisionContact> detectSelfCollisionContacts(
    SimulationModel* model,
    const std::vector<std::pair<unsigned int, unsigned int>>& pieceRanges,
    const RopeConfig& ropeCfg,
    const SelfCollisionConfig& collisionCfg)
{
    if (model == nullptr)
        return {};
    return detectSelfCollisionContactsFromParticles(
        model->getParticles(), pieceRanges, ropeCfg, collisionCfg);
}

std::vector<SelfCollisionContact> detectSelfCollisionContactsForPieces(
    const std::vector<std::vector<Vector3r>>& pieces,
    const RopeConfig& ropeCfg,
    const SelfCollisionConfig& collisionCfg)
{
    ParticleData pd;
    std::vector<std::pair<unsigned int, unsigned int>> pieceRanges;

    for (const std::vector<Vector3r>& piece : pieces)
    {
        if (piece.size() < 2)
            continue;

        const unsigned int start = pd.getNumberOfParticles();
        for (const Vector3r& point : piece)
            pd.addVertex(point);
        pieceRanges.push_back(
            { start, start + static_cast<unsigned int>(piece.size()) - 1 });
    }

    return detectSelfCollisionContactsFromParticles(
        pd, pieceRanges, ropeCfg, collisionCfg);
}

void applySelfCollision(
    SimulationModel* model,
    Real dt,
    const std::vector<std::pair<unsigned int, unsigned int>>& pieceRanges,
    const RopeConfig& ropeCfg,
    const SelfCollisionConfig& collisionCfg,
    SelfCollisionDiagnostics* diagnostics)
{
    ParticleData& pd = model->getParticles();
    const unsigned int n = pd.getNumberOfParticles();
    if (n < 4)
        return;

    const Real minDist = static_cast<Real>(2.0) * collisionCfg.ropeRadius;
    const Real eps = static_cast<Real>(1e-8);
    const Real alpha = collisionCfg.compliance / (dt * dt);

    // These buffers implement a Jacobi-style contact pass: all endpoint
    // corrections are accumulated first and applied after the contacts are
    // processed. Repeating applySelfCollision several times gives iterative
    // convergence while reducing contact-order dependence inside one pass.
    static std::vector<Vector3r> posDelta;
    static std::vector<unsigned int> contactCount;
    static std::vector<VelocityPairCorrection> velocityCorrections;
    static std::unordered_map<std::uint64_t, Vector3r> previousContactNormals;

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
    std::unordered_map<std::uint64_t, Vector3r> currentContactNormals;
    currentContactNormals.reserve(n * 2);
    unsigned int normalFlipViolationCount = 0;

    const std::vector<CableEdgeRef> edgeRefs = collectCableEdges(pd, pieceRanges);
    const std::vector<std::pair<unsigned int, unsigned int>> candidates =
        collectSelfCollisionCandidatesSpatialHash(
            pd,
            edgeRefs,
            minDist,
            collisionCfg.hashCellSize);

    // Candidate pairs are broad-phase pairs after arc-length/adjacency
    // filtering. Penetration statistics are computed in this narrow-phase
    // loop, before any compliant, capped, averaged correction is applied.
    size_t candidatePairCount = 0;
    size_t penetratingContactCount = 0;
    Real maxPenetration = static_cast<Real>(0.0);
    Real sumPenetration = static_cast<Real>(0.0);
    size_t correctionEndpointCount = 0;
    size_t cappedCorrectionCount = 0;
    Real maxAppliedCorrection = static_cast<Real>(0.0);

    for (const auto& pair : candidates)
    {
        const CableEdgeRef& edgeA = edgeRefs[pair.first];
        const CableEdgeRef& edgeB = edgeRefs[pair.second];
        const unsigned int i = edgeA.startParticle;
        const unsigned int j = edgeB.startParticle;

        if (edgeA.pieceIndex == edgeB.pieceIndex)
        {
            const Real arcDist =
                static_cast<Real>(
                    std::abs(
                        static_cast<int>(edgeB.localEdgeIndex) -
                        static_cast<int>(edgeA.localEdgeIndex))) *
                ropeCfg.dx;
            if (arcDist < collisionCfg.excludeArcLength)
                continue;
        }

        if (i + 1 >= n || j + 1 >= n)
            continue;

        ++candidatePairCount;

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

        // Contact is represented at closest points c1 and c2, usually inside
        // the two segments. Parameters s and t become barycentric weights for
        // distributing the correction to the four endpoints.
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

        const Real penetration = -C;
        ++penetratingContactCount;
        sumPenetration += penetration;
        if (penetration > maxPenetration)
            maxPenetration = penetration;

        const std::uint64_t contactKey = makeEdgePairKey(i, j);
        const auto previousNormalIt = previousContactNormals.find(contactKey);
        if (previousNormalIt != previousContactNormals.end())
        {
            const Real normalDot = previousNormalIt->second.dot(nrm);
            if (normalDot < static_cast<Real>(0.0))
                ++normalFlipViolationCount;
        }
        currentContactNormals[contactKey] = nrm;

        // Endpoint weights for c1=(1-s)p0+s p1 and c2=(1-t)q0+t q1.
        const Real wP0 = static_cast<Real>(1.0) - s;
        const Real wP1 = s;
        const Real wQ0 = static_cast<Real>(1.0) - t;
        const Real wQ1 = t;

        const Vector3r gradP0 = -wP0 * nrm;
        const Vector3r gradP1 = -wP1 * nrm;
        const Vector3r gradQ0 = wQ0 * nrm;
        const Vector3r gradQ1 = wQ1 * nrm;

        const Real denom =
            gradP0.squaredNorm() +
            gradP1.squaredNorm() +
            gradQ0.squaredNorm() +
            gradQ1.squaredNorm() +
            alpha;
        if (denom <= eps)
            continue;

        // XPBD-style compliant projection for C = distance - 2r. The normal is
        // the constraint gradient direction; alpha softens the contact.
        const Real deltaLambda = -C / denom;

        const Vector3r rawDp0 = deltaLambda * gradP0;
        const Vector3r rawDp1 = deltaLambda * gradP1;
        const Vector3r rawDq0 = deltaLambda * gradQ0;
        const Vector3r rawDq1 = deltaLambda * gradQ1;
        Vector3r dp0 = clampVec(rawDp0, collisionCfg.maxCorrection);
        Vector3r dp1 = clampVec(rawDp1, collisionCfg.maxCorrection);
        Vector3r dq0 = clampVec(rawDq0, collisionCfg.maxCorrection);
        Vector3r dq1 = clampVec(rawDq1, collisionCfg.maxCorrection);

        correctionEndpointCount += 4;
        if (dp0.squaredNorm() + static_cast<Real>(1e-18) < rawDp0.squaredNorm())
            ++cappedCorrectionCount;
        if (dp1.squaredNorm() + static_cast<Real>(1e-18) < rawDp1.squaredNorm())
            ++cappedCorrectionCount;
        if (dq0.squaredNorm() + static_cast<Real>(1e-18) < rawDq0.squaredNorm())
            ++cappedCorrectionCount;
        if (dq1.squaredNorm() + static_cast<Real>(1e-18) < rawDq1.squaredNorm())
            ++cappedCorrectionCount;

        maxAppliedCorrection = std::max(maxAppliedCorrection, dp0.norm());
        maxAppliedCorrection = std::max(maxAppliedCorrection, dp1.norm());
        maxAppliedCorrection = std::max(maxAppliedCorrection, dq0.norm());
        maxAppliedCorrection = std::max(maxAppliedCorrection, dq1.norm());

        posDelta[i0] += dp0;
        contactCount[i0] += 1;
        posDelta[i1] += dp1;
        contactCount[i1] += 1;
        posDelta[j0] += dq0;
        contactCount[j0] += 1;
        posDelta[j1] += dq1;
        contactCount[j1] += 1;

        velocityCorrections.push_back({ i0, i1, j0, j1, wP0, wP1, wQ0, wQ1, nrm });
    }

    // Apply the averaged accumulated position corrections. The 0.5 scale keeps
    // a particle shared by many simultaneous contacts from overreacting.
    for (unsigned int i = 0; i < n; ++i)
    {
        if (contactCount[i] == 0)
            continue;

        const Real averageScale =
            static_cast<Real>(0.5) / static_cast<Real>(contactCount[i]);
        Vector3r corr = clampVec(averageScale * posDelta[i], collisionCfg.maxCorrection);
        pd.setPosition(i, pd.getPosition(i) + corr);
    }

    for (const auto& vc : velocityCorrections)
    {
        // Velocity damping is evaluated at the same segment contact points as
        // the position constraint, then distributed to endpoints with the same
        // barycentric weights.
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
            pd.getVelocity(vc.i0) -= vc.wP0 * impulse;
            pd.getVelocity(vc.i1) -= vc.wP1 * impulse;
            pd.getVelocity(vc.j0) += vc.wQ0 * impulse;
            pd.getVelocity(vc.j1) += vc.wQ1 * impulse;
        };

        if (vn < static_cast<Real>(0.0) && collisionCfg.normalDamping > static_cast<Real>(0.0))
            applyImpulse(-collisionCfg.normalDamping * vn * vc.nrm);

        // Tangential damping removes sliding energy along the contact plane but
        // does not act like Coulomb friction; it is a stabilizing velocity blend.
        const Vector3r vt = relV - relV.dot(vc.nrm) * vc.nrm;
        const Real vtNorm = vt.norm();
        if (vtNorm > eps)
        {
            if (collisionCfg.tangentialDamping > static_cast<Real>(0.0))
                applyImpulse(-collisionCfg.tangentialDamping * vt);
        }
    }

    // Residual penetration after all position and velocity corrections of this
    // iteration, measured on the current positions. This is the quantity that
    // bounds what could still be missed by a single discrete pass.
    size_t residualCandidateCount = 0;
    size_t residualPenetratingContactCount = 0;
    Real residualMaxPenetration = static_cast<Real>(0.0);
    Real residualSumPenetration = static_cast<Real>(0.0);
    if (diagnostics != nullptr)
    {
        for (const auto& pair : candidates)
        {
            const CableEdgeRef& edgeA = edgeRefs[pair.first];
            const CableEdgeRef& edgeB = edgeRefs[pair.second];
            if (edgeA.pieceIndex == edgeB.pieceIndex)
            {
                const Real arcDist =
                    static_cast<Real>(
                        std::abs(
                            static_cast<int>(edgeB.localEdgeIndex) -
                            static_cast<int>(edgeA.localEdgeIndex))) *
                    ropeCfg.dx;
                if (arcDist < collisionCfg.excludeArcLength)
                    continue;
            }

            const unsigned int i0 = edgeA.startParticle;
            const unsigned int j0 = edgeB.startParticle;
            if (i0 + 1 >= n || j0 + 1 >= n)
                continue;

            ++residualCandidateCount;

            Real s, t;
            Vector3r c1, c2;
            closestPtSegmentSegment(
                pd.getPosition(i0), pd.getPosition(i0 + 1),
                pd.getPosition(j0), pd.getPosition(j0 + 1),
                s, t, c1, c2);

            const Real residualPenetration = minDist - (c2 - c1).norm();
            if (residualPenetration > static_cast<Real>(0.0))
            {
                ++residualPenetratingContactCount;
                residualSumPenetration += residualPenetration;
                if (residualPenetration > residualMaxPenetration)
                    residualMaxPenetration = residualPenetration;
            }
        }

        diagnostics->candidatePairCount = candidatePairCount;
        diagnostics->penetratingContactCount = penetratingContactCount;
        diagnostics->maxPenetration = maxPenetration;
        diagnostics->sumPenetration = sumPenetration;
        diagnostics->correctionEndpointCount = correctionEndpointCount;
        diagnostics->cappedCorrectionCount = cappedCorrectionCount;
        diagnostics->maxAppliedCorrection = maxAppliedCorrection;
        diagnostics->normalFlipViolationCount = normalFlipViolationCount;
        diagnostics->residualCandidateCount = residualCandidateCount;
        diagnostics->residualPenetratingContactCount = residualPenetratingContactCount;
        diagnostics->residualMaxPenetration = residualMaxPenetration;
        diagnostics->residualSumPenetration = residualSumPenetration;
    }

    if (normalFlipViolationCount > 0)
    {
        std::cout << "[self-collision] possible side-switch normal flips: "
                  << normalFlipViolationCount << std::endl;
    }

    previousContactNormals.swap(currentContactNormals);
}
