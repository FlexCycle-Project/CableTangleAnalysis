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

Real clamp01(Real x)
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

static Vector3r clampVec(const Vector3r& v, Real maxNorm)
{
    const Real len = v.norm();
    if (len > maxNorm && len > static_cast<Real>(1e-12))
        return v * (maxNorm / len);
    return v;
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
    ParticleData& pd = model->getParticles();
    const unsigned int n = pd.getNumberOfParticles();
    if (oldPositions.size() != n || dt <= static_cast<Real>(1e-12))
        return;

    blend = std::max(static_cast<Real>(0.0), std::min(static_cast<Real>(1.0), blend));

    for (unsigned int i = 0; i < n; ++i)
    {
        if (pd.getMass(i) <= 0.0)
            continue;

        const Vector3r vPos = (pd.getPosition(i) - oldPositions[i]) / dt;
        Vector3r& v = pd.getVelocity(i);
        v = (static_cast<Real>(1.0) - blend) * v + blend * vPos;
    }
}

static std::vector<CableEdgeRef> collectCableEdges(const ParticleData& pd)
{
    std::vector<CableEdgeRef> edges;
    if (g_cablePieceRanges.empty())
        return edges;

    for (size_t pieceIndex = 0; pieceIndex < g_cablePieceRanges.size(); ++pieceIndex)
    {
        const std::pair<unsigned int, unsigned int>& range = g_cablePieceRanges[pieceIndex];
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

void applyPullForces(SimulationModel* model, Real dt, const PhysicsConfig& physicsCfg)
{
    ParticleData& particles = model->getParticles();

    auto pullOneParticle = [&](unsigned int particleIndex, const Vector3r& direction)
    {
        if (particleIndex >= particles.getNumberOfParticles())
            return;
        if (particles.getMass(particleIndex) <= 0.0)
            return;

        Vector3r& velocity = particles.getVelocity(particleIndex);
        Vector3r acceleration =
            physicsCfg.pullSpring * direction -
            physicsCfg.pullDamping * velocity;

        const Real accelNorm = acceleration.norm();
        if (accelNorm > physicsCfg.pullMaxAccel)
            acceleration *= (physicsCfg.pullMaxAccel / accelNorm);

        velocity += dt * acceleration;
    };

    for (const std::pair<unsigned int, unsigned int>& range : g_cablePieceRanges)
    {
        const unsigned int firstParticle = range.first;
        const unsigned int lastParticle = range.second;
        if (firstParticle >= particles.getNumberOfParticles() ||
            lastParticle >= particles.getNumberOfParticles() ||
            firstParticle >= lastParticle)
        {
            continue;
        }

        Vector3r startDirection = particles.getPosition(firstParticle) -
            particles.getPosition(firstParticle + 1);
        Vector3r endDirection = particles.getPosition(lastParticle) -
            particles.getPosition(lastParticle - 1);

        if (startDirection.squaredNorm() > static_cast<Real>(1e-12))
            startDirection.normalize();
        else
            startDirection.setZero();

        if (endDirection.squaredNorm() > static_cast<Real>(1e-12))
            endDirection.normalize();
        else
            endDirection.setZero();

        pullOneParticle(firstParticle, startDirection);
        pullOneParticle(lastParticle, endDirection);
    }
}

void applyGroundConstraint(Real groundY)
{
    SimulationModel* model = Simulation::getCurrent()->getModel();
    ParticleData& particles = model->getParticles();

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

std::vector<SelfCollisionContact> detectSelfCollisionContacts(
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
    const std::vector<CableEdgeRef> edgeRefs = collectCableEdges(pd);
    const std::vector<std::pair<unsigned int, unsigned int>> candidates =
        collectSelfCollisionCandidatesSpatialHash(
            pd,
            edgeRefs,
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

void applySelfCollisionXPBD(
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

    const std::vector<CableEdgeRef> edgeRefs = collectCableEdges(pd);
    const std::vector<std::pair<unsigned int, unsigned int>> candidates =
        collectSelfCollisionCandidatesSpatialHash(
            pd,
            edgeRefs,
            minDist,
            collisionCfg.hashCellSize);

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
