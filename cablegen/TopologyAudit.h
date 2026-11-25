#ifndef MY_CABLE_DEMO_TOPOLOGY_AUDIT_H
#define MY_CABLE_DEMO_TOPOLOGY_AUDIT_H

#include "Common/Common.h"

#include <limits>
#include <string>
#include <vector>

namespace TopologyAudit
{
    // A crossing of the XZ projection. heightDelta is yA - yB at the
    // projected crossing point. Because the two points have the same x and z,
    // |heightDelta| is also their 3D centre-line distance at the crossing.
    struct Crossing
    {
        size_t pieceA = 0;
        int edgeA = -1;
        size_t pieceB = 0;
        int edgeB = -1;
        Real x = static_cast<Real>(0.0);
        Real z = static_cast<Real>(0.0);
        Real heightDelta = static_cast<Real>(0.0);
        bool firstOverSecond = false;
    };

    struct Settings
    {
        // A centre-line touch/pass-through threshold. This is deliberately
        // much smaller than the cable radius; it is a topology threshold.
        Real contactEpsilon = static_cast<Real>(1e-4);
        // Only used to report surface overlap, not to declare a topology
        // violation.
        Real surfaceRadius = static_cast<Real>(0.08);
        // XZ broad-phase cell size. Edges are short compared with this.
        Real cellSize = static_cast<Real>(1.0);
        Real parallelEpsilon = static_cast<Real>(1e-10);
        Real endpointEpsilon = static_cast<Real>(1e-8);
    };

    struct Violation
    {
        std::string type;
        unsigned int step = 0;
        Crossing previous;
        Crossing current;
        bool hasPrevious = false;
        bool hasCurrent = false;
        Real previousHeight = static_cast<Real>(0.0);
        Real currentHeight = static_cast<Real>(0.0);
    };

    struct StepResult
    {
        size_t crossingCount = 0;
        size_t centerlineContactCount = 0;
        size_t surfaceOverlapCount = 0;
        Real minAbsHeight = std::numeric_limits<Real>::infinity();
        std::vector<Violation> violations;
    };

    // Computes all XZ projection crossings for non-adjacent edges.
    std::vector<Crossing> computeCrossings(
        const std::vector<std::vector<Vector3r>>& pieces,
        const Settings& settings);

    // Stateful auditor. begin() establishes the reference crossing diagram;
    // audit() compares the current diagram with the previous one.
    class Auditor
    {
    public:
        void reset();
        void begin(
            const std::vector<std::vector<Vector3r>>& pieces,
            const Settings& settings);
        StepResult audit(
            unsigned int step,
            const std::vector<std::vector<Vector3r>>& pieces);

        bool hasReference() const { return hasPrevious_; }
        const Settings& settings() const { return settings_; }

    private:
        Settings settings_;
        std::vector<Crossing> previous_;
        bool hasPrevious_ = false;
    };
}

#endif
