#ifndef MY_CABLE_DEMO_TANGLEMENT_ANALYZER_H
#define MY_CABLE_DEMO_TANGLEMENT_ANALYZER_H

#include "Common/Common.h"
#include "CablePhysics.h"
#include "KymoKnotAnalyzer.h"

#include <string>
#include <utility>
#include <vector>

class TanglementAnalyzer
{
public:
    struct LocalKnotRegion
    {
        int start = -1;
        int end = -1;
        size_t contactCount = 0;
    };

    struct CrossingEvent
    {
        Real order = static_cast<Real>(0.0);
        size_t pieceIndex = 0;
        int firstEdge = -1;
        int secondEdge = -1;
        Real heightDelta = static_cast<Real>(0.0);
        bool firstOverSecond = false;
    };

    struct HoopCandidate
    {
        int startEdge = -1;
        int endEdge = -1;
        Real order = static_cast<Real>(0.0);
        int pairSumScore = 0;
    };

    struct TanglementMetricsSnapshot
    {
        size_t pieceCount = 0;
        size_t totalPointCount = 0;
        size_t selfCollisionCount = 0;
        Real knotSpanFraction = static_cast<Real>(0.0);
        int gravityCrossings = 0;
        size_t contactParticipationCount = 0;
    };

    struct TanglementAnalysisResult
    {
        std::vector<CrossingEvent> allCrossingEvents;
        std::vector<std::vector<LocalKnotRegion>> knotRegions;
        std::vector<PBD::KymoKnotAnalysisResult> kymoResults;
        std::vector<std::vector<int>> crossingEdges;
        std::vector<std::vector<int>> crossingCenterEdges;
        std::vector<std::vector<int>> crossingMiddleEdges;
        std::vector<std::vector<HoopCandidate>> scoredHoops;
        std::vector<std::vector<int>> problematicHoopMiddleEdges;
        std::vector<std::vector<int>> targetHoopMiddleEdges;
        TanglementMetricsSnapshot metrics;
    };

    static TanglementAnalysisResult analyzeTanglement(
        const std::vector<std::vector<Vector3r>>& ropePieces,
        const std::vector<SelfCollisionContact>& contacts,
        PBD::KymoKnotSearchMode searchMode);

    static void printTanglementReport(
        const std::vector<std::vector<Vector3r>>& ropePieces,
        const TanglementAnalysisResult& result);

};

#endif
