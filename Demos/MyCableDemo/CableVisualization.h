#ifndef MY_CABLE_DEMO_CABLE_VISUALIZATION_H
#define MY_CABLE_DEMO_CABLE_VISUALIZATION_H

#include "Common/Common.h"
#include "CablePhysics.h"
#include "TanglementAnalyzer.h"

#include <cstddef>
#include <limits>
#include <vector>

namespace CableVisualization
{
    struct HighlightSelection
    {
        bool selected = false;
        size_t pieceIndex = std::numeric_limits<size_t>::max();
        size_t edgeIndex = std::numeric_limits<size_t>::max();
    };

    struct HoopMarkerSelection
    {
        bool selected = false;
        size_t pieceIndex = std::numeric_limits<size_t>::max();
        size_t hoopIndex = std::numeric_limits<size_t>::max();
    };

    struct CablePickResult
    {
        size_t pieceIndex = std::numeric_limits<size_t>::max();
        size_t edgeIndex = std::numeric_limits<size_t>::max();
        Real rayDepthT = static_cast<Real>(0.0);
        Real edgeT = static_cast<Real>(0.0);
        Vector3r grabbedPoint = Vector3r::Zero();
    };

    struct VisualizationState
    {
        std::vector<std::vector<TanglementAnalyzer::LocalKnotRegion>> knotRegions;
        std::vector<bool> unknotPieces;
        std::vector<std::vector<int>> crossingEdges;
        std::vector<std::vector<int>> crossingCenterEdges;
        std::vector<std::vector<int>> crossingMiddleEdges;
        std::vector<std::vector<TanglementAnalyzer::HoopCandidate>> scoredHoops;
        std::vector<std::vector<int>> problematicHoopMiddleEdges;
        std::vector<std::vector<int>> targetHoopMiddleEdges;

        void clear();
    };

    VisualizationState buildVisualizationState(
        const TanglementAnalyzer::TanglementAnalysisResult& tanglement);

    void renderCable(
        PBD::SimulationModel* model,
        const RopeConfig& ropeCfg,
        const HighlightSelection& selection,
        const HoopMarkerSelection& hoopSelection,
        const VisualizationState& state);

    bool buildPickRay(
        const ::Vector2i& screenPos,
        Vector3r& rayStart,
        Vector3r& rayEnd);

    bool selectCableEdgeAtScreenPosition(
        const ::Vector2i& screenPos,
        const std::vector<std::vector<Vector3r>>& ropePieces,
        Real renderRadius,
        CablePickResult& outPick);

    bool selectHoopMarkerAtScreenPosition(
        const ::Vector2i& screenPos,
        const std::vector<std::vector<Vector3r>>& ropePieces,
        Real renderRadius,
        const VisualizationState& state,
        HoopMarkerSelection& outPick);
}

#endif
