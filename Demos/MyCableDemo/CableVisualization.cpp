#include "CableVisualization.h"

#include "Demos/Visualization/MiniGL.h"
#include "Simulation/SimulationModel.h"

#include <algorithm>
#include <cmath>

namespace CableVisualization
{
    namespace
    {
        static bool buildPickRayInternal(
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

        static float kBlueColor[4] = {0.0f, 0.0f, 0.8f, 1.0f};
        static float kGreenColor[4] = {0.20f, 0.85f, 0.25f, 1.0f};
        static float kRedColor[4] = {0.98f, 0.25f, 0.18f, 1.0f};
          static float kYellowColor[4] = {1.0f, 0.95f, 0.20f, 1.0f};
          static float kCyanColor[4] = {0.20f, 0.95f, 1.00f, 1.0f};
          static float kGrayColor[4] = {0.35f, 0.35f, 0.35f, 1.0f};
        static float kBlackColor[4] = {0.0f, 0.0f, 0.0f, 1.0f};
        static float kWhiteColor[4] = {1.0f, 1.0f, 1.0f, 1.0f};
        static float kPurpleColor[4] = {0.72f, 0.34f, 0.95f, 1.0f};

        static void setGrayscaleColor(Real t, float outColor[4])
        {
            const Real clampedT = std::max(static_cast<Real>(0.0), std::min(static_cast<Real>(1.0), t));
            const float grayscale = static_cast<float>(1.0 - clampedT);
            outColor[0] = grayscale;
            outColor[1] = grayscale;
            outColor[2] = grayscale;
            outColor[3] = 1.0f;
        }

        static float* grayscaleColorFromValue(Real t)
        {
            static float color[4];
            setGrayscaleColor(t, color);
            return color;
        }

        static bool isLocalIndexInKnotRegion(
            int localIndex,
            const TanglementAnalyzer::LocalKnotRegion& region)
        {
            return localIndex >= region.start && localIndex <= region.end;
        }

        static bool isLocalIndexHighlighted(
            size_t pieceIndex,
            int localIndex,
            const VisualizationState& state)
        {
            if (pieceIndex >= state.knotRegions.size())
                return false;

            const std::vector<TanglementAnalyzer::LocalKnotRegion>& regions =
                state.knotRegions[pieceIndex];
            for (const TanglementAnalyzer::LocalKnotRegion& region : regions)
            {
                if (isLocalIndexInKnotRegion(localIndex, region))
                    return true;
            }
            return false;
        }

        static bool isLocalEdgeHighlighted(
            size_t pieceIndex,
            int localEdgeIndex,
            const VisualizationState& state)
        {
            return isLocalIndexHighlighted(pieceIndex, localEdgeIndex, state) ||
                   isLocalIndexHighlighted(pieceIndex, localEdgeIndex + 1, state);
        }

        static bool isLocalEdgeInHighlight(
            const std::vector<std::vector<int>>& edgesPerPiece,
            size_t pieceIndex,
            int localEdgeIndex)
        {
            if (pieceIndex >= edgesPerPiece.size())
                return false;

            const std::vector<int>& edges = edgesPerPiece[pieceIndex];
            return std::find(edges.begin(), edges.end(), localEdgeIndex) != edges.end();
        }

    }

    bool buildPickRay(
        const ::Vector2i& screenPos,
        Vector3r& rayStart,
        Vector3r& rayEnd)
    {
        return buildPickRayInternal(screenPos, rayStart, rayEnd);
    }

        void VisualizationState::clear()
        {
            knotRegions.clear();
            unknotPieces.clear();
            crossingEdges.clear();
            crossingCenterEdges.clear();
            crossingMiddleEdges.clear();
            scoredHoops.clear();
            problematicHoopMiddleEdges.clear();
            targetHoopMiddleEdges.clear();
            segmentScores.clear();
    }

    VisualizationState buildVisualizationState(
        const TanglementAnalyzer::TanglementAnalysisResult& tanglement)
    {
        VisualizationState state;
        state.knotRegions = tanglement.knotRegions;
        state.unknotPieces.resize(tanglement.kymoResults.size(), false);
        state.crossingMiddleEdges = tanglement.crossingMiddleEdges;
        state.scoredHoops = tanglement.scoredHoops;
        state.problematicHoopMiddleEdges = tanglement.problematicHoopMiddleEdges;
        state.targetHoopMiddleEdges = tanglement.targetHoopMiddleEdges;
        state.segmentScores = tanglement.segmentScores;

        for (size_t pieceIndex = 0; pieceIndex < tanglement.kymoResults.size(); ++pieceIndex)
        {
            const PBD::KymoKnotAnalysisResult& knotResult = tanglement.kymoResults[pieceIndex];
            if (PBD::KymoKnotAnalyzer::isUnknotResult(knotResult))
            {
                state.unknotPieces[pieceIndex] = true;
            }
        }
        return state;
    }

    void renderCable(
        PBD::SimulationModel* model,
        const RopeConfig& ropeCfg,
        const HighlightSelection& selection,
        const HoopMarkerSelection& hoopSelection,
        const VisualizationState& state,
        bool showSegmentScores)
    {
        if (model == nullptr)
            return;

        const PBD::ParticleData& particles = model->getParticles();
        std::vector<Vector3r> targetMarkers;
        for (unsigned int m = 0; m < model->getLineModels().size(); ++m)
        {
            PBD::LineModel* lineModel = model->getLineModels()[m];
            if (lineModel == nullptr)
                continue;

            const unsigned int pointOffset = lineModel->getIndexOffset();
            const bool pieceIsUnknot =
                m < state.unknotPieces.size() && state.unknotPieces[m];

            for (unsigned int e = 0; e < lineModel->getEdges().size(); ++e)
            {
                const unsigned int i0 = lineModel->getEdges()[e].m_vert[0] + pointOffset;
                const unsigned int i1 = lineModel->getEdges()[e].m_vert[1] + pointOffset;

                const Vector3r& p0 = particles.getPosition(i0);
                const Vector3r& p1 = particles.getPosition(i1);

                const int localEdgeIndex = static_cast<int>(e);
                const bool highlightKnot = isLocalEdgeHighlighted(m, localEdgeIndex, state);
                bool highlightSelectedHoop = false;
                bool highlightGrayHoop = false;
                int grayLeftStart = -1;
                int grayLeftEnd = -1;
                int grayRightStart = -1;
                int grayRightEnd = -1;
                if (hoopSelection.selected &&
                    m == hoopSelection.pieceIndex &&
                    m < state.scoredHoops.size() &&
                    hoopSelection.hoopIndex < state.scoredHoops[m].size())
                {
                    const TanglementAnalyzer::HoopCandidate& selectedHoop =
                        state.scoredHoops[m][hoopSelection.hoopIndex];
                    const int selectedStart = std::min(selectedHoop.startEdge, selectedHoop.endEdge);
                    const int selectedEnd = std::max(selectedHoop.startEdge, selectedHoop.endEdge);

                    for (size_t hoopIndex = 0; hoopIndex < state.scoredHoops[m].size(); ++hoopIndex)
                    {
                        if (hoopIndex == hoopSelection.hoopIndex)
                            continue;

                        const TanglementAnalyzer::HoopCandidate& hoop = state.scoredHoops[m][hoopIndex];
                        if (hoop.pairSumScore == 0)
                            continue;

                        const int hoopStart = std::min(hoop.startEdge, hoop.endEdge);
                        const int hoopEnd = std::max(hoop.startEdge, hoop.endEdge);
                        const bool overlapsSelected =
                            hoopStart <= selectedEnd && hoopEnd >= selectedStart;
                        if (!overlapsSelected)
                            continue;

                        if (hoopStart < selectedStart &&
                            (grayLeftStart < 0 || hoopStart < grayLeftStart))
                        {
                            grayLeftStart = hoopStart;
                            grayLeftEnd = std::min(hoopEnd, selectedStart - 1);
                        }

                        if (hoopEnd > selectedEnd &&
                            (grayRightEnd < 0 || hoopEnd > grayRightEnd))
                        {
                            grayRightStart = std::max(hoopStart, selectedEnd + 1);
                            grayRightEnd = hoopEnd;
                        }
                    }

                    if (localEdgeIndex >= selectedStart && localEdgeIndex <= selectedEnd)
                    {
                        highlightSelectedHoop = true;
                    }
                    else if ((grayLeftStart >= 0 &&
                              localEdgeIndex >= grayLeftStart &&
                              localEdgeIndex <= grayLeftEnd) ||
                             (grayRightStart >= 0 &&
                              localEdgeIndex >= grayRightStart &&
                              localEdgeIndex <= grayRightEnd))
                    {
                        highlightGrayHoop = true;
                    }
                }
                const bool highlightTarget =
                    isLocalEdgeInHighlight(state.targetHoopMiddleEdges, m, localEdgeIndex);
                const bool highlightSelected =
                    selection.selected &&
                    m == selection.pieceIndex &&
                    e == selection.edgeIndex;

                float* activeColor = kBlueColor;
                float* sphereColor = kBlueColor;
                if (pieceIsUnknot)
                {
                    activeColor = kGreenColor;
                    sphereColor = kGreenColor;
                }
                // if (highlightKnot)
                // {
                //     activeColor = kRedColor;
                //     sphereColor = kRedColor;
                // }
                // if (highlightGrayHoop)
                // {
                //     activeColor = kGrayColor;
                //     sphereColor = kGrayColor;
                // }
                // if (highlightSelectedHoop)
                // {
                //     activeColor = kWhiteColor;
                //     sphereColor = kWhiteColor;
                // }
                if (highlightSelected)
                {
                    activeColor = kWhiteColor;
                    sphereColor = kWhiteColor;
                }

                if (hoopSelection.selected &&
                    m == hoopSelection.pieceIndex &&
                    m < state.scoredHoops.size() &&
                    hoopSelection.hoopIndex < state.scoredHoops[m].size())
                {
                    const TanglementAnalyzer::HoopCandidate& selectedHoop =
                        state.scoredHoops[m][hoopSelection.hoopIndex];
                    const int selectedHoopStart = std::max(0, std::min(selectedHoop.startEdge, selectedHoop.endEdge));
                    const int selectedHoopEnd = std::min(
                        static_cast<int>(lineModel->getEdges().size()) - 1,
                        std::max(selectedHoop.startEdge, selectedHoop.endEdge));
                    if (localEdgeIndex >= selectedHoopStart &&
                        localEdgeIndex <= selectedHoopEnd)
                    {
                        activeColor = kWhiteColor;
                        sphereColor = kWhiteColor;
                    }
                }

                float* scoreColor = nullptr;
                if (showSegmentScores &&
                    m < state.segmentScores.size() &&
                    static_cast<size_t>(localEdgeIndex) < state.segmentScores[m].size())
                {
                    const std::vector<Real>& pieceSegmentScores = state.segmentScores[m];
                    const auto minmaxScores = std::minmax_element(
                        pieceSegmentScores.begin(),
                        pieceSegmentScores.end());
                    const Real minScore = (minmaxScores.first != pieceSegmentScores.end())
                        ? *minmaxScores.first
                        : static_cast<Real>(0.0);
                    const Real maxScore = (minmaxScores.second != pieceSegmentScores.end())
                        ? *minmaxScores.second
                        : static_cast<Real>(0.0);
                    const Real score = pieceSegmentScores[static_cast<size_t>(localEdgeIndex)];
                    const Real normalizedScore =
                        std::abs(maxScore - minScore) <= static_cast<Real>(1e-8)
                            ? static_cast<Real>(0.5)
                            : (score - minScore) / (maxScore - minScore);
                    scoreColor = grayscaleColorFromValue(normalizedScore);
                    activeColor = scoreColor;
                }

                PBD::MiniGL::drawSphere(
                    p0,
                    0.03f,
                    sphereColor);
                if (e == lineModel->getEdges().size() - 1)
                {
                    PBD::MiniGL::drawSphere(
                        p1,
                        0.03f,
                        sphereColor);
                }

                PBD::MiniGL::drawCylinder(
                    p0,
                    p1,
                    activeColor,
                    static_cast<float>(ropeCfg.renderRadius));

                if (m < state.scoredHoops.size())
                {
                    for (size_t hoopIndex = 0; hoopIndex < state.scoredHoops[m].size(); ++hoopIndex)
                    {
                        const TanglementAnalyzer::HoopCandidate& hoop = state.scoredHoops[m][hoopIndex];
                        int hoopMidEdge = static_cast<int>(
                            std::round(
                                static_cast<Real>(0.5) *
                                (static_cast<Real>(hoop.startEdge) +
                                 static_cast<Real>(hoop.endEdge))));
                        hoopMidEdge = std::max(0, std::min(hoopMidEdge, static_cast<int>(lineModel->getEdges().size()) - 1));
                        if (localEdgeIndex != hoopMidEdge)
                            continue;

                        float* markerColor = hoop.pairSumScore > 0 ? kRedColor : kCyanColor;
                        const bool markerSelected =
                            hoopSelection.selected &&
                            m == hoopSelection.pieceIndex &&
                            hoopSelection.hoopIndex == hoopIndex;
                        if (markerSelected)
                            markerColor = kWhiteColor;

                        const GLboolean depthWasEnabled = glIsEnabled(GL_DEPTH_TEST);
                        if (depthWasEnabled)
                            glDisable(GL_DEPTH_TEST);
                        PBD::MiniGL::drawSphere(
                            static_cast<Real>(0.5) * (p0 + p1),
                            markerSelected ? 0.065f : 0.045f,
                            markerColor);
                        if (depthWasEnabled)
                            glEnable(GL_DEPTH_TEST);
                    }
                }

                // if (highlightTarget)
                // {
                //     const Vector3r midpoint = static_cast<Real>(0.5) * (p0 + p1);
                //     targetMarkers.push_back(midpoint + Vector3r(0.0, 0.015, 0.0));
                // }
            }
        }

        for (const Vector3r& markerPos : targetMarkers)
        {
            const GLboolean depthWasEnabled = glIsEnabled(GL_DEPTH_TEST);
            if (depthWasEnabled)
                glDisable(GL_DEPTH_TEST);
            PBD::MiniGL::drawSphere(
                markerPos,
                0.03f,
                kRedColor);
            if (depthWasEnabled)
                glEnable(GL_DEPTH_TEST);
        }
    }

    bool selectCableEdgeAtScreenPosition(
        const ::Vector2i& screenPos,
        const std::vector<std::vector<Vector3r>>& ropePieces,
        Real renderRadius,
        CablePickResult& outPick)
    {
        outPick = CablePickResult();
        if (ropePieces.empty())
            return false;

        Vector3r rayStart(0.0, 0.0, 0.0);
        Vector3r rayEnd(0.0, 0.0, 0.0);
        if (!buildPickRayInternal(screenPos, rayStart, rayEnd))
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
            static_cast<Real>(2.0) * renderRadius);
        if (bestPiece == std::numeric_limits<size_t>::max() || bestDistance > dragThreshold)
            return false;

        outPick.pieceIndex = bestPiece;
        outPick.edgeIndex = bestEdge;
        outPick.rayDepthT = bestRayT;
        outPick.edgeT = bestEdgeT;
        outPick.grabbedPoint = bestGrabbedPoint;
        return true;
    }

    bool selectHoopMarkerAtScreenPosition(
        const ::Vector2i& screenPos,
        const std::vector<std::vector<Vector3r>>& ropePieces,
        Real renderRadius,
        const VisualizationState& state,
        HoopMarkerSelection& outPick)
    {
        outPick = HoopMarkerSelection();
        if (ropePieces.empty())
            return false;

        Vector3r rayStart(0.0, 0.0, 0.0);
        Vector3r rayEnd(0.0, 0.0, 0.0);
        if (!buildPickRayInternal(screenPos, rayStart, rayEnd))
            return false;

        auto distanceRayPoint = [](const Vector3r& a, const Vector3r& b, const Vector3r& p) -> Real
        {
            const Vector3r ab = b - a;
            const Real abLen2 = ab.squaredNorm();
            if (abLen2 <= static_cast<Real>(1e-12))
                return (p - a).norm();

            Real t = (p - a).dot(ab) / abLen2;
            t = std::max(static_cast<Real>(0.0), std::min(static_cast<Real>(1.0), t));
            const Vector3r closest = a + t * ab;
            return (closest - p).norm();
        };

        size_t bestPiece = std::numeric_limits<size_t>::max();
        size_t bestHoopIndex = std::numeric_limits<size_t>::max();
        Real bestDistance = std::numeric_limits<Real>::max();

        for (size_t pieceIndex = 0; pieceIndex < ropePieces.size(); ++pieceIndex)
        {
            if (pieceIndex >= state.scoredHoops.size())
                continue;

            const std::vector<Vector3r>& piece = ropePieces[pieceIndex];
            for (size_t hoopIndex = 0; hoopIndex < state.scoredHoops[pieceIndex].size(); ++hoopIndex)
            {
                const TanglementAnalyzer::HoopCandidate& hoop = state.scoredHoops[pieceIndex][hoopIndex];
                const int hoopStart = std::min(hoop.startEdge, hoop.endEdge);
                const int hoopEnd = std::max(hoop.startEdge, hoop.endEdge);
                const int middleEdge = static_cast<int>(
                    std::round(
                        static_cast<Real>(0.5) *
                        (static_cast<Real>(hoop.startEdge) +
                         static_cast<Real>(hoop.endEdge))));
                if (hoopStart < 0 || hoopEnd + 1 >= static_cast<int>(piece.size()) ||
                    middleEdge < 0 || middleEdge + 1 >= static_cast<int>(piece.size()))
                    continue;

                const Vector3r markerPos =
                    static_cast<Real>(0.5) *
                    (piece[static_cast<size_t>(middleEdge)] +
                     piece[static_cast<size_t>(middleEdge + 1)]) +
                    Vector3r(0.0, 0.015, 0.0);
                const Real distance = distanceRayPoint(rayStart, rayEnd, markerPos);
                if (distance < bestDistance)
                {
                    bestDistance = distance;
                    bestPiece = pieceIndex;
                    bestHoopIndex = hoopIndex;
                }
            }
        }

        const Real dragThreshold = std::max(static_cast<Real>(0.10), static_cast<Real>(2.0) * renderRadius);
        if (bestPiece == std::numeric_limits<size_t>::max() || bestDistance > dragThreshold)
            return false;

        outPick.selected = true;
        outPick.pieceIndex = bestPiece;
        outPick.hoopIndex = bestHoopIndex;
        return true;
    }

}
