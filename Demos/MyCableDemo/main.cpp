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
#include <unistd.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
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
static void printSelectedHoopScoreIfChanged();

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
static void restoreDemoWindow();
static bool loadKnotSeedPresetFromFile(const std::string& selector);
static bool selectKnotSeedPresetByType(const std::string& knotType);
static std::vector<std::string> detectCurrentKnotTypes(SimulationModel* model);
static void saveCurrentKnotSeedPreset();
static void evaluateCurrentCableTanglement(bool saveSeedPreset);
static void evaluateAndRecordCatalogSample(bool printReport);
static void evaluateAndRecordCutHeuristicSample();
static void evaluateAndRecordCutsToZeroSample();
static void appendCatalogRows(
    const std::vector<std::vector<Vector3r>>& ropePieces,
    const std::vector<SelfCollisionContact>& contacts,
    const TanglementAnalyzer::TanglementAnalysisResult& tanglement);
static void prepareCatalogSeedForNextRun();
static bool storeKnotSeedPresetEntry(
    const std::string& knotType,
    unsigned int seed,
    bool appendToFile);
static bool isNontrivialKnotType(const std::string& knotType);
static std::vector<std::string> detectKnotTypesForPieces(
    const std::vector<std::vector<Vector3r>>& ropePieces,
    const std::vector<SelfCollisionContact>& contacts);
static std::string resolveDefaultKnotSeedBankPath();
static std::string resolveDefaultKnotSeedCatalogPath();
static std::string resolveDefaultCutEvalOutputPath();
static std::string resolveDefaultCutsToZeroEvalOutputPath();
static size_t countCurvatureSamples(const std::vector<Vector3r>& points);
static Real computeTotalCurvature(const std::vector<Vector3r>& points);
static Real computeCurvatureOutsideKnotRegions(
    const std::vector<Vector3r>& points,
    const PBD::KymoKnotAnalysisResult& knotResult,
    size_t* sampleCount = nullptr);
static void generateNextCatalogSeed();
static bool loadCutEvalCatalog();
static bool prepareCutEvalOutputFile();
static bool prepareCutsToZeroEvalOutputFile();
static void prepareCutEvalSeed(size_t index);
static void generateNextCutEvalSeed();

static DemoBase* base = nullptr;
static AppConfig g_cfg;
static unsigned int g_restCreepCounter = 0;
static std::vector<PBD::KymoKnotAnalysisResult> g_cachedCableKymoKnotResults;
static bool g_hasCableKymoKnotResults = false;
static CableVisualization::VisualizationState g_visualizationState;
static PBD::KymoKnotSearchMode g_kymoKnotSearchMode = PBD::KymoKnotSearchMode::TopDown;
static std::vector<std::vector<Vector3r>> g_cablePieces;
static std::vector<std::pair<unsigned int, unsigned int>> g_cablePieceRanges;
static std::map<std::string, std::vector<unsigned int>> g_knotSeedPresets;
static std::string g_knotSeedBankPath;
static std::string g_knotSeedCatalogPath;
static bool g_useKnotSeedPreset = false;
static std::string g_selectedKnotType;
static unsigned int g_selectedKnotSeed = 0;
static bool g_useDirectRandomSeed = false;
static unsigned int g_directRandomSeed = 0;
static bool g_catalogRunEnabled = false;
static unsigned int g_catalogSamplesWritten = 0;
static unsigned int g_catalogMaxSamples = 0;
static unsigned int g_catalogNextSeed = 0;
static unsigned int g_catalogStepsPerFrame = 1;
static bool g_catalogRenderEnabled = false;
static std::mt19937 g_catalogRng;
static bool g_evalRunEnabled = false;
static std::string g_evalOutputPath;
static unsigned int g_evalRandomTrials = 10;
static bool g_evalCutsToZeroEnabled = false;
static std::string g_evalCutsToZeroOutputPath;
static unsigned int g_evalCutsToZeroMaxCuts = 10;
static unsigned int g_evalCutsToZeroRandomTrials = 20;
static size_t g_evalCurrentIndex = 0;
static size_t g_evalSamplesWritten = 0;
struct CutEvalSeedRecord
{
    std::string knotType;
    unsigned int seed = 0;

    CutEvalSeedRecord() = default;
    CutEvalSeedRecord(const std::string& type, unsigned int s)
        : knotType(type), seed(s) {}
};
static std::vector<CutEvalSeedRecord> g_evalSeeds;
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
static CableVisualization::HoopMarkerSelection g_lastPrintedHoopSelection;
static bool g_autoKnotEvaluationDone = false;
static const Real g_mouseDragSpring = static_cast<Real>(1800.0);
static const Real g_mouseDragDamping = static_cast<Real>(24.0);
static const Real g_mouseDragMaxAccel = static_cast<Real>(3000.0);

static void printSelectedHoopScoreIfChanged()
{
    if (!g_hoopMarkerSelection.selected)
    {
        g_lastPrintedHoopSelection = CableVisualization::HoopMarkerSelection();
        return;
    }

    if (g_lastPrintedHoopSelection.selected &&
        g_lastPrintedHoopSelection.pieceIndex == g_hoopMarkerSelection.pieceIndex &&
        g_lastPrintedHoopSelection.hoopIndex == g_hoopMarkerSelection.hoopIndex)
    {
        return;
    }

    if (g_hoopMarkerSelection.pieceIndex >= g_visualizationState.scoredHoops.size() ||
        g_hoopMarkerSelection.hoopIndex >= g_visualizationState.scoredHoops[g_hoopMarkerSelection.pieceIndex].size())
    {
        return;
    }

    const TanglementAnalyzer::HoopCandidate& hoop =
        g_visualizationState.scoredHoops[g_hoopMarkerSelection.pieceIndex][g_hoopMarkerSelection.hoopIndex];
    std::cout << "Selected hoop score: piece " << g_hoopMarkerSelection.pieceIndex
              << ", hoop " << g_hoopMarkerSelection.hoopIndex
              << ", score=" << hoop.pairSumScore << std::endl;
    TanglementAnalyzer::printHoopScoreHitLists(
        collectCurrentCablePieces(Simulation::getCurrent()->getModel()),
        g_hoopMarkerSelection.pieceIndex,
        hoop);
    g_lastPrintedHoopSelection = g_hoopMarkerSelection;
}

static std::string resolveDefaultKnotSeedBankPath()
{
    char exePath[4096] = {0};
    const ssize_t len = readlink("/proc/self/exe", exePath, sizeof(exePath) - 1);
    if (len <= 0)
        return "Demos/MyCableDemo/knot_seed_bank.txt";

    exePath[len] = '\0';
    std::string rawPath(exePath);
    const size_t lastSlash = rawPath.find_last_of('/');
    if (lastSlash == std::string::npos)
        return "Demos/MyCableDemo/knot_seed_bank.txt";

    const std::string exeDir = rawPath.substr(0, lastSlash);
    const std::string candidate = exeDir + "/../Demos/MyCableDemo/knot_seed_bank.txt";

    char resolved[4096] = {0};
    if (realpath(candidate.c_str(), resolved) == nullptr)
        return candidate;

    return std::string(resolved);
}

static std::string resolveDefaultKnotSeedCatalogPath()
{
    std::string catalogPath = resolveDefaultKnotSeedBankPath();
    const size_t slashPos = catalogPath.find_last_of("\\/");
    if (slashPos == std::string::npos)
        return "knot_seed_catalog.tsv";

    return catalogPath.substr(0, slashPos + 1) + "knot_seed_catalog.tsv";
}

static std::string resolveDefaultCutEvalOutputPath()
{
    std::string catalogPath = resolveDefaultKnotSeedCatalogPath();
    const size_t slashPos = catalogPath.find_last_of("\\/");
    if (slashPos == std::string::npos)
        return "knot_cut_eval.tsv";

    return catalogPath.substr(0, slashPos + 1) + "knot_cut_eval.tsv";
}

static std::string resolveDefaultCutsToZeroEvalOutputPath()
{
    std::string catalogPath = resolveDefaultKnotSeedCatalogPath();
    const size_t slashPos = catalogPath.find_last_of("\\/");
    if (slashPos == std::string::npos)
        return "knot_cuts_to_zero.tsv";

    return catalogPath.substr(0, slashPos + 1) + "knot_cuts_to_zero.tsv";
}

static std::vector<std::string> splitTabSeparatedLine(const std::string& line)
{
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string field;
    while (std::getline(ss, field, '\t'))
        fields.push_back(field);
    return fields;
}

static bool loadCutEvalCatalog()
{
    g_evalSeeds.clear();

    std::ifstream in(g_knotSeedCatalogPath);
    if (!in)
    {
        std::cout << "[cut-eval] could not open input catalog "
                  << g_knotSeedCatalogPath << std::endl;
        return false;
    }

    std::string line;
    if (!std::getline(in, line))
    {
        std::cout << "[cut-eval] input catalog is empty." << std::endl;
        return false;
    }

    const std::vector<std::string> header = splitTabSeparatedLine(line);
    int knotTypeColumn = -1;
    int seedColumn = -1;
    for (size_t i = 0; i < header.size(); ++i)
    {
        if (header[i] == "knotType")
            knotTypeColumn = static_cast<int>(i);
        else if (header[i] == "seed")
            seedColumn = static_cast<int>(i);
    }

    if (knotTypeColumn < 0 || seedColumn < 0)
    {
        std::cout << "[cut-eval] input catalog must contain knotType and seed columns." << std::endl;
        return false;
    }

    std::set<unsigned int> seenSeeds;
    while (std::getline(in, line))
    {
        if (line.empty())
            continue;

        const std::vector<std::string> fields = splitTabSeparatedLine(line);
        const int requiredColumn = std::max(knotTypeColumn, seedColumn);
        if (static_cast<int>(fields.size()) <= requiredColumn)
            continue;

        try
        {
            const unsigned int seed = static_cast<unsigned int>(
                std::stoul(fields[static_cast<size_t>(seedColumn)]));
            if (seed == 0 || seenSeeds.find(seed) != seenSeeds.end())
                continue;

            seenSeeds.insert(seed);
            g_evalSeeds.push_back(
                CutEvalSeedRecord{
                    fields[static_cast<size_t>(knotTypeColumn)],
                    seed});
        }
        catch (...)
        {
        }
    }

    std::cout << "[cut-eval] loaded " << g_evalSeeds.size()
              << " unique seed(s) from " << g_knotSeedCatalogPath << std::endl;
    return !g_evalSeeds.empty();
}

static bool prepareCutEvalOutputFile()
{
    if (g_evalOutputPath.empty())
        g_evalOutputPath = resolveDefaultCutEvalOutputPath();

    std::ofstream out(g_evalOutputPath, std::ios::trunc);
    if (!out)
    {
        std::cout << "[cut-eval] could not open output file "
                  << g_evalOutputPath << std::endl;
        return false;
    }

    out << "seed\tknotType\tT_before\tdelta_heuristic\tdelta_random_mean"
        << "\tdelta_random_std\theuristic_hoop_score\tskipped"
        << "\tdelta_random_hoop_mean\tdelta_random_hoop_std\tdelta_biggest_hoop"
        << "\tbiggest_hoop_span\tbiggest_hoop_score\thoop_count\n";
    return true;
}

static bool prepareCutsToZeroEvalOutputFile()
{
    if (g_evalCutsToZeroOutputPath.empty())
        g_evalCutsToZeroOutputPath = resolveDefaultCutsToZeroEvalOutputPath();

    std::ofstream out(g_evalCutsToZeroOutputPath, std::ios::trunc);
    if (!out)
    {
        std::cout << "[cuts-to-zero] could not open output file "
                  << g_evalCutsToZeroOutputPath << std::endl;
        return false;
    }

    out << "seed\tknotType\tT_before\theuristic_cuts\theuristic_resolved"
        << "\trandom_cuts_mean\trandom_cuts_std\trandom_resolved_rate"
        << "\trandom_hoop_cuts_mean\trandom_hoop_cuts_std\trandom_hoop_resolved_rate"
        << "\tbiggest_hoop_cuts\tbiggest_hoop_resolved\thoop_count\n";
    return true;
}

static bool loadKnotSeedPresetFile()
{
    g_knotSeedPresets.clear();
    if (g_knotSeedBankPath.empty())
        return false;

    std::ifstream in(g_knotSeedBankPath);
    if (!in)
    {
        std::cout << "[knot-seeds] could not open " << g_knotSeedBankPath
                  << ", using current random-cable settings." << std::endl;
        return false;
    }

    std::string line;
    while (std::getline(in, line))
    {
        const size_t hashPos = line.find('#');
        if (hashPos != std::string::npos)
            line = line.substr(0, hashPos);

        std::istringstream iss(line);
        std::string knotType;
        unsigned int seed = 0;
        if (!(iss >> knotType >> seed))
            continue;

        g_knotSeedPresets[knotType].push_back(seed);
    }

    size_t presetCount = 0;
    for (const auto& entry : g_knotSeedPresets)
        presetCount += entry.second.size();

    std::cout << "[knot-seeds] loaded " << presetCount
              << " preset(s) from " << g_knotSeedBankPath << std::endl;
    return presetCount > 0;
}

static bool storeKnotSeedPresetEntry(
    const std::string& knotType,
    unsigned int seed,
    bool appendToFile)
{
    if (knotType.empty())
        return false;

    std::vector<unsigned int>& seeds = g_knotSeedPresets[knotType];
    if (std::find(seeds.begin(), seeds.end(), seed) != seeds.end())
        return false;

    seeds.push_back(seed);
    if (appendToFile)
    {
        if (g_knotSeedBankPath.empty())
            return true;

        std::ofstream out(g_knotSeedBankPath, std::ios::app);
        if (out)
            out << knotType << ' ' << seed << '\n';
    }

    return true;
}

static bool isNontrivialKnotType(const std::string& knotType)
{
    return !knotType.empty() &&
        knotType != "UN" &&
        knotType != "unknot" &&
        knotType != "UNKNOT" &&
        knotType != "0_1" &&
        knotType != "none";
}

static void restoreDemoWindow()
{
    GLFWwindow* window = MiniGL::getWindow();
    if (window == nullptr)
        return;

    glfwShowWindow(window);
    if (glfwGetWindowAttrib(window, GLFW_ICONIFIED))
        glfwRestoreWindow(window);
    glfwSetWindowSize(window, 1280, 1024);
    glfwSetWindowPos(window, 160, 120);
    glfwFocusWindow(window);
}

static bool loadKnotSeedPresetFromFile(const std::string& selector)
{
    std::string knotType = selector;
    size_t occurrenceIndex = 0;
    const size_t colonPos = selector.find(':');
    if (colonPos != std::string::npos)
    {
        knotType = selector.substr(0, colonPos);
        try
        {
            const int parsedIndex = std::stoi(selector.substr(colonPos + 1));
            occurrenceIndex = (parsedIndex > 0) ? static_cast<size_t>(parsedIndex - 1) : 0;
        }
        catch (...)
        {
            occurrenceIndex = 0;
        }
    }

    const auto found = g_knotSeedPresets.find(knotType);
    if (found == g_knotSeedPresets.end() || found->second.empty())
    {
        std::cout << "[knot-seeds] no preset found for '" << selector << "'" << std::endl;
        return false;
    }

    const std::vector<unsigned int>& seeds = found->second;
    if (occurrenceIndex >= seeds.size())
    {
        std::cout << "[knot-seeds] '" << knotType << "' only has "
                  << seeds.size() << " seed(s)." << std::endl;
        return false;
    }

    g_useKnotSeedPreset = true;
    g_selectedKnotType = knotType;
    g_selectedKnotSeed = seeds[occurrenceIndex];
    g_cfg.randomCable.enabled = true;
    g_cfg.randomCable.seed = g_selectedKnotSeed;
    std::cout << "[knot-seeds] selected " << g_selectedKnotType
              << " seed=" << g_selectedKnotSeed
              << " (index " << (occurrenceIndex + 1) << ')'
              << std::endl;
    return true;
}

static bool selectKnotSeedPresetByType(const std::string& knotType)
{
    const auto found = g_knotSeedPresets.find(knotType);
    if (found == g_knotSeedPresets.end() || found->second.empty())
    {
        std::cout << "[knot-seeds] no preset found for knot type '" << knotType << "'" << std::endl;
        return false;
    }

    return loadKnotSeedPresetFromFile(knotType);
}

static std::vector<std::string> detectKnotTypesForPieces(
    const std::vector<std::vector<Vector3r>>& ropePieces,
    const std::vector<SelfCollisionContact>& contacts)
{
    const TanglementAnalyzer::TanglementAnalysisResult tanglement =
        TanglementAnalyzer::analyzeTanglement(
            ropePieces,
            contacts,
            g_cfg.rope,
            g_cfg.selfCollision,
            g_kymoKnotSearchMode);

    std::vector<std::string> knotTypes;
    for (const PBD::KymoKnotAnalysisResult& knotResult : tanglement.kymoResults)
    {
        if (!knotResult.available || !knotResult.success || knotResult.segments.empty())
            continue;

        const std::string& knotType = knotResult.segments.front().knotType;
        if (isNontrivialKnotType(knotType) &&
            std::find(knotTypes.begin(), knotTypes.end(), knotType) == knotTypes.end())
        {
            knotTypes.push_back(knotType);
        }
    }

    return knotTypes;
}

static std::vector<std::string> detectCurrentKnotTypes(SimulationModel* model)
{
    if (model == nullptr)
        return {};

    const std::vector<std::vector<Vector3r>> ropePieces = collectCurrentCablePieces(model);
    const std::vector<SelfCollisionContact> contacts =
        detectSelfCollisionContacts(
            model,
            g_cablePieceRanges,
            g_cfg.rope,
            g_cfg.selfCollision);
    return detectKnotTypesForPieces(ropePieces, contacts);
}

static void saveCurrentKnotSeedPreset()
{
    SimulationModel* model = Simulation::getCurrent()->getModel();
    if (model == nullptr)
    {
        std::cout << "[knot-seeds] no simulation model available." << std::endl;
        return;
    }

    const std::vector<std::string> knotTypes = detectCurrentKnotTypes(model);
    if (knotTypes.empty())
    {
        std::cout << "[knot-seeds] no knot type detected for the current cable." << std::endl;
        return;
    }

    if (g_knotSeedBankPath.empty())
    {
        std::cout << "[knot-seeds] no seed bank file path configured." << std::endl;
        return;
    }

    bool wroteAny = false;
    for (const std::string& knotType : knotTypes)
    {
        if (storeKnotSeedPresetEntry(knotType, g_cfg.randomCable.seed, true))
        {
            wroteAny = true;
            std::cout << "[knot-seeds] saved " << knotType
                      << " seed=" << g_cfg.randomCable.seed
                      << " to " << g_knotSeedBankPath << std::endl;
        }
    }

    if (!wroteAny)
        std::cout << "[knot-seeds] the current seed was already catalogued." << std::endl;
}

static size_t countCurvatureSamples(const std::vector<Vector3r>& points)
{
    return points.size() < 3 ? 0 : points.size() - 2;
}

static Real computeTotalCurvature(const std::vector<Vector3r>& points)
{
    if (points.size() < 3)
        return static_cast<Real>(0.0);

    Real total = static_cast<Real>(0.0);
    for (size_t i = 1; i + 1 < points.size(); ++i)
    {
        const Vector3r a = points[i] - points[i - 1];
        const Vector3r b = points[i + 1] - points[i];
        const Real aNorm = a.norm();
        const Real bNorm = b.norm();
        if (aNorm <= static_cast<Real>(1e-8) || bNorm <= static_cast<Real>(1e-8))
            continue;

        const Real c = std::max(
            static_cast<Real>(-1.0),
            std::min(static_cast<Real>(1.0), a.dot(b) / (aNorm * bNorm)));
        total += std::acos(c);
    }

    return total;
}

static Real computeCurvatureOutsideKnotRegions(
    const std::vector<Vector3r>& points,
    const PBD::KymoKnotAnalysisResult& knotResult,
    size_t* sampleCount)
{
    if (sampleCount)
        *sampleCount = 0;

    if (points.size() < 3)
        return static_cast<Real>(0.0);

    std::vector<bool> insideKnot(points.size(), false);
    for (const PBD::KymoKnotSegment& segment : knotResult.segments)
    {
        const int start = std::max(segment.start, 0);
        const int end = std::min(segment.end, static_cast<int>(points.size()) - 1);
        for (int i = start; i <= end; ++i)
            insideKnot[static_cast<size_t>(i)] = true;
    }

    Real total = static_cast<Real>(0.0);
    for (size_t i = 1; i + 1 < points.size(); ++i)
    {
        if (insideKnot[i])
            continue;

        if (sampleCount)
            ++(*sampleCount);

        const Vector3r a = points[i] - points[i - 1];
        const Vector3r b = points[i + 1] - points[i];
        const Real aNorm = a.norm();
        const Real bNorm = b.norm();
        if (aNorm <= static_cast<Real>(1e-8) || bNorm <= static_cast<Real>(1e-8))
            continue;

        const Real c = std::max(
            static_cast<Real>(-1.0),
            std::min(static_cast<Real>(1.0), a.dot(b) / (aNorm * bNorm)));
        total += std::acos(c);
    }

    return total;
}

static void appendCatalogRows(
    const std::vector<std::vector<Vector3r>>& ropePieces,
    const std::vector<SelfCollisionContact>& contacts,
    const TanglementAnalyzer::TanglementAnalysisResult& tanglement)
{
    (void)contacts;

    if (g_knotSeedCatalogPath.empty())
        return;

    bool needsHeader = true;
    bool writeLegacyClusterColumn = false;
    bool writeCurvatureDensityColumns = true;
    {
        std::ifstream in(g_knotSeedCatalogPath);
        needsHeader = !in.good() || in.peek() == std::ifstream::traits_type::eof();
        if (!needsHeader)
        {
            std::string headerLine;
            std::getline(in, headerLine);
            writeLegacyClusterColumn =
                headerLine.find("\tclusterPenalty\t") != std::string::npos;
            writeCurvatureDensityColumns =
                headerLine.find("\ttotalCurvatureDensity\t") != std::string::npos;
        }
    }

    std::ofstream out(g_knotSeedCatalogPath, std::ios::app);
    if (!out)
    {
        std::cout << "[knot-catalog] could not open " << g_knotSeedCatalogPath << std::endl;
        return;
    }

    if (needsHeader)
    {
        out << "knotType\tseed\tpieceIndex\tknotCount\textraKnotCount\tpieceCount\tpointCount\tsegmentCount"
            << "\ttotalHoopCount\tthreadedHoopCount\ttotalThreadingCount\tmaxHoopScore\taverageHoopScore"
            << "\tselfCollisionCount\tclosePairCount\tsamePieceClosePairCount\tcrossPieceClosePairCount"
            << "\tcontactParticipationCount\tminPairDistance\tknotSpanFraction\tproximityScore"
            << "\tcrossPieceProximityScore\tsamePieceProximityScore"
            << "\ttotalCurvatureDensity\toutsideKnotCurvatureDensity\tentanglementScore\n";
    }

    const TanglementAnalyzer::ThreadingMetricsSnapshot threading =
        TanglementAnalyzer::computeThreadingMetrics(tanglement);
    const TanglementAnalyzer::TanglementMetricsSnapshot& metrics = tanglement.metrics;
    const size_t totalHoopCount = tanglement.allCrossingEvents.size();

    Real totalCurvature = static_cast<Real>(0.0);
    size_t totalCurvatureSampleCount = 0;
    for (const std::vector<Vector3r>& piece : ropePieces)
    {
        totalCurvature += computeTotalCurvature(piece);
        totalCurvatureSampleCount += countCurvatureSamples(piece);
    }
    const Real totalCurvatureDensity = totalCurvatureSampleCount == 0
        ? static_cast<Real>(0.0)
        : totalCurvature / static_cast<Real>(totalCurvatureSampleCount);

    std::vector<std::pair<std::string, size_t>> records;
    Real outsideKnotCurvature = totalCurvature;
    size_t outsideKnotCurvatureSampleCount = totalCurvatureSampleCount;
    size_t detectedKnotCount = 0;
    bool sawUnknot = false;

    for (size_t pieceIndex = 0; pieceIndex < tanglement.kymoResults.size(); ++pieceIndex)
    {
        const PBD::KymoKnotAnalysisResult& knotResult = tanglement.kymoResults[pieceIndex];
        if (!knotResult.available || !knotResult.success || knotResult.segments.empty())
            continue;

        if (PBD::KymoKnotAnalyzer::isUnknotResult(knotResult))
            sawUnknot = true;

        for (const PBD::KymoKnotSegment& segment : knotResult.segments)
        {
            if (!isNontrivialKnotType(segment.knotType))
                continue;

            ++detectedKnotCount;
            if (pieceIndex < ropePieces.size())
            {
                size_t outsideSamples = 0;
                outsideKnotCurvature =
                    computeCurvatureOutsideKnotRegions(ropePieces[pieceIndex], knotResult, &outsideSamples);
                outsideKnotCurvatureSampleCount = outsideSamples;
            }

            const std::string knotType = segment.knotType;
            records.emplace_back(knotType, pieceIndex);
        }
    }

    const Real outsideKnotCurvatureDensity = outsideKnotCurvatureSampleCount == 0
        ? static_cast<Real>(0.0)
        : outsideKnotCurvature / static_cast<Real>(outsideKnotCurvatureSampleCount);

    if (records.empty())
        records.emplace_back(sawUnknot ? "UN" : "none", std::numeric_limits<size_t>::max());

    const size_t extraKnotCount = detectedKnotCount > 0 ? detectedKnotCount - 1 : 0;
    const Real finiteMinDistance =
        std::isfinite(metrics.minimumPairDistance) ? metrics.minimumPairDistance : static_cast<Real>(-1.0);

    // Catalog-level entanglement score: each term is min-max normalized using
    // catalog-derived constants so raw scale differences do not dominate.
    // Normalization ranges are fixed from knot_seed_catalog_v3.tsv (403 entries).
    // Values outside the catalog range are clamped to [0, 1].
    auto normClamp = [](Real v, Real mn, Real mx) -> Real {
        if (mx <= mn) return static_cast<Real>(0);
        return std::max(static_cast<Real>(0), std::min(static_cast<Real>(1), (v - mn) / (mx - mn)));
    };
    const Real nThreading   = normClamp(static_cast<Real>(threading.totalThreadingCount),           static_cast<Real>(0),     static_cast<Real>(110));
    const Real nClosePair   = normClamp(static_cast<Real>(metrics.closePairCount),                  static_cast<Real>(914),   static_cast<Real>(3186));
    const Real nContactPart = normClamp(static_cast<Real>(metrics.contactParticipationCount),       static_cast<Real>(4),     static_cast<Real>(137));
    const Real nKnotSpan    = normClamp(metrics.knotSpanFraction,                                   static_cast<Real>(0),     static_cast<Real>(0.898));
    const Real nOutsideCurv = normClamp(outsideKnotCurvatureDensity,                                static_cast<Real>(0.0),   static_cast<Real>(0.10));
    const Real nProximity   = normClamp(metrics.proximityScore,                                     static_cast<Real>(0.016), static_cast<Real>(0.044));
    // extraKnotCount kept unnormalized: it is 0 for all catalog entries but carries a
    // large penalty (5.0) when multiple simultaneous knots are detected.
    const Real catalogEntanglementScore =
        static_cast<Real>(3.00) * nThreading   +
        static_cast<Real>(0.24) * nClosePair   +
        static_cast<Real>(0.17) * nContactPart +
        static_cast<Real>(0.56) * nKnotSpan    +
        static_cast<Real>(0.12) * nOutsideCurv +
        static_cast<Real>(5.00) * static_cast<Real>(extraKnotCount) +
        static_cast<Real>(0.42) * nProximity;

    out << std::fixed << std::setprecision(6);
    for (const auto& record : records)
    {
        out << record.first << '\t'
            << g_cfg.randomCable.seed << '\t'
            << (record.second == std::numeric_limits<size_t>::max() ? -1 : static_cast<int>(record.second)) << '\t'
            << detectedKnotCount << '\t'
            << extraKnotCount << '\t'
            << metrics.pieceCount << '\t'
            << metrics.totalPointCount << '\t'
            << metrics.totalSegmentCount << '\t'
            << totalHoopCount << '\t'
            << threading.threadedHoopCount << '\t'
            << threading.totalThreadingCount << '\t'
            << threading.maxHoopScore << '\t'
            << threading.averageHoopScore << '\t'
            << metrics.selfCollisionCount << '\t'
            << metrics.closePairCount << '\t'
            << metrics.samePieceClosePairCount << '\t'
            << metrics.crossPieceClosePairCount << '\t'
            << metrics.contactParticipationCount << '\t'
            << finiteMinDistance << '\t'
            << metrics.knotSpanFraction << '\t'
            << metrics.proximityScore << '\t'
            << metrics.crossPieceProximityScore << '\t'
            << metrics.samePieceProximityScore << '\t';
        if (writeLegacyClusterColumn)
            out << static_cast<Real>(0.0) << '\t';
        out << (writeCurvatureDensityColumns ? totalCurvatureDensity : totalCurvature) << '\t'
            << (writeCurvatureDensityColumns ? outsideKnotCurvatureDensity : outsideKnotCurvature) << '\t'
            << catalogEntanglementScore << '\n';
    }

    std::cout << "[knot-catalog] wrote " << records.size()
              << " row(s) for seed=" << g_cfg.randomCable.seed
              << " to " << g_knotSeedCatalogPath << std::endl;
}

static void evaluateCurrentCableTanglement(bool saveSeedPreset)
{
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
            g_cfg.rope,
            g_cfg.selfCollision,
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
    g_lastPrintedHoopSelection = CableVisualization::HoopMarkerSelection();

    TanglementAnalyzer::printTanglementReport(ropePieces, tanglement);
    if (saveSeedPreset)
        saveCurrentKnotSeedPreset();
}

static void evaluateAndRecordCatalogSample(bool printReport)
{
    SimulationModel* model = Simulation::getCurrent()->getModel();
    if (model == nullptr)
        return;

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
            g_cfg.rope,
            g_cfg.selfCollision,
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
    g_lastPrintedHoopSelection = CableVisualization::HoopMarkerSelection();

    if (printReport)
        TanglementAnalyzer::printTanglementReport(ropePieces, tanglement);

    std::vector<std::string> knotTypes;
    for (const PBD::KymoKnotAnalysisResult& knotResult : tanglement.kymoResults)
    {
        if (!knotResult.available || !knotResult.success || knotResult.segments.empty())
            continue;

        for (const PBD::KymoKnotSegment& segment : knotResult.segments)
        {
            if (isNontrivialKnotType(segment.knotType) &&
                std::find(knotTypes.begin(), knotTypes.end(), segment.knotType) == knotTypes.end())
            {
                knotTypes.push_back(segment.knotType);
            }
        }
    }

    for (const std::string& knotType : knotTypes)
        storeKnotSeedPresetEntry(knotType, g_cfg.randomCable.seed, true);

    appendCatalogRows(ropePieces, contacts, tanglement);
    ++g_catalogSamplesWritten;
}

namespace
{
    struct VirtualCutTarget
    {
        size_t pieceIndex = std::numeric_limits<size_t>::max();
        size_t edgeIndex = std::numeric_limits<size_t>::max();
        int hoopScore = 0;
        int hoopSpan = 0;
        bool valid = false;
    };

    static bool makeHoopMidpointCutTarget(
        size_t pieceIndex,
        const std::vector<Vector3r>& piece,
        const TanglementAnalyzer::HoopCandidate& hoop,
        VirtualCutTarget& target)
    {
        if (piece.size() < 4 || hoop.pairSumScore <= 0)
            return false;

        const size_t maxEdge = piece.size() - 2;
        const int rawMidpoint = (hoop.startEdge + hoop.endEdge) / 2;
        const size_t edgeIndex = static_cast<size_t>(
            std::max(1, std::min(rawMidpoint, static_cast<int>(maxEdge - 1))));
        if (edgeIndex + 2 >= piece.size())
            return false;

        target.pieceIndex = pieceIndex;
        target.edgeIndex = edgeIndex;
        target.hoopScore = hoop.pairSumScore;
        target.hoopSpan = hoop.endEdge - hoop.startEdge;
        target.valid = true;
        return true;
    }

    static std::vector<VirtualCutTarget> collectHoopCutTargets(
        const TanglementAnalyzer::TanglementAnalysisResult& tanglement,
        const std::vector<std::vector<Vector3r>>& ropePieces)
    {
        std::vector<VirtualCutTarget> targets;
        for (size_t pieceIndex = 0; pieceIndex < tanglement.scoredHoops.size(); ++pieceIndex)
        {
            if (pieceIndex >= ropePieces.size())
                continue;

            const std::vector<Vector3r>& piece = ropePieces[pieceIndex];
            for (const TanglementAnalyzer::HoopCandidate& hoop : tanglement.scoredHoops[pieceIndex])
            {
                VirtualCutTarget target;
                if (makeHoopMidpointCutTarget(pieceIndex, piece, hoop, target))
                    targets.push_back(target);
            }
        }
        return targets;
    }

    static VirtualCutTarget findBestHoopMidpointCut(
        const TanglementAnalyzer::TanglementAnalysisResult& tanglement,
        const std::vector<std::vector<Vector3r>>& ropePieces)
    {
        VirtualCutTarget best;
        for (const VirtualCutTarget& candidate : collectHoopCutTargets(tanglement, ropePieces))
        {
            if (!best.valid || candidate.hoopScore > best.hoopScore)
                best = candidate;
        }
        return best;
    }

    static VirtualCutTarget findBiggestHoopMidpointCut(
        const TanglementAnalyzer::TanglementAnalysisResult& tanglement,
        const std::vector<std::vector<Vector3r>>& ropePieces)
    {
        VirtualCutTarget best;
        for (const VirtualCutTarget& candidate : collectHoopCutTargets(tanglement, ropePieces))
        {
            if (!best.valid ||
                candidate.hoopSpan > best.hoopSpan ||
                (candidate.hoopSpan == best.hoopSpan && candidate.hoopScore > best.hoopScore))
            {
                best = candidate;
            }
        }
        return best;
    }

    static bool selectRandomHoopCutTarget(
        const TanglementAnalyzer::TanglementAnalysisResult& tanglement,
        const std::vector<std::vector<Vector3r>>& ropePieces,
        std::mt19937& rng,
        VirtualCutTarget& target)
    {
        const std::vector<VirtualCutTarget> targets =
            collectHoopCutTargets(tanglement, ropePieces);
        if (targets.empty())
            return false;

        std::uniform_int_distribution<size_t> dist(0, targets.size() - 1);
        target = targets[dist(rng)];
        return target.valid;
    }

    static std::vector<VirtualCutTarget> collectValidRandomCutTargets(
        const std::vector<std::vector<Vector3r>>& ropePieces,
        Real excludeArcLength)
    {
        std::vector<VirtualCutTarget> targets;
        for (size_t pieceIndex = 0; pieceIndex < ropePieces.size(); ++pieceIndex)
        {
            const std::vector<Vector3r>& piece = ropePieces[pieceIndex];
            if (piece.size() < 4)
                continue;

            std::vector<Real> cumulative(piece.size(), static_cast<Real>(0.0));
            for (size_t i = 1; i < piece.size(); ++i)
                cumulative[i] = cumulative[i - 1] + (piece[i] - piece[i - 1]).norm();

            const Real totalLength = cumulative.back();
            for (size_t edgeIndex = 1; edgeIndex + 2 < piece.size(); ++edgeIndex)
            {
                const Real edgeMidArc =
                    static_cast<Real>(0.5) * (cumulative[edgeIndex] + cumulative[edgeIndex + 1]);
                if (edgeMidArc < excludeArcLength ||
                    totalLength - edgeMidArc < excludeArcLength)
                {
                    continue;
                }

                VirtualCutTarget target;
                target.pieceIndex = pieceIndex;
                target.edgeIndex = edgeIndex;
                target.valid = true;
                targets.push_back(target);
            }
        }
        return targets;
    }

    static int computeThreadingCountGeometryOnly(
        const std::vector<std::vector<Vector3r>>& ropePieces)
    {
        const TanglementAnalyzer::TanglementAnalysisResult tanglement =
            TanglementAnalyzer::analyzeTanglementGeometryOnly(
                ropePieces,
                {},
                g_cfg.rope,
                g_cfg.selfCollision);
        return TanglementAnalyzer::computeThreadingMetrics(tanglement).totalThreadingCount;
    }

    static int computeVirtualCutDeltaThreading(
        const std::vector<std::vector<Vector3r>>& ropePieces,
        int beforeThreading,
        const VirtualCutTarget& target)
    {
        if (!target.valid)
            return 0;

        std::vector<std::vector<Vector3r>> splitPieces = ropePieces;
        if (!splitCablePieceAtEdge(splitPieces, target.pieceIndex, target.edgeIndex))
            return 0;

        const int afterThreading = computeThreadingCountGeometryOnly(splitPieces);
        return beforeThreading - afterThreading;
    }

    struct CutsToZeroOutcome
    {
        unsigned int cuts = 0;
        bool resolved = false;
    };

    static CutsToZeroOutcome runHeuristicCutsToZero(
        const std::vector<std::vector<Vector3r>>& startPieces,
        unsigned int maxCuts)
    {
        std::vector<std::vector<Vector3r>> pieces = startPieces;
        CutsToZeroOutcome outcome;

        while (outcome.cuts < maxCuts)
        {
            const TanglementAnalyzer::TanglementAnalysisResult tanglement =
                TanglementAnalyzer::analyzeTanglementGeometryOnly(
                    pieces,
                    {},
                    g_cfg.rope,
                    g_cfg.selfCollision);
            const int threading =
                TanglementAnalyzer::computeThreadingMetrics(tanglement).totalThreadingCount;
            if (threading == 0)
            {
                outcome.resolved = true;
                return outcome;
            }

            const VirtualCutTarget target = findBestHoopMidpointCut(tanglement, pieces);
            if (!target.valid || !splitCablePieceAtEdge(pieces, target.pieceIndex, target.edgeIndex))
                break;

            ++outcome.cuts;
        }

        outcome.resolved = computeThreadingCountGeometryOnly(pieces) == 0;
        if (!outcome.resolved)
            outcome.cuts = maxCuts;
        return outcome;
    }

    static CutsToZeroOutcome runBiggestHoopCutsToZero(
        const std::vector<std::vector<Vector3r>>& startPieces,
        unsigned int maxCuts)
    {
        std::vector<std::vector<Vector3r>> pieces = startPieces;
        CutsToZeroOutcome outcome;

        while (outcome.cuts < maxCuts)
        {
            const TanglementAnalyzer::TanglementAnalysisResult tanglement =
                TanglementAnalyzer::analyzeTanglementGeometryOnly(
                    pieces,
                    {},
                    g_cfg.rope,
                    g_cfg.selfCollision);
            const int threading =
                TanglementAnalyzer::computeThreadingMetrics(tanglement).totalThreadingCount;
            if (threading == 0)
            {
                outcome.resolved = true;
                return outcome;
            }

            const VirtualCutTarget target = findBiggestHoopMidpointCut(tanglement, pieces);
            if (!target.valid || !splitCablePieceAtEdge(pieces, target.pieceIndex, target.edgeIndex))
                break;

            ++outcome.cuts;
        }

        outcome.resolved = computeThreadingCountGeometryOnly(pieces) == 0;
        if (!outcome.resolved)
            outcome.cuts = maxCuts;
        return outcome;
    }

    static CutsToZeroOutcome runRandomHoopCutsToZero(
        const std::vector<std::vector<Vector3r>>& startPieces,
        unsigned int maxCuts,
        std::mt19937& rng)
    {
        std::vector<std::vector<Vector3r>> pieces = startPieces;
        CutsToZeroOutcome outcome;

        while (outcome.cuts < maxCuts)
        {
            const TanglementAnalyzer::TanglementAnalysisResult tanglement =
                TanglementAnalyzer::analyzeTanglementGeometryOnly(
                    pieces,
                    {},
                    g_cfg.rope,
                    g_cfg.selfCollision);
            const int threading =
                TanglementAnalyzer::computeThreadingMetrics(tanglement).totalThreadingCount;
            if (threading == 0)
            {
                outcome.resolved = true;
                return outcome;
            }

            VirtualCutTarget target;
            if (!selectRandomHoopCutTarget(tanglement, pieces, rng, target) ||
                !splitCablePieceAtEdge(pieces, target.pieceIndex, target.edgeIndex))
            {
                break;
            }

            ++outcome.cuts;
        }

        outcome.resolved = computeThreadingCountGeometryOnly(pieces) == 0;
        if (!outcome.resolved)
            outcome.cuts = maxCuts;
        return outcome;
    }

    static CutsToZeroOutcome runRandomCutsToZero(
        const std::vector<std::vector<Vector3r>>& startPieces,
        unsigned int maxCuts,
        std::mt19937& rng)
    {
        std::vector<std::vector<Vector3r>> pieces = startPieces;
        CutsToZeroOutcome outcome;

        while (outcome.cuts < maxCuts)
        {
            const int threading = computeThreadingCountGeometryOnly(pieces);
            if (threading == 0)
            {
                outcome.resolved = true;
                return outcome;
            }

            const std::vector<VirtualCutTarget> targets =
                collectValidRandomCutTargets(pieces, g_cfg.selfCollision.excludeArcLength);
            if (targets.empty())
                break;

            std::uniform_int_distribution<size_t> dist(0, targets.size() - 1);
            const VirtualCutTarget& target = targets[dist(rng)];
            if (!splitCablePieceAtEdge(pieces, target.pieceIndex, target.edgeIndex))
                break;

            ++outcome.cuts;
        }

        outcome.resolved = computeThreadingCountGeometryOnly(pieces) == 0;
        if (!outcome.resolved)
            outcome.cuts = maxCuts;
        return outcome;
    }
}

static void evaluateAndRecordCutHeuristicSample()
{
    if (g_evalCurrentIndex >= g_evalSeeds.size())
        return;

    SimulationModel* model = Simulation::getCurrent()->getModel();
    if (model == nullptr)
        return;

    const CutEvalSeedRecord& record = g_evalSeeds[g_evalCurrentIndex];
    const std::vector<std::vector<Vector3r>> ropePieces = collectCurrentCablePieces(model);
    const TanglementAnalyzer::TanglementAnalysisResult tanglement =
        TanglementAnalyzer::analyzeTanglementGeometryOnly(
            ropePieces,
            {},
            g_cfg.rope,
            g_cfg.selfCollision);
    const TanglementAnalyzer::ThreadingMetricsSnapshot threading =
        TanglementAnalyzer::computeThreadingMetrics(tanglement);
    const int beforeThreading = threading.totalThreadingCount;
    const std::vector<VirtualCutTarget> hoopTargets =
        collectHoopCutTargets(tanglement, ropePieces);
    const int hoopCount = static_cast<int>(hoopTargets.size());

    int skipped = 0;
    int heuristicDelta = 0;
    int heuristicHoopScore = 0;
    if (beforeThreading > 0)
    {
        const VirtualCutTarget heuristicTarget = findBestHoopMidpointCut(tanglement, ropePieces);
        if (heuristicTarget.valid)
        {
            heuristicHoopScore = heuristicTarget.hoopScore;
            heuristicDelta = computeVirtualCutDeltaThreading(
                ropePieces,
                beforeThreading,
                heuristicTarget);
        }
        else
        {
            skipped = 1;
        }
    }

    int biggestHoopDelta = 0;
    int biggestHoopSpan = 0;
    int biggestHoopScore = 0;
    if (beforeThreading > 0)
    {
        const VirtualCutTarget biggestTarget =
            findBiggestHoopMidpointCut(tanglement, ropePieces);
        if (biggestTarget.valid)
        {
            biggestHoopDelta = computeVirtualCutDeltaThreading(
                ropePieces,
                beforeThreading,
                biggestTarget);
            biggestHoopSpan = biggestTarget.hoopSpan;
            biggestHoopScore = biggestTarget.hoopScore;
        }
    }

    Real randomMean = static_cast<Real>(0.0);
    Real randomStd = static_cast<Real>(0.0);
    if (beforeThreading > 0)
    {
        const std::vector<VirtualCutTarget> randomTargets =
            collectValidRandomCutTargets(ropePieces, g_cfg.selfCollision.excludeArcLength);
        if (!randomTargets.empty() && g_evalRandomTrials > 0)
        {
            std::mt19937 rng(record.seed ^ 0x9E3779B9u);
            std::uniform_int_distribution<size_t> dist(0, randomTargets.size() - 1);
            std::vector<Real> randomDeltas;
            randomDeltas.reserve(g_evalRandomTrials);

            for (unsigned int trial = 0; trial < g_evalRandomTrials; ++trial)
            {
                const VirtualCutTarget& randomTarget = randomTargets[dist(rng)];
                const int delta = computeVirtualCutDeltaThreading(
                    ropePieces,
                    beforeThreading,
                    randomTarget);
                randomDeltas.push_back(static_cast<Real>(delta));
                randomMean += static_cast<Real>(delta);
            }

            randomMean /= static_cast<Real>(randomDeltas.size());
            Real variance = static_cast<Real>(0.0);
            for (Real delta : randomDeltas)
            {
                const Real d = delta - randomMean;
                variance += d * d;
            }
            variance /= static_cast<Real>(randomDeltas.size());
            randomStd = std::sqrt(variance);
        }
    }

    Real randomHoopMean = static_cast<Real>(0.0);
    Real randomHoopStd = static_cast<Real>(0.0);
    if (beforeThreading > 0 && !hoopTargets.empty() && g_evalRandomTrials > 0)
    {
        std::mt19937 rng(record.seed ^ 0xA5A5A5A5u);
        std::uniform_int_distribution<size_t> dist(0, hoopTargets.size() - 1);
        std::vector<Real> randomHoopDeltas;
        randomHoopDeltas.reserve(g_evalRandomTrials);

        for (unsigned int trial = 0; trial < g_evalRandomTrials; ++trial)
        {
            const VirtualCutTarget& randomHoopTarget = hoopTargets[dist(rng)];
            const int delta = computeVirtualCutDeltaThreading(
                ropePieces,
                beforeThreading,
                randomHoopTarget);
            randomHoopDeltas.push_back(static_cast<Real>(delta));
            randomHoopMean += static_cast<Real>(delta);
        }

        randomHoopMean /= static_cast<Real>(randomHoopDeltas.size());
        Real variance = static_cast<Real>(0.0);
        for (Real delta : randomHoopDeltas)
        {
            const Real d = delta - randomHoopMean;
            variance += d * d;
        }
        variance /= static_cast<Real>(randomHoopDeltas.size());
        randomHoopStd = std::sqrt(variance);
    }

    std::ofstream out(g_evalOutputPath, std::ios::app);
    if (!out)
    {
        std::cout << "[cut-eval] could not append to " << g_evalOutputPath << std::endl;
        g_evalRunEnabled = false;
        return;
    }

    out << std::fixed << std::setprecision(6)
        << record.seed << '\t'
        << record.knotType << '\t'
        << beforeThreading << '\t'
        << heuristicDelta << '\t'
        << randomMean << '\t'
        << randomStd << '\t'
        << heuristicHoopScore << '\t'
        << skipped << '\t'
        << randomHoopMean << '\t'
        << randomHoopStd << '\t'
        << biggestHoopDelta << '\t'
        << biggestHoopSpan << '\t'
        << biggestHoopScore << '\t'
        << hoopCount << '\n';

    ++g_evalSamplesWritten;
    std::cout << "[cut-eval] seed=" << record.seed
              << " knotType=" << record.knotType
              << " T_before=" << beforeThreading
              << " heuristicDelta=" << heuristicDelta
              << " randomMean=" << randomMean
              << " randomHoopMean=" << randomHoopMean
              << " biggestHoopDelta=" << biggestHoopDelta
              << " hoops=" << hoopCount
              << " (" << (g_evalCurrentIndex + 1) << "/" << g_evalSeeds.size() << ")"
              << std::endl;
}

static void evaluateAndRecordCutsToZeroSample()
{
    if (g_evalCurrentIndex >= g_evalSeeds.size())
        return;

    SimulationModel* model = Simulation::getCurrent()->getModel();
    if (model == nullptr)
        return;

    const CutEvalSeedRecord& record = g_evalSeeds[g_evalCurrentIndex];
    const std::vector<std::vector<Vector3r>> ropePieces = collectCurrentCablePieces(model);
    const TanglementAnalyzer::TanglementAnalysisResult initialTanglement =
        TanglementAnalyzer::analyzeTanglementGeometryOnly(
            ropePieces,
            {},
            g_cfg.rope,
            g_cfg.selfCollision);
    const int beforeThreading =
        TanglementAnalyzer::computeThreadingMetrics(initialTanglement).totalThreadingCount;
    const int hoopCount =
        static_cast<int>(collectHoopCutTargets(initialTanglement, ropePieces).size());

    CutsToZeroOutcome heuristicOutcome;
    CutsToZeroOutcome biggestHoopOutcome;
    Real randomMean = static_cast<Real>(0.0);
    Real randomStd = static_cast<Real>(0.0);
    Real randomResolvedRate = static_cast<Real>(1.0);
    Real randomHoopMean = static_cast<Real>(0.0);
    Real randomHoopStd = static_cast<Real>(0.0);
    Real randomHoopResolvedRate = static_cast<Real>(1.0);

    if (beforeThreading == 0)
    {
        heuristicOutcome.resolved = true;
        biggestHoopOutcome.resolved = true;
    }
    else
    {
        heuristicOutcome = runHeuristicCutsToZero(
            ropePieces,
            g_evalCutsToZeroMaxCuts);
        biggestHoopOutcome = runBiggestHoopCutsToZero(
            ropePieces,
            g_evalCutsToZeroMaxCuts);

        randomResolvedRate = static_cast<Real>(0.0);
        std::mt19937 edgeRng(record.seed ^ 0x9E3779B9u);
        std::vector<Real> randomCuts;
        randomCuts.reserve(g_evalCutsToZeroRandomTrials);

        for (unsigned int trial = 0; trial < g_evalCutsToZeroRandomTrials; ++trial)
        {
            const CutsToZeroOutcome randomOutcome =
                runRandomCutsToZero(
                    ropePieces,
                    g_evalCutsToZeroMaxCuts,
                    edgeRng);
            randomCuts.push_back(static_cast<Real>(randomOutcome.cuts));
            randomMean += static_cast<Real>(randomOutcome.cuts);
            if (randomOutcome.resolved)
                randomResolvedRate += static_cast<Real>(1.0);
        }

        if (!randomCuts.empty())
        {
            randomMean /= static_cast<Real>(randomCuts.size());
            randomResolvedRate /= static_cast<Real>(randomCuts.size());
            Real variance = static_cast<Real>(0.0);
            for (Real cuts : randomCuts)
            {
                const Real d = cuts - randomMean;
                variance += d * d;
            }
            variance /= static_cast<Real>(randomCuts.size());
            randomStd = std::sqrt(variance);
        }

        randomHoopResolvedRate = static_cast<Real>(0.0);
        std::mt19937 hoopRng(record.seed ^ 0xA5A5A5A5u);
        std::vector<Real> randomHoopCuts;
        randomHoopCuts.reserve(g_evalCutsToZeroRandomTrials);

        for (unsigned int trial = 0; trial < g_evalCutsToZeroRandomTrials; ++trial)
        {
            const CutsToZeroOutcome randomHoopOutcome =
                runRandomHoopCutsToZero(
                    ropePieces,
                    g_evalCutsToZeroMaxCuts,
                    hoopRng);
            randomHoopCuts.push_back(static_cast<Real>(randomHoopOutcome.cuts));
            randomHoopMean += static_cast<Real>(randomHoopOutcome.cuts);
            if (randomHoopOutcome.resolved)
                randomHoopResolvedRate += static_cast<Real>(1.0);
        }

        if (!randomHoopCuts.empty())
        {
            randomHoopMean /= static_cast<Real>(randomHoopCuts.size());
            randomHoopResolvedRate /= static_cast<Real>(randomHoopCuts.size());
            Real variance = static_cast<Real>(0.0);
            for (Real cuts : randomHoopCuts)
            {
                const Real d = cuts - randomHoopMean;
                variance += d * d;
            }
            variance /= static_cast<Real>(randomHoopCuts.size());
            randomHoopStd = std::sqrt(variance);
        }
    }

    std::ofstream out(g_evalCutsToZeroOutputPath, std::ios::app);
    if (!out)
    {
        std::cout << "[cuts-to-zero] could not append to "
                  << g_evalCutsToZeroOutputPath << std::endl;
        g_evalCutsToZeroEnabled = false;
        return;
    }

    out << std::fixed << std::setprecision(6)
        << record.seed << '\t'
        << record.knotType << '\t'
        << beforeThreading << '\t'
        << heuristicOutcome.cuts << '\t'
        << (heuristicOutcome.resolved ? 1 : 0) << '\t'
        << randomMean << '\t'
        << randomStd << '\t'
        << randomResolvedRate << '\t'
        << randomHoopMean << '\t'
        << randomHoopStd << '\t'
        << randomHoopResolvedRate << '\t'
        << biggestHoopOutcome.cuts << '\t'
        << (biggestHoopOutcome.resolved ? 1 : 0) << '\t'
        << hoopCount << '\n';

    ++g_evalSamplesWritten;
    std::cout << "[cuts-to-zero] seed=" << record.seed
              << " knotType=" << record.knotType
              << " T_before=" << beforeThreading
              << " heuristicCuts=" << heuristicOutcome.cuts
              << " heuristicResolved=" << (heuristicOutcome.resolved ? 1 : 0)
              << " randomMean=" << randomMean
              << " randomResolvedRate=" << randomResolvedRate
              << " randomHoopMean=" << randomHoopMean
              << " randomHoopResolvedRate=" << randomHoopResolvedRate
              << " biggestHoopCuts=" << biggestHoopOutcome.cuts
              << " biggestHoopResolved=" << (biggestHoopOutcome.resolved ? 1 : 0)
              << " hoops=" << hoopCount
              << " (" << (g_evalCurrentIndex + 1) << "/" << g_evalSeeds.size() << ")"
              << std::endl;
}

static void prepareCatalogSeedForNextRun()
{
    unsigned int nextSeed = 0;
    if (g_catalogNextSeed != 0)
    {
        nextSeed = g_catalogNextSeed++;
        if (g_catalogNextSeed == 0)
            g_catalogNextSeed = 1;
    }
    else
    {
        do
        {
            nextSeed = g_catalogRng();
        }
        while (nextSeed == 0);
    }

    g_useDirectRandomSeed = false;
    g_useKnotSeedPreset = false;
    g_selectedKnotType.clear();
    g_cfg.randomCable.enabled = true;
    g_cfg.randomCable.seed = nextSeed;
}

static void prepareCutEvalSeed(size_t index)
{
    if (index >= g_evalSeeds.size())
        return;

    g_evalCurrentIndex = index;
    g_useDirectRandomSeed = false;
    g_useKnotSeedPreset = false;
    g_selectedKnotType.clear();
    g_cfg.randomCable.enabled = true;
    g_cfg.randomCable.seed = g_evalSeeds[g_evalCurrentIndex].seed;
}

static void generateNextCatalogSeed()
{
    prepareCatalogSeedForNextRun();
    g_cablePieces.clear();
    std::cout << "[knot-catalog] next seed=" << g_cfg.randomCable.seed << std::endl;
    reset();
}

static void generateNextCutEvalSeed()
{
    prepareCutEvalSeed(g_evalCurrentIndex);
    g_cablePieces.clear();
    std::cout << "[cut-eval] next seed=" << g_cfg.randomCable.seed << std::endl;
    reset();
}

int main(int argc, char** argv)
{
    base = new DemoBase();

    std::vector<char*> frameworkArgs;
    frameworkArgs.push_back(argv[0]);
    for (int i = 1; i < argc; ++i)
    {
        const std::string arg(argv[i]);
        if (arg == "--no-cache")
            frameworkArgs.push_back(argv[i]);
    }
    base->init(static_cast<int>(frameworkArgs.size()), frameworkArgs.data(), "Random Cable Mess Demo");
    g_knotSeedBankPath = resolveDefaultKnotSeedBankPath();
    g_knotSeedCatalogPath = resolveDefaultKnotSeedCatalogPath();
    g_evalOutputPath = resolveDefaultCutEvalOutputPath();
    g_evalCutsToZeroOutputPath = resolveDefaultCutsToZeroEvalOutputPath();
    std::cout << "[knot-seeds] bank path: " << g_knotSeedBankPath << std::endl;
    std::cout << "[knot-catalog] catalog path: " << g_knotSeedCatalogPath << std::endl;

    std::string knotSelector;
    bool hasDirectSeed = false;
    unsigned int directSeed = 0;
    std::random_device catalogSeedSource;
    g_catalogRng.seed(catalogSeedSource());

    for (int i = 1; i < argc; ++i)
    {
        const std::string arg(argv[i]);
        if (arg.rfind("--seed-bank=", 0) == 0)
        {
            g_knotSeedBankPath = arg.substr(12);
        }
        else if (arg.rfind("--catalog-file=", 0) == 0)
        {
            g_knotSeedCatalogPath = arg.substr(15);
        }
        else if (arg == "--eval-run")
        {
            g_evalRunEnabled = true;
        }
        else if (arg == "--eval-cuts-to-zero")
        {
            g_evalCutsToZeroEnabled = true;
        }
        else if (arg.rfind("--eval-output=", 0) == 0)
        {
            g_evalOutputPath = arg.substr(14);
        }
        else if (arg.rfind("--eval-ctz-output=", 0) == 0)
        {
            g_evalCutsToZeroOutputPath = arg.substr(18);
        }
        else if (arg.rfind("--eval-random-trials=", 0) == 0)
        {
            try
            {
                g_evalRandomTrials = std::max(
                    1u,
                    static_cast<unsigned int>(std::stoul(arg.substr(21))));
            }
            catch (...)
            {
                g_evalRandomTrials = 10;
            }
        }
        else if (arg.rfind("--eval-ctz-max-cuts=", 0) == 0)
        {
            try
            {
                g_evalCutsToZeroMaxCuts = std::max(
                    1u,
                    static_cast<unsigned int>(std::stoul(arg.substr(20))));
            }
            catch (...)
            {
                g_evalCutsToZeroMaxCuts = 10;
            }
        }
        else if (arg.rfind("--eval-ctz-random-trials=", 0) == 0)
        {
            try
            {
                g_evalCutsToZeroRandomTrials = std::max(
                    1u,
                    static_cast<unsigned int>(std::stoul(arg.substr(25))));
            }
            catch (...)
            {
                g_evalCutsToZeroRandomTrials = 20;
            }
        }
        else if (arg == "--catalog-run")
        {
            g_catalogRunEnabled = true;
        }
        else if (arg == "--catalog-render")
        {
            g_catalogRenderEnabled = true;
        }
        else if (arg.rfind("--catalog-steps-per-frame=", 0) == 0)
        {
            try
            {
                g_catalogStepsPerFrame =
                    std::max(1u, static_cast<unsigned int>(std::stoul(arg.substr(26))));
                g_catalogRunEnabled = true;
            }
            catch (...)
            {
                g_catalogStepsPerFrame = 1;
            }
        }
        else if (arg.rfind("--catalog-max=", 0) == 0)
        {
            try
            {
                g_catalogMaxSamples = static_cast<unsigned int>(std::stoul(arg.substr(14)));
                g_catalogRunEnabled = true;
            }
            catch (...)
            {
                g_catalogMaxSamples = 0;
            }
        }
        else if (arg.rfind("--catalog-start=", 0) == 0)
        {
            try
            {
                g_catalogNextSeed = static_cast<unsigned int>(std::stoul(arg.substr(16)));
                if (g_catalogNextSeed == 0)
                    g_catalogNextSeed = 1;
                g_catalogRunEnabled = true;
            }
            catch (...)
            {
                g_catalogNextSeed = 0;
            }
        }
        else if (arg.rfind("--seed=", 0) == 0)
        {
            try
            {
                directSeed = static_cast<unsigned int>(std::stoul(arg.substr(7)));
                hasDirectSeed = true;
            }
            catch (...)
            {
                hasDirectSeed = false;
            }
        }
        else if (arg.rfind("--random-seed=", 0) == 0)
        {
            try
            {
                directSeed = static_cast<unsigned int>(std::stoul(arg.substr(14)));
                hasDirectSeed = true;
            }
            catch (...)
            {
                hasDirectSeed = false;
            }
        }
        else if (arg.rfind("--knot-type=", 0) == 0)
        {
            knotSelector = arg.substr(12);
        }
        else if (arg.rfind("--knot=", 0) == 0)
        {
            knotSelector = arg.substr(7);
        }
        else if (arg == "--knot" && i + 1 < argc)
        {
            knotSelector = argv[++i];
        }
    }

    loadKnotSeedPresetFile();

    if (hasDirectSeed)
    {
        if (directSeed == 0)
            directSeed = 1;
        g_useDirectRandomSeed = true;
        g_directRandomSeed = directSeed;
        g_cfg.randomCable.enabled = true;
        g_cfg.randomCable.seed = g_directRandomSeed;
        g_useKnotSeedPreset = false;
        g_selectedKnotType.clear();
        std::cout << "[knot-seeds] direct seed selected: " << g_directRandomSeed << std::endl;
    }
    else if (!knotSelector.empty())
    {
        loadKnotSeedPresetFromFile(knotSelector);
    }

    if (g_evalRunEnabled)
    {
        g_catalogRunEnabled = false;
        g_catalogRenderEnabled = false;
        if (!loadCutEvalCatalog() || !prepareCutEvalOutputFile())
        {
            g_evalRunEnabled = false;
        }
        else
        {
            prepareCutEvalSeed(0);
            std::cout << "[cut-eval] running unattended cut evaluation"
                      << " seed=" << g_cfg.randomCable.seed
                      << " samples=" << g_evalSeeds.size()
                      << " randomTrials=" << g_evalRandomTrials
                      << " output=" << g_evalOutputPath << std::endl;
        }
    }

    if (g_evalCutsToZeroEnabled)
    {
        g_catalogRunEnabled = false;
        g_evalRunEnabled = false;
        g_catalogRenderEnabled = false;
        if (!loadCutEvalCatalog() || !prepareCutsToZeroEvalOutputFile())
        {
            g_evalCutsToZeroEnabled = false;
        }
        else
        {
            prepareCutEvalSeed(0);
            std::cout << "[cuts-to-zero] running unattended cuts-to-zero evaluation"
                      << " seed=" << g_cfg.randomCable.seed
                      << " samples=" << g_evalSeeds.size()
                      << " maxCuts=" << g_evalCutsToZeroMaxCuts
                      << " randomTrials=" << g_evalCutsToZeroRandomTrials
                      << " output=" << g_evalCutsToZeroOutputPath << std::endl;
        }
    }

    if (g_catalogRunEnabled)
    {
        prepareCatalogSeedForNextRun();
        std::cout << "[knot-catalog] running unattended catalogue mode"
                  << " seed=" << g_cfg.randomCable.seed;
        if (g_catalogMaxSamples > 0)
            std::cout << " maxSamples=" << g_catalogMaxSamples;
        std::cout << " stepsPerFrame=" << g_catalogStepsPerFrame
                  << " render=" << (g_catalogRenderEnabled ? "on" : "off")
                  << " file=" << g_knotSeedCatalogPath << std::endl;
    }

    SimulationModel* model = new SimulationModel();
    model->init();
    Simulation::getCurrent()->setModel(model);

    buildModel();
    captureGravityReference();
    applyGravityState();
    restoreDemoWindow();

    // Keyboard input for this demo lives here.
    MiniGL::setClientIdleFunc(timeStep);
    MiniGL::setClientSceneFunc(render);
    MiniGL::addMousePressFunc(handleMousePress);
    MiniGL::addMouseMoveFunc(handleMouseMove);
    MiniGL::addKeyFunc('r', reset);
    MiniGL::addKeyFunc('t', []() { evaluateCurrentCableTanglement(false); });
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
    // MiniGL::addKeyFunc('P', []() {
    //     MiniGL::toggleMonochromeMode();
    // });
    // MiniGL::addKeyFunc('p', []() {
    //     MiniGL::toggleMonochromeMode();
    // });
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
    g_lastPrintedHoopSelection = CableVisualization::HoopMarkerSelection();
    g_cachedCableKymoKnotResults.clear();
    g_hasCableKymoKnotResults = false;
    g_visualizationState.clear();
    g_cablePieces.clear();
    g_autoKnotEvaluationDone = false;
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
    g_lastPrintedHoopSelection = CableVisualization::HoopMarkerSelection();
    g_autoKnotEvaluationDone = false;
    g_useMathematicalKnot = false;
    g_cfg.randomCable.enabled = true;
    g_cablePieces.clear();

    if (g_useDirectRandomSeed)
    {
        g_cfg.randomCable.seed = g_directRandomSeed;
    }
    else if (g_useKnotSeedPreset)
    {
        g_cfg.randomCable.seed = g_selectedKnotSeed;
    }
    else
    {
        std::random_device rd;
        unsigned int newSeed = rd();
        if (newSeed == 0)
            newSeed = 1;
        g_cfg.randomCable.seed = newSeed;
    }

    std::cout << "Regenerating random cable with seed " << g_cfg.randomCable.seed << std::endl;
    reset();
}

static void regenerateMathematicalKnot()
{
    g_mouseDrag = MouseDragState();
    g_hoopMarkerSelection = CableVisualization::HoopMarkerSelection();
    g_lastPrintedHoopSelection = CableVisualization::HoopMarkerSelection();
    g_autoKnotEvaluationDone = false;
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
    g_lastPrintedHoopSelection = CableVisualization::HoopMarkerSelection();
    g_autoKnotEvaluationDone = false;
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
    if (g_useKnotSeedPreset)
    {
        ropePoints = CableGeneration::generateRandomCableMessPolyline(
            g_cfg.layout,
            g_cfg.rope,
            g_cfg.randomCable,
            g_cfg.groundY);
    }
    else if (g_useMathematicalKnot)
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
    {
        return false;
    }

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
    if (!g_catalogRunEnabled &&
        !g_evalRunEnabled &&
        !g_evalCutsToZeroEnabled &&
        base->getValue<bool>(DemoBase::PAUSE))
    {
        return;
    }

    SimulationModel* model = Simulation::getCurrent()->getModel();
    ParticleData& particles = model->getParticles();

    const unsigned int stepsPerFrame =
        (g_catalogRunEnabled || g_evalRunEnabled || g_evalCutsToZeroEnabled)
            ? g_catalogStepsPerFrame
            : base->getValue<unsigned int>(DemoBase::NUM_STEPS_PER_RENDER);
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
        if (!g_autoKnotEvaluationDone &&
            g_cfg.restCreep.enabled &&
            g_restCreepCounter >= g_cfg.restCreep.activationDelaySteps)
        {
            g_autoKnotEvaluationDone = true;
            if (g_catalogRunEnabled)
            {
                evaluateAndRecordCatalogSample(false);
                if (g_catalogMaxSamples > 0 && g_catalogSamplesWritten >= g_catalogMaxSamples)
                {
                    g_catalogRunEnabled = false;
                    std::cout << "[knot-catalog] reached maxSamples="
                              << g_catalogSamplesWritten
                              << ", output=" << g_knotSeedCatalogPath << std::endl;
                    MiniGL::leaveMainLoop();
                    return;
                }

                generateNextCatalogSeed();
                return;
            }
            else if (g_evalRunEnabled)
            {
                evaluateAndRecordCutHeuristicSample();
                ++g_evalCurrentIndex;
                if (g_evalCurrentIndex >= g_evalSeeds.size())
                {
                    g_evalRunEnabled = false;
                    std::cout << "[cut-eval] finished " << g_evalSamplesWritten
                              << " sample(s), output=" << g_evalOutputPath << std::endl;
                    MiniGL::leaveMainLoop();
                    return;
                }

                generateNextCutEvalSeed();
                return;
            }
            else if (g_evalCutsToZeroEnabled)
            {
                evaluateAndRecordCutsToZeroSample();
                ++g_evalCurrentIndex;
                if (g_evalCurrentIndex >= g_evalSeeds.size())
                {
                    g_evalCutsToZeroEnabled = false;
                    std::cout << "[cuts-to-zero] finished " << g_evalSamplesWritten
                              << " sample(s), output=" << g_evalCutsToZeroOutputPath << std::endl;
                    MiniGL::leaveMainLoop();
                    return;
                }

                generateNextCutEvalSeed();
                return;
            }
            else
            {
                if (!g_useDirectRandomSeed)
                    evaluateCurrentCableTanglement(true);
            }
        }
    }
}

static void render()
{
    if ((g_catalogRunEnabled && !g_catalogRenderEnabled) ||
        g_evalRunEnabled ||
        g_evalCutsToZeroEnabled)
    {
        return;
    }

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
        g_visualizationState,
        false);

    printSelectedHoopScoreIfChanged();
}
