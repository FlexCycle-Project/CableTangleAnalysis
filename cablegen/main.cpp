#include "Common/Common.h"
#include "Simulation/TimeManager.h"
#include "Simulation/SimulationModel.h"
#include "Simulation/Constraints.h"
#include "Simulation/Simulation.h"
#include "CablePhysics.h"
#include "CableGeneration.h"
#include "DeterministicRandom.h"
#include "TanglementAnalyzer.h"
#include "TopologyAudit.h"
#include "KymoKnotAnalyzer.h"
#include "Utils/Logger.h"
#include "Utils/Timing.h"

#include <Eigen/Dense>
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
#include <unordered_map>
#include <utility>
#include <filesystem>
#include <chrono>

#if defined(_DEBUG) && !defined(EIGEN_ALIGN)
#define new DEBUG_NEW
#endif

using namespace PBD;
using namespace Eigen;
using namespace std;

INIT_LOGGING
INIT_TIMING

static void buildModel();
static void reset();
static void timeStep();

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
static void evaluateAndRecordCatalogSample(bool printReport);
static void evaluateAndRecordCutHeuristicSample();
static void evaluateAndRecordCutsToZeroSample();
static void appendCatalogRows(
    const std::vector<std::vector<Vector3r>>& ropePieces,
    const std::vector<SelfCollisionContact>& contacts,
    const TanglementAnalyzer::TanglementAnalysisResult& tanglement);
static bool storeKnotSeedPresetEntry(
    const std::string& knotType,
    unsigned int seed,
    bool appendToFile);
static bool isNontrivialKnotType(const std::string& knotType);
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
static bool loadCutEvalCatalog();
static bool prepareCutEvalOutputFile();
static bool prepareCutsToZeroEvalOutputFile();
static void prepareCutEvalSeed(size_t index);
static void generateNextCutEvalSeed();
static std::string resolveDefaultTopologyAuditPrefix();
static bool prepareTopologyAuditOutputFiles();
static void beginTopologyAuditPhase(
    const std::string& phase,
    const std::vector<std::vector<Vector3r>>& pieces);
static bool auditTopologyStep(
    const std::vector<std::vector<Vector3r>>& pieces);
static bool applyTopologyAuditCut();
static void finishTopologyAuditSeed();
static void finishTopologyAuditRun();
static void writeTopologyAuditSummary();

static AppConfig g_cfg;
static unsigned int g_restCreepCounter = 0;
static std::vector<PBD::KymoKnotAnalysisResult> g_cachedCableKymoKnotResults;
static bool g_hasCableKymoKnotResults = false;
static PBD::KymoKnotSearchMode g_kymoKnotSearchMode = PBD::KymoKnotSearchMode::TopDown;
static std::vector<std::vector<Vector3r>> g_cablePieces;
static std::vector<std::pair<unsigned int, unsigned int>> g_cablePieceRanges;
static std::map<std::string, std::vector<unsigned int>> g_knotSeedPresets;
static std::string g_knotSeedBankPath;
static std::string g_knotSeedCatalogPath;
static bool g_catalogRunEnabled = false;
static unsigned int g_catalogSamplesWritten = 0;
static unsigned int g_catalogStepsPerFrame = 1;
static bool g_evalRunEnabled = false;
static std::string g_evalOutputPath;
static unsigned int g_evalRandomTrials = 10;
static bool g_evalOracleEnabled = false;
static bool g_evalCutsToZeroEnabled = false;
static std::string g_evalCutsToZeroOutputPath;
static std::string g_evalStageMetricsOutputPath;
static std::string g_evalCutsToZeroStageMetricsOutputPath;
static unsigned int g_evalCutsToZeroMaxCuts = 10;
static unsigned int g_evalCutsToZeroRandomTrials = 20;
static bool g_evalCutsToZeroOracleEnabled = false;
static unsigned long long g_evalCutsToZeroOracleMaxNodes = 5000000ULL;
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

struct TopologyAuditPhaseSummary
{
    std::string phase;
    unsigned int steps = 0;
    size_t crossings = 0;
    size_t violations = 0;
    size_t centerlineContacts = 0;
    size_t surfaceOverlaps = 0;
    Real minAbsHeight = std::numeric_limits<Real>::infinity();
    bool hasFirstViolation = false;
    TopologyAudit::Violation firstViolation;
    bool cutApplied = false;
    size_t cutPiece = 0;
    int cutEdge = -1;

    // Contact/solver diagnostics accumulated over all audited substeps of the
    // phase. "maxPenetration" and "maxGroundPenetration" are pre-correction /
    // pre-clamp maxima (what the solver had to resolve). "residual*" fields are
    // measured after the final self-collision iteration of each substep, before
    // the ground clamp, and bound what remained unresolved in that discrete
    // state.
    size_t candidatePairCount = 0;
    size_t penetratingContactCount = 0;
    Real maxPenetration = static_cast<Real>(0.0);
    Real sumPenetration = static_cast<Real>(0.0);
    size_t correctionEndpointCount = 0;
    size_t cappedCorrectionCount = 0;
    Real maxAppliedCorrection = static_cast<Real>(0.0);
    size_t normalFlipViolationCount = 0;
    size_t residualCandidateCount = 0;
    size_t residualPenetratingContactCount = 0;
    Real residualMaxPenetration = static_cast<Real>(0.0);
    Real residualSumPenetration = static_cast<Real>(0.0);
    Real maxGroundPenetration = static_cast<Real>(0.0);
};

static bool g_topologyAuditRunEnabled = false;
static bool g_topologyAuditActive = false;
static bool g_topologyAuditTraceEnabled = false;
static unsigned int g_topologyAuditSteps = 1000;
static unsigned int g_topologyAuditMaxCuts = 0;
static unsigned int g_topologyAuditMaxSeeds = 0;
static unsigned int g_topologyAuditCutsApplied = 0;
static unsigned int g_topologyAuditStepInPhase = 0;
static std::string g_topologyAuditOutputPrefix;
static std::ofstream g_topologyAuditSeedsOut;
static std::ofstream g_topologyAuditEventsOut;
static std::ofstream g_topologyAuditTraceOut;
static std::ofstream g_topologyAuditSummaryOut;
static TopologyAudit::Auditor g_topologyAuditor;
static TopologyAudit::Settings g_topologyAuditSettings;
static TopologyAuditPhaseSummary g_topologyPhase;
static size_t g_topologyAuditTotalSeeds = 0;
static size_t g_topologyAuditTotalPhases = 0;
static size_t g_topologyAuditTotalSteps = 0;
static size_t g_topologyAuditTotalViolations = 0;
static size_t g_topologyAuditTotalCenterlineContacts = 0;
static size_t g_topologyAuditTotalSurfaceOverlaps = 0;
static size_t g_topologyAuditSeedsWithViolations = 0;
static bool g_topologyAuditCurrentSeedHadViolation = false;
static size_t g_topologyAuditTotalCandidatePairs = 0;
static size_t g_topologyAuditTotalPenetratingContacts = 0;
static Real g_topologyAuditMaxPenetration = static_cast<Real>(0.0);
static Real g_topologyAuditSumPenetration = static_cast<Real>(0.0);
static size_t g_topologyAuditTotalCappedCorrections = 0;
static size_t g_topologyAuditTotalCorrectionEndpoints = 0;
static Real g_topologyAuditMaxAppliedCorrection = static_cast<Real>(0.0);
static size_t g_topologyAuditTotalNormalFlips = 0;
static size_t g_topologyAuditTotalResidualCandidates = 0;
static size_t g_topologyAuditTotalResidualPenetratingContacts = 0;
static Real g_topologyAuditResidualMaxPenetration = static_cast<Real>(0.0);
static Real g_topologyAuditResidualSumPenetration = static_cast<Real>(0.0);
static Real g_topologyAuditMaxGroundPenetration = static_cast<Real>(0.0);
static SelfCollisionDiagnostics g_lastSubstepCollision;
static Real g_lastSubstepGroundPenetration = static_cast<Real>(0.0);

static Vector3r g_gravityReference = Vector3r(0.0, -9.81, 0.0);
static bool g_autoKnotEvaluationDone = false;
static bool g_stopRequested = false;
static void captureGravityReference();
static void applyGravityState();

static std::string resolveDefaultKnotSeedBankPath()
{
    return "knot_seed_bank.txt";
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

static std::string deriveStageMetricsPath(const std::string& evalPath)
{
    const size_t slashPos = evalPath.find_last_of("\\/");
    const size_t dotPos = evalPath.find_last_of('.');
    if (dotPos == std::string::npos ||
        (slashPos != std::string::npos && dotPos < slashPos))
    {
        return evalPath + "_stage_metrics.tsv";
    }
    return evalPath.substr(0, dotPos) + "_stage_metrics.tsv";
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

static bool prepareStageMetricsOutputFile(
    std::string& path,
    const std::string& fallbackPath)
{
    if (path.empty())
        path = deriveStageMetricsPath(fallbackPath);

    std::ofstream out(path, std::ios::trunc);
    if (!out)
    {
        std::cout << "[stage-metrics] could not open output file "
                  << path << std::endl;
        return false;
    }

    out << "seed\tknotType\tstrategy\ttrial\tstage\tcuts\tthreading\tresolved"
        << "\tpieceCount\tpointCount\tsegmentCount\tselfCollisionCount"
        << "\tclosePairCount\ttotalCurvature\toutsideKnotCurvature\n";
    return true;
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
        << "\tbiggest_hoop_span\tbiggest_hoop_score\thoop_count"
        << "\tdelta_smallest_hoop\tsmallest_hoop_span\tsmallest_hoop_score"
        << "\tdelta_max_curvature\tmax_curvature"
        << "\tdelta_oracle\toracle_evaluated\toracle_piece_index\toracle_edge_index\n";

    return prepareStageMetricsOutputFile(
        g_evalStageMetricsOutputPath,
        g_evalOutputPath);
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
        << "\tbiggest_hoop_cuts\tbiggest_hoop_resolved\thoop_count"
        << "\tsmallest_hoop_cuts\tsmallest_hoop_resolved"
        << "\tmax_curvature_cuts\tmax_curvature_resolved"
        << "\toracle_cuts\toracle_resolved\toracle_exhaustive"
        << "\toracle_nodes\toracle_sequence\n";

    return prepareStageMetricsOutputFile(
        g_evalCutsToZeroStageMetricsOutputPath,
        g_evalCutsToZeroOutputPath);
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
    bool writeCurvatureRawColumns = true;
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
            writeCurvatureRawColumns =
                headerLine.find("\ttotalCurvature\t") != std::string::npos;
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
            << "\ttotalCurvature\toutsideKnotCurvature"
            << "\ttotalCurvatureDensity\toutsideKnotCurvatureDensity\n";
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

    // The entanglement/complexity score is derived from this generated catalog
    // by the Julia analysis; the catalog stores raw metrics only.

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
            << metrics.samePieceProximityScore;
        if (writeLegacyClusterColumn)
            out << '\t' << static_cast<Real>(0.0);
        if (writeCurvatureRawColumns)
            out << '\t' << totalCurvature << '\t' << outsideKnotCurvature;
        if (writeCurvatureDensityColumns)
            out << '\t' << totalCurvatureDensity << '\t' << outsideKnotCurvatureDensity;
        out << '\n';
    }

    std::cout << "[knot-catalog] wrote " << records.size()
              << " row(s) for seed=" << g_cfg.randomCable.seed
              << " to " << g_knotSeedCatalogPath << std::endl;
}

struct CutStageLogContext
{
    std::string eval;
    unsigned int seed = 0;
    std::string knotType;
    std::string strategy;
    unsigned int trial = 0;
    bool enabled = false;
};

struct CutStageMetrics
{
    int threading = 0;
    bool resolved = false;
    size_t pieceCount = 0;
    size_t pointCount = 0;
    size_t segmentCount = 0;
    size_t selfCollisionCount = 0;
    size_t closePairCount = 0;
    Real totalCurvature = static_cast<Real>(0.0);
    Real outsideKnotCurvature = static_cast<Real>(0.0);
};

static CutStageMetrics computeCutStageMetrics(
    const std::vector<std::vector<Vector3r>>& ropePieces)
{
    CutStageMetrics stage;
    const std::vector<SelfCollisionContact> contacts =
        detectSelfCollisionContactsForPieces(
            ropePieces,
            g_cfg.rope,
            g_cfg.selfCollision);
    const TanglementAnalyzer::TanglementAnalysisResult tanglement =
        TanglementAnalyzer::analyzeTanglement(
            ropePieces,
            contacts,
            g_cfg.rope,
            g_cfg.selfCollision,
            g_kymoKnotSearchMode);
    const TanglementAnalyzer::TanglementMetricsSnapshot& metrics = tanglement.metrics;

    stage.threading =
        TanglementAnalyzer::computeThreadingMetrics(tanglement).totalThreadingCount;
    stage.resolved = stage.threading == 0;
    stage.pieceCount = metrics.pieceCount;
    stage.pointCount = metrics.totalPointCount;
    stage.segmentCount = metrics.totalSegmentCount;
    stage.selfCollisionCount = contacts.size();
    stage.closePairCount = metrics.closePairCount;

    Real totalCurvature = static_cast<Real>(0.0);
    for (const std::vector<Vector3r>& piece : ropePieces)
        totalCurvature += computeTotalCurvature(piece);
    stage.totalCurvature = totalCurvature;

    Real outsideKnotCurvature = totalCurvature;
    for (size_t pieceIndex = 0; pieceIndex < tanglement.kymoResults.size(); ++pieceIndex)
    {
        if (pieceIndex >= ropePieces.size())
            continue;

        const PBD::KymoKnotAnalysisResult& knotResult = tanglement.kymoResults[pieceIndex];
        if (!knotResult.available || !knotResult.success || knotResult.segments.empty())
            continue;

        bool hasNontrivialKnot = false;
        for (const PBD::KymoKnotSegment& segment : knotResult.segments)
        {
            if (isNontrivialKnotType(segment.knotType))
            {
                hasNontrivialKnot = true;
                break;
            }
        }
        if (!hasNontrivialKnot)
            continue;

        outsideKnotCurvature = computeCurvatureOutsideKnotRegions(
            ropePieces[pieceIndex],
            knotResult,
            nullptr);
    }
    stage.outsideKnotCurvature = outsideKnotCurvature;

    return stage;
}

static void appendCutStageMetricsRow(
    const CutStageLogContext& context,
    unsigned int stage,
    unsigned int cuts,
    const CutStageMetrics& metrics)
{
    const std::string& path = context.eval == "cuts_to_zero"
        ? g_evalCutsToZeroStageMetricsOutputPath
        : g_evalStageMetricsOutputPath;
    if (!context.enabled || path.empty())
        return;

    std::ofstream out(path, std::ios::app);
    if (!out)
    {
        std::cout << "[stage-metrics] could not append to " << path << std::endl;
        return;
    }

    out << std::fixed << std::setprecision(6)
        << context.seed << '\t'
        << context.knotType << '\t'
        << context.strategy << '\t'
        << context.trial << '\t'
        << stage << '\t'
        << cuts << '\t'
        << metrics.threading << '\t'
        << (metrics.resolved ? 1 : 0) << '\t'
        << metrics.pieceCount << '\t'
        << metrics.pointCount << '\t'
        << metrics.segmentCount << '\t'
        << metrics.selfCollisionCount << '\t'
        << metrics.closePairCount << '\t'
        << metrics.totalCurvature << '\t'
        << metrics.outsideKnotCurvature << '\n';
}

static void recordCutStage(
    const CutStageLogContext& context,
    unsigned int stage,
    unsigned int cuts,
    const std::vector<std::vector<Vector3r>>& ropePieces)
{
    if (!context.enabled)
        return;
    appendCutStageMetricsRow(
        context,
        stage,
        cuts,
        computeCutStageMetrics(ropePieces));
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
        Real curvature = static_cast<Real>(0.0);
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

    static VirtualCutTarget findSmallestHoopMidpointCut(
        const TanglementAnalyzer::TanglementAnalysisResult& tanglement,
        const std::vector<std::vector<Vector3r>>& ropePieces)
    {
        VirtualCutTarget best;
        for (const VirtualCutTarget& candidate : collectHoopCutTargets(tanglement, ropePieces))
        {
            if (!best.valid ||
                candidate.hoopSpan < best.hoopSpan ||
                (candidate.hoopSpan == best.hoopSpan && candidate.hoopScore < best.hoopScore))
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

        target = targets[DeterministicRandom::uniformIndex(rng, targets.size())];
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

    static Real computeVertexCurvature(
        const std::vector<Vector3r>& piece,
        size_t vertexIndex)
    {
        if (vertexIndex == 0 || vertexIndex + 1 >= piece.size())
            return static_cast<Real>(0.0);

        const Vector3r a = piece[vertexIndex] - piece[vertexIndex - 1];
        const Vector3r b = piece[vertexIndex + 1] - piece[vertexIndex];
        const Real aNorm = a.norm();
        const Real bNorm = b.norm();
        if (aNorm <= static_cast<Real>(1e-8) || bNorm <= static_cast<Real>(1e-8))
            return static_cast<Real>(0.0);

        const Real c = std::max(
            static_cast<Real>(-1.0),
            std::min(static_cast<Real>(1.0), a.dot(b) / (aNorm * bNorm)));
        return std::acos(c);
    }

    static Real computeEdgeCurvature(
        const std::vector<Vector3r>& piece,
        size_t edgeIndex)
    {
        return std::max(
            computeVertexCurvature(piece, edgeIndex),
            computeVertexCurvature(piece, edgeIndex + 1));
    }

    static VirtualCutTarget findMaxCurvatureCutTarget(
        const std::vector<std::vector<Vector3r>>& ropePieces,
        Real excludeArcLength)
    {
        VirtualCutTarget best;
        for (VirtualCutTarget target : collectValidRandomCutTargets(ropePieces, excludeArcLength))
        {
            target.curvature = computeEdgeCurvature(
                ropePieces[target.pieceIndex],
                target.edgeIndex);
            if (!best.valid || target.curvature > best.curvature)
                best = target;
        }
        return best;
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

    static int computeAndRecordVirtualCutDeltaThreading(
        const std::vector<std::vector<Vector3r>>& ropePieces,
        int beforeThreading,
        const VirtualCutTarget& target,
        const CutStageLogContext& context)
    {
        if (!target.valid)
            return 0;

        std::vector<std::vector<Vector3r>> splitPieces = ropePieces;
        if (!splitCablePieceAtEdge(
                splitPieces,
                target.pieceIndex,
                target.edgeIndex))
        {
            return 0;
        }

        const int afterThreading = computeThreadingCountGeometryOnly(splitPieces);
        recordCutStage(context, 1, 1, splitPieces);
        return beforeThreading - afterThreading;
    }

    struct SingleCutOracleResult
    {
        VirtualCutTarget target;
        int delta = 0;
        bool evaluated = false;
    };

    static SingleCutOracleResult findBestSingleCutOracle(
        const std::vector<std::vector<Vector3r>>& ropePieces,
        int beforeThreading,
        Real excludeArcLength)
    {
        SingleCutOracleResult result;
        if (beforeThreading <= 0)
            return result;

        const std::vector<VirtualCutTarget> targets =
            collectValidRandomCutTargets(ropePieces, excludeArcLength);
        result.evaluated = !targets.empty();

        for (const VirtualCutTarget& target : targets)
        {
            const int delta = computeVirtualCutDeltaThreading(
                ropePieces,
                beforeThreading,
                target);
            if (!result.target.valid || delta > result.delta)
            {
                result.target = target;
                result.delta = delta;
            }
        }
        return result;
    }

    struct CutsToZeroOutcome
    {
        unsigned int cuts = 0;
        bool resolved = false;
    };

    static CutsToZeroOutcome runHeuristicCutsToZero(
        const std::vector<std::vector<Vector3r>>& startPieces,
        unsigned int maxCuts,
        const CutStageLogContext& stageContext)
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
            recordCutStage(stageContext, outcome.cuts, outcome.cuts, pieces);
        }

        outcome.resolved = computeThreadingCountGeometryOnly(pieces) == 0;
        if (!outcome.resolved)
            outcome.cuts = maxCuts;
        return outcome;
    }

    static CutsToZeroOutcome runBiggestHoopCutsToZero(
        const std::vector<std::vector<Vector3r>>& startPieces,
        unsigned int maxCuts,
        const CutStageLogContext& stageContext)
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
            recordCutStage(stageContext, outcome.cuts, outcome.cuts, pieces);
        }

        outcome.resolved = computeThreadingCountGeometryOnly(pieces) == 0;
        if (!outcome.resolved)
            outcome.cuts = maxCuts;
        return outcome;
    }

    static CutsToZeroOutcome runSmallestHoopCutsToZero(
        const std::vector<std::vector<Vector3r>>& startPieces,
        unsigned int maxCuts,
        const CutStageLogContext& stageContext)
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

            const VirtualCutTarget target = findSmallestHoopMidpointCut(tanglement, pieces);
            if (!target.valid || !splitCablePieceAtEdge(pieces, target.pieceIndex, target.edgeIndex))
                break;

            ++outcome.cuts;
            recordCutStage(stageContext, outcome.cuts, outcome.cuts, pieces);
        }

        outcome.resolved = computeThreadingCountGeometryOnly(pieces) == 0;
        if (!outcome.resolved)
            outcome.cuts = maxCuts;
        return outcome;
    }

    static CutsToZeroOutcome runMaxCurvatureCutsToZero(
        const std::vector<std::vector<Vector3r>>& startPieces,
        unsigned int maxCuts,
        const CutStageLogContext& stageContext)
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

            const VirtualCutTarget target = findMaxCurvatureCutTarget(
                pieces,
                g_cfg.selfCollision.excludeArcLength);
            if (!target.valid || !splitCablePieceAtEdge(pieces, target.pieceIndex, target.edgeIndex))
                break;

            ++outcome.cuts;
            recordCutStage(stageContext, outcome.cuts, outcome.cuts, pieces);
        }

        outcome.resolved = computeThreadingCountGeometryOnly(pieces) == 0;
        if (!outcome.resolved)
            outcome.cuts = maxCuts;
        return outcome;
    }

    static CutsToZeroOutcome runRandomHoopCutsToZero(
        const std::vector<std::vector<Vector3r>>& startPieces,
        unsigned int maxCuts,
        std::mt19937& rng,
        const CutStageLogContext& stageContext)
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
            recordCutStage(stageContext, outcome.cuts, outcome.cuts, pieces);
        }

        outcome.resolved = computeThreadingCountGeometryOnly(pieces) == 0;
        if (!outcome.resolved)
            outcome.cuts = maxCuts;
        return outcome;
    }

    static CutsToZeroOutcome runRandomCutsToZero(
        const std::vector<std::vector<Vector3r>>& startPieces,
        unsigned int maxCuts,
        std::mt19937& rng,
        const CutStageLogContext& stageContext)
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

            const VirtualCutTarget& target =
                targets[DeterministicRandom::uniformIndex(rng, targets.size())];
            if (!splitCablePieceAtEdge(pieces, target.pieceIndex, target.edgeIndex))
                break;

            ++outcome.cuts;
            recordCutStage(stageContext, outcome.cuts, outcome.cuts, pieces);
        }

        outcome.resolved = computeThreadingCountGeometryOnly(pieces) == 0;
        if (!outcome.resolved)
            outcome.cuts = maxCuts;
        return outcome;
    }

    // ── Complete hoop-midpoint oracle ────────────────────────────────────
    //
    // The oracle evaluates only hoop-midpoint cuts, but it looks ahead until
    // it can prove the smallest number of cuts that reaches threading == 0.
    // Cuts are identified by their original (piece, edge) coordinates so the
    // search state is independent of cut order.

    struct OracleCut
    {
        size_t pieceIndex = 0;
        size_t edgeIndex = 0;
    };

    struct OracleState
    {
        std::vector<std::vector<Vector3r>> pieces;
        std::vector<size_t> originalPiece;
        std::vector<size_t> originalOffset;
        std::vector<OracleCut> cuts;
        int threading = 0;
    };

    struct OracleCandidate
    {
        OracleCut cut;
        OracleState child;
        int threadingAfter = 0;
        int delta = 0;
        int hoopScore = 0;
    };

    struct OraclePathResult
    {
        std::vector<OracleCut> path;
        unsigned int cuts = 0;
        bool resolved = false;
        bool exhaustive = false;
        unsigned long long nodes = 0;
    };

    static bool oracleCutLess(const OracleCut& a, const OracleCut& b)
    {
        if (a.pieceIndex != b.pieceIndex)
            return a.pieceIndex < b.pieceIndex;
        return a.edgeIndex < b.edgeIndex;
    }

    static std::string oracleStateKey(const OracleState& state)
    {
        std::string key;
        key.reserve(state.cuts.size() * sizeof(uint64_t));
        for (const OracleCut& cut : state.cuts)
        {
            const uint64_t packed =
                (static_cast<uint64_t>(cut.pieceIndex) << 32) |
                static_cast<uint64_t>(cut.edgeIndex);
            key.append(
                reinterpret_cast<const char*>(&packed),
                sizeof(packed));
        }
        return key;
    }

    static OracleState makeOracleState(
        const std::vector<std::vector<Vector3r>>& startPieces)
    {
        OracleState state;
        state.pieces = startPieces;
        state.originalPiece.resize(startPieces.size());
        state.originalOffset.assign(startPieces.size(), 0);
        for (size_t i = 0; i < startPieces.size(); ++i)
            state.originalPiece[i] = i;
        state.threading = computeThreadingCountGeometryOnly(state.pieces);
        return state;
    }

    static bool applyOracleCut(
        const OracleState& state,
        size_t pieceIndex,
        size_t edgeIndex,
        OracleState& child,
        OracleCut& cut)
    {
        if (pieceIndex >= state.pieces.size())
            return false;
        const std::vector<Vector3r>& piece = state.pieces[pieceIndex];
        if (piece.size() < 4 || edgeIndex < 1 || edgeIndex + 2 >= piece.size())
            return false;

        child = state;
        std::vector<Vector3r>& target = child.pieces[pieceIndex];
        std::vector<Vector3r> left(
            target.begin(),
            target.begin() + static_cast<std::ptrdiff_t>(edgeIndex + 1));
        std::vector<Vector3r> right(
            target.begin() + static_cast<std::ptrdiff_t>(edgeIndex + 1),
            target.end());
        if (left.size() < 2 || right.size() < 2)
            return false;

        target = std::move(left);
        child.pieces.insert(
            child.pieces.begin() + static_cast<std::ptrdiff_t>(pieceIndex + 1),
            std::move(right));

        const size_t originPiece = state.originalPiece[pieceIndex];
        const size_t originOffset = state.originalOffset[pieceIndex];
        child.originalPiece.insert(
            child.originalPiece.begin() + static_cast<std::ptrdiff_t>(pieceIndex + 1),
            originPiece);
        child.originalOffset.insert(
            child.originalOffset.begin() + static_cast<std::ptrdiff_t>(pieceIndex + 1),
            originOffset + edgeIndex + 1);

        cut.pieceIndex = originPiece;
        cut.edgeIndex = originOffset + edgeIndex;
        const auto insertAt = std::lower_bound(
            child.cuts.begin(),
            child.cuts.end(),
            cut,
            oracleCutLess);
        child.cuts.insert(insertAt, cut);
        child.threading = computeThreadingCountGeometryOnly(child.pieces);
        return true;
    }

    static bool applyOracleCutById(
        const OracleState& state,
        const OracleCut& cut,
        OracleState& child)
    {
        for (size_t pieceIndex = 0; pieceIndex < state.pieces.size(); ++pieceIndex)
        {
            if (state.originalPiece[pieceIndex] != cut.pieceIndex)
                continue;

            const size_t offset = state.originalOffset[pieceIndex];
            if (cut.edgeIndex < offset)
                continue;
            if (state.pieces[pieceIndex].size() < 2)
                continue;

            const size_t localEdge = cut.edgeIndex - offset;
            if (localEdge >= state.pieces[pieceIndex].size() - 1)
                continue;

            OracleCut applied;
            return applyOracleCut(state, pieceIndex, localEdge, child, applied);
        }
        return false;
    }

    static void sortOracleCandidates(
        std::vector<OracleCandidate>& candidates,
        bool byHoopScore)
    {
        if (byHoopScore)
        {
            std::sort(
                candidates.begin(),
                candidates.end(),
                [](const OracleCandidate& a, const OracleCandidate& b)
                {
                    if (a.hoopScore != b.hoopScore)
                        return a.hoopScore > b.hoopScore;
                    if (a.delta != b.delta)
                        return a.delta > b.delta;
                    if (a.threadingAfter != b.threadingAfter)
                        return a.threadingAfter < b.threadingAfter;
                    return oracleCutLess(a.cut, b.cut);
                });
            return;
        }

        std::sort(
            candidates.begin(),
            candidates.end(),
            [](const OracleCandidate& a, const OracleCandidate& b)
            {
                if (a.delta != b.delta)
                    return a.delta > b.delta;
                if (a.threadingAfter != b.threadingAfter)
                    return a.threadingAfter < b.threadingAfter;
                if (a.hoopScore != b.hoopScore)
                    return a.hoopScore > b.hoopScore;
                return oracleCutLess(a.cut, b.cut);
            });
    }

    static std::vector<OracleCandidate> collectOracleCandidates(
        const OracleState& state)
    {
        std::vector<OracleCandidate> candidates;
        if (state.threading == 0)
            return candidates;

        const TanglementAnalyzer::TanglementAnalysisResult tanglement =
            TanglementAnalyzer::analyzeTanglementGeometryOnly(
                state.pieces,
                {},
                g_cfg.rope,
                g_cfg.selfCollision);
        const std::vector<VirtualCutTarget> targets =
            collectHoopCutTargets(tanglement, state.pieces);

        std::set<std::pair<size_t, size_t>> seen;
        for (const VirtualCutTarget& target : targets)
        {
            if (!target.valid)
                continue;
            const std::pair<size_t, size_t> key(
                target.pieceIndex,
                target.edgeIndex);
            if (!seen.insert(key).second)
                continue;

            OracleState child;
            OracleCut cut;
            if (!applyOracleCut(
                    state,
                    target.pieceIndex,
                    target.edgeIndex,
                    child,
                    cut))
            {
                continue;
            }

            OracleCandidate candidate;
            candidate.cut = cut;
            candidate.threadingAfter = child.threading;
            candidate.delta = state.threading - child.threading;
            candidate.hoopScore = target.hoopScore;
            candidate.child = std::move(child);
            candidates.push_back(std::move(candidate));
        }

        sortOracleCandidates(candidates, false);
        return candidates;
    }

    static OraclePathResult runOracleGreedy(
        const OracleState& start,
        unsigned int maxCuts,
        bool byHoopScore)
    {
        OraclePathResult result;
        OracleState state = start;

        while (state.threading > 0 && result.path.size() < maxCuts)
        {
            std::vector<OracleCandidate> candidates =
                collectOracleCandidates(state);
            if (candidates.empty())
                break;

            if (byHoopScore)
                sortOracleCandidates(candidates, true);

            const OracleCandidate& best = candidates.front();
            result.path.push_back(best.cut);
            state = std::move(best.child);
        }

        result.resolved = state.threading == 0;
        result.cuts = static_cast<unsigned int>(result.path.size());
        if (!result.resolved)
            result.cuts = maxCuts;
        return result;
    }

    struct OracleSearchBudget
    {
        unsigned long long maxNodes = 0;
        unsigned long long nodes = 0;
        bool aborted = false;
    };

    static bool oracleDepthLimitedSearch(
        const OracleState& state,
        int remaining,
        std::unordered_map<std::string, int>& visited,
        std::vector<OracleCut>& path,
        OracleSearchBudget& budget)
    {
        if (state.threading == 0)
            return true;
        if (remaining == 0)
            return false;

        const std::string key = oracleStateKey(state);
        const auto existing = visited.find(key);
        if (existing != visited.end() && existing->second >= remaining)
            return false;
        visited[key] = remaining;

        ++budget.nodes;
        if (budget.maxNodes > 0 && budget.nodes > budget.maxNodes)
        {
            budget.aborted = true;
            return false;
        }

        std::vector<OracleCandidate> candidates =
            collectOracleCandidates(state);
        for (const OracleCandidate& candidate : candidates)
        {
            path.push_back(candidate.cut);
            if (oracleDepthLimitedSearch(
                    candidate.child,
                    remaining - 1,
                    visited,
                    path,
                    budget))
            {
                return true;
            }
            path.pop_back();

            if (budget.aborted)
                return false;
        }
        return false;
    }

    static OraclePathResult runOracleExact(
        const OracleState& start,
        const OraclePathResult& fallback,
        unsigned int maxCuts,
        unsigned long long maxNodes)
    {
        OraclePathResult result;
        result.path = fallback.path;
        result.cuts = fallback.cuts;
        result.resolved = fallback.resolved;
        result.exhaustive = false;

        if (start.threading == 0)
        {
            result.path.clear();
            result.cuts = 0;
            result.resolved = true;
            result.exhaustive = true;
            return result;
        }

        const int upper = fallback.resolved
            ? static_cast<int>(fallback.cuts)
            : static_cast<int>(maxCuts) + 1;
        OracleSearchBudget budget;
        budget.maxNodes = maxNodes;
        bool allLowerExhausted = true;

        for (int limit = 1;
             limit < upper && limit <= static_cast<int>(maxCuts);
             ++limit)
        {
            std::unordered_map<std::string, int> visited;
            std::vector<OracleCut> path;
            if (oracleDepthLimitedSearch(
                    start,
                    limit,
                    visited,
                    path,
                    budget))
            {
                result.path = path;
                result.cuts = static_cast<unsigned int>(limit);
                result.resolved = true;
                result.exhaustive = allLowerExhausted;
                result.nodes = budget.nodes;
                return result;
            }

            if (budget.aborted)
            {
                allLowerExhausted = false;
                break;
            }
        }

        result.nodes = budget.nodes;
        result.exhaustive = allLowerExhausted && !budget.aborted;
        if (!fallback.resolved)
        {
            result.cuts = maxCuts;
            result.resolved = false;
        }
        return result;
    }

    static OraclePathResult runCompleteOracleCutsToZero(
        const std::vector<std::vector<Vector3r>>& startPieces,
        unsigned int maxCuts,
        unsigned long long maxNodes)
    {
        const OracleState start = makeOracleState(startPieces);
        OraclePathResult deltaGreedy = runOracleGreedy(start, maxCuts, false);
        OraclePathResult hoopGreedy = runOracleGreedy(start, maxCuts, true);

        OraclePathResult fallback = deltaGreedy;
        if (hoopGreedy.resolved &&
            (!fallback.resolved || hoopGreedy.cuts < fallback.cuts))
        {
            fallback = hoopGreedy;
        }

        return runOracleExact(start, fallback, maxCuts, maxNodes);
    }

    static bool replayOraclePath(
        const std::vector<std::vector<Vector3r>>& startPieces,
        const std::vector<OracleCut>& path,
        const CutStageLogContext& stageContext,
        unsigned int& finalCuts,
        bool& finalResolved)
    {
        OracleState state = makeOracleState(startPieces);
        recordCutStage(stageContext, 0, 0, state.pieces);

        unsigned int stage = 0;
        for (const OracleCut& cut : path)
        {
            OracleState child;
            if (!applyOracleCutById(state, cut, child))
                return false;
            state = std::move(child);
            ++stage;
            recordCutStage(stageContext, stage, stage, state.pieces);
        }

        finalCuts = stage;
        finalResolved = state.threading == 0;
        return true;
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

    auto makeStageContext = [&](const std::string& strategy, unsigned int trial)
    {
        CutStageLogContext context;
        context.eval = "cut_eval";
        context.seed = record.seed;
        context.knotType = record.knotType;
        context.strategy = strategy;
        context.trial = trial;
        context.enabled = !g_evalStageMetricsOutputPath.empty();
        return context;
    };

    if (!g_evalStageMetricsOutputPath.empty())
        recordCutStage(makeStageContext("initial", 0), 0, 0, ropePieces);

    int skipped = 0;
    int heuristicDelta = 0;
    int heuristicHoopScore = 0;
    if (beforeThreading > 0)
    {
        const VirtualCutTarget heuristicTarget = findBestHoopMidpointCut(tanglement, ropePieces);
        if (heuristicTarget.valid)
        {
            heuristicHoopScore = heuristicTarget.hoopScore;
            heuristicDelta = computeAndRecordVirtualCutDeltaThreading(
                ropePieces,
                beforeThreading,
                heuristicTarget,
                makeStageContext("heuristic", 0));
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
            biggestHoopDelta = computeAndRecordVirtualCutDeltaThreading(
                ropePieces,
                beforeThreading,
                biggestTarget,
                makeStageContext("biggest_hoop", 0));
            biggestHoopSpan = biggestTarget.hoopSpan;
            biggestHoopScore = biggestTarget.hoopScore;
        }
    }

    int smallestHoopDelta = 0;
    int smallestHoopSpan = 0;
    int smallestHoopScore = 0;
    if (beforeThreading > 0)
    {
        const VirtualCutTarget smallestTarget =
            findSmallestHoopMidpointCut(tanglement, ropePieces);
        if (smallestTarget.valid)
        {
            smallestHoopDelta = computeAndRecordVirtualCutDeltaThreading(
                ropePieces,
                beforeThreading,
                smallestTarget,
                makeStageContext("smallest_hoop", 0));
            smallestHoopSpan = smallestTarget.hoopSpan;
            smallestHoopScore = smallestTarget.hoopScore;
        }
    }

    int maxCurvatureDelta = 0;
    Real maxCurvature = static_cast<Real>(0.0);
    if (beforeThreading > 0)
    {
        const VirtualCutTarget curvatureTarget = findMaxCurvatureCutTarget(
            ropePieces,
            g_cfg.selfCollision.excludeArcLength);
        if (curvatureTarget.valid)
        {
            maxCurvatureDelta = computeAndRecordVirtualCutDeltaThreading(
                ropePieces,
                beforeThreading,
                curvatureTarget,
                makeStageContext("max_curvature", 0));
            maxCurvature = curvatureTarget.curvature;
        }
    }

    SingleCutOracleResult oracle;
    if (g_evalOracleEnabled && beforeThreading > 0)
    {
        oracle = findBestSingleCutOracle(
            ropePieces,
            beforeThreading,
            g_cfg.selfCollision.excludeArcLength);
        if (oracle.target.valid)
        {
            computeAndRecordVirtualCutDeltaThreading(
                ropePieces,
                beforeThreading,
                oracle.target,
                makeStageContext("oracle", 0));
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
            std::vector<Real> randomDeltas;
            randomDeltas.reserve(g_evalRandomTrials);

            for (unsigned int trial = 0; trial < g_evalRandomTrials; ++trial)
            {
                const VirtualCutTarget& randomTarget =
                    randomTargets[DeterministicRandom::uniformIndex(rng, randomTargets.size())];
                const int delta = computeAndRecordVirtualCutDeltaThreading(
                    ropePieces,
                    beforeThreading,
                    randomTarget,
                    makeStageContext("random_edge", trial + 1));
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
        std::vector<Real> randomHoopDeltas;
        randomHoopDeltas.reserve(g_evalRandomTrials);

        for (unsigned int trial = 0; trial < g_evalRandomTrials; ++trial)
        {
            const VirtualCutTarget& randomHoopTarget =
                hoopTargets[DeterministicRandom::uniformIndex(rng, hoopTargets.size())];
            const int delta = computeAndRecordVirtualCutDeltaThreading(
                ropePieces,
                beforeThreading,
                randomHoopTarget,
                makeStageContext("random_hoop", trial + 1));
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
        << hoopCount << '\t'
        << smallestHoopDelta << '\t'
        << smallestHoopSpan << '\t'
        << smallestHoopScore << '\t'
        << maxCurvatureDelta << '\t'
        << maxCurvature << '\t'
        << oracle.delta << '\t'
        << (oracle.evaluated ? 1 : 0) << '\t'
        << (oracle.target.valid ? static_cast<int>(oracle.target.pieceIndex) : -1) << '\t'
        << (oracle.target.valid ? static_cast<int>(oracle.target.edgeIndex) : -1) << '\n';

    ++g_evalSamplesWritten;
    std::cout << "[cut-eval] seed=" << record.seed
              << " knotType=" << record.knotType
              << " T_before=" << beforeThreading
              << " heuristicDelta=" << heuristicDelta
              << " randomMean=" << randomMean
              << " randomHoopMean=" << randomHoopMean
              << " biggestHoopDelta=" << biggestHoopDelta
              << " smallestHoopDelta=" << smallestHoopDelta
              << " maxCurvatureDelta=" << maxCurvatureDelta
              << " oracleDelta=" << oracle.delta
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

    auto makeStageContext = [&](const std::string& strategy, unsigned int trial)
    {
        CutStageLogContext context;
        context.eval = "cuts_to_zero";
        context.seed = record.seed;
        context.knotType = record.knotType;
        context.strategy = strategy;
        context.trial = trial;
        context.enabled = !g_evalCutsToZeroStageMetricsOutputPath.empty();
        return context;
    };

    if (!g_evalCutsToZeroStageMetricsOutputPath.empty())
        recordCutStage(makeStageContext("initial", 0), 0, 0, ropePieces);

    CutsToZeroOutcome heuristicOutcome;
    CutsToZeroOutcome biggestHoopOutcome;
    CutsToZeroOutcome smallestHoopOutcome;
    CutsToZeroOutcome maxCurvatureOutcome;
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
        smallestHoopOutcome.resolved = true;
        maxCurvatureOutcome.resolved = true;
    }
    else
    {
        heuristicOutcome = runHeuristicCutsToZero(
            ropePieces,
            g_evalCutsToZeroMaxCuts,
            makeStageContext("heuristic", 0));
        biggestHoopOutcome = runBiggestHoopCutsToZero(
            ropePieces,
            g_evalCutsToZeroMaxCuts,
            makeStageContext("biggest_hoop", 0));
        smallestHoopOutcome = runSmallestHoopCutsToZero(
            ropePieces,
            g_evalCutsToZeroMaxCuts,
            makeStageContext("smallest_hoop", 0));
        maxCurvatureOutcome = runMaxCurvatureCutsToZero(
            ropePieces,
            g_evalCutsToZeroMaxCuts,
            makeStageContext("max_curvature", 0));

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
                    edgeRng,
                    makeStageContext("random_edge", trial + 1));
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
                    hoopRng,
                    makeStageContext("random_hoop", trial + 1));
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

    // Complete (multi-cut) oracle.  This is the exhaustive best-path search,
    // as opposed to the single-cut oracle recorded in cut_eval.tsv.
    const bool haveOracle = g_evalCutsToZeroOracleEnabled;
    unsigned int oracleCuts = 0;
    bool oracleResolved = false;
    bool oracleExhaustive = false;
    unsigned long long oracleNodes = 0;
    std::string oracleSequence;

    if (haveOracle)
    {
        OraclePathResult oracleResult;
        if (beforeThreading == 0)
        {
            oracleResult.resolved = true;
            oracleResult.exhaustive = true;
        }
        else
        {
            oracleResult = runCompleteOracleCutsToZero(
                ropePieces,
                g_evalCutsToZeroMaxCuts,
                g_evalCutsToZeroOracleMaxNodes);
        }

        unsigned int replayCuts = 0;
        bool replayResolved = false;
        const bool replayOk = replayOraclePath(
            ropePieces,
            oracleResult.path,
            makeStageContext("oracle", 0),
            replayCuts,
            replayResolved);

        oracleNodes = oracleResult.nodes;
        oracleCuts = oracleResult.resolved ? replayCuts : g_evalCutsToZeroMaxCuts;
        oracleResolved = oracleResult.resolved && replayResolved;
        oracleExhaustive = oracleResult.exhaustive && replayOk;

        if (!replayOk)
        {
            oracleResolved = false;
            oracleCuts = g_evalCutsToZeroMaxCuts;
        }
        else if (oracleResult.resolved && replayCuts != oracleResult.cuts)
        {
            std::cout << "[ctz-oracle] warning: search found "
                      << oracleResult.cuts << " cuts but replay used "
                      << replayCuts << " cuts for seed=" << record.seed
                      << std::endl;
        }

        for (size_t i = 0; i < oracleResult.path.size(); ++i)
        {
            if (i > 0)
                oracleSequence += ';';
            oracleSequence += std::to_string(oracleResult.path[i].pieceIndex);
            oracleSequence += ':';
            oracleSequence += std::to_string(oracleResult.path[i].edgeIndex);
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
        << hoopCount << '\t'
        << smallestHoopOutcome.cuts << '\t'
        << (smallestHoopOutcome.resolved ? 1 : 0) << '\t'
        << maxCurvatureOutcome.cuts << '\t'
        << (maxCurvatureOutcome.resolved ? 1 : 0) << '\t';
    if (haveOracle)
    {
        out << oracleCuts << '\t'
            << (oracleResolved ? 1 : 0) << '\t'
            << (oracleExhaustive ? 1 : 0) << '\t'
            << oracleNodes << '\t'
            << oracleSequence;
    }
    else
    {
        out << "\t\t\t\t";
    }
    out << '\n';

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
              << " smallestHoopCuts=" << smallestHoopOutcome.cuts
              << " smallestHoopResolved=" << (smallestHoopOutcome.resolved ? 1 : 0)
              << " maxCurvatureCuts=" << maxCurvatureOutcome.cuts
              << " maxCurvatureResolved=" << (maxCurvatureOutcome.resolved ? 1 : 0)
              << " hoops=" << hoopCount;
    if (haveOracle)
    {
        std::cout << " oracleCuts=" << oracleCuts
                  << " oracleResolved=" << (oracleResolved ? 1 : 0)
                  << " oracleExhaustive=" << (oracleExhaustive ? 1 : 0)
                  << " oracleNodes=" << oracleNodes;
    }
    std::cout << " (" << (g_evalCurrentIndex + 1) << "/" << g_evalSeeds.size() << ")"
              << std::endl;
}

static void prepareCutEvalSeed(size_t index)
{
    if (index >= g_evalSeeds.size())
        return;

    g_evalCurrentIndex = index;
    g_cfg.randomCable.enabled = true;
    g_cfg.randomCable.seed = g_evalSeeds[g_evalCurrentIndex].seed;
}

static void generateNextCutEvalSeed()
{
    prepareCutEvalSeed(g_evalCurrentIndex);
    g_cablePieces.clear();
    std::cout << "[cut-eval] next seed=" << g_cfg.randomCable.seed << std::endl;
    reset();
}

static std::string resolveDefaultTopologyAuditPrefix()
{
    const std::string catalogPath = resolveDefaultKnotSeedCatalogPath();
    const size_t slashPos = catalogPath.find_last_of("\\/");
    if (slashPos == std::string::npos)
        return "topology_audit";

    return catalogPath.substr(0, slashPos + 1) + "topology_audit";
}

static bool prepareTopologyAuditOutputFiles()
{
    if (g_topologyAuditOutputPrefix.empty())
        g_topologyAuditOutputPrefix = resolveDefaultTopologyAuditPrefix();

    const std::string seedsPath = g_topologyAuditOutputPrefix + "_seeds.tsv";
    const std::string eventsPath = g_topologyAuditOutputPrefix + "_events.tsv";
    const std::string summaryPath = g_topologyAuditOutputPrefix + "_summary.txt";

    g_topologyAuditSeedsOut.open(seedsPath, std::ios::out | std::ios::trunc);
    g_topologyAuditEventsOut.open(eventsPath, std::ios::out | std::ios::trunc);
    g_topologyAuditSummaryOut.open(summaryPath, std::ios::out | std::ios::trunc);

    if (!g_topologyAuditSeedsOut || !g_topologyAuditEventsOut || !g_topologyAuditSummaryOut)
    {
        std::cout << "[topology-audit] could not open output files under prefix "
                  << g_topologyAuditOutputPrefix << std::endl;
        return false;
    }

    g_topologyAuditSeedsOut
        << "seed\tknotType\tphase\tsteps\tcutApplied\tcutPiece\tcutEdge"
        << "\tviolations\tcenterlineContacts\tsurfaceOverlaps\tcrossingsPerStep"
        << "\tminAbsHeight"
        << "\tcandidatePairs\tpenetratingContacts\tmaxPenetration\tmeanPenetration"
        << "\tresidualCandidatePairs\tresidualPenetratingContacts"
        << "\tresidualMaxPenetration\tresidualMeanPenetration"
        << "\tmaxGroundPenetration\tcappedCorrections\tcorrectionEndpoints"
        << "\tmaxAppliedCorrection\tnormalFlips"
        << "\tfirstViolationStep\tfirstViolationType"
        << "\tfirstPieceA\tfirstEdgeA\tfirstPieceB\tfirstEdgeB"
        << "\tfirstPreviousHeight\tfirstCurrentHeight\n";

    g_topologyAuditEventsOut
        << "seed\tknotType\tphase\tstep\ttype"
        << "\tpieceA\tedgeA\tpieceB\tedgeB\tx\tz"
        << "\tpreviousX\tpreviousZ"
        << "\tpreviousHeight\tcurrentHeight"
        << "\tcandidatePairs\tpenetratingContacts\tmaxPenetration"
        << "\tresidualPenetratingContacts\tresidualMaxPenetration"
        << "\tmaxGroundPenetration\tnormalFlips\n";

    if (g_topologyAuditTraceEnabled)
    {
        const std::string tracePath = g_topologyAuditOutputPrefix + "_trace.tsv";
        g_topologyAuditTraceOut.open(tracePath, std::ios::out | std::ios::trunc);
        if (!g_topologyAuditTraceOut)
        {
            std::cout << "[topology-audit] could not open trace file "
                      << tracePath << std::endl;
            return false;
        }

        g_topologyAuditTraceOut
            << "seed\tknotType\tphase\tstep\tcrossings"
            << "\tcenterlineContacts\tsurfaceOverlaps\tminAbsHeight\tviolations"
            << "\tcandidatePairs\tpenetratingContacts\tmaxPenetration\tmeanPenetration"
            << "\tresidualCandidatePairs\tresidualPenetratingContacts"
            << "\tresidualMaxPenetration\tresidualMeanPenetration"
            << "\tmaxGroundPenetration\tcappedCorrections\tcorrectionEndpoints"
            << "\tmaxAppliedCorrection\tnormalFlips\n";
    }

    std::cout << "[topology-audit] writing " << seedsPath
              << ", " << eventsPath << ", " << summaryPath;
    if (g_topologyAuditTraceEnabled)
        std::cout << ", " << g_topologyAuditOutputPrefix << "_trace.tsv";
    std::cout << std::endl;
    return true;
}

static void writeTopologyAuditEvent(
    const CutEvalSeedRecord& record,
    const std::string& phase,
    unsigned int step,
    const TopologyAudit::Violation& violation)
{
    if (!g_topologyAuditEventsOut)
        return;

    const TopologyAudit::Crossing& crossing =
        violation.hasCurrent ? violation.current : violation.previous;

    g_topologyAuditEventsOut << std::fixed << std::setprecision(9)
        << record.seed << '\t'
        << record.knotType << '\t'
        << phase << '\t'
        << step << '\t'
        << violation.type << '\t'
        << crossing.pieceA << '\t'
        << crossing.edgeA << '\t'
        << crossing.pieceB << '\t'
        << crossing.edgeB << '\t'
        << crossing.x << '\t'
        << crossing.z << '\t'
        << (violation.hasPrevious ? violation.previous.x : static_cast<Real>(-1.0)) << '\t'
        << (violation.hasPrevious ? violation.previous.z : static_cast<Real>(-1.0)) << '\t'
        << violation.previousHeight << '\t'
        << violation.currentHeight << '\t'
        << g_lastSubstepCollision.candidatePairCount << '\t'
        << g_lastSubstepCollision.penetratingContactCount << '\t'
        << g_lastSubstepCollision.maxPenetration << '\t'
        << g_lastSubstepCollision.residualPenetratingContactCount << '\t'
        << g_lastSubstepCollision.residualMaxPenetration << '\t'
        << g_lastSubstepGroundPenetration << '\t'
        << g_lastSubstepCollision.normalFlipViolationCount << '\n';
}

static void beginTopologyAuditPhase(
    const std::string& phase,
    const std::vector<std::vector<Vector3r>>& pieces)
{
    if (phase == "pre_cut")
    {
        g_topologyAuditCurrentSeedHadViolation = false;
        ++g_topologyAuditTotalSeeds;
    }

    g_topologyPhase = TopologyAuditPhaseSummary();
    g_topologyPhase.phase = phase;
    g_topologyAuditStepInPhase = 0;
    g_topologyAuditor.begin(pieces, g_topologyAuditSettings);
    g_topologyAuditActive = true;

    const CutEvalSeedRecord& record = g_evalSeeds[g_evalCurrentIndex];
    std::cout << "[topology-audit] seed=" << record.seed
              << " knotType=" << record.knotType
              << " phase=" << phase
              << " steps=" << g_topologyAuditSteps
              << " cutStep=" << g_topologyAuditCutsApplied
              << std::endl;
}

static void finishTopologyAuditPhaseRow(const CutEvalSeedRecord& record)
{
    if (!g_topologyAuditSeedsOut)
        return;

    const Real crossingsPerStep =
        g_topologyPhase.steps == 0
            ? static_cast<Real>(0.0)
            : static_cast<Real>(g_topologyPhase.crossings) /
                  static_cast<Real>(g_topologyPhase.steps);

    const Real meanPenetration =
        (g_topologyPhase.penetratingContactCount > 0)
            ? g_topologyPhase.sumPenetration /
                  static_cast<Real>(g_topologyPhase.penetratingContactCount)
            : static_cast<Real>(0.0);
    const Real residualMeanPenetration =
        (g_topologyPhase.residualPenetratingContactCount > 0)
            ? g_topologyPhase.residualSumPenetration /
                  static_cast<Real>(g_topologyPhase.residualPenetratingContactCount)
            : static_cast<Real>(0.0);

    g_topologyAuditSeedsOut << std::fixed << std::setprecision(9)
        << record.seed << '\t'
        << record.knotType << '\t'
        << g_topologyPhase.phase << '\t'
        << g_topologyPhase.steps << '\t'
        << (g_topologyPhase.cutApplied ? 1 : 0) << '\t'
        << static_cast<long long>(g_topologyPhase.cutPiece) << '\t'
        << g_topologyPhase.cutEdge << '\t'
        << g_topologyPhase.violations << '\t'
        << g_topologyPhase.centerlineContacts << '\t'
        << g_topologyPhase.surfaceOverlaps << '\t'
        << crossingsPerStep << '\t'
        << (std::isfinite(g_topologyPhase.minAbsHeight)
                ? g_topologyPhase.minAbsHeight
                : static_cast<Real>(-1.0))
        << '\t'
        << g_topologyPhase.candidatePairCount << '\t'
        << g_topologyPhase.penetratingContactCount << '\t'
        << g_topologyPhase.maxPenetration << '\t'
        << meanPenetration << '\t'
        << g_topologyPhase.residualCandidateCount << '\t'
        << g_topologyPhase.residualPenetratingContactCount << '\t'
        << g_topologyPhase.residualMaxPenetration << '\t'
        << residualMeanPenetration << '\t'
        << g_topologyPhase.maxGroundPenetration << '\t'
        << g_topologyPhase.cappedCorrectionCount << '\t'
        << g_topologyPhase.correctionEndpointCount << '\t'
        << g_topologyPhase.maxAppliedCorrection << '\t'
        << g_topologyPhase.normalFlipViolationCount
        << '\t';

    if (g_topologyPhase.hasFirstViolation)
    {
        const TopologyAudit::Violation& first = g_topologyPhase.firstViolation;
        const TopologyAudit::Crossing& crossing =
            first.hasCurrent ? first.current : first.previous;
        g_topologyAuditSeedsOut
            << first.step << '\t'
            << first.type << '\t'
            << crossing.pieceA << '\t'
            << crossing.edgeA << '\t'
            << crossing.pieceB << '\t'
            << crossing.edgeB << '\t'
            << first.previousHeight << '\t'
            << first.currentHeight;
    }
    else
    {
        g_topologyAuditSeedsOut
            << -1 << '\t'
            << "none" << '\t'
            << -1 << '\t'
            << -1 << '\t'
            << -1 << '\t'
            << -1 << '\t'
            << 0.0 << '\t'
            << 0.0;
    }

    g_topologyAuditSeedsOut << '\n';
}

static bool auditTopologyStep(
    const std::vector<std::vector<Vector3r>>& pieces)
{
    if (!g_topologyAuditActive)
        return false;

    const CutEvalSeedRecord& record = g_evalSeeds[g_evalCurrentIndex];
    const TopologyAudit::StepResult result =
        g_topologyAuditor.audit(g_topologyAuditStepInPhase, pieces);

    ++g_topologyPhase.steps;
    g_topologyPhase.crossings += result.crossingCount;
    g_topologyPhase.violations += result.violations.size();
    g_topologyPhase.centerlineContacts += result.centerlineContactCount;
    g_topologyPhase.surfaceOverlaps += result.surfaceOverlapCount;
    if (std::isfinite(result.minAbsHeight))
        g_topologyPhase.minAbsHeight =
            std::min(g_topologyPhase.minAbsHeight, result.minAbsHeight);

    // Solver/contact diagnostics captured for the substep that just executed.
    const SelfCollisionDiagnostics& diag = g_lastSubstepCollision;
    g_topologyPhase.candidatePairCount += diag.candidatePairCount;
    g_topologyPhase.penetratingContactCount += diag.penetratingContactCount;
    g_topologyPhase.maxPenetration = std::max(g_topologyPhase.maxPenetration, diag.maxPenetration);
    g_topologyPhase.sumPenetration += diag.sumPenetration;
    g_topologyPhase.correctionEndpointCount += diag.correctionEndpointCount;
    g_topologyPhase.cappedCorrectionCount += diag.cappedCorrectionCount;
    g_topologyPhase.maxAppliedCorrection =
        std::max(g_topologyPhase.maxAppliedCorrection, diag.maxAppliedCorrection);
    g_topologyPhase.normalFlipViolationCount += diag.normalFlipViolationCount;
    g_topologyPhase.residualCandidateCount += diag.residualCandidateCount;
    g_topologyPhase.residualPenetratingContactCount += diag.residualPenetratingContactCount;
    g_topologyPhase.residualMaxPenetration =
        std::max(g_topologyPhase.residualMaxPenetration, diag.residualMaxPenetration);
    g_topologyPhase.residualSumPenetration += diag.residualSumPenetration;
    g_topologyPhase.maxGroundPenetration =
        std::max(g_topologyPhase.maxGroundPenetration, g_lastSubstepGroundPenetration);

    ++g_topologyAuditTotalSteps;
    g_topologyAuditTotalViolations += result.violations.size();
    g_topologyAuditTotalCenterlineContacts += result.centerlineContactCount;
    g_topologyAuditTotalSurfaceOverlaps += result.surfaceOverlapCount;
    g_topologyAuditTotalCandidatePairs += diag.candidatePairCount;
    g_topologyAuditTotalPenetratingContacts += diag.penetratingContactCount;
    g_topologyAuditMaxPenetration = std::max(g_topologyAuditMaxPenetration, diag.maxPenetration);
    g_topologyAuditSumPenetration += diag.sumPenetration;
    g_topologyAuditTotalCappedCorrections += diag.cappedCorrectionCount;
    g_topologyAuditTotalCorrectionEndpoints += diag.correctionEndpointCount;
    g_topologyAuditMaxAppliedCorrection =
        std::max(g_topologyAuditMaxAppliedCorrection, diag.maxAppliedCorrection);
    g_topologyAuditTotalNormalFlips += diag.normalFlipViolationCount;
    g_topologyAuditTotalResidualCandidates += diag.residualCandidateCount;
    g_topologyAuditTotalResidualPenetratingContacts += diag.residualPenetratingContactCount;
    g_topologyAuditResidualMaxPenetration =
        std::max(g_topologyAuditResidualMaxPenetration, diag.residualMaxPenetration);
    g_topologyAuditResidualSumPenetration += diag.residualSumPenetration;
    g_topologyAuditMaxGroundPenetration =
        std::max(g_topologyAuditMaxGroundPenetration, g_lastSubstepGroundPenetration);

    if (!result.violations.empty())
    {
        g_topologyAuditCurrentSeedHadViolation = true;
        if (!g_topologyPhase.hasFirstViolation)
        {
            g_topologyPhase.hasFirstViolation = true;
            g_topologyPhase.firstViolation = result.violations.front();
        }

        for (const TopologyAudit::Violation& violation : result.violations)
            writeTopologyAuditEvent(record, g_topologyPhase.phase, g_topologyAuditStepInPhase, violation);
    }

    if (g_topologyAuditTraceOut)
    {
        const Real traceMeanPen = (diag.penetratingContactCount > 0)
            ? diag.sumPenetration / static_cast<Real>(diag.penetratingContactCount)
            : static_cast<Real>(0.0);
        const Real traceResidualMeanPen = (diag.residualPenetratingContactCount > 0)
            ? diag.residualSumPenetration / static_cast<Real>(diag.residualPenetratingContactCount)
            : static_cast<Real>(0.0);
        g_topologyAuditTraceOut << std::fixed << std::setprecision(9)
            << record.seed << '\t'
            << record.knotType << '\t'
            << g_topologyPhase.phase << '\t'
            << g_topologyAuditStepInPhase << '\t'
            << result.crossingCount << '\t'
            << result.centerlineContactCount << '\t'
            << result.surfaceOverlapCount << '\t'
            << (std::isfinite(result.minAbsHeight)
                    ? result.minAbsHeight
                    : static_cast<Real>(-1.0))
            << '\t'
            << result.violations.size() << '\t'
            << diag.candidatePairCount << '\t'
            << diag.penetratingContactCount << '\t'
            << diag.maxPenetration << '\t'
            << traceMeanPen << '\t'
            << diag.residualCandidateCount << '\t'
            << diag.residualPenetratingContactCount << '\t'
            << diag.residualMaxPenetration << '\t'
            << traceResidualMeanPen << '\t'
            << g_lastSubstepGroundPenetration << '\t'
            << diag.cappedCorrectionCount << '\t'
            << diag.correctionEndpointCount << '\t'
            << diag.maxAppliedCorrection << '\t'
            << diag.normalFlipViolationCount << '\n';
    }

    ++g_topologyAuditStepInPhase;
    if (g_topologyAuditStepInPhase < g_topologyAuditSteps)
        return false;

    g_topologyAuditActive = false;
    finishTopologyAuditPhaseRow(record);
    ++g_topologyAuditTotalPhases;

    std::cout << "[topology-audit] seed=" << record.seed
              << " phase=" << g_topologyPhase.phase
              << " steps=" << g_topologyPhase.steps
              << " violations=" << g_topologyPhase.violations
              << " contacts=" << g_topologyPhase.centerlineContacts
              << " minAbsHeight=" << g_topologyPhase.minAbsHeight
              << " maxPenetration=" << g_topologyPhase.maxPenetration
              << " residualMaxPenetration=" << g_topologyPhase.residualMaxPenetration
              << std::endl;
    return true;
}

static bool applyTopologyAuditCut()
{
    SimulationModel* model = Simulation::getCurrent()->getModel();
    if (model == nullptr)
        return false;

    std::vector<std::vector<Vector3r>> pieces = collectCurrentCablePieces(model);
    if (pieces.empty())
        return false;

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
        std::cout << "[topology-audit] seed=" << g_evalSeeds[g_evalCurrentIndex].seed
                  << " has no remaining threading; no further cut applied" << std::endl;
        return false;
    }

    const VirtualCutTarget target = findBestHoopMidpointCut(tanglement, pieces);
    if (!target.valid ||
        !splitCablePieceAtEdge(pieces, target.pieceIndex, target.edgeIndex))
    {
        std::cout << "[topology-audit] no valid heuristic cut for seed="
                  << g_evalSeeds[g_evalCurrentIndex].seed
                  << " (threading=" << threading << ")" << std::endl;
        return false;
    }

    g_cablePieces = std::move(pieces);
    g_cachedCableKymoKnotResults.clear();
    g_hasCableKymoKnotResults = false;

    ++g_topologyAuditCutsApplied;
    std::cout << "[topology-audit] cut seed=" << g_evalSeeds[g_evalCurrentIndex].seed
              << " piece=" << target.pieceIndex
              << " edge=" << target.edgeIndex
              << " cutStep=" << g_topologyAuditCutsApplied << std::endl;

    model->cleanup();
    buildModel();

    const std::vector<std::vector<Vector3r>> afterCut =
        collectCurrentCablePieces(model);
    const std::string phase =
        "post_cut_" + std::to_string(g_topologyAuditCutsApplied);
    beginTopologyAuditPhase(phase, afterCut);
    g_topologyPhase.cutApplied = true;
    g_topologyPhase.cutPiece = target.pieceIndex;
    g_topologyPhase.cutEdge = target.edgeIndex;
    return true;
}

static void writeTopologyAuditSummary()
{
    if (!g_topologyAuditSummaryOut)
        return;

    const Real meanPenetration =
        (g_topologyAuditTotalPenetratingContacts > 0)
            ? g_topologyAuditSumPenetration /
                  static_cast<Real>(g_topologyAuditTotalPenetratingContacts)
            : static_cast<Real>(0.0);
    const Real residualMeanPenetration =
        (g_topologyAuditTotalResidualPenetratingContacts > 0)
            ? g_topologyAuditResidualSumPenetration /
                  static_cast<Real>(g_topologyAuditTotalResidualPenetratingContacts)
            : static_cast<Real>(0.0);

    g_topologyAuditSummaryOut << std::setprecision(17)
        << "seeds=" << g_topologyAuditTotalSeeds << '\n'
        << "phases=" << g_topologyAuditTotalPhases << '\n'
        << "steps=" << g_topologyAuditTotalSteps << '\n'
        << "violations=" << g_topologyAuditTotalViolations << '\n'
        << "centerline_contacts=" << g_topologyAuditTotalCenterlineContacts << '\n'
        << "surface_overlaps=" << g_topologyAuditTotalSurfaceOverlaps << '\n'
        << "seeds_with_violations=" << g_topologyAuditSeedsWithViolations << '\n'
        << "candidate_pairs=" << g_topologyAuditTotalCandidatePairs << '\n'
        << "penetrating_contacts=" << g_topologyAuditTotalPenetratingContacts << '\n'
        << "max_penetration=" << g_topologyAuditMaxPenetration << '\n'
        << "mean_penetration=" << meanPenetration << '\n'
        << "residual_candidate_pairs=" << g_topologyAuditTotalResidualCandidates << '\n'
        << "residual_penetrating_contacts=" << g_topologyAuditTotalResidualPenetratingContacts << '\n'
        << "residual_max_penetration=" << g_topologyAuditResidualMaxPenetration << '\n'
        << "residual_mean_penetration=" << residualMeanPenetration << '\n'
        << "max_ground_penetration=" << g_topologyAuditMaxGroundPenetration << '\n'
        << "capped_corrections=" << g_topologyAuditTotalCappedCorrections << '\n'
        << "correction_endpoints=" << g_topologyAuditTotalCorrectionEndpoints << '\n'
        << "max_applied_correction=" << g_topologyAuditMaxAppliedCorrection << '\n'
        << "normal_flips=" << g_topologyAuditTotalNormalFlips << '\n'
        << "collision_iterations=" << g_cfg.selfCollision.iterations << '\n'
        << "collision_compliance=" << static_cast<double>(g_cfg.selfCollision.compliance) << '\n'
        << "collision_max_correction=" << static_cast<double>(g_cfg.selfCollision.maxCorrection) << '\n'
        << "normal_damping=" << static_cast<double>(g_cfg.selfCollision.normalDamping) << '\n'
        << "tangential_damping=" << static_cast<double>(g_cfg.selfCollision.tangentialDamping) << '\n'
        << "hash_cell_size=" << static_cast<double>(g_cfg.selfCollision.hashCellSize) << '\n'
        << "dt=" << static_cast<double>(g_cfg.physics.timeStep) << '\n'
        << "creep=" << static_cast<double>(g_cfg.restCreep.creep) << '\n';
}

static void finishTopologyAuditRun()
{
    g_topologyAuditActive = false;
    g_topologyAuditRunEnabled = false;
    writeTopologyAuditSummary();
    std::cout << "[topology-audit] finished seeds=" << g_topologyAuditTotalSeeds
              << " phases=" << g_topologyAuditTotalPhases
              << " steps=" << g_topologyAuditTotalSteps
              << " violations=" << g_topologyAuditTotalViolations
              << " seedsWithViolations=" << g_topologyAuditSeedsWithViolations
              << " maxPenetration=" << g_topologyAuditMaxPenetration
              << " residualMaxPenetration=" << g_topologyAuditResidualMaxPenetration
              << " maxGroundPenetration=" << g_topologyAuditMaxGroundPenetration
              << std::endl;
    g_stopRequested = true;
}

static void finishTopologyAuditSeed()
{
    g_topologyAuditActive = false;

    if (g_topologyAuditCurrentSeedHadViolation)
        ++g_topologyAuditSeedsWithViolations;

    ++g_evalCurrentIndex;
    g_topologyAuditCutsApplied = 0;

    if (g_evalCurrentIndex >= g_evalSeeds.size())
    {
        finishTopologyAuditRun();
        return;
    }

    std::cout << "[topology-audit] next seed=" << g_evalSeeds[g_evalCurrentIndex].seed
              << " (" << (g_evalCurrentIndex + 1) << "/" << g_evalSeeds.size() << ")"
              << std::endl;
    generateNextCutEvalSeed();
}

int main(int argc, char** argv)
{
    unsigned int seed = 0;
    bool haveSeed = false;
    std::string outDir;
    unsigned int randomTrials = 10;
    unsigned int ctzMaxCuts = 10;
    unsigned int ctzRandomTrials = 20;
    unsigned int auditSteps = 300;
    unsigned int auditCuts = 3;
    bool runOracle = true;
    bool runAudit = true;
    bool auditTrace = false;

    g_cfg.restCreep.enabled = true;
    g_cfg.restCreep.creep = static_cast<Real>(0.02);
    g_cfg.restCreep.evaluationDelaySteps = 300;
    g_cfg.restCreep.activationThreshold = static_cast<Real>(1e-4);
    g_topologyAuditSettings = TopologyAudit::Settings();

    for (int i = 1; i < argc; ++i)
    {
        const std::string arg(argv[i]);
        if (arg.rfind("--seed=", 0) == 0)
        {
            seed = static_cast<unsigned int>(std::stoul(arg.substr(7)));
            haveSeed = seed != 0;
        }
        else if (arg.rfind("--out=", 0) == 0)
        {
            outDir = arg.substr(6);
        }
        else if (arg.rfind("--creep=", 0) == 0)
        {
            g_cfg.restCreep.creep = static_cast<Real>(std::max(
                static_cast<double>(0.0),
                std::min(static_cast<double>(1.0), std::stod(arg.substr(8)))));
        }
        else if (arg.rfind("--creep-settle-steps=", 0) == 0)
        {
            g_cfg.restCreep.evaluationDelaySteps =
                static_cast<unsigned int>(std::stoul(arg.substr(21)));
        }
        else if (arg.rfind("--creep-threshold=", 0) == 0)
        {
            g_cfg.restCreep.activationThreshold = static_cast<Real>(std::stod(arg.substr(18)));
        }
        else if (arg.rfind("--dt=", 0) == 0)
        {
            g_cfg.physics.timeStep = static_cast<Real>(std::stod(arg.substr(5)));
        }
        else if (arg.rfind("--random-trials=", 0) == 0)
        {
            randomTrials = static_cast<unsigned int>(std::stoul(arg.substr(16)));
        }
        else if (arg.rfind("--ctz-max-cuts=", 0) == 0)
        {
            ctzMaxCuts = static_cast<unsigned int>(std::stoul(arg.substr(15)));
        }
        else if (arg.rfind("--ctz-random-trials=", 0) == 0)
        {
            ctzRandomTrials = static_cast<unsigned int>(std::stoul(arg.substr(20)));
        }
        else if (arg.rfind("--ctz-oracle-max-nodes=", 0) == 0)
        {
            g_evalCutsToZeroOracleMaxNodes = std::stoull(arg.substr(23));
        }
        else if (arg.rfind("--audit-steps=", 0) == 0)
        {
            auditSteps = static_cast<unsigned int>(std::stoul(arg.substr(14)));
        }
        else if (arg.rfind("--audit-cuts=", 0) == 0)
        {
            auditCuts = static_cast<unsigned int>(std::stoul(arg.substr(13)));
        }
        else if (arg.rfind("--audit-contact-eps=", 0) == 0)
        {
            g_topologyAuditSettings.contactEpsilon = static_cast<Real>(std::stod(arg.substr(20)));
        }
        else if (arg == "--no-oracle")
        {
            runOracle = false;
        }
        else if (arg == "--no-audit")
        {
            runAudit = false;
        }
        else if (arg == "--audit-trace")
        {
            auditTrace = true;
        }
        else if (arg == "--help" || arg == "-h")
        {
            std::cout
                << "usage: cablegen --seed=N [--out=DIR]\n"
                << "  [--creep=V] [--creep-settle-steps=N] [--creep-threshold=V]\n"
                << "  [--dt=V] [--random-trials=N]\n"
                << "  [--ctz-max-cuts=N] [--ctz-random-trials=N]\n"
                << "  [--ctz-oracle-max-nodes=N]\n"
                << "  [--audit-steps=N] [--audit-cuts=N] [--audit-contact-eps=V]\n"
                << "  [--no-oracle] [--no-audit] [--audit-trace]\n";
            return 0;
        }
        else
        {
            std::cerr << "unknown argument: " << arg << "\n";
            return 2;
        }
    }

    if (!haveSeed)
    {
        std::cerr << "error: --seed=N (non-zero) is required\n";
        return 2;
    }

    const auto startedAt = std::chrono::system_clock::now();

    if (outDir.empty())
        outDir = "results/seeds/seed_" + std::to_string(seed);
    std::filesystem::create_directories(outDir);

    const std::string catalogPath = outDir + "/catalog.tsv";
    const std::string cutEvalPath = outDir + "/cut_eval.tsv";
    const std::string cutEvalStagePath = outDir + "/cut_eval_stage_metrics.tsv";
    const std::string ctzPath = outDir + "/cuts_to_zero.tsv";
    const std::string ctzStagePath = outDir + "/cuts_to_zero_stage_metrics.tsv";
    const std::string auditPrefix = outDir + "/topology_audit";

    g_cfg.randomCable.enabled = true;
    g_cfg.randomCable.seed = seed;
    KymoKnotAnalyzer::setRandomSeed(seed);
    g_knotSeedBankPath = outDir + "/seed_bank.txt";
    g_knotSeedCatalogPath = catalogPath;
    {
        std::ofstream truncate(catalogPath, std::ios::trunc);
    }

    SimulationModel* model = new SimulationModel();
    model->init();
    Simulation::getCurrent()->setModel(model);
    TimeManager::getCurrent()->setTimeStepSize(g_cfg.physics.timeStep);
    captureGravityReference();
    applyGravityState();

    auto resetStage = [&]()
    {
        g_stopRequested = false;
        g_autoKnotEvaluationDone = false;
        reset();
    };

    auto runSettleStage = [&]()
    {
        size_t safety = 0;
        while (!g_autoKnotEvaluationDone && !g_stopRequested && safety++ < 10000000)
            timeStep();
        g_stopRequested = false;
    };

    bool ok = true;

    // 1. Catalog / tanglement metrics.
    g_catalogRunEnabled = true;
    g_catalogStepsPerFrame = 1;
    g_catalogSamplesWritten = 0;
    resetStage();
    std::cout << "[cablegen] catalog seed=" << seed << " output=" << catalogPath << std::endl;
    runSettleStage();
    g_catalogRunEnabled = false;
    ok = ok && g_catalogSamplesWritten >= 1 && std::filesystem::exists(catalogPath);

    if (ok && !loadCutEvalCatalog())
        ok = false;
    if (!ok)
    {
        std::cerr << "[cablegen] catalog stage failed for seed " << seed << std::endl;
    }

    // 2. Single-cut evaluation (with exhaustive oracle).
    if (ok)
    {
        g_evalOutputPath = cutEvalPath;
        g_evalStageMetricsOutputPath = cutEvalStagePath;
        g_evalRunEnabled = true;
        g_evalOracleEnabled = runOracle;
        g_evalRandomTrials = randomTrials;
        g_evalSamplesWritten = 0;
        if (!prepareCutEvalOutputFile())
        {
            ok = false;
        }
        else
        {
            prepareCutEvalSeed(0);
            resetStage();
            std::cout << "[cablegen] cut-eval seed=" << seed << " output=" << cutEvalPath << std::endl;
            runSettleStage();
            g_evalRunEnabled = false;
            ok = ok && g_evalCurrentIndex >= g_evalSeeds.size()
                && std::filesystem::exists(cutEvalPath);
        }
    }

    // 3. Cuts-to-zero.
    if (ok)
    {
        g_evalCutsToZeroOutputPath = ctzPath;
        g_evalCutsToZeroStageMetricsOutputPath = ctzStagePath;
        g_evalCutsToZeroEnabled = true;
        g_evalCutsToZeroOracleEnabled = runOracle;
        g_evalCutsToZeroMaxCuts = ctzMaxCuts;
        g_evalCutsToZeroRandomTrials = ctzRandomTrials;
        g_evalSamplesWritten = 0;
        if (!prepareCutsToZeroEvalOutputFile())
        {
            ok = false;
        }
        else
        {
            prepareCutEvalSeed(0);
            resetStage();
            std::cout << "[cablegen] cuts-to-zero seed=" << seed << " output=" << ctzPath << std::endl;
            runSettleStage();
            g_evalCutsToZeroEnabled = false;
            g_evalCutsToZeroOracleEnabled = false;
            ok = ok && g_evalCurrentIndex >= g_evalSeeds.size()
                && std::filesystem::exists(ctzPath);
        }
    }

    // 4. Topology audit.
    if (ok && runAudit)
    {
        g_topologyAuditOutputPrefix = auditPrefix;
        g_topologyAuditSteps = auditSteps;
        g_topologyAuditMaxCuts = auditCuts;
        g_topologyAuditMaxSeeds = 0;
        g_topologyAuditTraceEnabled = auditTrace;
        g_topologyAuditRunEnabled = true;
        g_topologyAuditActive = false;
        g_topologyAuditCutsApplied = 0;
        g_topologyAuditTotalSeeds = 0;
        g_topologyAuditTotalPhases = 0;
        g_topologyAuditTotalSteps = 0;
        g_topologyAuditTotalViolations = 0;
        g_topologyAuditTotalCenterlineContacts = 0;
        g_topologyAuditTotalSurfaceOverlaps = 0;
        g_topologyAuditSeedsWithViolations = 0;
        g_topologyAuditCurrentSeedHadViolation = false;
        g_topologyAuditTotalCandidatePairs = 0;
        g_topologyAuditTotalPenetratingContacts = 0;
        g_topologyAuditMaxPenetration = static_cast<Real>(0.0);
        g_topologyAuditSumPenetration = static_cast<Real>(0.0);
        g_topologyAuditTotalCappedCorrections = 0;
        g_topologyAuditTotalCorrectionEndpoints = 0;
        g_topologyAuditMaxAppliedCorrection = static_cast<Real>(0.0);
        g_topologyAuditTotalNormalFlips = 0;
        g_topologyAuditTotalResidualCandidates = 0;
        g_topologyAuditTotalResidualPenetratingContacts = 0;
        g_topologyAuditResidualMaxPenetration = static_cast<Real>(0.0);
        g_topologyAuditResidualSumPenetration = static_cast<Real>(0.0);
        g_topologyAuditMaxGroundPenetration = static_cast<Real>(0.0);
        g_lastSubstepCollision = SelfCollisionDiagnostics();
        g_lastSubstepGroundPenetration = static_cast<Real>(0.0);
        g_topologyAuditor.reset();

        if (!prepareTopologyAuditOutputFiles())
        {
            ok = false;
        }
        else
        {
            prepareCutEvalSeed(0);
            resetStage();
            std::cout << "[cablegen] topology-audit seed=" << seed << " output=" << auditPrefix << std::endl;
            size_t safety = 0;
            while (!g_stopRequested && safety++ < 100000000)
                timeStep();
            g_topologyAuditRunEnabled = false;
            g_topologyAuditSeedsOut.flush();
            g_topologyAuditEventsOut.flush();
            g_topologyAuditSummaryOut.flush();
            if (g_topologyAuditTraceOut)
                g_topologyAuditTraceOut.flush();
            ok = ok && g_topologyAuditTotalSeeds >= 1
                && std::filesystem::exists(auditPrefix + "_seeds.tsv");
        }
    }

    std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
    std::time_t nowTime = std::chrono::system_clock::to_time_t(now);
    std::time_t startTime = std::chrono::system_clock::to_time_t(startedAt);
    std::ofstream meta(outDir + "/meta.tsv", std::ios::trunc);
    meta << "key\tvalue\n"
         << "version\t" << CABLEGEN_VERSION << "\n"
         << "seed\t" << seed << "\n"
         << "creep\t" << g_cfg.restCreep.creep << "\n"
         << "creep_settle_steps\t" << g_cfg.restCreep.evaluationDelaySteps << "\n"
         << "creep_threshold\t" << g_cfg.restCreep.activationThreshold << "\n"
         << "dt\t" << g_cfg.physics.timeStep << "\n"
         << "random_trials\t" << randomTrials << "\n"
         << "ctz_max_cuts\t" << ctzMaxCuts << "\n"
         << "ctz_random_trials\t" << ctzRandomTrials << "\n"
         << "ctz_oracle_max_nodes\t" << g_evalCutsToZeroOracleMaxNodes << "\n"
         << "audit_steps\t" << auditSteps << "\n"
         << "audit_cuts\t" << auditCuts << "\n"
         << "audit_contact_eps\t" << g_topologyAuditSettings.contactEpsilon << "\n"
         << "audit_trace\t" << (auditTrace ? 1 : 0) << "\n"
         << "oracle\t" << (runOracle ? 1 : 0) << "\n"
         << "audit\t" << (runAudit ? 1 : 0) << "\n"
         << "status\t" << (ok ? "complete" : "failed") << "\n"
         << "started_at\t"
         << std::put_time(std::localtime(&startTime), "%Y-%m-%dT%H:%M:%S") << "\n"
         << "finished_at\t"
         << std::put_time(std::localtime(&nowTime), "%Y-%m-%dT%H:%M:%S") << "\n";

    delete Simulation::getCurrent();
    delete model;
    return ok ? 0 : 1;
}


static void reset()
{
    g_cachedCableKymoKnotResults.clear();
    g_hasCableKymoKnotResults = false;
    g_cablePieces.clear();
    g_autoKnotEvaluationDone = false;
    g_restCreepCounter = 0;
    Simulation::getCurrent()->reset();

    Simulation::getCurrent()->getModel()->cleanup();
    buildModel();
    applyGravityState();
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

    Simulation::getCurrent()->setVecValue<Real>(
        Simulation::GRAVITATION, g_gravityReference.data());
}

static void initializeCablePieces()
{
    g_cablePieces.clear();
    std::vector<Vector3r> ropePoints =
        CableGeneration::generateRandomCableMessPolyline(
            g_cfg.layout,
            g_cfg.rope,
            g_cfg.randomCable,
            g_cfg.groundY);
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

    SimulationModel* model = Simulation::getCurrent()->getModel();
    ParticleData& particles = model->getParticles();

    const unsigned int stepsPerFrame = g_catalogStepsPerFrame;
    const Real dt = TimeManager::getCurrent()->getTimeStepSize();

    for (unsigned int step = 0; step < stepsPerFrame; ++step)
    {
        Simulation::getCurrent()->getTimeStep()->step(*model);

        // Per-substep contact/solver diagnostics. Penetration and normal-flip
        // counts accumulate across all collision iterations; residual fields are
        // measured after the final iteration, before the ground clamp.
        SelfCollisionDiagnostics substepCollision;
        Real substepGroundPenetration = static_cast<Real>(0.0);

        if (g_cfg.selfCollision.enabled)
        {
            std::vector<Vector3r> oldPositions(particles.getNumberOfParticles());
            for (unsigned int i = 0; i < particles.getNumberOfParticles(); ++i)
                oldPositions[i] = particles.getPosition(i);

            for (unsigned int k = 0; k < g_cfg.selfCollision.iterations; ++k)
            {
                SelfCollisionDiagnostics iterDiag;
                applySelfCollision(
                    model,
                    dt,
                    g_cablePieceRanges,
                    g_cfg.rope,
                    g_cfg.selfCollision,
                    &iterDiag);
                substepCollision.accumulate(iterDiag);
            }

            substepGroundPenetration =
                applyGroundConstraintToModel(model, g_cfg.groundY);

            recomputeVelocitiesFromPositions(
                model,
                oldPositions,
                dt,
                g_cfg.physics.velocityRecomputeBlend);
        }
        else
        {
            substepGroundPenetration =
                applyGroundConstraintToModel(model, g_cfg.groundY);
        }

        applyBendTwistRestStateCreep(model, g_cfg.restCreep, g_restCreepCounter);

        g_lastSubstepCollision = substepCollision;
        g_lastSubstepGroundPenetration = substepGroundPenetration;
        if (!g_autoKnotEvaluationDone &&
            g_cfg.restCreep.enabled &&
            g_restCreepCounter >=
                g_cfg.restCreep.activationDelaySteps + g_cfg.restCreep.evaluationDelaySteps)
        {
            g_autoKnotEvaluationDone = true;
            if (g_topologyAuditRunEnabled)
            {
                const std::vector<std::vector<Vector3r>> settledPieces =
                    collectCurrentCablePieces(model);
                beginTopologyAuditPhase("pre_cut", settledPieces);
            }
            else if (g_catalogRunEnabled)
            {
                evaluateAndRecordCatalogSample(false);
                g_catalogRunEnabled = false;
                std::cout << "[knot-catalog] wrote " << g_catalogSamplesWritten
                          << " sample(s), output=" << g_knotSeedCatalogPath << std::endl;
                g_stopRequested = true;
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
                    g_stopRequested = true;
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
                    g_stopRequested = true;
                    return;
                }

                generateNextCutEvalSeed();
                return;
            }
        }

        if (g_topologyAuditActive)
        {
            const std::vector<std::vector<Vector3r>> currentPieces =
                collectCurrentCablePieces(model);
            if (auditTopologyStep(currentPieces))
            {
                if (g_topologyAuditCutsApplied < g_topologyAuditMaxCuts)
                {
                    if (!applyTopologyAuditCut())
                        finishTopologyAuditSeed();
                }
                else
                {
                    finishTopologyAuditSeed();
                }
                return;
            }
        }
    }
}

