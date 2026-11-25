#include "KymoKnotAnalyzer.h"

#include "KymoKnotNativeBridge.h"

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

namespace PBD
{
	namespace
	{
		unsigned int g_kymoKnotSeed = 1;

		static std::string describeSegments(
			const std::vector<KymoKnotSegment>& segments,
			const std::vector<Vector3r>* points)
		{
			if (segments.empty())
				return "unrecognized";

			KymoKnotAnalysisResult tempResult;
			tempResult.available = true;
			tempResult.success = true;
			tempResult.segments = segments;
			if (KymoKnotAnalyzer::isUnknotResult(tempResult))
				return "unknot";

			if (points == nullptr || points->size() < 2)
				return "KymoKnot";

			const std::vector<Vector3r>& cablePoints = *points;
			const size_t pieceEdges = cablePoints.size() - 1;
			size_t pieceCoverage = 0;
			size_t pieceLongest = 0;
			std::ostringstream out;
			out << "KymoKnot ";
			for (size_t i = 0; i < segments.size(); ++i)
			{
				const KymoKnotSegment& segment = segments[i];
				const size_t segmentLength = static_cast<size_t>(std::max(0, segment.length));
				pieceCoverage += std::min(segmentLength, pieceEdges);
				pieceLongest = std::max(pieceLongest, segmentLength);
				if (i > 0)
					out << ", ";
				out << segment.knotType
				    << "[" << segment.start
				    << "-" << segment.end
				    << ":" << segment.length << "]";
			}
			const Real spanFraction = static_cast<Real>(pieceLongest) /
				static_cast<Real>(std::max<size_t>(1, pieceEdges));
			const Real coverageFraction = static_cast<Real>(pieceCoverage) /
				static_cast<Real>(std::max<size_t>(1, pieceEdges));
			out << " | span=" << spanFraction
			    << " coverage=" << coverageFraction;
			return out.str();
		}
	}

	void KymoKnotAnalyzer::setRandomSeed(unsigned int seed)
	{
		g_kymoKnotSeed = (seed == 0) ? 1u : seed;
	}

	KymoKnotAnalysisResult KymoKnotAnalyzer::analyze(
		const CableGeometry& cable,
		KymoKnotSearchMode searchMode)
	{
		return analyzeLinear(cable.points, searchMode);
	}

	KymoKnotAnalysisResult KymoKnotAnalyzer::analyzeLinear(
		const std::vector<Vector3r>& points,
		KymoKnotSearchMode searchMode)
	{
		KymoKnotAnalysisResult result;
		if (points.size() < 2)
		{
			result.message = "Need at least two cable points.";
			return result;
		}

		std::vector<double> xyz;
		xyz.reserve(points.size() * 3);
		for (const Vector3r& point : points)
		{
			xyz.push_back(static_cast<double>(point.x()));
			xyz.push_back(static_cast<double>(point.y()));
			xyz.push_back(static_cast<double>(point.z()));
		}

		const int mode =
			(searchMode == KymoKnotSearchMode::BottomUp) ? 1 : 0;

		KymoKnotNativeResult nativeResult;
		nativeResult.available = 0;
		nativeResult.success = 0;
		nativeResult.count = 0;
		nativeResult.segments = nullptr;

		const int produced = kymoknot_native_analyze_linear(
			xyz.data(),
			static_cast<int>(xyz.size() / 3),
			mode,
			g_kymoKnotSeed,
			&nativeResult);

		if (!produced)
		{
			result.message = "KymoKnot native analysis failed.";
			kymoknot_native_free_result(&nativeResult);
			return result;
		}

		result.available = nativeResult.available != 0;
		result.success = nativeResult.success != 0;
		result.segments.reserve(static_cast<size_t>(nativeResult.count));
		for (int i = 0; i < nativeResult.count; ++i)
		{
			const KymoKnotNativeSegment& nativeSegment = nativeResult.segments[i];
			KymoKnotSegment segment;
			segment.knotType = nativeSegment.knotType;
			segment.start = nativeSegment.start;
			segment.end = nativeSegment.end;
			segment.length = nativeSegment.length;
			result.segments.push_back(segment);
		}

		if (!result.success)
			result.message = "KymoKnot completed, but no segments were parsed.";

		kymoknot_native_free_result(&nativeResult);
		return result;
	}

	bool KymoKnotAnalyzer::isUnknotResult(const KymoKnotAnalysisResult& result)
	{
		if (!result.available || !result.success || result.segments.size() != 1)
			return false;

		const KymoKnotSegment& segment = result.segments.front();
		return segment.knotType == "UN" ||
			(segment.start < 0 && segment.end < 0 && segment.length < 0);
	}

	std::string KymoKnotAnalyzer::describeAnalysisResult(
		const KymoKnotAnalysisResult& result,
		const std::vector<Vector3r>* points)
	{
		if (!result.available)
			return result.message;
		if (result.segments.empty())
			return "unrecognized";
		if (isUnknotResult(result))
			return "unknot";
		return describeSegments(result.segments, points);
	}
}
