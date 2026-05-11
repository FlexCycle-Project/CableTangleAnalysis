#ifndef MY_CABLE_DEMO_KYMOKNOT_ANALYZER_H
#define MY_CABLE_DEMO_KYMOKNOT_ANALYZER_H

#include "Common/Common.h"

#include <string>
#include <vector>

namespace PBD
{
	enum class KymoKnotSearchMode
	{
		TopDown,
		BottomUp
	};

	struct KymoKnotSegment
	{
		std::string knotType;
		int start = -1;
		int end = -1;
		int length = -1;
	};

	struct KymoKnotAnalysisResult
	{
		bool available = false;
		bool success = false;
		std::string message;
		std::vector<KymoKnotSegment> segments;
	};

	struct CableGeometry
	{
		std::vector<Vector3r> points;
	};

	class KymoKnotAnalyzer
	{
	public:
		static KymoKnotAnalysisResult analyze(
			const CableGeometry& cable,
			KymoKnotSearchMode searchMode = KymoKnotSearchMode::TopDown);

		static KymoKnotAnalysisResult analyzeLinear(
			const std::vector<Vector3r>& points,
			KymoKnotSearchMode searchMode = KymoKnotSearchMode::TopDown);

		static bool isUnknotResult(const KymoKnotAnalysisResult& result);

		static std::string describeAnalysisResult(
			const KymoKnotAnalysisResult& result,
			const std::vector<Vector3r>* points = nullptr);
	};
}

#endif
