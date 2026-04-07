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

	class KymoKnotAnalyzer
	{
	public:
		static KymoKnotAnalysisResult analyzeLinear(
			const std::vector<Vector3r>& points,
			const std::string& bridgeRoot,
			KymoKnotSearchMode searchMode = KymoKnotSearchMode::TopDown);
	};
}

#endif
