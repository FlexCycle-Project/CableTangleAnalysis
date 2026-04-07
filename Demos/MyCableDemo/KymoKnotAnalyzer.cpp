#include "KymoKnotAnalyzer.h"

#include "Utils/FileSystem.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>

#ifdef _WIN32
#include <windows.h>
#endif

namespace PBD
{
	namespace
	{
		static std::string escapeForSingleQuotedBash(const std::string& value)
		{
			std::string escaped;
			escaped.reserve(value.size() + 8);
			for (char c : value)
			{
				if (c == '\'')
					escaped += "'\"'\"'";
				else
					escaped += c;
			}
			return escaped;
		}

		static std::string toWslPath(const std::string& path)
		{
			std::string normalized = Utilities::FileSystem::normalizePath(path);
			if (normalized.size() >= 2 && normalized[1] == ':')
			{
				const char drive = static_cast<char>(std::tolower(static_cast<unsigned char>(normalized[0])));
				return std::string("/mnt/") + drive + normalized.substr(2);
			}
			return normalized;
		}

		static bool isIntegerToken(const std::string& token)
		{
			if (token.empty())
				return false;

			size_t index = 0;
			if (token[0] == '-' || token[0] == '+')
				index = 1;
			if (index >= token.size())
				return false;

			for (; index < token.size(); ++index)
			{
				if (!std::isdigit(static_cast<unsigned char>(token[index])))
					return false;
			}
			return true;
		}

		static std::string getKymoKnotRoot()
		{
			const char* env = std::getenv("KYMOKNOT_ROOT");
			if (env != nullptr && env[0] != '\0')
				return Utilities::FileSystem::normalizePath(env);

			return "C:/Users/alpse/Desktop/KymoKnot-clone";
		}

		static bool writeLinearSnapshot(
			const std::vector<Vector3r>& points,
			const std::string& fileName)
		{
			std::ofstream out(fileName.c_str(), std::ios::out | std::ios::trunc);
			if (!out.is_open())
				return false;

			out << points.size() << "\n";
			out << std::fixed << std::setprecision(9);
			for (const Vector3r& p : points)
				out << p.x() << " " << p.y() << " " << p.z() << "\n";
			return true;
		}

		static std::vector<KymoKnotSegment> parseKymoKnotOutput(const std::string& fileName)
		{
			std::ifstream in(fileName.c_str());
			std::vector<KymoKnotSegment> segments;
			if (!in.is_open())
				return segments;

			std::string line;
			while (std::getline(in, line))
			{
				if (line.empty() || line[0] == '#')
					continue;

				std::istringstream iss(line);
				std::vector<std::string> tokens;
				std::string token;
				while (iss >> token)
					tokens.push_back(token);

				for (size_t i = 0; i + 3 < tokens.size(); ++i)
				{
					const std::string& knotToken = tokens[i];
					if (knotToken != "UN" && knotToken.find('_') == std::string::npos)
						continue;
					if (!isIntegerToken(tokens[i + 1]) ||
						!isIntegerToken(tokens[i + 2]) ||
						!isIntegerToken(tokens[i + 3]))
					{
						continue;
					}

					KymoKnotSegment segment;
					segment.knotType = knotToken;
					segment.start = std::atoi(tokens[i + 1].c_str());
					segment.end = std::atoi(tokens[i + 2].c_str());
					segment.length = std::atoi(tokens[i + 3].c_str());
					segments.push_back(segment);
				}
			}

			return segments;
		}

		static std::string buildRunCommand(
			const std::string& rootDir,
			const std::string& inputFile,
			const std::string& stdoutFile,
			const std::string& stderrFile,
			KymoKnotSearchMode searchMode)
		{
			const char* searchFlag =
				(searchMode == KymoKnotSearchMode::BottomUp) ? "-b" : "-t";
#ifdef _WIN32
			const std::string wslRoot = toWslPath(rootDir);
			const std::string wslInput = toWslPath(inputFile);
			const std::string wslStdout = toWslPath(stdoutFile);
			const std::string wslStderr = toWslPath(stderrFile);
			return
				"wsl.exe -d Ubuntu -- bash -lc \"cd '" + escapeForSingleQuotedBash(wslRoot) +
				"' && ./KymoKnot_linear.x " + searchFlag + " '" + escapeForSingleQuotedBash(wslInput) +
				"' > '" + escapeForSingleQuotedBash(wslStdout) +
				"' 2> '" + escapeForSingleQuotedBash(wslStderr) + "'\"";
#else
			return
				"cd '" + escapeForSingleQuotedBash(rootDir) +
				"' && ./KymoKnot_linear.x " + searchFlag + " '" + escapeForSingleQuotedBash(inputFile) +
				"' > '" + escapeForSingleQuotedBash(stdoutFile) +
				"' 2> '" + escapeForSingleQuotedBash(stderrFile) + "'";
#endif
		}

		static std::string getSearchOutputPrefix(KymoKnotSearchMode searchMode)
		{
			return (searchMode == KymoKnotSearchMode::BottomUp) ? "BU__" : "TD__";
		}
	}

	KymoKnotAnalysisResult KymoKnotAnalyzer::analyzeLinear(
		const std::vector<Vector3r>& points,
		const std::string& bridgeRoot,
		KymoKnotSearchMode searchMode)
	{
		KymoKnotAnalysisResult result;

		if (points.size() < 2)
		{
			result.message = "Need at least two cable points.";
			return result;
		}

		const std::string rootDir = getKymoKnotRoot();
		const std::string bridgeDir = Utilities::FileSystem::normalizePath(
			bridgeRoot.empty() ? (rootDir + "/kymoknot_bridge") : bridgeRoot);
		Utilities::FileSystem::makeDirs(bridgeDir);

		static unsigned int runCounter = 0;
		++runCounter;

		const std::string inputFile = bridgeDir + "/linear_input_" + std::to_string(runCounter) + ".dat";
		const std::string stdoutFile = bridgeDir + "/linear_run_" + std::to_string(runCounter) + ".out";
		const std::string stderrFile = bridgeDir + "/linear_run_" + std::to_string(runCounter) + ".err";
		const std::string outputFile =
			bridgeDir + "/" + getSearchOutputPrefix(searchMode) +
			"linear_input_" + std::to_string(runCounter) + ".dat";

		if (!writeLinearSnapshot(points, inputFile))
		{
			result.message = "Failed to write KymoKnot input snapshot.";
			return result;
		}

		const std::string command =
			buildRunCommand(rootDir, inputFile, stdoutFile, stderrFile, searchMode);
		const int exitCode = std::system(command.c_str());
		if (exitCode != 0 && !Utilities::FileSystem::isFile(outputFile))
		{
			result.message = "KymoKnot command failed.";
			return result;
		}

		result.segments = parseKymoKnotOutput(outputFile);
		result.available = true;
		result.success = !result.segments.empty();
		if (!result.success)
			result.message = "KymoKnot completed, but no segments were parsed.";
		return result;
	}
}
