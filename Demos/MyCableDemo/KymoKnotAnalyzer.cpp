#include "KymoKnotAnalyzer.h"

#include "Utils/FileSystem.h"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <chrono>

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

			const std::string bundledRoot = Utilities::FileSystem::normalizePath(
				Utilities::FileSystem::getProgramPath() + "/../KymoKnot-clone");
			if (Utilities::FileSystem::fileExists(bundledRoot + "/KymoKnot_linear.x"))
				return bundledRoot;

			return Utilities::FileSystem::normalizePath("KymoKnot-clone");
		}

		static std::string makeTempWorkDirectory()
		{
			const auto tick = std::chrono::high_resolution_clock::now().time_since_epoch().count();
			static std::atomic<unsigned int> counter{0};
			const unsigned int id = ++counter;
			std::string tempRoot;
#ifdef _WIN32
			char tempBuffer[MAX_PATH] = {0};
			const DWORD length = GetTempPathA(MAX_PATH, tempBuffer);
			if (length > 0 && length < MAX_PATH)
				tempRoot = tempBuffer;
#endif
			if (tempRoot.empty())
				tempRoot = "C:/Temp";

			const std::string workDir = Utilities::FileSystem::normalizePath(
				tempRoot + "/MyCableDemo/KymoKnot/run_" +
				std::to_string(tick) + "_" + std::to_string(id));
			Utilities::FileSystem::makeDirs(workDir);
			return workDir;
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
					if (isIntegerToken(tokens[i]))
						continue;
					if (!isIntegerToken(tokens[i + 1]) ||
						!isIntegerToken(tokens[i + 2]) ||
						!isIntegerToken(tokens[i + 3]))
					{
						continue;
					}

					KymoKnotSegment segment;
					segment.knotType = tokens[i];
					segment.start = std::atoi(tokens[i + 1].c_str());
					segment.end = std::atoi(tokens[i + 2].c_str());
					segment.length = std::atoi(tokens[i + 3].c_str());
					segments.push_back(segment);
					break;
				}
			}

			return segments;
		}

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

		static KymoKnotAnalysisResult analyzeWithWorkDir(
			const std::vector<Vector3r>& points,
			const std::string& bridgeDir,
			KymoKnotSearchMode searchMode)
		{
			KymoKnotAnalysisResult result;
			if (points.size() < 2)
			{
				result.message = "Need at least two cable points.";
				return result;
			}

			const std::string rootDir = getKymoKnotRoot();
			Utilities::FileSystem::makeDirs(bridgeDir);

			static unsigned int runCounter = 0;
			++runCounter;

			const std::string inputFile = bridgeDir + "/linear_input_" + std::to_string(runCounter) + ".dat";
			const std::string stdoutFile = bridgeDir + "/linear_run_" + std::to_string(runCounter) + ".out";
			const std::string stderrFile = bridgeDir + "/linear_run_" + std::to_string(runCounter) + ".err";
			const std::string outputFile =
				bridgeDir + "/" + getSearchOutputPrefix(searchMode) +
				"linear_input_" + std::to_string(runCounter) + ".dat";

			struct TempCleanup
			{
				std::string directory;
				std::string inputFile;
				std::string stdoutFile;
				std::string stderrFile;
				std::string outputFile;
				~TempCleanup()
				{
					if (directory.empty())
						return;
#ifdef _WIN32
					if (!inputFile.empty())
						DeleteFileA(inputFile.c_str());
					if (!stdoutFile.empty())
						DeleteFileA(stdoutFile.c_str());
					if (!stderrFile.empty())
						DeleteFileA(stderrFile.c_str());
					if (!outputFile.empty())
						DeleteFileA(outputFile.c_str());
					RemoveDirectoryA(directory.c_str());
#endif
				}
			} cleanup{bridgeDir};
			cleanup.inputFile = inputFile;
			cleanup.stdoutFile = stdoutFile;
			cleanup.stderrFile = stderrFile;
			cleanup.outputFile = outputFile;

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

	KymoKnotAnalysisResult KymoKnotAnalyzer::analyze(
		const CableGeometry& cable,
		KymoKnotSearchMode searchMode)
	{
		const std::string bridgeDir = makeTempWorkDirectory();
		return analyzeWithWorkDir(cable.points, bridgeDir, searchMode);
	}

	KymoKnotAnalysisResult KymoKnotAnalyzer::analyzeLinear(
		const std::vector<Vector3r>& points,
		KymoKnotSearchMode searchMode)
	{
		CableGeometry cable;
		cable.points = points;
		return analyze(cable, searchMode);
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
