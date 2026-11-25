#pragma once

#include "CablePhysics.h"

#include <string>
#include <vector>

namespace CableGeneration
{
std::vector<Vector3r> generateStraightCablePolyline(
    const CableLayoutConfig& layoutCfg,
    const RopeConfig& ropeCfg);

std::vector<Vector3r> generateRandomCableMessPolyline(
    const CableLayoutConfig& layoutCfg,
    const RopeConfig& ropeCfg,
    const RandomCableConfig& randomCfg,
    Real groundY);

std::vector<Vector3r> generateMathematicalKnotPolyline(
    const std::string& knotName);

void smoothPolyline(
    std::vector<Vector3r>& points,
    int iterations,
    Real alpha);

std::vector<Vector3r> resamplePolylineByLength(
    const std::vector<Vector3r>& polyline,
    int wantedCount);
}
