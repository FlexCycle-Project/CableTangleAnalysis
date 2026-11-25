#pragma once

#include "CablePhysics.h"

#include <string>
#include <vector>

namespace CableGeneration
{

std::vector<Vector3r> generateRandomCableMessPolyline(
    const CableLayoutConfig& layoutCfg,
    const RopeConfig& ropeCfg,
    const RandomCableConfig& randomCfg,
    Real groundY);


void smoothPolyline(
    std::vector<Vector3r>& points,
    int iterations,
    Real alpha);

std::vector<Vector3r> resamplePolylineByLength(
    const std::vector<Vector3r>& polyline,
    int wantedCount);
}
