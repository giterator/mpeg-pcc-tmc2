#pragma once
#include"PCCMath.h"
#include "PCCCommon.h"
#include "PCCGroupOfFrames.h"
#include"PCCEncoderParameters.h"
#include"PCCPointSet.h"
#include<unordered_set>
using namespace pcc;
using namespace std;

void downscale_pc( PCCGroupOfFrames& sources, int frameCount, int downscale_factor );