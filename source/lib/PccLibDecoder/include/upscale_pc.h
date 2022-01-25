#pragma once
#include "PCCMath.h"
#include "PCCCommon.h"
#include "PCCGroupOfFrames.h"
#include "PCCDecoderParameters.h"
#include "PCCPointSet.h"
#include <unordered_set>
using namespace pcc;
using namespace std;

void upscale_pcs( PCCGroupOfFrames& reconstructs, int frameCount, int downscale_factor );

void upscale_pc( PCCPointSet3& reconstruct, int downscale_factor );