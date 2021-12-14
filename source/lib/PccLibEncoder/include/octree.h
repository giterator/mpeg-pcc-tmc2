#pragma once
#include"PCCMath.h"
#include "PCCCommon.h"
#include "PCCGroupOfFrames.h"
#include"PCCEncoderParameters.h"
#include"PCCPointSet.h"
#include<unordered_set>
using namespace pcc;

int num_points_in_box( std::vector<PCCPoint3D>&, std::vector<PCCPoint3D>&, int, int, int, int, int, int );
void octree_recurse_decomp( std::vector<PCCPoint3D>&,
                            std::vector<std::vector<Range>>&,
                            int,
                            int,
                            int,
                            int,
                            int,
                            int,
                            int,
                            int );
void octree_decomp( const PCCPointSet3&              ,
                    std::vector<std::vector<Range>>& ,
                    PCCEncoderParameters&            ,
                    std::vector<int>&                ,
                    std::vector<int>&                ,
                    std::vector<int>&                ,
                    std::vector<int>&                ,
                    std::vector<int>&                ,
                    std::vector<int>&                );
std::vector<std::vector<Range>> get_octree_decomp_chunks( const PCCPointSet3&, PCCEncoderParameters& );

std::pair<PCCPoint3D, PCCColor3B> get_centroid( const PCCPointSet3&, std::vector<int> );

void threeDD_voxel_grid_filter( PCCGroupOfFrames& , int , PCCEncoderParameters& );

