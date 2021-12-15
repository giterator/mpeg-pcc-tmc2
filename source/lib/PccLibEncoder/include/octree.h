#pragma once
#include"PCCMath.h"
#include "PCCCommon.h"
#include "PCCGroupOfFrames.h"
#include"PCCEncoderParameters.h"
#include"PCCPointSet.h"
#include<unordered_set>
using namespace pcc;
using namespace std;

int                             num_points_in_box( std::vector<pair<PCCPoint3D, PCCColor3B>>&,
                                                   std::vector<pair<PCCPoint3D, PCCColor3B>>&,
                                                   int,
                                                   int,
                                                   int,
                                                   int,
                                                   int,
                                                   int );
void                        octree_recurse_decomp( std::vector<pair<PCCPoint3D, PCCColor3B>>&,
                            std::vector<std::vector<Range>>&,
                            std::vector<pair<PCCPoint3D, PCCColor3B>>&,
                            int,
                            int,
                            int,
                            int,
                            int,
                            int,
                            int,
                            int,
                            bool );
void octree_decomp( const PCCPointSet3&              ,
                    std::vector<std::vector<Range>>& ,
                    PCCEncoderParameters&            ,
                    std::vector<int>&                ,
                    std::vector<int>&                ,
                    std::vector<int>&                ,
                    std::vector<int>&                ,
                    std::vector<int>&                ,
                    std::vector<int>&                );

std::vector<std::pair<PCCPoint3D, PCCColor3B>> get_octree_decomp_centroids( const PCCPointSet3&, PCCEncoderParameters& );

std::pair<PCCPoint3D, PCCColor3B> get_centroid( std::vector<pair<PCCPoint3D, PCCColor3B>>& );

void threeDD_voxel_grid_filter( PCCGroupOfFrames& , int , PCCEncoderParameters& );

