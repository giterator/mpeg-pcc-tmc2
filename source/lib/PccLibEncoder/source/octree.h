/**
 * @file octree.cpp
 * @author Pranav Venkatram (e0552200@u.nus.edu)
 * @brief 
 * @date 2021-12-07
 * 
 * @copyright Copyright (c) 2021
 * 
 */



#include "PCCCommon.h"
#include<unordered_set>
using namespace pcc;

//overwrite params_.numROIs_ since used for patch generation and packing

int num_points_in_box( std::unordered_set<PCCPoint3D>& us_points_global,
                       std::unordered_set<PCCPoint3D>& us_points_local,
                       int x_min,
                       int x_max,
                       int y_min,
                       int y_max,
                       int z_min,
                       int z_max ) {
  int count      = 0;

  auto point = us_points_global.begin();
  while (point != us_points_global.end()) {
    auto x = (*point)[0];
    auto y = (*point)[1];
    auto z = (*point)[2];

    if ( x_min <= x && x <= x_max && y_min <= y && y <= y_max && z_min <= z && z <= z_max ) {
        count++; 
        us_points_local.insert( *point );
        point = us_points_global.erase( point );
    }
    else {
      point++;
    }
  }

  /*
  for ( auto point : us_points ) {
    auto x = point[0];
    auto y = point[1];
    auto z = point[2];

    if ( x_min <= x && x <= x_max && y_min <= y && y <= y_max && z_min <= z && z <= z_max ) { count++; }
  }*/

  return count;
}


void octree_recurse_decomp( std::unordered_set<PCCPoint3D>&  us_points,
                            std::vector<std::vector<Range>>& chunks,
                            int                              num_points,
                            int                              x_min,
                            int                              x_max,
                            int                              y_min,
                            int                              y_max,
                            int                              z_min,
                            int                              z_max,
                            int                              maxPointsPerVoxel ) {
  if ( num_points > maxPointsPerVoxel ) {
    int bisect_x = ( x_min + x_max ) / 2;
    int bisect_y = ( y_min + y_max ) / 2;
    int bisect_z = ( z_min + z_max ) / 2;

    std::unordered_set<PCCPoint3D> us_points_front_bottom_left;
    int front_bottom_left_count = 
        num_points_in_box( us_points, us_points_front_bottom_left, x_min, bisect_x, y_min, bisect_y, z_min, bisect_z );
    octree_recurse_decomp( us_points, chunks, front_bottom_left_count, x_min, bisect_x, y_min, bisect_y, z_min, bisect_z,
                           maxPointsPerVoxel );

    std::unordered_set<PCCPoint3D> us_points_front_top_left;
    int front_top_left_count = 
        num_points_in_box( us_points, us_points_front_top_left, x_min, bisect_x, bisect_y + 1, y_max, z_min, bisect_z );
    octree_recurse_decomp( us_points, chunks, front_top_left_count, x_min, bisect_x, bisect_y + 1, y_max, z_min,
                           bisect_z,
                           maxPointsPerVoxel );

    std::unordered_set<PCCPoint3D> us_points_front_bottom_right;
    int front_bottom_right_count = 
        num_points_in_box( us_points, us_points_front_bottom_right, bisect_x + 1, x_max, y_min,
                                                      bisect_y, z_min, bisect_z );
    octree_recurse_decomp( us_points, chunks, front_bottom_right_count, bisect_x + 1, x_max, y_min, bisect_y, z_min,
                           bisect_z, maxPointsPerVoxel );

    std::unordered_set<PCCPoint3D> us_points_front_top_right;
    int front_top_right_count = 
        num_points_in_box( us_points, us_points_front_top_right, bisect_x + 1, x_max,
                                                   bisect_y + 1, y_max, z_min, bisect_z );
    octree_recurse_decomp( us_points, chunks, front_top_right_count, bisect_x + 1, x_max, bisect_y + 1, y_max, z_min,
                           bisect_z, maxPointsPerVoxel );

    std::unordered_set<PCCPoint3D> us_points_back_bottom_left;
    int back_bottom_left_count = 
        num_points_in_box( us_points, us_points_back_bottom_left, x_min, bisect_x, y_min, bisect_y, bisect_z + 1, z_max );
    octree_recurse_decomp( us_points, chunks, back_bottom_left_count, x_min, bisect_x, y_min, bisect_y, bisect_z + 1,
                           z_max, maxPointsPerVoxel );

    std::unordered_set<PCCPoint3D> us_points_back_top_left;
    int back_top_left_count = 
        num_points_in_box( us_points, us_points_back_top_left, x_min, bisect_x, bisect_y + 1, y_max,
                                                 bisect_z + 1, z_max );
    octree_recurse_decomp( us_points, chunks, back_top_left_count, x_min, bisect_x, bisect_y + 1, y_max, bisect_z + 1,
                           z_max, maxPointsPerVoxel );

    std::unordered_set<PCCPoint3D> us_points_back_bottom_right;
    int back_bottom_right_count = 
        num_points_in_box( us_points, us_points_back_bottom_right, bisect_x + 1, x_max, y_min,
                                                     bisect_y, bisect_z + 1, z_max );
    octree_recurse_decomp( us_points, chunks, back_bottom_right_count, bisect_x + 1, x_max, y_min, bisect_y,
                           bisect_z + 1,
                           z_max, maxPointsPerVoxel );

    std::unordered_set<PCCPoint3D> us_points_back_top_right;
    int back_top_right_count = 
        num_points_in_box( us_points, us_points_back_top_right, bisect_x + 1, x_max, bisect_y + 1,
                                                  y_max, bisect_z, z_max );
    octree_recurse_decomp( us_points, chunks, back_top_right_count, bisect_x + 1, x_max, bisect_y + 1, y_max, bisect_z,
                           z_max, maxPointsPerVoxel );

  } else {
    Range rangeX = Range( x_min, x_max );
    Range rangeY = Range( y_min, y_max );
    Range rangeZ = Range( z_min, z_max );

    chunks.push_back( std::vector<Range>( { rangeX, rangeY, rangeZ } ) );
  }
}

/* Perform octree decomposition of PC. 
      1. Compute bounding box for PC
      2. Divide PC into ROIs using octree
      3. Overwrite:
        userParams.numROIs

        params.roiBoundingBoxMinX
        params.roiBoundingBoxMaxX
        params.roiBoundingBoxMinY
        params.roiBoundingBoxMaxY
        params.roiBoundingBoxMinZ
        params.roiBoundingBoxMaxZ

        chunks where each chunk is a vector of 3 Ranges (1 per axis). Each Range is a pair of values (min, max) corresponding to the 
        bounding box's min and max coordinates along a particular axis (same as the params.roiBoundingBox.... values)
*/
void octree_decomp(const PCCPointSet3&                 points,
                    std::vector<std::vector<Range>>&    chunks,
                    PCCEncoderParameters&               userParams,
                    std::vector<int>&                 roiBoundingBoxMinX,
                    std::vector<int>&                 roiBoundingBoxMaxX,
                    std::vector<int>&                 roiBoundingBoxMinY,
                    std::vector<int>&                 roiBoundingBoxMaxY,
                    std::vector<int>&                 roiBoundingBoxMinZ,
                    std::vector<int>&                 roiBoundingBoxMaxZ ) {
    
    //Bounding Box for PC
    int x_min;
    int x_max;
    int y_min;
    int y_max;
    int z_min;
    int z_max;
    x_min = y_min = z_min = ( std::numeric_limits<int>::max )();
    x_max = y_max = z_max = ( std::numeric_limits<int>::min )();

    int num_points = points.getPointCount();
    for ( int i = 0; i < num_points; ++i ) {
        auto x = points[i][0];
        auto y = points[i][1];
        auto z = points[i][2];
        x_min = ( x < x_min ) ? x : x_min;
        y_min = ( y < y_min ) ? y : y_min;
        z_min = ( z < z_min ) ? z : z_min;
        x_max = ( x > x_max ) ? x : x_max;
        y_max = ( y > y_max ) ? y : y_max;
        z_max = ( z > z_max ) ? z : z_max;
    }

    std::unordered_set<PCCPoint3D> us_points;
    for ( int i = 0; i < points.getPointCount(); i++) us_points.insert( points[i] );

    //Divide PC into voxels
    octree_recurse_decomp( us_points, chunks, num_points, x_min, x_max, y_min, y_max, z_min, z_max,
                           userParams.maxPointsPerVoxelOctree );

    //OVERWRITE user defined numROIs
    userParams.numROIs_ = chunks.size();

    //OVERWRITE user defined bounding boxes for each ROI with that of chunks i.e. in new impl, chunk and ROI are treated the same.
    roiBoundingBoxMinX.clear();
    roiBoundingBoxMaxX.clear();
    roiBoundingBoxMinY.clear();
    roiBoundingBoxMaxY.clear();
    roiBoundingBoxMinZ.clear();
    roiBoundingBoxMaxZ.clear();

    for (int i = 0; i < chunks.size(); i++) { 
        roiBoundingBoxMinX.push_back( chunks[i][0].first); //min x
        roiBoundingBoxMaxX.push_back( chunks[i][0].second ); //max x

        roiBoundingBoxMinY.push_back( chunks[i][1].first ); //min y
        roiBoundingBoxMaxY.push_back( chunks[i][1].second ); //max y

        roiBoundingBoxMinZ.push_back( chunks[i][2].first ); //min z
        roiBoundingBoxMaxZ.push_back( chunks[i][2].second ); //max z
    }

}
