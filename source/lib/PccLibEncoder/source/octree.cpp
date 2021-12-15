/**
 * @file octree.cpp
 * @author Pranav Venkatram (e0552200@u.nus.edu)
 * @brief 
 * @date 2021-12-07
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "octree.h"
using namespace std;
//overwrite params_.numROIs_ since used for patch generation and packing

int num_points_in_box( std::vector<pair<PCCPoint3D, PCCColor3B>>& points_global,
                       std::vector<pair<PCCPoint3D, PCCColor3B>>& points_local,
                       int x_min,
                       int x_max,
                       int y_min,
                       int y_max,
                       int z_min,
                       int z_max ) {
  int count      = 0;
  
  for ( int i = 0; i < points_global.size() && !points_global.empty(); i++ ) {
    auto x = points_global[i].first[0];
    auto y = points_global[i].first[1];
    auto z = points_global[i].first[2];

    if ( x_min <= x && x <= x_max && y_min <= y && y <= y_max && z_min <= z && z <= z_max ) {
      count++;
      points_local.push_back( points_global[i] );
      
      std::swap( points_global[i], *( points_global.end() - 1 ) );
      points_global.pop_back();
      i--;

      
    }
  }

  return count;
}


void octree_recurse_decomp( std::vector<pair<PCCPoint3D, PCCColor3B>>& points,
                            std::vector<std::vector<Range>>& chunks,
                            std::vector<pair<PCCPoint3D, PCCColor3B>>& centroids,
                            int                              num_points,
                            int                              x_min,
                            int                              x_max,
                            int                              y_min,
                            int                              y_max,
                            int                              z_min,
                            int                              z_max,
                            int                              maxPointsPerVoxel,
                            bool                             threeDD) {

  if ( num_points > maxPointsPerVoxel ) {
    int bisect_x = ( x_min + x_max ) / 2;
    int bisect_y = ( y_min + y_max ) / 2;
    int bisect_z = ( z_min + z_max ) / 2;

    std::vector<pair<PCCPoint3D, PCCColor3B>> points_front_bottom_left;
    int front_bottom_left_count = 
        num_points_in_box( points, points_front_bottom_left, x_min, bisect_x, y_min, bisect_y, z_min, bisect_z );
    octree_recurse_decomp( points_front_bottom_left, chunks, centroids, front_bottom_left_count, x_min, bisect_x, y_min,
                           bisect_y, z_min, bisect_z, maxPointsPerVoxel, threeDD );

    std::vector<pair<PCCPoint3D, PCCColor3B>> points_front_top_left;
    int front_top_left_count = 
        num_points_in_box( points, points_front_top_left, x_min, bisect_x, bisect_y + 1, y_max, z_min, bisect_z );
    octree_recurse_decomp( points_front_top_left, chunks, centroids, front_top_left_count, x_min, bisect_x,
                           bisect_y + 1, y_max, z_min, bisect_z, maxPointsPerVoxel, threeDD );

    std::vector<pair<PCCPoint3D, PCCColor3B>> points_front_bottom_right;
    int front_bottom_right_count = 
        num_points_in_box( points, points_front_bottom_right, bisect_x + 1, x_max, y_min,
                                                      bisect_y, z_min, bisect_z );
    octree_recurse_decomp( points_front_bottom_right, chunks, centroids, front_bottom_right_count, bisect_x + 1, x_max,
                           y_min, bisect_y, z_min, bisect_z, maxPointsPerVoxel, threeDD );

   std::vector<pair<PCCPoint3D, PCCColor3B>> points_front_top_right;
    int front_top_right_count = 
        num_points_in_box( points, points_front_top_right, bisect_x + 1, x_max,
                                                   bisect_y + 1, y_max, z_min, bisect_z );
   octree_recurse_decomp( points_front_top_right, chunks, centroids, front_top_right_count, bisect_x + 1, x_max,
                          bisect_y + 1, y_max, z_min, bisect_z, maxPointsPerVoxel, threeDD );

    std::vector<pair<PCCPoint3D, PCCColor3B>> points_back_bottom_left;
    int back_bottom_left_count = 
        num_points_in_box( points, points_back_bottom_left, x_min, bisect_x, y_min, bisect_y, bisect_z + 1, z_max );
    octree_recurse_decomp( points_back_bottom_left, chunks, centroids, back_bottom_left_count, x_min, bisect_x, y_min,
                           bisect_y, bisect_z + 1, z_max, maxPointsPerVoxel, threeDD );

    std::vector<pair<PCCPoint3D, PCCColor3B>> points_back_top_left;
    int back_top_left_count = 
        num_points_in_box( points, points_back_top_left, x_min, bisect_x, bisect_y + 1, y_max,
                                                 bisect_z + 1, z_max );
    octree_recurse_decomp( points_back_top_left, chunks, centroids, back_top_left_count, x_min, bisect_x, bisect_y + 1,
                           y_max, bisect_z + 1, z_max, maxPointsPerVoxel, threeDD );

    std::vector<pair<PCCPoint3D, PCCColor3B>> points_back_bottom_right;
    int back_bottom_right_count = 
        num_points_in_box( points, points_back_bottom_right, bisect_x + 1, x_max, y_min,
                                                     bisect_y, bisect_z + 1, z_max );
    octree_recurse_decomp( points_back_bottom_right, chunks, centroids, back_bottom_right_count, bisect_x + 1, x_max,
                           y_min, bisect_y, bisect_z + 1, z_max, maxPointsPerVoxel, threeDD );

    std::vector<pair<PCCPoint3D, PCCColor3B>> points_back_top_right;
    int back_top_right_count = 
        num_points_in_box( points, points_back_top_right, bisect_x + 1, x_max, bisect_y + 1,
                                                  y_max, bisect_z, z_max );
    octree_recurse_decomp( points_back_top_right, chunks, centroids, back_top_right_count, bisect_x + 1, x_max,
                           bisect_y + 1, y_max, bisect_z, z_max, maxPointsPerVoxel, threeDD );

  } else {
    Range rangeX = Range( x_min, x_max );
    Range rangeY = Range( y_min, y_max );
    Range rangeZ = Range( z_min, z_max );

    chunks.push_back( std::vector<Range>( { rangeX, rangeY, rangeZ } ) );

    if ( threeDD ) {
      centroids.push_back( get_centroid( points ) );
    }
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
                    std::vector<int>&                 roiBoundingBoxMaxZ) {
    
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

    std::vector<pair<PCCPoint3D, PCCColor3B>> vec_points;
    for ( int i = 0; i < points.getPointCount(); i++ )
      vec_points.push_back( make_pair( points[i], points.getColor( i ) ) );

    //Divide PC into voxels
    std::vector<pair<PCCPoint3D, PCCColor3B>> centroids; //DUMMY VAR, ONLY USED WHEN --threeDD
    octree_recurse_decomp( vec_points, chunks, centroids, num_points, x_min, x_max, y_min, y_max, z_min, z_max,
                           userParams.maxPointsPerVoxelOctree, userParams.threeDD);

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

//////////////////3DD using voxel grid filter///////////////////////////////////////////
std::vector<pair<PCCPoint3D, PCCColor3B>>
    get_octree_decomp_centroids( const PCCPointSet3&points, PCCEncoderParameters& userParams ) {
  // Bounding Box for PC
  int x_min;
  int x_max;
  int y_min;
  int y_max;
  int z_min;
  int z_max;
  x_min = y_min = z_min = ( std::numeric_limits<int>::max )();
  x_max = y_max = z_max = ( std::numeric_limits<int>::min )();
  int num_points        = points.getPointCount();
  for ( int i = 0; i < num_points; ++i ) {
    auto x = points[i][0];
    auto y = points[i][1];
    auto z = points[i][2];
    x_min  = ( x < x_min ) ? x : x_min;
    y_min  = ( y < y_min ) ? y : y_min;
    z_min  = ( z < z_min ) ? z : z_min;
    x_max  = ( x > x_max ) ? x : x_max;
    y_max  = ( y > y_max ) ? y : y_max;
    z_max  = ( z > z_max ) ? z : z_max;
  }
  std::vector<pair<PCCPoint3D, PCCColor3B>> vec_points;
  for ( int i = 0; i < points.getPointCount(); i++ ) vec_points.push_back( make_pair(points[i], points.getColor(i)) );
  std::vector<std::vector<Range>> chunks;

  // Divide PC into voxels, compute centroids
  std::vector<pair<PCCPoint3D, PCCColor3B>> centroids;
  octree_recurse_decomp( vec_points, chunks, centroids, num_points, x_min, x_max, y_min, y_max, z_min, z_max,
                         userParams.threeDDPointsPerVoxel, userParams.threeDD );
  return centroids;
}


std::pair<PCCPoint3D, PCCColor3B> get_centroid( std::vector<pair<PCCPoint3D, PCCColor3B>>& points) {
  int x_sum = 0;
  int    y_sum = 0;
  int    z_sum = 0;
  int    r_sum = 0;
  int    g_sum = 0;
  int    b_sum = 0;
  for ( int i = 0; i < points.size(); i++ ) {
    x_sum += points[i].first[0];
    y_sum += points[i].first[1];
    z_sum += points[i].first[2];
    r_sum += points[i].second[0];
    g_sum += points[i].second[1];
    b_sum += points[i].second[2];
  }
  x_sum = points.size() > 0 ? x_sum / points.size() : 0;
  y_sum = points.size() > 0 ? y_sum / points.size() : 0;
  z_sum = points.size() > 0 ? z_sum / points.size() : 0;
  r_sum = points.size() > 0 ? r_sum / points.size() : 0;
  g_sum = points.size() > 0 ? g_sum / points.size() : 0;
  b_sum = points.size() > 0 ? b_sum / points.size() : 0;

  //cout << "calc done, need to consturct pair" << endl;
  PCCPoint3D avg_geometry( (int)x_sum, (int)y_sum, (int)z_sum );
  //cout << "geom done, need to consturct pair" << endl;
  PCCColor3B avg_colour( (int)r_sum, (int)g_sum, (int)b_sum );
  //cout << "colour done, need to consturct pair" << endl;
  return std::make_pair( avg_geometry, avg_colour );
}

// averages geometry and colour of points in the same voxel
void threeDD_voxel_grid_filter( PCCGroupOfFrames& sources, int frameCount, PCCEncoderParameters& userParams ) {
  for ( int i = 0; i < frameCount; i++ ) {
    //std::vector<PCCBox3D>           boundingBoxes;
    std::vector<pair<PCCPoint3D, PCCColor3B>> centroids = get_octree_decomp_centroids( sources[i], userParams );
    //cout << "NUMBER OF CENTROIDS:      " << centroids.size() << endl;
    /*for ( auto& chunk : chunks ) {
      PCCBox3D box;
      box.min_[0] = chunk[0].first;
      box.max_[0] = chunk[0].second;
      box.min_[1] = chunk[1].first;
      box.max_[1] = chunk[1].second;
      box.min_[2] = chunk[2].first;
      box.max_[2] = chunk[2].second;
      boundingBoxes.push_back( box );
    }
    cout << "GEN BOUNDING BOXS DONE" << endl;
    std::vector<std::vector<int>> indexes;
    indexes.resize( boundingBoxes.size() );
    for ( int j = 0; j < sources[i].getPointCount(); j++ ) {
      for ( int k = 0; k < boundingBoxes.size(); k++ ) {
        if ( boundingBoxes[k].fullyContains( sources[i][j] ) ) { indexes[k].push_back( j ); }
      }
    }
    cout << "GET INDEXES FOR EACH BOUNDING BOX DONE" << endl;*/
    PCCPointSet3 downsampled_points;
    for ( auto &point : centroids ) {
      downsampled_points.addPoint( point.first, point.second );  // geometry and attribute of centroid
    }
    sources[i] = downsampled_points;
    //cout << "REPLACE OLD FRAME DONE" << endl;
  }
}