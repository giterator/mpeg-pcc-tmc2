#include "upscale_pc.h"
#include "PCCKdTree.h"
#include <climits>

using namespace std;

void upscale_pcs( PCCGroupOfFrames& reconstructs, int frameCount, int upscale_factor ) {
	
	for ( int i = 0; i < frameCount; i++ ) {
    /*vector<PCCColor3B> cols  = reconstructs[i].getColors();
    vector<PCCPoint3D> coors = reconstructs[i].getPositions();

    reconstructs[i].resize( 0 );

    for ( int j = 0; j < coors.size(); j++ ) {
      coors[j] *= upscale_factor;

      reconstructs[i].addPoint( coors[j], cols[j] );
    }*/
        upscale_pc( reconstructs[i], upscale_factor );
    }
}


void upscale_pc( PCCPointSet3& reconstruct, int upscale_factor ) {
  //rescale geometry
  for ( int i = 0; i < reconstruct.getPointCount(); i++ ) { reconstruct[i] *= upscale_factor; }

  // rescaling geometry complete, need to interpolate points
  /*
    for each point in reconstruct:
        add to map of processed points
        query neighbors  wihtin neighbour radius
        for each neighbour:
            if (neighbour not in processed points i.e. if map[x][y][z].size() == 0)
                add to map of processed points
                interpolate (2 * upscale_factor) - 1 equidistant points along the direction vector (neighbour - point)
                    colour of interpolated point = weighted average of point & neighbour's colours = (point's colour *
    dist bet. point & interpolated/dist bet. point & neighbour) + (neighbours colour * dist bet. neighbour &
    interpolated/dist bet. point & neighbour)
                    [(point's colour * dist bet. point & interpolated) + (neighbours colour * dist bet. neighbour &
    interpolated)] * (dist bet. point & neighbour)
                    => add these points to new PCCPointSet3 interpolated
    */

  if (upscale_factor > 1) {
      //perform 3D interpolation
    double neighbour_radius = pow(upscale_factor, 2.0) * 9.0; // 9.0 is maxAllowedDist2RawPointsDetection_ in encoding using ctc-common config

    auto kd_tree = PCCKdTree( reconstruct );
    std::map<int, std::map<int, std::map<int, std::vector<size_t>>>> map;

    PCCPointSet3 interpolated;

    for ( int i = 0; i < reconstruct.getPointCount(); i++ ) { 
        int x = reconstruct[i][0];
        int y = reconstruct[i][1];
        int z = reconstruct[i][2];
        map[x][y][z].push_back( i );

        PCCNNResult neighbours;
        kd_tree.searchRadius( reconstruct[i], INT_MAX, neighbour_radius, neighbours );

        for (int j = 0; j < neighbours.size(); j++) {
          int index = neighbours.indices( j );
          double dist    = neighbours.dist( j );

          int x_n     = reconstruct[index][0];
          int y_n     = reconstruct[index][1];
          int z_n     = reconstruct[index][2];

          if ( map[x_n][y_n][z_n].size() == 0 ) { 
              map[x_n][y_n][z_n].push_back( index );
              auto np = reconstruct[index] - reconstruct[i];  // PCCPoint3D(x_n - x, y_n - y, z_n - z);

              for ( int k = 1; k <= pow( upscale_factor, 2.0 ) -1; k++) {
                  //geometry & colour of interpolated points
                  double geom_x = (double)x + ( (double)np[0] * (double)k / pow( upscale_factor, 2.0 ) );
                  double geom_y = (double)y + ( (double)np[1] * (double)k / pow( upscale_factor, 2.0 ) );
                  double geom_z = (double)z + ( (double)np[2] * (double)k / pow( upscale_factor, 2.0 ) );

                  PCCPoint3D geometry = PCCPoint3D(geom_x, geom_y, geom_z );   // reconstruct[i] + (np * (double)k / ( 2.0 * upscale_factor ));

                  PCCColor3B col_point = reconstruct.getColor(i);
                  PCCColor3B col_neighbour = reconstruct.getColor( index );

                  double pi_dist = dist * (double)k / ( pow( upscale_factor, 2.0 ) );
                  double ni_dist = dist - pi_dist;

                  double col_r = ( col_point[0] * ni_dist + col_neighbour[0] * pi_dist ) * ( 1 / dist );
                  double col_g = ( col_point[1] * ni_dist + col_neighbour[1] * pi_dist ) * ( 1 / dist );
                  double col_b = ( col_point[2] * ni_dist + col_neighbour[2] * pi_dist ) * ( 1 / dist );
                  PCCColor3B colour = PCCColor3B( col_r, col_g, col_b ); //PCCColor3B colour = (col_point * pi_dist + col_neighbour * ni_dist) * ( 1 / dist );

                  interpolated.addPoint( geometry, colour );
              }
          }
        }
    }

    
    reconstruct.appendPointSet( interpolated );

    cout << "3D interpolated point count: " << interpolated.getPointCount() << endl;
    cout << "Total point count after 3D interpolation: " << reconstruct.getPointCount() << endl;

  }
}