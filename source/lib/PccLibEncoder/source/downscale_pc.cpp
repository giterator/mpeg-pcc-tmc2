#include"downscale_pc.h"

using namespace std;

void downscale_pc( PCCGroupOfFrames& sources, int frameCount, int downscale_factor ) {
  for ( int i = 0; i < frameCount; i++ ) {
    vector<PCCColor3B> cols  = sources[i].getColors();
    vector<PCCPoint3D> coors = sources[i].getPositions();

    sources[i].resize( 0 );

    for ( int j = 0; j < coors.size(); j++ ) {
      //coors[j] /= downscale_factor; //change to rounding logic - down & up -> see shift
      float exact_x = coors[j][0] / (float)downscale_factor;
      float exact_y = coors[j][1] / (float)downscale_factor;
      float exact_z = coors[j][2] / (float)downscale_factor;

      coors[j][0] = round(exact_x);
      coors[j][1] = round( exact_y );
      coors[j][2] = round( exact_z );

      sources[i].addPoint( coors[j], cols[j] );
    }
  }
}

