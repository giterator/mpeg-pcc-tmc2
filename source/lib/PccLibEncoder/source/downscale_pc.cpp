#include"downscale_pc.h"

using namespace std;

void downscale_pc( PCCGroupOfFrames& sources, int frameCount, int downscale_factor ) {
  for ( int i = 0; i < frameCount; i++ ) {
    vector<PCCColor3B> cols  = sources[i].getColors();
    vector<PCCPoint3D> coors = sources[i].getPositions();

    sources[i].resize( 0 );

    for ( int j = 0; j < coors.size(); j++ ) {
      coors[j] /= downscale_factor;

      sources[i].addPoint( coors[j], cols[j] );
    }
  }
}

