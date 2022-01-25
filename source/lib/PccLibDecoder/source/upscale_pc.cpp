#include "upscale_pc.h"

using namespace std;

void upscale_pcs( PCCGroupOfFrames& reconstructs, int frameCount, int upscale_factor ) {
	
	for ( int i = 0; i < frameCount; i++ ) {
    vector<PCCColor3B> cols  = reconstructs[i].getColors();
    vector<PCCPoint3D> coors = reconstructs[i].getPositions();

    reconstructs[i].resize( 0 );

    for ( int j = 0; j < coors.size(); j++ ) {
      coors[j] *= upscale_factor;

      reconstructs[i].addPoint( coors[j], cols[j] );
    }
  }
}


void upscale_pc( PCCPointSet3& reconstruct, int upscale_factor ) {
  for ( int i = 0; i < reconstruct.getPointCount(); i++ ) { reconstruct[i] *= upscale_factor; }
}