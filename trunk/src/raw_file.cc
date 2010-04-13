/**
 *
 *    Render Volume GPU
 *
 *  File: raw_file.cc
 *
 *  Authors:
 *    Andre Maximo
 *    Ricardo Marroquim
 *
 *  Last Update: May 04, 2006
 *
 */

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

typedef unsigned int uint;
typedef unsigned short ushort;

static bool readRaw(const char* fn, volume& obj, const uint rows, const uint cols, const uint ns);

bool readRaw(const char* fn, volume& obj,
	     const uint rows, const uint cols, const uint ns)
{
  uint dimZ = (ns)/SLICE_OFFSET;
  ushort data;
  string filename;
  ifstream slice_file;

  //create buffers
  uint dimX =  (rows - (CULL_OFFSET_X_END + CULL_OFFSET_X_BEGIN)) / SAMPLE_OFFSET;  
  uint dimY =  (cols - (CULL_OFFSET_Y_END + CULL_OFFSET_Y_BEGIN)) / SAMPLE_OFFSET;

  uint numTets = (dimX-1)*(dimY-1)*(dimZ-1)*6;
  uint numVerts = dimX*dimY*dimZ;

  obj.createScalarBuffer(rows, cols, ns);

  for (uint z = 0; z < dimZ; ++z) {

    stringstream ssfn;
    ssfn << fn << "." << (z*SLICE_OFFSET + 1);
    filename = ssfn.str();

    slice_file.open(filename.c_str(), ios::binary);

    if(slice_file.fail())
      return false;

    for (uint y = 0; y < cols; ++y) {
      for (uint x = 0; x < rows; ++x) {
	char buf[2];
	slice_file.read(buf, 2);

	data = (ushort)((unsigned char)buf[0] + 256*(unsigned char)buf[1]);
// 	if ((x%SAMPLE_OFFSET == 0) && (y%SAMPLE_OFFSET == 0))
// 	  if ((x >= CULL_OFFSET_X_BEGIN) && (x <= rows - CULL_OFFSET_X_END - SAMPLE_OFFSET) &&
// 	      (y >= CULL_OFFSET_Y_BEGIN) && (y <= cols - CULL_OFFSET_Y_END - SAMPLE_OFFSET))
	    obj.addScalar(x, y, z, data)
      }
    }

    slice_file.close();
  }

  obj.Normalize();
    
  return true;
}
