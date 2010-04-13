/**
 *
 *    Render Volume GPU
 *
 *  File: volume.h
 *
 *  Authors:
 *    Andre Maximo
 *    Ricardo Marroquim
 *
 *  Last Update: Mar 03, 2006
 *
 */

#ifndef NO_NVIDIA
#include "glslKernel.h"
#endif

extern "C" {
#include <GL/glut.h>
#include <float.h>
}

#include <vector>
#include "transferFunction.h"



#define MINORTHOSIZE -1.0
#define MAXORTHOSIZE 1.0

/// Object file format type
enum fileType { BIN8, BIN16BE, BIN16LE,
		MULTIRAW8, MULTIRAW16BE, MULTIRAW16LE,
		SINGLERAW8, SINGLERAW16, PNM8 };


typedef unsigned int uint;

using namespace std;

typedef struct _pairTet {

  GLuint id;
  GLfloat cZ;

} pairTet;

//#define TEX_FORMAT GL_TEXTURE_2D
//#define TEX_TYPE GL_RGBA
//#define TEX_FORMAT GL_TEXTURE_RECTANGLE_NV
//#define TEX_TYPE GL_FLOAT_RGBA32_NV
#define TEX_FORMAT GL_TEXTURE_2D
#define TEX_TYPE GL_RGBA32F_ARB

/// Extern global variable declarations from main

extern bool debug_shaders, debug_cout, debug_setup, integrating, sorting;
extern GLint modelWinWidth, modelWinHeight;

///Order of vertices to extract 5 tetrahedra from one hexahedra
const uint tetVertsIds[5][4] = {
  {0, 2, 3, 6},
  {3, 5, 6, 7},
  {0, 4, 5, 6},
  {0, 1, 3, 5},
  {0, 3, 5, 6}};


class volume
{
 public:
  volume();
  ~volume() {

  }

  /// Function to read volume files
  bool readMedBIN(const char* filename, const uint& _rows,
		  const uint& _cols, const uint& _slices,
		  const fileType& ft);


  void setSelectionBox(GLfloat minX, GLfloat maxX, GLfloat minY, GLfloat maxY, int);
  void setSelectionBox(GLint, GLint, int);

  void CreateTextures();
  void CreateInputTextures(void);
  void CreateShaders(void);
  void createTransferFunctionTex(void);
  void createArrays(void);

  void createScalarBuffer(const uint&, const uint&, const uint&);

  void addScalar(const uint&, const uint&, const uint&, const GLfloat&);

  void normalize();

  void draw(void);
  void drawBBox(bool select = false);
		
  void drawGrid(int);

  uint getNumVerts(void) { return numVerts; }
  uint getNumTets(void) { return numTets; }
  uint getPsiGamaTableSize(void) { return preIntTexSize; }
  uint getExpTexSize(void) { return preIntTexSize; }

  uint getBBoxMinX(void) { return selectionBox[0][0]; }
  uint getBBoxMinY(void) { return selectionBox[0][1]; }
  uint getBBoxMaxX(void) { return selectionBox[1][0]; }
  uint getBBoxMaxY(void) { return selectionBox[1][1]; }

  uint getRows(void) { return rows; }
  uint getCols(void) { return cols; }
  uint getSlices(void) { return slices; }

  uint getCurNumTets(void) { return (selectionBox[1][0] - selectionBox[0][0] )
      * (selectionBox[1][1] - selectionBox[0][1] )
      * (selectionBox[1][2] - selectionBox[0][2] ) * 5; }

  void resetSelectionBox(void);

  void projectBasisHexahedra(void);

  ///updates the hexa mask for discarding opacity zero cells
  void updateHexaMask(void);

  ///Single hexahedra volume for testing.
  void createTest(void);

  transferFunction tf;
  
  //Modelview, ModelViewProjection and Projection matrix
  double mv[16], mvpj[16], pj[16];

 private:

  GLfloat *tetrahedralBuffer, *positionBuffer;
		
#ifndef NO_NVIDIA

  glslKernel *shaders_1st, *shaders_2nd; // currently used shaders
  glslKernel *shaders_1st_with_sort, *shaders_1st_no_sort; // shader with/no sort
  glslKernel *shaders_2nd_with_int, *shaders_2nd_no_int; // shader with/no integrate
  GLuint frameBuffer, tetOutputTex0, tetOutputTex1; // FBO
  GLuint vertexPosTex, tetrahedralTex, orderTableTex; // Frag 1
  GLuint tfTex, expTex, psiGammaTableTex; // Frag 2

#else

  void compute_1st_shader_on_cpu(void);
  void compute_2nd_shader_on_cpu(void);

#endif
		
  ///Buffer containing all scalars read from RAW file [x][y][slice].
  unsigned short int *scalarBuffer;

  ///Buffer with mask for not rendering hexas with zero opacity
  unsigned char *hexaMaskBuffer;

  ///Buffer with transfer function
  float *tfTexBuffer;

  ///Number of rows/cols/slices on RAW data file.
  uint rows, cols, slices;

  ///Scaling factors to normalize dimensions.
  GLfloat scaleX, scaleY, scaleZ;

  ///Parameters for normalizing scalars
  GLfloat scaleScalar;
  unsigned short int minScalar;

  ///Selection area
  uint selectionBox[2][3];
  
  ///Basis Hexahedra vertices
  double basisVerts[8][4];

  ///Basis Hexahedra projected vertices
  GLfloat basisProjVerts[8*4];

  ///Basis Hexahedra thick vertices
  GLfloat basisThickVerts[ 5 * 4 ];

  ///Basis Hexahedra number of triangles per tetrahedra
  uint basisCount[5];

  ///Basis Hexahedra index to classification table per tetrahedra
  uint basisClassificationTableRow[5];

  ///Basis Hexahedra intersection parameters
  double basisIntersectionParams[5][2];

  ///Basis Hexahedra thickness per tetrahedra
  double basisThickness[5];

  ///Number of tetrahedra and vertices
  uint numHexs, numTets, numVerts;
  uint vertTexSize, tetTexSize, expTexSize, preIntTexSize;
	
  ///VERTEX ARRAY data structures
  GLint vertex_array[90];
  GLfloat normal_array[90], color_array[90];

  ///Data structures for glMultiDrawElements
  GLuint **indices;
  GLint *count;
  GLvoid **ids;

  void sortHexahedron(int *);
  void sort(double mask[3], int max, int& beg, int& end, int& step);
	
  void computeThickScalars(const uint& pos, const uint& tetId, GLfloat* scalars, GLfloat& sf, GLfloat& sb);

  ///updates the scale factor for the xyz axis
  void updateScaleFactors(void);

  ///Integrate ray
  TColor integrateColor(double sf, double sb, double l);
  //Computes a color given a scalar value and a thickness
  TColor computeColor(double sf, double sb, double thickness);
  ///Returns the scalar position on the buffer relative to vertex 0 of a hexahedral.
  uint hexahedralVertexWalk(uint vertId);	
  ///Reset Basis Hexahedra to its intial state
  void resetBasisHexahedra(void);
  ///Project a tetrahedron to screen space and compute its thickness
  void projectTetrahedra(uint tetId);
  ///Draw one tetrahedron
  void drawTetrahedron(uint, int, int, int, uint);
};
