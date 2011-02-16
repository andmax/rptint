/**
 *
 *    Render Volume GPU (regular)
 *
 *  File: volume.cc
 *
 *  Authors:
 *    Andre Maximo
 *    Ricardo Marroquim
 *
 *  Last Update: Oct 01, 2006
 *
 */

#define GL_GLEXT_PROTOTYPES

#include "volume.h"

extern "C" {
#include <math.h>
#include <GL/glut.h>
#ifndef NO_NVIDIA
#include <GL/glext.h>
#endif
}

#include <iostream>
#include <iomanip>

#include "tables.h"
#include <cstdlib>

/// Shaders CPU version

#include "../psiGammaTable512.h"

#ifndef NO_NVIDIA

/// Static local function
static void errcheck(const char * place);

/// Error check
/// @arg place were this functions was called

void errcheck(const char * place)
{
  static GLenum errCode;
  const GLubyte *errString;
    
  if ((errCode = glGetError()) != GL_NO_ERROR)
    {
      errString = gluErrorString(errCode);
      fprintf (stderr, "%s : OpenGL Error: %s\n", place, errString);
      exit(1);
    }
}

#endif

volume::volume()
{
  preIntTexSize = PSI_GAMMA_SIZE_FRONT;
}


void volume::updateHexaMask(void)
{
  uint pos = 0;
  float alpha = 0; 
  GLfloat scalar = 0;

  //draw all tetrahedra
  for(uint z = 0; z < slices-1; ++z)
    for(uint y = 0; y < cols-1; ++y)
      for(uint x = 0; x < rows-1; ++x)
	{
	  //position of first vertex of hexahedral on scalar buffer
		pos = x + y*(rows-1) + z*(rows-1)*(cols-1);
	  alpha = 0;

	  for (int i = 0; i < 8; ++i)
	    {
	      scalar = scalarBuffer[ pos + hexahedralVertexWalk( i )];
	      scalar = (scalar - minScalar) * scaleScalar;
	      alpha += tfTexBuffer[(uint)(scalar*255)*4 + 3];
	    }

	  if (alpha > 0.0)
	    hexaMaskBuffer[pos] = 1;
	  else
	    hexaMaskBuffer[pos] = 0;
	}
}

/// Create Transfer Function Texture
/// Texture:  { Transfer function (tf[i].r, tf[i].g, tf[i].b, tf[i].s) }

void volume::createTransferFunctionTex(void)
{
#ifndef NO_NVIDIA

  for (uint i = 0; i < 256; ++i)
    {
      for (uint j = 0; j < 4; ++j)
	//tfTexBuffer[i*4 + j] = tf[i][j];
	tfTexBuffer[i*4 + j] = tf.getColor(i).getValue(j);
    }

  glActiveTexture(GL_TEXTURE5);
  glBindTexture(GL_TEXTURE_1D, tfTex);
  glTexSubImage1D(GL_TEXTURE_1D, 0, 0, 256, GL_RGBA, GL_FLOAT, tfTexBuffer);
  
  shaders_2nd_with_int->use();
  shaders_2nd_with_int->set_uniform("tfTex", 5);
  shaders_2nd_with_int->set_uniform("brightness", tf.getBrightness());
  shaders_2nd_with_int->use(0);
  shaders_2nd_no_int->use();
  shaders_2nd_no_int->set_uniform("tfTex", 5);
  shaders_2nd_no_int->set_uniform("brightness", tf.getBrightness());
  shaders_2nd_no_int->use(0);

  errcheck("tf rebuild");
#endif
}

/// Create Input Texture
/// Texture:  { TF }
/// Texture:  { Exponential }
/// Texture:  { PsigammaTable }

void volume::CreateInputTextures(void)
{
#ifndef NO_NVIDIA

  tfTexBuffer = new float[256*4];

  for (int i = 0; i < 256; ++i)
    {
      for (int j = 0; j < 4; ++j)
	tfTexBuffer[i*4 + j] = tf[i][j];
    } 

  //Generate Transfer Function texture
  glGenTextures(1, &tfTex);
  glActiveTexture(GL_TEXTURE5);
  glBindTexture(GL_TEXTURE_1D, tfTex);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA, 256, 0, GL_RGBA, GL_FLOAT, tfTexBuffer);

  //Generate Exponential texture

  // This defines the exponential 1d texture size

#define EXP_SIZE 512

  expTexSize = EXP_SIZE;

  GLfloat *expTexBuffer;
  expTexBuffer = new GLfloat[expTexSize];

  for (int u = 0; u < (int)expTexSize; ++u)
    expTexBuffer[u] = (GLfloat)exp( -u / (GLfloat)(expTexSize - 1) );

  cout << "*** Exponential Texture Size : " << setw(10) << expTexSize << " ***" << endl;

  glGenTextures(1, &expTex);
  glActiveTexture(GL_TEXTURE6);
  glBindTexture(GL_TEXTURE_1D, expTex);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  //   glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  //   glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  //   glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  //   glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA, EXP_SIZE, 0, GL_ALPHA, GL_FLOAT, expTexBuffer);

  delete expTexBuffer;

  //Generate PsiGammaTable texture

  // This header includes the psiGammaTable matrix variable for the texture,
  // the matrix is a local variable that is cleaned in the end of this function.

#include "../psiGammaTable512.h"
    
  preIntTexSize = PSI_GAMMA_SIZE_FRONT; // always 2D quad texture
    
  cout << "*** psiGammaTex Texture Size : " << setw(10) << preIntTexSize << " ***" << endl;

  glGenTextures(1, &psiGammaTableTex);
  glActiveTexture(GL_TEXTURE7);
  glBindTexture(TEX_FORMAT, psiGammaTableTex);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, preIntTexSize, preIntTexSize,
   	       0, GL_ALPHA, GL_FLOAT, psiGammaTable);

  glTexParameteri(TEX_FORMAT, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(TEX_FORMAT, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glEnable(GL_TEXTURE_2D);

  errcheck("texCreation");
#else

  preIntTexSize = PSI_GAMMA_SIZE_FRONT; // always 2D quad texture

  if (debug_cout)
    cout << "*** psiGamma Table Size      : " << setw(10) << preIntTexSize << " ***" << endl;

#endif
}

///  Creates the glsl shaders.
///  * Second
///     * Vertex Shader: simply switch between w = 0 for the intersection
///                      vertex (already ModelViewProjected) and w = 1 for
///                      all others (doing the ModelViewProjection).
///     * Fragment Shader: look up the transfer function color by the scalar
///                        and compute the exponential alpha

/// Create Arrays
/// vertex_array, color_array, normal_array
/// indices, ids, count -- data structure for glMultiDrawElements

void volume::createArrays(void)
{

  indices = new GLuint*[5];

  for (uint i = 0; i < 5; ++i)
    indices[i] = new GLuint[6];

//   if (count)
//     delete count;

  count = new GLint[5];

  ids = new GLvoid*[5];

  for (uint i = 0; i < 5; ++i)
    { 
      ids[i] = (GLvoid*)indices[i];
      count[i] = 6;
    }

#ifndef NO_NVIDIA
  glEnableClientState(GL_COLOR_ARRAY);
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glColorPointer(3, GL_FLOAT, 0, &color_array[0]);
  glVertexPointer(3, GL_INT, 0, &vertex_array[0]);
  glNormalPointer(GL_FLOAT, 0, &normal_array[0]);
#endif

}

/// Create Textures
/// Input and Output Textures

void volume::CreateTextures()
{
  if (debug_cout) {
    cout << "*********************************************" << endl;
    cout << "*** Number of Vertices       : " << setw(10) << numVerts << " ***" << endl;
    cout << "*** Number of Tetrahedra     : " << setw(10) <<  numTets << " ***" << endl;
    cout << "*** Number of Hexahedra      : " << setw(10) <<  numHexs << " ***" << endl;
    cout << "*** Volume rows              : " << setw(10) <<     rows << " ***" << endl;
    cout << "*** Volume cols              : " << setw(10) <<     cols << " ***" << endl;
    cout << "*** Volume slices            : " << setw(10) <<   slices << " ***" << endl;
  }   

  CreateInputTextures();

}

void volume::CreateShaders()
{
#ifndef NO_NVIDIA
  /************ Install Shaders ******************************/
  // 2nd Shaders

  shaders_2nd_with_int = new glslKernel();
  shaders_2nd_with_int->vertex_source("shader/vert.shader");
  shaders_2nd_with_int->fragment_source("shader/frag_with_int.shader");
  shaders_2nd_with_int->install(debug_shaders);
  shaders_2nd_no_int = new glslKernel();
  shaders_2nd_no_int->vertex_source("shader/vert.shader"); 
  shaders_2nd_no_int->fragment_source("shader/frag_no_int.shader");
  shaders_2nd_no_int->install(debug_shaders);
  /*
  float xOffset = ((rows-1)/2);
  float yOffset = ((cols-1)/2);
  float zOffset = ((slices-1)/2);
  */
  shaders_2nd = shaders_2nd_with_int;
  shaders_2nd->use();

  // shaders_2nd->set_uniform("offset",xOffset, yOffset, zOffset);
  shaders_2nd->set_uniform("scale", (float)scaleX, (float)scaleY, (float)scaleZ);
  shaders_2nd->set_uniform("basisThickVerts", basisThickVerts, (GLuint)4, (GLsizei)5);
  shaders_2nd->set_uniform("basisProjVerts", basisProjVerts, (GLuint)4, (GLsizei)8);

  shaders_2nd->set_uniform("brightness", tf.getBrightness());
  shaders_2nd->set_uniform("tfTex", 5);
  shaders_2nd->set_uniform("expTex", 6);
  shaders_2nd->set_uniform("psiGammaTableTex", 7);
  shaders_2nd->set_uniform("preIntTexSize", (GLfloat)preIntTexSize);
  shaders_2nd->use(0);

  shaders_2nd = shaders_2nd_no_int;
  shaders_2nd->use();
  //  shaders_2nd->set_uniform("offset",xOffset, yOffset, zOffset);
  shaders_2nd->set_uniform("scale", (float)scaleX, (float)scaleY, (float)scaleZ);
  shaders_2nd->set_uniform("basisThickVerts", basisThickVerts, (GLuint)4, (GLsizei)5);
  shaders_2nd->set_uniform("basisProjVerts", basisProjVerts, (GLuint)4, (GLsizei)8);
  shaders_2nd->set_uniform("brightness", tf.getBrightness());
  shaders_2nd->set_uniform("tfTex", 5);
  shaders_2nd->set_uniform("expTex", 6);
  shaders_2nd->use(0);

  if (integrating) shaders_2nd = shaders_2nd_with_int;
  else shaders_2nd = shaders_2nd_no_int;

  errcheck("createShaders");

#endif
}

/**
 * Find the corresponding scalar value of a cube vertex on the scalar buffer.
 * Vertices id order:
 * Hexahedral front and back faces
 *
 *  2----3    6----7
 *  |    |    |    |
 *  |    |    |    |
 *  0----1    4----5
 *
 * The return value is the displacement on the buffer relative from vertex 0,
 * that must be applied to reach another vertex.
 * For example, to reach vertex 2, must add 'row + 1', where row is the number
 * of vertices on the x axis.
 * @param vertId Vertex to be retrieved.
 * @return Displacement.
 **/
uint volume::hexahedralVertexWalk(uint vertId)
{
  switch (vertId)
    {
    case 0:
      return 0;
    case 1:
      return 1;
    case 2:
      return rows;
    case 3:
      return rows + 1;
    case 4:
      return rows * cols;
    case 5:
      return (rows * cols) + 1;
    case 6:
      return (rows * cols) + rows;
    case 7:
      return (rows * cols) + rows + 1;
    default:
      return 0;
    }
}

/// Return the thick scalars of each tetrahedron

void volume::computeThickScalars(const uint& pos, const uint& tetId, GLfloat* scalars, GLfloat& sf, GLfloat& sb)
{
  //row of projection class on classification table
  int tableRow = basisClassificationTableRow[ tetId ];

  //compute front and back scalars depending on projection class
  GLfloat sb_tmp;
  switch (basisCount[tetId])
    {
    case 6:
      sf = scalars[ order_table[tableRow][0] ] + 
	basisIntersectionParams[tetId][0]*
	(scalars[ order_table[tableRow][2] ] - scalars[ order_table[tableRow][0] ]);

      sb = scalars[ order_table[tableRow][1] ] + 
	basisIntersectionParams[tetId][1]*
	(scalars[ order_table[tableRow][3] ] - scalars[ order_table[tableRow][1] ]);
      break;
    case 5:
      sf = scalars[ order_table[tableRow][2] ];
      sb_tmp = scalars[ order_table[tableRow][1] ] + 
	basisIntersectionParams[tetId][1]*
	(scalars[ order_table[tableRow][3] ] - scalars[ order_table[tableRow][1] ]);
      sb = scalars[ order_table[tableRow][0] ] + 
	basisIntersectionParams[tetId][0]*
	(sb_tmp - scalars[ order_table[tableRow][0] ]);
      break;
    case 4:
      sf = scalars[ order_table[tableRow][2] ];
      sb = scalars[ order_table[tableRow][1] ] + 
	basisIntersectionParams[tetId][1]*
	(scalars[ order_table[tableRow][3] ] - scalars[ order_table[tableRow][1] ]);
      break;
    case 3:
    default:
      sf = scalars[ order_table[tableRow][0] ];
      sb = scalars[ order_table[tableRow][1] ];
      break;
    }
}

/**
 * Computes a color given a scalar value.
 * @param Given scalar value.
 * @param thickness Thickness value to compute alpha.
 * @return Color.
 **/
TColor volume::computeColor(double sf, double sb, double thickness)
{
  
  double s = (sf + sb) * 0.5;
  TColor color = tf[ (int) (s*255) ];  
  //  color[3] = thickness * color[3];
  color[3] = 1.0 - exp( -thickness * color[3] ); 
  for (int j = 0; j < 3; ++j)
    color[j] *= color[3];

  return color;
}


/**
 * Integrate ray.
 * @param sf Front scalar.
 * @param sb Back scalar.
 * @param l Thickness.
 * @return Pixel color.
 **/
TColor volume::integrateColor(double sf, double sb, double l)
{
  TColor color;

  TColor colorFront = tf[ (int) (sf * 255) ];
  TColor colorBack = tf[ (int) (sb * 255) ];

  double tau[2];

  tau[0] = colorFront.a * l;
  tau[1] = colorBack.a * l;

  double zeta = exp( -(tau[0] * 0.5 + tau[1] * 0.5) );

  double gamma[2];

  gamma[0] = tau[0] / (1.0 + tau[0]);
  gamma[1] = tau[1] / (1.0 + tau[1]);
	
  double psi = psiGammaTable[(uint)(gamma[0] + 0.5/(double)preIntTexSize)*preIntTexSize]
    [(uint)(gamma[1] + 0.5/(double)preIntTexSize)*preIntTexSize];

  color[0] = colorFront[0]*(1.0 - psi) + colorBack[0]*(psi - zeta);
  color[1] = colorFront[1]*(1.0 - psi) + colorBack[1]*(psi - zeta);
  color[2] = colorFront[2]*(1.0 - psi) + colorBack[2]*(psi - zeta);
  color[3] = 1.0 - zeta;

  return color;
}

/// Sort regular volume
void volume::sort(double mask[3], int max, int& beg, int& end, int& step)
{
  double v1[4]  = { 1.0*mask[0], 1.0*mask[1], 1.0*mask[2], 1.0};
  double v2[4]  = {-1.0*mask[0], -1.0*mask[1], -1.0*mask[2], 1.0};

  double res1[4] = { 0.0, 0.0, 0.0, 0.0};
  double res2[4] = { 0.0, 0.0, 0.0, 0.0};

  for (uint r = 0; r < 4; ++r)
    for (uint c = 0; c < 4; ++c) {
      res1[r] += v1[c] * mv[c*4 + r];
      res2[r] += v2[c] * mv[c*4 + r];
    }

  if (res1[2] < res2[2]) { step = +1; beg = 0; end = max; }
  else { step = -1; beg = max-1; end = -1; }
}

/**
 * Sort the five tetrahedra inside the one hexahedron using
 * the modelview matrix.
 * @param Reference to an array containing the correct order.
**/
void volume::sortHexahedron(int *tetOrder)
{
  double hexaVerts[4][4] = {
	  {-0.5, 0.5, -0.5, 1.0},  //vertex 2 of basis hexahedron
	  {0.5, 0.5, 0.5, 1.0},  //vertex 7 of basis hexahedron
	  {-0.5, -0.5, 0.5, 1.0},  //vertex 4 of basis hexahedron
	  {0.5, -0.5, -0.5, 1.0}}; //vertex 1 of basis hexahedron

  double rotatedHexaVerts[4][4] = {
	  {0, 0, 0, 1},
	  {0, 0, 0, 1},
	  {0, 0, 0, 1},
	  {0, 0, 0, 1}};

  //rotate basis hexahedron vertices
  for (uint i = 0; i < 4; ++i)
    for (uint r = 0; r < 4; ++r)
      for (uint c = 0; c < 4; ++c)
	rotatedHexaVerts[i][r] += hexaVerts[i][c] * mv[c*4 + r];
  
  //sort rotated vertices
  for (int i = 0; i < 4; ++i)
    tetOrder[i] = i;

  double minZ;
  for (int i = 0; i < 3; ++i)
    {
      minZ = rotatedHexaVerts[i][2];
      for (int j = i+1; j < 4; ++j)
	{
	  if (rotatedHexaVerts[j][2] > minZ)
	    {
	      minZ = rotatedHexaVerts[j][2];
	      int tmp = tetOrder[i];
	      tetOrder[i] = tetOrder[j];
	      tetOrder[j] = tmp;
	    }
	}
    }
    
  tetOrder[4] = tetOrder[3];
  tetOrder[3] = tetOrder[2];
  tetOrder[2] = 4;
}


/**
 * Draw a single tetrahedron.
 * Tranlates the basis projected tetrahedra and draw it.
 * Compute the vertices coordinates using the indices of the
 * scalar matrix.
 * @param Indices of vertices in scalarBuffer matrix.
 * @param Index of the correponding tetrahedron of the basis hexahedra.
 **/
void volume::drawTetrahedron(uint tetId, int x, int y, int z, uint pos)
{
  double thickness = basisThickness[tetId];
 
  if (thickness == 0.0)
    return;

  //multiply thickness by brightness
  thickness *= tf.getBrightness();

  //row of projection class on classification table
  int tableRow = basisClassificationTableRow[tetId];

#ifdef NO_NVIDIA

  double xOffset = ((rows-1)/2);
  double yOffset = ((cols-1)/2);
  double zOffset = ((slices-1)/2);

  //compute rotated translation
  //centralize current coordinates on origin
  double curr[4] = {(double)x - xOffset,
		    (double)y - yOffset,
		    (double)z - zOffset, 1.0};
  
  //translation vector
  double transl[4] = {0.0, 0.0, 0.0, 0.0};
  
  /// Apply ModelViewProjection Matrix mvpj on translation vector
  for (uint r = 0; r < 4; ++r)
    for (uint c = 0; c < 4; ++c)
      transl[r] += curr[c] * mvpj[c*4 + r];
  
  transl[0] *= scaleX;
  transl[1] *= scaleY;
  transl[2] *= scaleZ;

  double verts[4][3];

  //compute the vertices coordinates according to its position in the basis projected hexahedral
  for (int i = 0; i < 4; ++i)
    {
      uint vId = tetVertsIds[tetId][i]; //id of vertex in the basis hexahedral
      verts[i][0] = basisProjVerts[vId*4+0] + transl[0];
      verts[i][1] = basisProjVerts[vId*4+1] + transl[1];
      verts[i][2] = basisProjVerts[vId*4+2] + transl[2];
    }


#endif

  //retrieve scalar values for the four tetrahedron vertices
  GLfloat scalars[4];
  for (uint i = 0; i < 4; ++i)
    scalars[i] = scalarBuffer[pos + hexahedralVertexWalk( tetVertsIds[tetId][i] )];

  //compute front and back scalars depending on projection class
  double sf, sb, sb_tmp;
  switch (basisCount[tetId])
    {
    case 6:
      sf = scalars[ order_table[tableRow][0] ] + 
	basisIntersectionParams[tetId][0]*
	(scalars[ order_table[tableRow][2] ] - scalars[ order_table[tableRow][0] ]);

      sb = scalars[ order_table[tableRow][1] ] + 
	basisIntersectionParams[tetId][1]*
	(scalars[ order_table[tableRow][3] ] - scalars[ order_table[tableRow][1] ]);
      break;
    case 5:
      sf = scalars[ order_table[tableRow][2] ];
      sb_tmp = scalars[ order_table[tableRow][1] ] + 
	basisIntersectionParams[tetId][1]*
	(scalars[ order_table[tableRow][3] ] - scalars[ order_table[tableRow][1] ]);
      sb = scalars[ order_table[tableRow][0] ] + 
	basisIntersectionParams[tetId][0]*
	(sb_tmp - scalars[ order_table[tableRow][0] ]);
      break;
    case 4:
      sf = scalars[ order_table[tableRow][2] ];
      sb = scalars[ order_table[tableRow][1] ] + 
	basisIntersectionParams[tetId][1]*
	(scalars[ order_table[tableRow][3] ] - scalars[ order_table[tableRow][1] ]);
      break;
    case 3:
    default:
      sf = scalars[ order_table[tableRow][0] ];
      sb = scalars[ order_table[tableRow][1] ];
      break;
    }
 
  //   //compute thick color
  TColor color;

#ifndef NO_NVIDIA

  color[0] = sf; color[1] = sb; color[2] = thickness; color[3] = 0.0;

  if (color[2] == 0) //totally transparent tetrahedra(no need to render)
    return;
#else

  if (integrating)
    color = integrateColor(sf, sb, thickness);
  else
    color = computeColor(sf, sb, thickness);

  if (color[3] == 0) //totally transparent tetrahedra(no need to render)
    return;

#endif

  //draw the triangle fan of the projected tetrahedron
  glBegin(GL_TRIANGLE_FAN);

  //draw thick vertex first (center of the triangle fan)
  glColor4f(color[0], color[1], color[2], color[3]); //use thick color for thick vertex

  int orderId = triangle_fan_order_table[ tableRow ][0]; //first vertex in fan_order_table

#ifndef NO_NVIDIA

  uint vId;

  if (orderId == -1) //count = 6, use computed intersection point as thick vertex
    vId = 8;
  else
    vId = tetVertsIds[tetId][orderId]; //id of vertex in the basis hexahedral

  glNormal3f(tetId, vId, 1);
  glVertex3i(x, y, z);

#else

  if (orderId == -1) //count = 6, use computed intersection point as thick vertex
    glVertex4d(basisThickVerts[tetId*4+0] + transl[0],
	       basisThickVerts[tetId*4+1] + transl[1],
	       0.0, 1.0);
  else

    glVertex4d(verts[ orderId ][0],
	       verts[ orderId ][1],
	       verts[ orderId ][2], 1.0);
#endif

  thickness = 0;
  //draw remaining of the fan
  for (uint i = 1; i < basisCount[tetId]; ++i)
    {
      orderId = triangle_fan_order_table[ tableRow ][i];

      double s = scalars[ orderId ]; //use actual vertex color for remaining vertices

#ifndef NO_NVIDIA

      color[0] = s; color[1] = s; color[2] = thickness; color[3] = 0.0;

      vId = tetVertsIds[tetId][orderId]; //id of vertex in the basis hexahedral

      glNormal3f(tetId, vId, 1);
      glColor4f(color[0], color[1], color[2], color[3]);
      glVertex3i(x, y, z);

#else
      if (integrating)
	color = integrateColor(s, s, thickness);
      else
	color = computeColor(s, s, thickness);

      glVertex4d(verts[ orderId ][0], verts[ orderId ][1],
		 verts[ orderId ][2], 1.0);

#endif

    }
  glEnd();

}

/**
 * Draw the volume by dividing each hexahedron in 6 tetrahedron.
 **/
void volume::draw()
{
  glDisable(GL_CULL_FACE);
  glDisable(GL_NORMALIZE);

  glEnable(GL_BLEND);

  glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);

  glClear(GL_COLOR_BUFFER_BIT);

#ifndef NO_NVIDIA

  if (integrating) shaders_2nd = shaders_2nd_with_int;
  else shaders_2nd = shaders_2nd_no_int;

  shaders_2nd->use();

  shaders_2nd->set_uniform("basisThickVerts", basisThickVerts, (GLuint)4, (GLsizei)5);
  shaders_2nd->set_uniform("basisProjVerts", basisProjVerts, (GLuint)4, (GLsizei)8);

#endif

  int xbeg, xend, stepx, ybeg, yend, stepy, zbeg, zend, stepz;
  double maskx[3] = {-1.0, 0.0, 0.0};
  double masky[3] = {0.0, -1.0, 0.0};
  double maskz[3] = {0.0, 0.0, -1.0};

  sort(maskx, rows-1, xbeg, xend, stepx);
  sort(masky, cols-1, ybeg, yend, stepy);
  sort(maskz, slices-1, zbeg, zend, stepz);

  int tetOrder[5];

  if (sorting)
    sortHexahedron(&tetOrder[0]);

  if (stepx > 0)
    {
      xbeg = selectionBox[0][0];
      xend = selectionBox[1][0];
    }
  else
    {
      xbeg = selectionBox[1][0]-1;
      xend = selectionBox[0][0]-1;
    }
  if (stepy > 0)
    {
      ybeg = selectionBox[0][1];
      yend = selectionBox[1][1];
    }
  else
    {
      ybeg = selectionBox[1][1]-1;
      yend = selectionBox[0][1]-1;
    }
  if (stepz > 0)
    {
      zbeg = selectionBox[0][2];
      zend = selectionBox[1][2];
    }
  else
    {
      zbeg = selectionBox[1][2]-1;
      zend = selectionBox[0][2]-1;
    }

  GLfloat xOffset = (selectionBox[0][0]+selectionBox[1][0])/2.0;
  GLfloat yOffset = (selectionBox[0][1]+selectionBox[1][1])/2.0;
  GLfloat zOffset = (selectionBox[0][2]+selectionBox[1][2])/2.0;

  int tableRow, orderId, vId, tetId;

  uint k=0;
   // cout << "--- draw --- " << endl;
  //compute normal and indices arrays
  for (uint i=0; i < 5; i++) 
    {    
      tetId = tetOrder[ i ];
      tableRow = basisClassificationTableRow[ tetId ];

      for (uint j=0; j < basisCount[ tetId ]; j++)
	{
	  orderId = triangle_fan_order_table[ tableRow ][ j ];

	  if (orderId == -1)
	    vId = 8;
	  else
	    vId = tetVertsIds[ tetId ][ orderId ];

	  normal_array[k + 0] = (GLfloat)tetId;
	  normal_array[k + 1] = (GLfloat)vId;
	  normal_array[k + 2] = 1.0; // z values aren't used
	  indices[i][j] = k/3;
	  k += 3;
	}
      count[i] = basisCount[ tetId ];
    }

  uint array_size = k;

  GLfloat s, thickness, tetScalars[4], scalars[8], sf, sb;
  uint pos;

  glTranslatef(-xOffset, -yOffset, -zOffset);

  //numTets = (rows-1)*(cols-1)*(slices-1)*5;
  numTets = 0;

//   double totalSetup = 0, totalDraw = 0, totalFetch = 0;
//   int sta_sh=0, end_sh=0;

  //draw all tetrahedra
  for(GLint z = zbeg; z != zend; z += stepz)
    for(GLint y = ybeg; y != yend; y += stepy)
      for(GLint x = xbeg; x != xend; x += stepx)
	{	  
	  //position of first vertex of hexahedral on scalar buffer
		pos = x + y*(rows-1) + z*(rows-1)*(cols-1);

	  if (!hexaMaskBuffer[pos])
	    continue;

	  for (int i = 0; i < 4; ++i)
 	    scalars[i] = (scalarBuffer[ pos + hexahedralVertexWalk( i )] - minScalar)*scaleScalar;

	  scalars[4] = scalars[0];
	  scalars[5] = scalars[1];
	  scalars[6] = scalars[2];
	  scalars[7] = scalars[3];

	  numTets += 5;
	  //compute vertex array for this hexahedron
	  for(uint i=0; i < array_size; i += 3)
	    {
	      vertex_array[i + 0] = x;
	      vertex_array[i + 1] = y;
	      vertex_array[i + 2] = z;
	    }

	  k = 0;
	  //compute the color array
	  for (uint i = 0; i < 5; i++)
	    {
	      tetId = tetOrder[ i ];
	      tableRow = basisClassificationTableRow[ tetId ];
	      thickness = basisThickness[ tetId ] * tf.getBrightness();

	      //retrieve scalar values for the four tetrahedron vertices
	      for (uint l = 0; l < 4; l++)
		tetScalars[l] = scalars[ tetVertsIds[ tetId ][ l ] ];

	      computeThickScalars(pos, tetId, tetScalars, sf, sb);
	      color_array[k + 0] = sf;
	      color_array[k + 1] = sb;
	      color_array[k + 2] = thickness;

	      k += 3;

	      for (uint j = 1; j < basisCount[ tetId ]; j++)
		{
		  orderId = triangle_fan_order_table[ tableRow ][ j ];
		  s = tetScalars[ orderId ];

		  color_array[k + 0] = s;
		  color_array[k + 1] = s;
		  color_array[k + 2] = 0.0;
		  k += 3;
		}
	    }
	  
//	  if (!debug_setup)
	    //	    Draw 5 triangle fans
	  glMultiDrawElements(GL_TRIANGLE_FAN, count,
			      GL_UNSIGNED_INT, (const GLvoid**) ids, 5);

	  /*
	  for (uint i = 0; i < 5; ++i)
	    if (sorting)
	      drawTetrahedron(tetOrder[i], x, y, z, pos);
	    else
	      drawTetrahedron(i, x, y, z, pos);
	  */
	}
   // cout << "fetch time : " << totalFetch / 1000 << endl;
   // cout << "setup time : " << totalSetup / 1000 << endl;
   // cout << "draw time : " << totalDraw/ 1000 << endl;

#ifndef NO_NVIDIA

  shaders_2nd->use(0);

#endif

  glDisable(GL_BLEND);

  glDisable(GL_CULL_FACE);

}


/**
 * Computes the projection class and the thick vertex of a tetrahedron.
 * @param vertIds Ids of vertices in basis cubes completing one tetrahedron.
 * @param tetId Identification of the tetrahedron in basis hexahedra.
 **/
void volume::projectTetrahedra(uint tetId)
{
  double vertOrder[4][4];
  int tests[4];

  //define vectors with tetrahedron's projected vertices
  double vc1_0[3], vc2_0[3], vc3_0[3], vc1_2[3], vc1_3[3];
  for (int i = 0; i < 3; ++i)
    {
      vc1_0[i] = basisProjVerts[ tetVertsIds[tetId][1] *4 + i] - basisProjVerts[ tetVertsIds[tetId][0] *4 + i];
      vc2_0[i] = basisProjVerts[ tetVertsIds[tetId][2] *4 + i] - basisProjVerts[ tetVertsIds[tetId][0] *4 + i];
      vc3_0[i] = basisProjVerts[ tetVertsIds[tetId][3] *4 + i] - basisProjVerts[ tetVertsIds[tetId][0] *4 + i];
      vc1_2[i] = basisProjVerts[ tetVertsIds[tetId][1] *4 + i] - basisProjVerts[ tetVertsIds[tetId][2] *4 + i];
      vc1_3[i] = basisProjVerts[ tetVertsIds[tetId][1] *4 + i] - basisProjVerts[ tetVertsIds[tetId][3] *4 + i];      
    }

  //compute cross products between vectors 
  double cross[4];
  cross[0] = (vc1_0[0] * vc2_0[1]) - (vc1_0[1] * vc2_0[0]);
  cross[1] = (vc1_0[0] * vc3_0[1]) - (vc1_0[1] * vc3_0[0]);
  cross[2] = (vc2_0[0] * vc3_0[1]) - (vc2_0[1] * vc3_0[0]);
  cross[3] = (vc1_2[0] * vc1_3[1]) - (vc1_2[1] * vc1_3[0]);

  //count the number of tests with result = 1 (degenerated cases)
  int num_crosses0 = 0;
  for (uint i = 0; i < 4; ++i)
    {  
      tests[i] = (cross[i] < 0.0) ? 0 : ( (cross[i] > 0.0) ? 2 : 1 );      
      if (tests[i] == 1)
	++num_crosses0;    
    }

  if (num_crosses0 == 2)
    basisCount[tetId] = 3;
  else if (num_crosses0 == 1)
    basisCount[tetId] = 4;
  else
    basisCount[tetId] = 5;
  
  //determine row of classification table
  int id_order = tests[0] * 27 + tests[1] * 9 + tests[2] * 3 + tests[3] * 1;  

  //--- Order the vertices (order_vertices) ---
  for(uint i = 0; i < 4; ++i)
    {   
      int ord = order_table[id_order][i];
      for(uint j = 0; j < 4; ++j)
	vertOrder[i][j] = basisProjVerts[ tetVertsIds[tetId][ord] *4 + j];      
    }

  //--- Compute parameters (compute_params) ---
  double denominator, numeratorU1, numeratorU2;
  double paramU1, paramU2;

  if (basisCount[tetId] == 3)
    {
      paramU1 = 1.0;
      paramU2 = 1.0;
    }
  else if (basisCount[tetId] == 4)
    {
      //line intersection denominator between v0->v2 and v1->v3
      denominator = ((vertOrder[3][1] - vertOrder[1][1]) * (vertOrder[2][0] - vertOrder[0][0])) -
	((vertOrder[3][0] - vertOrder[1][0]) * (vertOrder[2][1] - vertOrder[0][1]));
	  
      //line defined by vector v1->v3
      numeratorU2 = ((vertOrder[2][0] - vertOrder[0][0]) * (vertOrder[0][1] - vertOrder[1][1])) -
	((vertOrder[2][1] - vertOrder[0][1]) * (vertOrder[0][0] - vertOrder[1][0]));
	  
      paramU1 = 1.0;
      paramU2 = numeratorU2 / denominator;
    }
  else
    {
      //line intersection denominator between v0->v2 and v1->v3
      denominator = ((vertOrder[3][1] - vertOrder[1][1]) * (vertOrder[2][0] - vertOrder[0][0])) -
	((vertOrder[3][0] - vertOrder[1][0]) * (vertOrder[2][1] - vertOrder[0][1]));
	  
      //line defined by vector v0->v2
      numeratorU1 = ((vertOrder[3][0] - vertOrder[1][0]) * (vertOrder[0][1] - vertOrder[1][1])) -
	((vertOrder[3][1] - vertOrder[1][1]) * (vertOrder[0][0] - vertOrder[1][0]));
	  
      //line defined by vector v1->v3
      numeratorU2 = ((vertOrder[2][0] - vertOrder[0][0]) * (vertOrder[0][1] - vertOrder[1][1])) -
	((vertOrder[2][1] - vertOrder[0][1]) * (vertOrder[0][0] - vertOrder[1][0]));
	  
      paramU1 = numeratorU1 / denominator;
      paramU2 = numeratorU2 / denominator;
    }

  //--- Compute intersection point (compute_intersection) ---
  double intersectionPoint[3];
  for (int i = 0; i < 3; ++i)
    intersectionPoint[i] = 0.0;

  if (basisCount[tetId] == 3)
    {
      basisThickness[tetId] = fabs(vertOrder[0][2] - vertOrder[1][2]);
    }
  else if (basisCount[tetId] == 4)
    {
      double zBackIntersection = vertOrder[1][2] + paramU2*(vertOrder[3][2] - vertOrder[1][2]);
      basisThickness[tetId] = fabs(vertOrder[2][2] - zBackIntersection);
    }
  else
    {
      //find z coordinate of back intersection point by interpolating original vertices (not projected)
      double zBackIntersection = vertOrder[1][2] + paramU2*(vertOrder[3][2] - vertOrder[1][2]);
	  
      //find ordered intersection point between the two ordered lines (basis graph)
      for(uint i = 0; i < 3; ++i)
	intersectionPoint[i] = (vertOrder[0][i] + paramU1*(vertOrder[2][i] - vertOrder[0][i]));
	  
      basisThickness[tetId] = fabs(intersectionPoint[2] - zBackIntersection);
    }

  //if Class 1 then paramU1 is greater than 1.0 (middle vertex inside projected triangle)
  //in this case the thickness value must be multiplied by the
  //ratio r=(|V0V4|/|V0I|), that is, r=1/paramU1,
  //where V0V4 is the distance from V0 to the middle vertex
  //and V0I is the distance from V0 to the intersection point
  if (paramU1 >= 1.0)
    {
      if (basisCount[tetId] == 5)
	{
	  basisThickness[tetId] /= paramU1;
	  paramU1 = 1.0 / paramU1;
	}
    }
  else // count == 6 ; thick vertex == intersection vertex
    {
      basisCount[tetId] = 6;
    }

  basisClassificationTableRow[tetId] = id_order;

  basisIntersectionParams[tetId][0] = paramU1;
  basisIntersectionParams[tetId][1] = paramU2;

  basisThickVerts[tetId*4 + 0] = intersectionPoint[0];
  basisThickVerts[tetId*4 + 1] = intersectionPoint[1];
  basisThickVerts[tetId*4 + 2] = 0.0;
  basisThickVerts[tetId*4 + 3] = 1.0;
}

/**
 * Computes the projection class for the six tetrahedra of
 * the basis hexahedra.
 * This projection classes are used for all remaining tetrahedra.
 **/
void volume::projectBasisHexahedra(void)
{
  resetBasisHexahedra();

  //--- Compute the ModelViewProjection Matrix ---
  glGetDoublev(GL_MODELVIEW_MATRIX, mv);
  glGetDoublev(GL_PROJECTION_MATRIX, pj);

  glMatrixMode(GL_MODELVIEW);
  //glMatrixMode(GL_PROJECTION);

  glPushMatrix();
  glMultMatrixd(pj);
  //glMultMatrixd(mv);
  glGetDoublev(GL_MODELVIEW_MATRIX, mvpj);
  //glGetDoublev(GL_PROJECTION_MATRIX, mvpj);
  glPopMatrix();

  //project the vertices of the basis hexahedra
  for (uint i = 0; i < 8; ++i)
    {          
      /// Apply ModelViewProjection Matrix mvpj
      for (uint r = 0; r < 4; ++r)
	for (uint c = 0; c < 4; ++c)
	  basisProjVerts[i*4 + r] += basisVerts[i][c] * mvpj[c*4 + r];

      /// Apply ModelView Matrix mv (z -> r=2)
      basisProjVerts[i*4 + 2] = 0.0;
      for (uint c = 0; c < 4; ++c)
	basisProjVerts[i*4 + 2] += basisVerts[i][c] * mv[c*4+2];

      /// Scale the vertex to the size of one hexahedra
      basisProjVerts[i*4 + 0] *= scaleX;
      basisProjVerts[i*4 + 1] *= scaleY;
      basisProjVerts[i*4 + 2] *= scaleZ;
    }

  for (uint i = 0; i < 5; ++i)
    projectTetrahedra(i);

}

/**
 * Create a buffer (3D matrix) to store all scalar values.
 * @param rows Number of rows per slice.
 * @param cols Number of columns per slice.
 * @param slices Number of slices.
 **/
void volume::createScalarBuffer(const uint& _rows, const uint& _cols, const uint& _slices)
{
  rows = _rows;
  cols = _cols;
  slices = _slices;

  numHexs = (rows-1)*(cols-1)*(slices-1);
  numTets = (rows-1)*(cols-1)*(slices-1)*5;
  numVerts = rows*cols*slices;

  hexaMaskBuffer = new unsigned char[numHexs];
  scalarBuffer = new unsigned short int[numVerts];

  //compute a scaling factor for all axis
  scaleX = 1/((double)rows-1);
  scaleY = 1/((double)cols-1);
  scaleZ = 1/((double)slices-1);

  double smax = scaleZ;
  if (scaleX < smax)
    smax = scaleX;
  if (scaleY < smax)
    smax = scaleY;

  scaleX = smax;
  scaleY = smax;
  scaleZ = smax;
}

/**
 * Sets the value for the scalar buffer at given position.
 * @param x X coordinate of buffer.
 * @param y Y coordinate of buffer.
 * @param z Z coordinate of buffer.
 * @param data Scalar value at given position.
 **/
void volume::addScalar(const uint& x, const uint& y, const uint& z, const float& data)
{
  scalarBuffer[x + y*rows + z*rows*cols] = (unsigned short int)data;
}

/**
 * Normalize all scalars between [0, 1]
 **/
void volume::normalize(void)
{
  unsigned short int data, min, max;

  min = scalarBuffer[0];
  max = scalarBuffer[0];
    
  for(uint x = 0; x < rows; ++x)
    for(uint y = 0; y < cols; ++y)
      for(uint z = 0; z < slices; ++z)
	{
	  data = scalarBuffer[x + y*rows + z*rows*cols];

	  if(min > data) min = data;
	  if(max < data) max = data;       
	}

  scaleScalar = (double)( 1.0 / (max - min) );
  minScalar = min;

//   for(uint z = 0; z < slices; ++z)
//     for(uint y = 0; y < cols; ++y)
//       for(uint x = 0; x < rows; ++x)
// 	{
// 	  uint pos = x + y*rows + z*rows*cols;
// 	  data = scalarBuffer[pos];
// 	  scalarBuffer[pos] = (data - min) * scale;
// 	}
}

void volume::resetSelectionBox(void) 
{
  selectionBox[0][0] = 0;
  selectionBox[1][0] = cols-1;
  selectionBox[0][1] = 0;
  selectionBox[1][1] = rows-1;
  selectionBox[0][2] = 0;
  selectionBox[1][2] = slices-1;
  updateScaleFactors();
}

/**
 * Return the basis hexahedra vertices to their original positions.
 **/
void volume::resetBasisHexahedra(void)
{
  basisVerts[0][0] =  0.0;
  basisVerts[0][1] =  0.0;
  basisVerts[0][2] =  0.0;
  basisVerts[0][3] =  1.0;

  basisVerts[1][0] =  1.0;
  basisVerts[1][1] =  0.0;
  basisVerts[1][2] =  0.0;
  basisVerts[1][3] =  1.0;

  basisVerts[2][0] =  0.0;
  basisVerts[2][1] =  1.0;
  basisVerts[2][2] =  0.0;
  basisVerts[2][3] =  1.0;

  basisVerts[3][0] =  1.0;
  basisVerts[3][1] =  1.0;
  basisVerts[3][2] =  0.0;
  basisVerts[3][3] =  1.0;

  basisVerts[4][0] =  0.0;
  basisVerts[4][1] =  0.0;
  basisVerts[4][2] =  1.0;
  basisVerts[4][3] =  1.0;

  basisVerts[5][0] =  1.0;
  basisVerts[5][1] =  0.0;
  basisVerts[5][2] =  1.0;
  basisVerts[5][3] =  1.0;

  basisVerts[6][0] =  0.0;
  basisVerts[6][1] =  1.0;
  basisVerts[6][2] =  1.0;
  basisVerts[6][3] =  1.0;

  basisVerts[7][0] =  1.0;
  basisVerts[7][1] =  1.0;
  basisVerts[7][2] =  1.0;
  basisVerts[7][3] =  1.0;

  for (int i = 0; i < 8; ++i)
    for (int j = 0; j < 4; ++j)
      basisProjVerts[i*4 + j] = 0.0;

  for (int i = 0; i < 5; ++i)
    {
      basisCount[i] = 0;
      basisThickness[i] = 0.0;
      basisClassificationTableRow[i] = 0;
      basisThickVerts[i*4 + 0] = 0.0;
      basisThickVerts[i*4 + 1] = 0.0;
      basisThickVerts[i*4 + 2] = 0.0;
      basisThickVerts[i*4 + 3] = 0.0;
    }
}

/**
 * Creates a test volume with one hexahedra.
 **/
void volume::createTest(void)
{
  rows = 2;
  cols = 2;
  slices = 2;

  scalarBuffer = new unsigned short int[rows * cols * slices];

  for (uint i = 0; i < rows * cols * slices; ++i)
    scalarBuffer[i] = 1;

  scalarBuffer[0] = 0;
  scalarBuffer[1] = 0;
  scalarBuffer[2] = 0;
  scalarBuffer[3] = 0;

  normalize();

  numTets = (rows-1)*(cols-1)*(slices-1)*5;

  //compute a scaling factor for all axis
  scaleX = 1/((double)rows-1);
  scaleY = 1/((double)cols-1);
  scaleZ = 1/((double)slices-1);

  double smax = scaleZ;
  if (scaleX < smax)
    smax = scaleX;
  if (scaleY < smax)
    smax = scaleY;

  scaleX = smax;
  scaleY = smax;
  scaleZ = smax;

  selectionBox[0][0] = 0;
  selectionBox[1][0] = cols-1;
  selectionBox[0][1] = 0;
  selectionBox[1][1] = rows-1;
  selectionBox[0][2] = 0;
  selectionBox[1][2] = slices-1;
}

void volume::updateScaleFactors(void)
{
  uint sizeX = selectionBox[1][0] - selectionBox[0][0];
  uint sizeY = selectionBox[1][1] - selectionBox[0][1];
  uint sizeZ = selectionBox[1][2] - selectionBox[0][2];

  //compute a scaling factor for all axis
  scaleX = 1/(double)sizeX;
  scaleY = 1/(double)sizeY;
  scaleZ = 1/(double)sizeZ;

  double smax = scaleZ;

  if (scaleX < smax)
    smax = scaleX;
  if (scaleY < smax)
    smax = scaleY;

  scaleX = smax;
  scaleY = smax;
  scaleZ = smax;

  shaders_2nd = shaders_2nd_with_int;
  shaders_2nd->use();
  shaders_2nd->set_uniform("scale", (float)scaleX, (float)scaleY, (float)scaleZ);
  shaders_2nd->use(0);

  shaders_2nd = shaders_2nd_no_int;
  shaders_2nd->use();
  shaders_2nd->set_uniform("scale", (float)scaleX, (float)scaleY, (float)scaleZ);
  shaders_2nd->use(0);
}


void volume::setSelectionBox(GLint pt0, GLint pt1, int focusFace)
{
  uint minX, minY, maxX, maxY;

  uint sizeX = selectionBox[1][0] - selectionBox[0][0];
  //  uint sizeY = selectionBox[1][1] - selectionBox[0][1];
  uint sizeZ = selectionBox[1][2] - selectionBox[0][2];

  uint size = sizeX;
  //lateral faces
  if ((focusFace == 2) || (focusFace == 3))
    size = sizeZ;

  minX = pt0 % size;
  minY = pt0 / size;
  maxX = pt1 % size;
  maxY = pt1 / size;

  if (minX > maxX)
    {
      uint tmp = minX;
      minX = maxX;
      maxX = tmp;
    }
  if (minY > maxY)
    {
      uint tmp = minY;
      minY = maxY;
      maxY = tmp;
    }

  //previous dimensions of the box, before this new cut
  //to add to the new cut
  uint offset[2];

  //depending on the selected face, the cut must be made on different
  //axis and directions of the selection box
  switch (focusFace)
    {
    case 0: //back face
      offset[0] = selectionBox[1][0];
      offset[1] = selectionBox[0][1];
      selectionBox[0][0] = offset[0] - maxX;
      selectionBox[0][1] = offset[1] + minY;
      selectionBox[1][0] = offset[0] - minX;
      selectionBox[1][1] = offset[1] + maxY;
      break;
    case 1: //front face
      offset[0] = selectionBox[0][0];
      offset[1] = selectionBox[0][1];
      selectionBox[0][0] = offset[0] + minX;
      selectionBox[0][1] = offset[1] + minY;
      selectionBox[1][0] = offset[0] + maxX;
      selectionBox[1][1] = offset[1] + maxY;
      break;
    case 2: //left face
      offset[0] = selectionBox[0][1];
      offset[1] = selectionBox[0][2];
      selectionBox[0][1] = offset[0] + minY;
      selectionBox[0][2] = offset[1] + minX;
      selectionBox[1][1] = offset[0] + maxY;
      selectionBox[1][2] = offset[1] + maxX;
      break;
    case 3: //right face
      offset[0] = selectionBox[0][1];
      offset[1] = selectionBox[1][2];
      selectionBox[0][1] = offset[0] + minY;
      selectionBox[0][2] = offset[1] - maxX;
      selectionBox[1][1] = offset[0] + maxY;
      selectionBox[1][2] = offset[1] - minX;
      break;
    case 4: //top face
      offset[0] = selectionBox[0][0];
      offset[1] = selectionBox[1][2];
      selectionBox[0][0] = offset[0] + minX;
      selectionBox[0][2] = offset[1] - maxY;
      selectionBox[1][0] = offset[0] + maxX;
      selectionBox[1][2] = offset[1] - minY;
      break;
    case 5: //bottom face
      offset[0] = selectionBox[0][0];
      offset[1] = selectionBox[0][2];
      selectionBox[0][0] = offset[0] + minX;
      selectionBox[0][2] = offset[1] + minY;
      selectionBox[1][0] = offset[0] + maxX;
      selectionBox[1][2] = offset[1] + maxY;
      break;
    }

  updateScaleFactors();
}
/**
 * Sets the selection box num_tets = vol->getCurNumTets()Xwhich is the drawable part of the volume.
 **/
void volume::setSelectionBox(GLfloat minX, GLfloat minY, GLfloat maxX, GLfloat maxY, int focusFace)
{
  if (minX > maxX)
    {
      GLfloat tmp = minX;
      minX = maxX;
      maxX = tmp;
    }

  if (minY > maxY)
    {
      GLfloat tmp = minY;
      minY = maxY;
      maxY = tmp;
    }

  if (minX < -0.5) minX = -0.5;
  if (minY < -0.5) minY = -0.5;
  if (maxX >  0.5) maxX =  0.5;
  if (maxY >  0.5) maxY =  0.5;

  if ((minX >= maxX) || (minY >= maxY))
    return;

  uint sizeX = selectionBox[1][0] - selectionBox[0][0];
  uint sizeY = selectionBox[1][1] - selectionBox[0][1];
  uint sizeZ = selectionBox[1][2] - selectionBox[0][2];

  minX += 0.5;
  minY += 0.5;
  maxX += 0.5;
  maxY += 0.5;

  double tmp;
  switch (focusFace)
    {
    case 0:
      tmp = minX;
      minX = 1.0 - maxX;
      maxX = 1.0 - tmp;
      selectionBox[0][0] = (uint)(minX * (sizeX));
      selectionBox[0][1] = (uint)(minY * (sizeY));
      selectionBox[1][0] = (uint)(maxX * (sizeY));
      selectionBox[1][1] = (uint)(maxY * (sizeX));
      break;
    case 1:
      selectionBox[0][0] = (uint)(minX * (sizeY));
      selectionBox[0][1] = (uint)(minY * (sizeX));
      selectionBox[1][0] = (uint)(maxX * (sizeY));
      selectionBox[1][1] = (uint)(maxY * (sizeX));
      break;
    case 2:
      selectionBox[0][1] = (uint)(minY * (sizeX));
      selectionBox[0][2] = (uint)(minX * (sizeZ));
      selectionBox[1][1] = (uint)(maxY * (sizeX));
      selectionBox[1][2] = (uint)(maxX * (sizeZ));
      break;
    case 3:
      tmp = minX;
      minX = 1.0 - maxX;
      maxX = 1.0 - tmp;
      selectionBox[0][1] = (uint)(minY * (sizeX));
      selectionBox[0][2] = (uint)(minX * (sizeZ));
      selectionBox[1][1] = (uint)(maxY * (sizeX));
      selectionBox[1][2] = (uint)(maxX * (sizeZ));
      break;
    case 4:
      tmp = minY;
      minY = 1.0 - maxY;
      maxY = 1.0 - tmp;
      selectionBox[0][0] = (uint)(minX * (sizeY));
      selectionBox[0][2] = (uint)(minY * (sizeZ));
      selectionBox[1][0] = (uint)(maxX * (sizeY));
      selectionBox[1][2] = (uint)(maxY * (sizeZ));
      break;
    case 5:
      selectionBox[0][0] = (uint)(minX * (sizeY));
      selectionBox[0][2] = (uint)(minY * (sizeZ));
      selectionBox[1][0] = (uint)(maxX * (sizeY));
      selectionBox[1][2] = (uint)(maxY * (sizeZ));
      break;
    }
  //cout << "min max : " << minX << " " << maxX << endl;
  updateScaleFactors();
}

/// Draw volume bounding box

void volume::drawBBox(bool select)
{
  uint currCols = selectionBox[1][0] - selectionBox[0][0];
  uint currRows = selectionBox[1][1] - selectionBox[0][1];
  uint currSlices = selectionBox[1][2] - selectionBox[0][2];
  uint dmax = currRows;

  if (currCols > currRows)
    dmax = currCols;
  if (currSlices > dmax)
    dmax = currSlices;

  GLfloat dv = 0.5 / (GLfloat)dmax;
  
  GLfloat minX = -1.0 * currCols * dv;
  GLfloat maxX = currCols * dv;
  GLfloat minY = -1.0 * currRows * dv;
  GLfloat maxY =  currRows * dv;
  GLfloat minZ = -1.0 * currSlices * dv;
  GLfloat maxZ =  currSlices * dv;

  if (select)
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  else
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  if (select) glPushName(0);
  glBegin(GL_QUADS);  
  //back face
  glVertex3f(minX, minY, minZ);
  glVertex3f(minX, maxY, minZ);
  glVertex3f(maxX, maxY, minZ);
  glVertex3f(maxX, minY, minZ);
  glEnd();
  if (select) glPopName();
 

  if (select) glPushName(1);
  glBegin(GL_QUADS);
  //front face
  glVertex3f(minX, minY, maxZ);
  glVertex3f(maxX, minY, maxZ);
  glVertex3f(maxX, maxY, maxZ);
  glVertex3f(minX, maxY, maxZ);
  glEnd();
  if (select) glPopName();

  if (select) glPushName(2);
  glBegin(GL_QUADS);
  //left face  cout << "pts : " << pt0 << " " << pt1 << endl;
  glVertex3f(minX, minY, maxZ);
  glVertex3f(minX, minY, minZ);
  glVertex3f(minX, maxY, minZ);
  glVertex3f(minX, maxY, maxZ);
  glEnd();
  if (select) glPopName();

  if (select) glPushName(3);
  glBegin(GL_QUADS);
  //right face
  glVertex3f(maxX, minY, maxZ);
  glVertex3f(maxX, maxY, maxZ);
  glVertex3f(maxX, maxY, minZ);
  glVertex3f(maxX, minY, minZ);
  glEnd();
  if (select) glPopName();

  if (select) glPushName(4);
  glBegin(GL_QUADS);
  //top face
  glVertex3f(minX, maxY, maxZ);
  glVertex3f(minX, maxY, minZ);
  glVertex3f(maxX, maxY, minZ);
  glVertex3f(maxX, maxY, maxZ);
  glEnd();
  if (select) glPopName();

  if (select) glPushName(5);
  glBegin(GL_QUADS);
  //bottom face
  glVertex3f(minX, minY, maxZ);
  glVertex3f(maxX, minY, maxZ);
  glVertex3f(maxX, minY, minZ);
  glVertex3f(minX, minY, minZ);
  glEnd();
  if (select) glPopName();

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

/// Draw volume bounding box

void volume::drawGrid(int f)
{
  //size of the selection box
  uint sizeX = selectionBox[1][0] - selectionBox[0][0];
  uint sizeY = selectionBox[1][1] - selectionBox[0][1];
  uint sizeZ = selectionBox[1][2] - selectionBox[0][2];
 
  //  GLfloat minX = -1.0 * sizeX * dv * 0.5;
  //  GLfloat maxX = sizeX * dv * 0.5;
  //  GLfloat minY = -1.0 * sizeY * dv * 0.5;
  //  GLfloat maxY =  sizeY * dv * 0.5;
  //  GLfloat minZ = -1.0 * sizeZ * dv * 0.5;
  //  GLfloat maxZ =  sizeZ * dv * 0.5;

  uint gridSizeX[6] = {sizeX, sizeX, sizeZ, sizeZ, sizeX, sizeX};
  uint gridSizeY[6] = {sizeY, sizeY, sizeY, sizeY, sizeZ, sizeZ};
  
  uint dmax = gridSizeX[f];
  if (gridSizeY[f] > gridSizeX[f])
    dmax = gridSizeY[f];

  GLfloat dv = 1.0 / (GLfloat)dmax;

  GLfloat minX = -1.0 * gridSizeX[f] * dv * 0.5;
  GLfloat minY = -1.0 * gridSizeY[f] * dv * 0.5;

  for (uint i = 0; i < gridSizeX[f]; ++i)
    {
      for (uint j = 0; j < gridSizeY[f]; ++j)
	{
	  glPushName(j * gridSizeX[f] + i);

	  glBegin(GL_QUADS);

	  glVertex3f(minX + dv *     i, minY + dv *     j, 1.0);
	  glVertex3f(minX + dv * (i+1), minY + dv *     j, 1.0);
	  glVertex3f(minX + dv * (i+1), minY + dv * (j+1), 1.0);
	  glVertex3f(minX + dv *     i, minY + dv * (j+1), 1.0);

	  glEnd();
	  
	  glPopName();
	}
    }
}

/// Read Medical BIN
/// @arg filename name of the file to be open
/// @arg rows number of rows of the binary file
/// @arg cols number of columns of the binary file
/// @arg nums number of slices of the binary file
/// @arg ft file type

bool volume::readMedBIN(const char* filename, const uint& _rows,
			const uint& _cols, const uint& _slices,
			const fileType& ft)
{
  uint nB = 1; //< number of Bytes for each scalar
  bool ffiles = false; //< true if it is fragmented files
  bool bigE = false; //< tells if the file is in big or little endian format
  bool pnm = false; //< true if it is a Netpbm file

  // File type selection
  if (ft == BIN8) ;
  if (ft == BIN16LE) nB = 2;
  if (ft == BIN16BE) { nB = 2; bigE = true; }
  if (ft == MULTIRAW16LE) { nB = 2; ffiles = true; }
  if (ft == MULTIRAW16BE) { nB = 2; ffiles = true; bigE = true; }
  if (ft == MULTIRAW8) { nB = 1; ffiles = true; bigE = true; }
  if (ft == SINGLERAW8) ;
  if (ft == SINGLERAW16) { nB = 2; bigE = true; }
  if (ft == PNM8) { ffiles = true; pnm = true; }

  createScalarBuffer(_rows, _cols, _slices);

  /// Get Offsets
  uint offset_x_begin = 0; uint offset_x_end = 0;
  uint offset_y_begin = 0; uint offset_y_end = 0;
  uint offset_z_begin = 0; uint offset_z_end = 0;
  uint xy_step = 1; uint z_step = 1;

  // File to be read
  ifstream volume_file;

  if (!ffiles) { // if its only one binary file
    volume_file.open(filename);
    if (volume_file.fail()) {
      cerr << "Can't open " << filename << " for reading." << endl;
      return false;
    }
  }

  for (uint z = offset_z_begin; z < (_slices - offset_z_end); z+=z_step) {

    if (ffiles) { // if its more than one binary file
      stringstream ss;
      if (pnm)
	ss << filename << "." << (z+1) << ".pnm";
      else
	ss << filename << "." << (z+1);
      volume_file.open(ss.str().c_str());
      if (volume_file.fail()) {
	cerr << "Can't open " << ss.str() << " for reading." << endl;
	return false;
      }
      if (pnm) { // jump the first 3 lines (header)
	char buf[256];
	volume_file.getline(buf, 256);
	volume_file.getline(buf, 256);
	volume_file.getline(buf, 256);
      }
    }

    for (uint y = offset_y_begin; y < (cols - offset_y_end); y+=xy_step) {

      for (uint x = offset_x_begin; x < (rows - offset_x_end); x+=xy_step) {

	// seek the file to the correct position
	if (ffiles) {
	  if (pnm)
	    volume_file.seekg( (uint)(y*rows*nB*3 + x*nB*3) );
	  else
	    volume_file.seekg( (uint)(y*rows*nB + x*nB) );
	}
	else
	  volume_file.seekg( (uint)(z*cols*rows*nB + y*rows*nB + x*nB) );

	char* buf;
	GLfloat data;
	buf = new char[nB];
	volume_file.read(buf, nB);

	if (nB > 1) {
	  if (bigE)
	    data = (GLfloat)((unsigned char)buf[1] + 256*(unsigned char)buf[0]);
	  else
	    data = (GLfloat)((unsigned char)buf[0] + 256*(unsigned char)buf[1]);
	}
	else
	  data = (GLfloat)((unsigned char)buf[0]);

	delete buf;

	addScalar(x, y, z, (float)data);
      }
    }
    if (ffiles) // if its more than one binary file
      volume_file.close();
  }

  if (!ffiles) // if its only one binary file
    volume_file.close();

  normalize();

  return true;
}
