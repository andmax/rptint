/*
 *    Render Volume GPU
 *
 *  File: cpu_frags.h
 *
 *  Authors:
 *    Andre Maximo
 *    Ricardo Marroquim
 *
 *  Last Update: Apr 24, 2006
 *
 */

/**
 * Note:
 * This header is used by volume.cc only in to compute the
 * same data (stored in outputBuffers) that is computed on
 * the first Fragment Shader, and to draw in the same way
 * as the second Fragment Shader.
 */

/// GLSL type definition

typedef GLfloat vec2 [2];
typedef GLfloat vec3 [3];
typedef GLfloat vec4 [4];

typedef GLint ivec4[4];

/// Compute the first shader on CPU

void volume::compute_1st_shader_on_cpu(void)
{
  GLfloat mv[16], mvpj[16], pj[16];

  //--- Compute the ModelViewProjection Matrix ---
  glGetFloatv(GL_MODELVIEW_MATRIX, mv);
  glGetFloatv(GL_PROJECTION_MATRIX, pj);

  glMatrixMode(GL_MODELVIEW);
  //glMatrixMode(GL_PROJECTION);
  glPushMatrix();

  glMultMatrixf(pj);
  //glMultMatrixf(mv);
  glGetFloatv(GL_MODELVIEW_MATRIX, mvpj);
  //glGetFloatv(GL_PROJECTION_MATRIX, mvpj);

  glPopMatrix();

  //--- For each tetrahedron do (fragment shader) ---
  for (uint idTet = 0; idTet < numTets; ++idTet)
    {
      //--- Global variables definitions ---
      vec4 vert_proj[4];
      vec4 vert_order[4];
      GLfloat scalar_order[4];
      GLfloat scalar_orig[4];
      ivec4 tests;
      GLfloat paramU1, paramU2;

      //--- (main) ---
      GLint id_order;
      vec3 intersectionPoint;
      GLfloat thickness;

      GLfloat centroidZ = 0.0;
      GLint count_tfan;

      GLfloat scalar_front = 0.0, scalar_back = 0.0;

      //--- Compute Vertex on Screen Space (vertex_data_retrieval) ---
      for (uint i = 0; i < 4; ++i) {

	vec4 v;

	GLuint vertId = (GLuint)tetrahedralBuffer[idTet*4 + i];

	for (uint j = 0; j < 4; ++j) {

	  v[j] = positionBuffer[vertId*4 + j];

	  vert_proj[i][j] = 0.0;

	}

	scalar_orig[i] = v[3];
	v[3] = 1.0;

	/// Apply ModelViewProjection Matrix mvpj
	for (uint r = 0; r < 4; ++r)
	  for (uint c = 0; c < 4; ++c)
	    vert_proj[i][r] += mvpj[r*4+c] * v[c];

	/// Apply ModelView Matrix mv (z -> r=2)
 	vert_proj[i][2] = 0.0;
	for (uint c = 0; c < 4; ++c)
	  vert_proj[i][2] += mv[2*4+c] * v[c];


	//if (sorting)
	  centroidZ += vert_proj[i][2];
	
	/// Update vertex array: in CPU we don't have
	/// to do the 2nd vertex shader
	for (uint j = 0; j < 4; ++j)
	  vertexArray[idTet*5*4 + (i+1)*4 + j] = vert_proj[i][j];

      }

      //if (sorting)
	centroidZ /= 4.0;
	
      //--- Compute vector crosses (pt_classification) ---
      vec4 cross;

      vec3 vc1_0, vc2_0, vc3_0, vc1_2, vc1_3;
      
      GLint num_crosses0 = 0;

      for (int i = 0; i < 3; ++i) {

	vc1_0[i] = vert_proj[1][i] - vert_proj[0][i];
	vc2_0[i] = vert_proj[2][i] - vert_proj[0][i];
	vc3_0[i] = vert_proj[3][i] - vert_proj[0][i];
	vc1_2[i] = vert_proj[1][i] - vert_proj[2][i];
	vc1_3[i] = vert_proj[1][i] - vert_proj[3][i];

      }

      cross[0] = (vc1_0[0] * vc2_0[1]) - (vc1_0[1] * vc2_0[0]);
      cross[1] = (vc1_0[0] * vc3_0[1]) - (vc1_0[1] * vc3_0[0]);
      cross[2] = (vc2_0[0] * vc3_0[1]) - (vc2_0[1] * vc3_0[0]);
      cross[3] = (vc1_2[0] * vc1_3[1]) - (vc1_2[1] * vc1_3[0]);

      for (uint i = 0; i < 4; ++i) {

	tests[i] = (cross[i] < 0.0) ? 0 : ( (cross[i] > 0.0) ? 2 : 1 );

	if (tests[i] == 1)
	  ++num_crosses0;

      }
	
      if (num_crosses0 == 2)
	count_tfan = 3;
      else if (num_crosses0 == 1)
	count_tfan = 4;
      else
	count_tfan = 5;
	
      id_order = tests[0] * 27 + tests[1] * 9 + tests[2] * 3 + tests[3] * 1;
	
      //--- Order the vertices (order_vertices) ---
      for(uint i = 0; i < 4; ++i) {

	GLuint ord = order_table[id_order][i];

	for(uint j = 0; j < 4; ++j)
	  vert_order[i][j] = vert_proj[ord][j];

	scalar_order[i] = scalar_orig[ord];

      }

      //--- Compute parameters (compute_params) ---
      GLfloat denominator, numeratorU1, numeratorU2;

      if (count_tfan == 3)
	{
	  paramU1 = 1.0;
	  paramU2 = 1.0;
	}
      else if (count_tfan == 4)
	{
	  //line intersection denominator between v0->v2 and v1->v3
	  denominator = ((vert_order[3][1] - vert_order[1][1]) * (vert_order[2][0] - vert_order[0][0])) -
	    ((vert_order[3][0] - vert_order[1][0]) * (vert_order[2][1] - vert_order[0][1]));
	  
	  //line defined by vector v1->v3
	  numeratorU2 = ((vert_order[2][0] - vert_order[0][0]) * (vert_order[0][1] - vert_order[1][1])) -
	    ((vert_order[2][1] - vert_order[0][1]) * (vert_order[0][0] - vert_order[1][0]));
	  
	  paramU1 = 1.0;
	  paramU2 = numeratorU2 / denominator;
	}
      else
	{
	  //line intersection denominator between v0->v2 and v1->v3
	  denominator = ((vert_order[3][1] - vert_order[1][1]) * (vert_order[2][0] - vert_order[0][0])) -
	    ((vert_order[3][0] - vert_order[1][0]) * (vert_order[2][1] - vert_order[0][1]));
	  
	  //line defined by vector v0->v2
	  numeratorU1 = ((vert_order[3][0] - vert_order[1][0]) * (vert_order[0][1] - vert_order[1][1])) -
	    ((vert_order[3][1] - vert_order[1][1]) * (vert_order[0][0] - vert_order[1][0]));
	  
	  //line defined by vector v1->v3
	  numeratorU2 = ((vert_order[2][0] - vert_order[0][0]) * (vert_order[0][1] - vert_order[1][1])) -
	    ((vert_order[2][1] - vert_order[0][1]) * (vert_order[0][0] - vert_order[1][0]));
	  
	  paramU1 = numeratorU1 / denominator;
	  paramU2 = numeratorU2 / denominator;
	}
	
      //--- Compute intersection point (compute_intersection) ---
      for (int i = 0; i < 3; ++i)
	intersectionPoint[i] = 0.0;

      if (count_tfan == 3)
	{
	  thickness = (vert_order[0][2] - vert_order[1][2]);
	}
      else if (count_tfan == 4)
	{
	  GLfloat zBackIntersection = vert_order[1][2] + paramU2*(vert_order[3][2] - vert_order[1][2]);
	  thickness = fabs(vert_order[2][2] - zBackIntersection);
	}
      else
	{
	  //find z coordinate of back intersection point by interpolating original vertices (not projected)
	  GLfloat zBackIntersection = vert_order[1][2] + paramU2*(vert_order[3][2] - vert_order[1][2]);
	  
	  //find ordered intersection point between the two ordered lines (basis graph)
	  for(uint i = 0; i < 3; ++i)
	    intersectionPoint[i] = (vert_order[0][i] + paramU1*(vert_order[2][i] - vert_order[0][i]));
	  
	  thickness = (intersectionPoint[2] - zBackIntersection);
	}
      
      //--- In main (outside the compute_intersection) ---

      //if Class 1 then paramU1 is greater than 1.0 (middle vertex inside projected triangle)
      //in this case the thickness value must be multiplied by the
      //ratio r=(|V0V4|/|V0I|), that is, r=1/paramU1,
      //where V0V4 is the distance from V0 to the middle vertex
      //and V0I is the distance from V0 to the intersection point
      if (paramU1 >= 1.0)
	{
	  if (count_tfan == 5) {
	    thickness /= paramU1;
	    paramU1 = 1.0 / paramU1;
	  }
	}
      else // count == 6 ; thick vertex == intersection vertex
	{
	  count_tfan = 6;
	}
	
      //--- Compute the scalars (compute_scalars) ---
	
      if (count_tfan == 6)
	{
	  scalar_front = scalar_order[0] + paramU1*(scalar_order[2] - scalar_order[0]);
	  scalar_back = scalar_order[1] + paramU2*(scalar_order[3] - scalar_order[1]);
	}
      else if (count_tfan == 5)
	{
	  scalar_front = scalar_order[2];
	  GLfloat sb_tmp = scalar_order[1] + paramU2*(scalar_order[3] - scalar_order[1]);
	  scalar_back = scalar_order[0] + (sb_tmp - scalar_order[0])*paramU1;
	}
      else if (count_tfan == 4)
	{
	  scalar_front = scalar_order[2];
	  scalar_back = scalar_order[1] + paramU2*(scalar_order[3] - scalar_order[1]);
	}
      else if (count_tfan == 3)
	{
	  scalar_front = scalar_order[0];
	  scalar_back = scalar_order[1];
	}
      
      if (thickness < 0.0)
	{
	  GLfloat tmp = scalar_back;
	  scalar_back = scalar_front;
	  scalar_front = tmp;  
	}
      
      //--- Write output data in buffers ---
	
      outputBuffer0[idTet*4 + 0] = intersectionPoint[0];
      outputBuffer0[idTet*4 + 1] = intersectionPoint[1];
      outputBuffer0[idTet*4 + 2] = centroidZ;
      outputBuffer0[idTet*4 + 3] = id_order;
      
      outputBuffer1[idTet*4 + 0] = scalar_front;
      outputBuffer1[idTet*4 + 1] = scalar_back;
      outputBuffer1[idTet*4 + 2] = fabs(thickness);
      outputBuffer1[idTet*4 + 3] = count_tfan;
    } // for numTets
}

/// Compute the second Shader on CPU

void volume::compute_2nd_shader_on_cpu(void)
{
  ///--- Erase the Modelview and Projection Matrix ---
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();

  //--- For each tetrahedron do (glMultiDrawElements) ---
  for (uint i = 0; i < numTets; ++i) {

    glBegin(GL_TRIANGLE_FAN);

    //--- For each triangle fan do (glMultiDrawElements) ---
    for (GLint j = 0; j < count[i]; ++j) {

      //--- Returned value ---
      TColor color;

      //--- Second fragment shader (main) ---
      GLfloat sf = colorArray[indices[i][j] * 4 + 1];
      GLfloat sb = colorArray[indices[i][j] * 4 + 2];
      GLfloat l  = colorArray[indices[i][j] * 4 + 3];

      if (integrating) {
	//--- With Integration ---
	TColor colorFront = tf[ (int) (sf * 255) ];
	TColor colorBack = tf[ (int) (sb * 255) ];

	vec2 tau;

	tau[0] = colorFront.a * l;
	tau[1] = colorBack.a * l;

	GLfloat zeta = exp( -(tau[0] * 0.5 + tau[1] * 0.5) );

	vec2 gamma;

	gamma[0] = tau[0] / (1.0 + tau[0]);
	gamma[1] = tau[1] / (1.0 + tau[1]);
	
	GLfloat psi = psiGammaTable[(uint)(gamma[0] + 0.5/(GLfloat)preIntTexSize)*preIntTexSize]
	  [(uint)(gamma[1] + 0.5/(GLfloat)preIntTexSize)*preIntTexSize];

	color[0] = colorFront[0]*(1.0 - psi) + colorBack[0]*(psi - zeta);
	color[1] = colorFront[1]*(1.0 - psi) + colorBack[1]*(psi - zeta);
	color[2] = colorFront[2]*(1.0 - psi) + colorBack[2]*(psi - zeta);
	color[3] = 1.0 - zeta;

      }
      else {
	//--- NO Integration ---
	GLfloat sT = (sf + sb) * 0.5; // average scalar
	color = tf[ (int) (sT*255) ];
	color[3] = 1.0 - exp( -l * color[3] );

	for (int j = 0; j < 3; ++j)
	  color[j] *= color[3];
      }
      /*
      /// Apply ModelViewProjection only in no-intersection vertices
      vec4 v;

      for (uint k = 0; k < 4; ++k)
	v[k] = 0.0;

      if (vertexArray[indices[i][j] + 3] == 0.0) {
	/// Do not multiply any Matrix
	for (uint k = 0; k < 4; ++k)
	  v[k] = vertexArray[indices[i][j] + k];
      }
      else {
	/// Multiply by ModelViewProjection Matrix mvpj
	for (uint r = 0; r < 4; ++r)
	  for (uint c = 0; c < 4; ++c)
	    v[r] += mvpj[r*4+c] * vertexArray[indices[i][j] + c];
      }
      */

      glColor4f(color[0], color[1], color[2], color[3]);
      //glVertex4f(v[0], v[1], v[2], 1.0);
      
      glVertex4f(vertexArray[indices[i][j] * 4 + 0],
		 vertexArray[indices[i][j] * 4 + 1],
		 vertexArray[indices[i][j] * 4 + 2],
		 1.0);

    } // j

    glEnd();

  } // i

  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}
