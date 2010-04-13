/**
 *
 *    Render Volume GPU
 *
 *  File: Trasfer Function
 *
 *  Authors:
 *    Andre Maximo
 *    Ricardo Marroquim
 *
 *  Created: May 5, 2006
 *
 */

#include "transferFunction.h"

/*
  Update the position and value of the picked point
  in the transfer function window
*/
void transferFunction::updateTF(int x, int y)
{
  if (picked_point != -1)
    {    
      if (picked_point == num_colors)
	{
	  GLfloat x_ortho = xOrtho(x);
	  if (x_ortho > 1.0) x_ortho = 1.0;
	  if (x_ortho < 0.0) x_ortho = 0.0;
	  brightness = (GLfloat)( x_ortho * BRIGHTNESS_MAX );
	}
      else
	{
	  GLfloat y_ortho = yOrtho(y);	 
	  if (y_ortho > 1.0) y_ortho = 1.0;
	  if (y_ortho < 0.0) y_ortho = 0.0;

	  int pId = getCurScalar();
	  colors[pId][3] = y_ortho;
	}
    }
}

bool transferFunction::writeTF(char* filename)
{

  cout << "Wrote transfer function to file" << endl;

  ofstream output(filename);

  if(output.fail())
    return false;

  output << "# Generated transfer function #" << endl;
  output << "%opacity" << endl << endl;
  output << "0   0.0" << endl;
  output << "255 0 0" << endl;
  output << "%voxel" << endl << endl;

  for(int i = 0; i < tf_size; i++)
    {
      output << i << " " << colors[i] << endl;
    }

  output.close();

  return true;
}

bool transferFunction::readTF(char* filename)
{
  int index;
  TColor c;
  char buffer[200];

  ifstream input(filename);

  if(input.fail())
    return false;

  //read the first 7 lines containing header and general opacity info
  for (int l = 0; l < 7; ++l)
    {
      input.getline(buffer, 200);
    }

  for(int i = 0; i < tf_size; i++)
    {
      input >> index >> c;
      colors[i][3] = c[3];
    }

  input.close();

  return true;
}

GLfloat transferFunction::gaussian(GLfloat x)
{
  GLfloat rho = 1.0;
  GLfloat amp = (GLfloat)( 1.0 / (sqrt(2 * 3.14 * rho)) );

  GLfloat exp = (GLfloat)pow(2.71, (-1*(x*x))/(rho*rho));
  GLfloat g = (GLfloat)( amp * exp );

  return g;
}

/*
  Computes the values of every transfer function point
*/

void transferFunction::computeTF(void)
{
  int xMax, xMin;
  GLfloat yMax, yMin, a, b;

  for (int layer = 0; layer < num_colors - 1; ++layer)
    {
		
      xMin = (int)(((GLfloat)layer)*layer_size - 1.0);
      xMax = (int)(((GLfloat)layer+1)*layer_size - 1.0);
      if (xMin == -1) xMin = 0;
      if (xMax > 253) xMax = 255;
      yMin = colors[xMin][3];
      yMax = colors[xMax][3];
      a = (yMax - yMin)/(GLfloat)(xMax - xMin);
      b = yMin - a*xMin;
      for (int i = xMin; i < xMax; ++i)
	{
	  colors[i][3] = a*i + b;
	  //colors[i][3] = gaussian(a*i +b);
	}
    }
	
}

void transferFunction::computeChromacity(void)
{
  int color_range = 0, offset = 0;
  double e = exp(1), sd = 0.0, su = 0.0;
  double sd1, sd2, su1, su2; //step down and step up values

  switch (colorCode) {
  case 0:
    num_colors = 16;
    layer_size = tf_size / (GLfloat) (num_colors-1);
    color_range = (int)tf_size;
    for (int i = 0; i <= color_range; ++i) {
      sd = (GLfloat)(1.0 - log(1.0 + (e - 1.0) * ((double)i / ((double)color_range - 1.0))));
      //(1, 1, 1) -> (0, 0, 0)
      colors[i][0] = sd;
      colors[i][1] = sd;
      colors[i][2] = sd;
    }
    break;
  case 1:
    num_colors = 16;
    layer_size = tf_size / (GLfloat) (num_colors-1);
    color_range = (int)tf_size;
    for (int i = 0; i <= color_range; ++i) {
      su = (GLfloat)log(1.0 + (e - 1.0) * ((double)i / ((double)color_range - 1.0)));
      //(0, 0, 0) -> (1, 1, 1)
      colors[i][0] = su;
      colors[i][1] = su;
      colors[i][2] = su;
    }
    break;
//   case 1:
//     for (int i = 0; i < tf_size; ++i) {
//       colors[i][0] = (GLfloat)log(1.0 + (e - 1.0) * (i / (GLfloat)(tf_size - 1.0)));
//       colors[i][1] = 0.0;
//       colors[i][2] = (GLfloat)(1.0 - log(1.0 + (e - 1.0) * (i / (GLfloat)(tf_size - 1.0))));
//     }
//     break;
  case 2:
    for (int i = 0; i < tf_size; ++i) {
      colors[i][0] = (GLfloat)(1.0 - log(1.0 + (e - 1.0) * (i / (GLfloat)(tf_size - 1.0))));
      colors[i][1] = 0.0;
      colors[i][2] = (GLfloat)log(1.0 + (e - 1.0) * (i / (GLfloat)(tf_size - 1.0)));
    }
    break;
  case 3:
    num_colors = 6;
    layer_size = tf_size / (GLfloat) (num_colors-1);
    color_range = (int)layer_size;
    for (int i = 0; i <= color_range; ++i) {
      su = (GLfloat)log(1.0 + (e - 1.0) * ((double)i / ((double)color_range - 1.0)));
      sd = (GLfloat)(1.0 - log(1.0 + (e - 1.0) * ((double)i / ((double)color_range - 1.0))));
      //(0, 1, 1) -> (0, 0, 1)
      offset = 0;
      colors[i + offset][0] = 0.0;
      colors[i + offset][1] = sd;
      colors[i + offset][2] = 1.0;
      //(0, 0, 1) -> (1, 0, 1)
      offset = color_range - 1;
      colors[i + offset][0] = su;
      colors[i + offset][1] = 0.0;
      colors[i + offset][2] = 1.0;      
      //(1, 0, 1) -> (0, 1, 0)
      offset = color_range*2 - 1;
      colors[i + offset][0] = sd;
      colors[i + offset][1] = su;
      colors[i + offset][2] = sd;
      //(0, 1, 0) -> (1, 1, 0)
      offset = color_range*3 - 1;
      colors[i + offset][0] = su;
      colors[i + offset][1] = 1.0;
      colors[i + offset][2] = 0.0;
      //(1, 1, 0) -> (1, 0, 0)
      offset = color_range*4 - 1;
      colors[i + offset][0] = 1.0;
      colors[i + offset][1] = sd;
      colors[i + offset][2] = 0.0;
    }
    colors[255][0] = 1.0;
    colors[255][1] = 0.0;
    colors[255][2] = 0.0;
    break;
  case 4:
    num_colors = 6;
    layer_size = tf_size / (GLfloat) (num_colors-1);
    color_range = (int)layer_size;
    for (int i = 0; i <= color_range; ++i) {
      su = (GLfloat)log(1.0 + (e - 1.0) * ((double)i / ((double)color_range - 1.0)));
      sd = (GLfloat)(1.0 - log(1.0 + (e - 1.0) * ((double)i / ((double)color_range - 1.0))));
      //(1, 0, 0) -> (1, 1, 0)
      offset = 0;
      colors[i + offset][0] = 1.0;
      colors[i + offset][1] = su;
      colors[i + offset][2] = 0.0;
      //(1, 1, 0) -> (0, 1, 0)
      offset = color_range*1 - 1;
      colors[i + offset][0] = sd;
      colors[i + offset][1] = 1.0;
      colors[i + offset][2] = 0.0;
      //(0, 1, 0) -> (1, 0, 1)
      offset = color_range*2 - 1;
      colors[i + offset][0] = su;
      colors[i + offset][1] = sd;
      colors[i + offset][2] = su;
      //(1, 0, 1) -> (0, 0, 1)
      offset = color_range*3 - 1;
      colors[i + offset][0] = sd;
      colors[i + offset][1] = 0.0;
      colors[i + offset][2] = 1.0;   
      //(0, 0, 1) -> (0, 1, 1)
      offset = color_range*4 - 1;
      colors[i + offset][0] = 0.0;
      colors[i + offset][1] = su;
      colors[i + offset][2] = 1.0;
    }
    colors[255][0] = 0.0;
    colors[255][1] = 1.0;
    colors[255][2] = 1.0;
    break;
  case 5:
    num_colors = 5;
    layer_size = tf_size / (GLfloat) (num_colors-1);
    color_range = (int)layer_size;
    for (int i = 0; i <= color_range; ++i) {
      su = (GLfloat)log(1.0 + (e - 1.0) * ((double)i / ((double)color_range - 1.0)));
      sd = (GLfloat)(1.0 - log(1.0 + (e - 1.0) * ((double)i / ((double)color_range - 1.0))));
      //(1, 0, 0) -> (1, 1, 0)
      offset = 0;
      colors[i + offset][0] = 1.0;
      colors[i + offset][1] = su;
      colors[i + offset][2] = 0.0;
      //(1, 1, 0) -> (0, 1, 0)
      offset = color_range - 1;
      colors[i + offset][0] = sd;
      colors[i + offset][1] = 1.0;
      colors[i + offset][2] = 0.0;
      //(0, 1, 0) -> (0, 1, 1)
      offset = color_range*2 - 1;
      colors[i + offset][0] = 0.0;
      colors[i + offset][1] = 1.0;
      colors[i + offset][2] = su;
      //(0, 1, 1) -> (0, 0, 1)
      offset = color_range*3 - 1;
      colors[i + offset][0] = 0.0;
      colors[i + offset][1] = sd;
      colors[i + offset][2] = 1.0;
    }
    colors[255][0] = 0.0;
    colors[255][1] = 0.0;
    colors[255][2] = 1.0;
    break;
  case 6:
    num_colors = 5;
    layer_size = tf_size / (GLfloat) (num_colors-1);
    color_range = (int)layer_size;
    for (int i = 0; i <= color_range; ++i) {
      su = (GLfloat)log(1.0 + (e - 1.0) * ((double)i / ((double)color_range - 1.0)));
      sd = (GLfloat)(1.0 - log(1.0 + (e - 1.0) * ((double)i / ((double)color_range - 1.0))));
      //(0, 0, 1) -> (0, 1, 1)
      offset = 0;
      colors[i + offset][0] = 0.0;
      colors[i + offset][1] = su;
      colors[i + offset][2] = 1.0;
      //(0, 1, 1) -> (0, 1, 0)
      offset = color_range - 1;
      colors[i + offset][0] = 0.0;
      colors[i + offset][1] = 1.0;
      colors[i + offset][2] = sd;
      //(0, 1, 0) -> (1, 1, 0)
      offset = color_range*2 - 1;
      colors[i + offset][0] = su;
      colors[i + offset][1] = 1.0;
      colors[i + offset][2] = 0.0;
      //(1, 1, 0) -> (1, 0, 0)
      offset = color_range*3 - 1;
      colors[i + offset][0] = 1.0;
      colors[i + offset][1] = sd;
      colors[i + offset][2] = 0.0;
    }
    colors[255][0] = 1.0;
    colors[255][1] = 0.0;
    colors[255][2] = 0.0;
    break;
 case 7:
   num_colors = 9;
   layer_size = tf_size / (GLfloat) (num_colors-1);
   color_range = (int)layer_size;
   for (int i = 0; i <= color_range; ++i) {
     su1 = (GLfloat)0.5*log(1.0 + (e - 1.0) * ((double)i / ((double)color_range - 1.0)));
     sd1 = (GLfloat)(0.5 - 0.5*log(1.0 + (e - 1.0) * ((double)i / ((double)color_range - 1.0))));
     su2 = (GLfloat)0.5 + 0.5*log(1.0 + (e - 1.0) * ((double)i / ((double)color_range - 1.0)));
     sd2 = (GLfloat)(1.0 - 0.5*log(1.0 + (e - 1.0) * ((double)i / ((double)color_range - 1.0))));
     //(1, 0, 0) -> (1, 0.5, 0)
     offset = 0;
     colors[i + offset][0] = 1.0;
     colors[i + offset][1] = su1;
     colors[i + offset][2] = 0.0;
     //(1, 0.5, 0) -> (1, 1, 0)
     offset = color_range - 1;
     colors[i + offset][0] = 1.0;
     colors[i + offset][1] = su2;
     colors[i + offset][2] = 0.0;
     //(1, 1, 0) -> (0.5, 1, 0)
     offset = color_range*2 - 1;
     colors[i + offset][0] = sd1;
     colors[i + offset][1] = 1.0;
     colors[i + offset][2] = 0.0;
     //(0.5, 1, 0) -> (0, 1, 0)
     offset = color_range*3 - 1;
     colors[i + offset][0] = sd2;
     colors[i + offset][1] = 1.0;
     colors[i + offset][2] = 0.0;
     //(0, 1, 0) -> (0, 1, 0.5)
     offset = color_range*4 - 1;
     colors[i + offset][0] = 0.0;
     colors[i + offset][1] = 1.0;
     colors[i + offset][2] = su1;
     //(0, 1, 0.5) -> (0, 1, 1)
     offset = color_range*5 - 1;
     colors[i + offset][0] = 0.0;
     colors[i + offset][1] = 1.0;
     colors[i + offset][2] = su2;
     //(0, 1, 1) -> (0, 0.5, 1)
     offset = color_range*6 - 1;
     colors[i + offset][0] = 0.0;
     colors[i + offset][1] = sd1;
     colors[i + offset][2] = 1.0;
     //(0, 0.5, 1) -> (0, 1, 1)
     offset = color_range*7 - 1;
     colors[i + offset][0] = 0.0;
     colors[i + offset][1] = sd2;
     colors[i + offset][2] = 1.0;
   }
   colors[255][0] = 0.0;
   colors[255][1] = 0.0;
   colors[255][2] = 1.0;
   break;
 case 8:
   for (int i = 0; i < tf_size; ++i) {
     colors[i][0] = (GLfloat)log(1.0 + (e - 1.0) * (i / (GLfloat)(tf_size - 1.0)));
     colors[i][1] = (GLfloat)(1.0 - log(1.0 + (e - 1.0) * (i / (GLfloat)(tf_size - 1.0))));
     colors[i][2] = 0.0;
   }
   break;
 case 9:
   num_colors = 6;
   layer_size = tf_size / (GLfloat) (num_colors-1);
   color_range = (int)layer_size;
   for (int i = 0; i <= color_range; ++i) {
     //(1, 0, 0)
     offset = 0;
     colors[i + offset][0] = 1.0;
     colors[i + offset][1] = 0.0;
     colors[i + offset][2] = 0.0;
     //(1, 1, 0)
     offset = color_range - 1;
     colors[i + offset][0] = 1.0;
     colors[i + offset][1] = 1.0;
     colors[i + offset][2] = 0.0;
     //(0, 1, 0)
     offset = color_range*2 - 1;
     colors[i + offset][0] = 0.0;
     colors[i + offset][1] = 1.0;
     colors[i + offset][2] = 0.0;
     //(0, 1, 1)
     offset = color_range*3 - 1;
     colors[i + offset][0] = 0.0;
     colors[i + offset][1] = 1.0;
     colors[i + offset][2] = 1.0;
     //(0, 0, 1)
     offset = color_range*4 - 1;
     colors[i + offset][0] = 0.0;
     colors[i + offset][1] = 0.0;
     colors[i + offset][2] = 1.0;
   }
   colors[255][0] = 0.0;
   colors[255][1] = 0.0;
   colors[255][2] = 1.0;
   break;
  }
}

/// Display
void transferFunction::draw(void)
{
  int pId = 0, pIdNext = 0;
  glBegin(GL_QUADS);
  for (int i = 0; i < num_colors - 1; ++i)
    {
      //pId = (int)(layer_size*(GLfloat)i) - 1;
      pId = getScalar(i);
      //pIdNext = (int)(layer_size*((GLfloat)(i+1))) - 1;
      pIdNext = getScalar(i+1);

      //if (pId == -1) pId = 0;

      glColor3f(colors[pId][0], colors[pId][1], colors[pId][2]);
      glVertex2d((layer_size*i)/(GLfloat)tf_size, 0.0);
      glVertex2d((layer_size*i)/(GLfloat)tf_size, colors[pId][3]);

      glColor3f(colors[pIdNext][0], colors[pIdNext][1], colors[pIdNext][2]);
      glVertex2d((layer_size*(i+1))/(GLfloat)tf_size, colors[pIdNext][3]);
      glVertex2d((layer_size*(i+1))/(GLfloat)tf_size, 0.0); 
    }
  glEnd();

  /// Draw Brightness Line

  glBegin(GL_LINES);
  glColor3f(1.0, 1.0, 1.0);
  glVertex2d(0.0, BRIGHTNESS_Y_POS);
  glColor3f(0.0, 0.0, 0.0);
  glVertex2d(1.0, BRIGHTNESS_Y_POS);
  glEnd();

  /// Draw Control Points

  glPointSize(5.0);
  glBegin(GL_POINTS);	
  for (int i = 0; i < num_colors; ++i)
    {
      glColor3f(0.0, 0.0, 0.0);
      if (i == picked_point)
	glColor3f(0.7, 0.3, 0.1);

      //pId = (int)(layer_size*(GLfloat)i) - 1;
      //if (pId == -1) pId = 0;
      pId = getScalar(i);
      glVertex2d((layer_size*i)/(GLfloat)tf_size, colors[pId][3]);
    }	

  glColor3f(0.0, 0.0, 0.0);
  glVertex2d(brightness/BRIGHTNESS_MAX, BRIGHTNESS_Y_POS); /// Brightness point
  glEnd();

  /// Draw Axis

  glColor3f(0.0, 0.0, 0.0);
  glBegin(GL_LINES);
  glVertex2d(0.0, 0.0);
  glVertex2d(0.0, 1.1);
  glVertex2d(0.0, 0.0);
  glVertex2d(1.1, 0.0);
	
  glVertex2d(0.0, 1.1);
  glVertex2d(0.02, 1.08);
  glVertex2d(0.0, 1.1);
  glVertex2d(-0.02, 1.08);
	
  glVertex2d(1.1, 0.0);
  glVertex2d(1.08, -0.02);
  glVertex2d(1.1, 0.0);
  glVertex2d(1.08, 0.02);
  glEnd();
	
  glEnable(GL_LINE_STIPPLE);
  glLineStipple(4, 0xAAAA);
  glBegin(GL_LINES);
  glVertex2d(0.0, 1.0);
  glVertex2d(1.0, 1.0);
  glVertex2d(1.0, 0.0);
  glVertex2d(1.0, 1.0);
  glEnd();
  glDisable(GL_LINE_STIPPLE);
}

/*
  Selects a transfer function control point on the screen
*/
void transferFunction::pick(int x, int y)
{
  int pId = 0;
  picked_point = -1;
  for (int i = 0; i < num_colors; ++i)
    {
      if ((x >= xScreen((layer_size*i)/(GLfloat)tf_size) - TF_PICK_BOX) &&
	  (x <= xScreen((layer_size*i)/(GLfloat)tf_size) + TF_PICK_BOX))
	{
	  //pId = (int)(layer_size*i) - 1;
	  //if (pId == -1) pId = 0;
	  pId = getScalar(i);

	  //if ((y >= yScreen(colors[pId][3]) - TF_PICK_BOX) && (y <= yScreen(colors[pId][3]) + TF_PICK_BOX))
	  if ((y >= yScreen(1.0) && (y <= yScreen(0.0)) ) )
	    {
	      picked_point = i;
	      return;
	    }
	}
    }
  if ((y >= yScreen(BRIGHTNESS_Y_POS) - TF_PICK_BOX) &&
      (y <= yScreen(BRIGHTNESS_Y_POS) + TF_PICK_BOX))
    {
      if ( (x >= xScreen((GLfloat)((brightness/BRIGHTNESS_MAX) - TF_PICK_BOX)) )
	   && (x <= xScreen((GLfloat)((brightness/BRIGHTNESS_MAX) + TF_PICK_BOX))) )
	{
	  picked_point = num_colors;
	  return;
	}
    }
}
