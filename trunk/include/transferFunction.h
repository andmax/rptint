/**
 *
 *    Render Volume GPU
 *
 *  File: transferFunction.h
 *
 *  Authors:
 *    Andre Maximo
 *    Ricardo Marroquim
 *
 *  Created: May 06, 2006
 *
 */
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

#include <GL/glut.h>

using namespace std;

typedef unsigned int uint;

/// Default values

#define INIT_TF_SIZE 256
#define INIT_NUM_COLORS 7
#define INIT_COLOR_CODE 7
#define INIT_BRIGHTNESS 1.0

#define TF_PICK_BOX 5
#define BRIGHTNESS_Y_POS -0.15
#define BRIGHTNESS_MAX 30.0

#define TF_ORTHO_MIN -0.2
#define TF_ORTHO_MAX 1.3
#define TF_ORTHO_SIZE (TF_ORTHO_MAX - TF_ORTHO_MIN)

/// Extern screen size

extern GLint tfWinWidth, tfWinHeight;

/// A simple Color Type for the double 4-array RGBA

struct TColor {
  /// Constructor
  TColor() : r(0.0), g(0.0), b(0.0), a(0.0) { }

  /// Constructor
  /// @arg RGBA _r _g _b _a color values
  TColor(GLfloat _r, GLfloat _g, GLfloat _b, GLfloat _a) : r(_r), g(_g), b(_b), a(_a) { }

  /// Destructor
  ~TColor() { }

  /// Operator to grant write/read access to the color values
  GLfloat& operator [](const uint i) {
    return (i==0) ? this->r : ( (i==1) ? this->g : ( (i==2) ? this->b : ( (i==3) ? this->a : (this->a) ) ) );
  }

  /// get color component
  GLfloat getValue(const uint i) {
    return (i==0) ? r : ( (i==1) ? g : ( (i==2) ? b : ( (i==3) ? a : (0.0) ) ) );
  }

  /// Operator to output the color values
  inline friend ostream& operator << (ostream& out, const TColor& c) {
    out << c.r << " " << c.g << " " << c.b << " " << c.a;
    return out;
  }

  /// Operator to input the color values
  inline friend istream& operator >> (istream& in, TColor& c) {
    in >> c.r >> c.g >> c.b >> c.a;
    return in;
  }

  GLfloat r, g, b, a; ///< RGBA color values
};

/// A class to handle the Transfer Function

class transferFunction
{
 public:
  /// Constructor
  transferFunction(GLfloat _bt = INIT_BRIGHTNESS, int _cc = INIT_COLOR_CODE,
		   int _nc = INIT_NUM_COLORS, int _size = INIT_TF_SIZE) :
    colorCode(_cc), picked_point(-1),  num_colors (_nc),
    tf_size(_size), brightness(_bt), layer_size(_size / (GLfloat)_nc)
  {
    colors = new TColor[(int)tf_size];
    computeColors();
  }

  /// Destructor
  ~transferFunction() { delete []colors; }

  void pick(int x, int y);
  void draw(void);    

  void computeTF(void);	
  bool writeTF(char* fn);
  bool readTF(char* fn);
  void updateTF(int x, int y);

  /// Operator to grant write/read access to the color values
  TColor& operator [](const uint i) { return this->colors[i]; }

  /// get color
  TColor getColor(const uint i) { return colors[i]; }

  int getScalar(int i) {
    int pId = (int)(layer_size*i) - 1;
    if (pId == -1) pId = 0;
    //if (pId == 254) pId = 255;
    return pId;    
  }

  void updateBrightness(GLfloat delta) {
    brightness += delta;
    if (brightness < 0.0) brightness = 0.0;
    if (brightness > 30.0) brightness = 30.0;
  }

  inline int getCurScalar(void) { return getScalar(picked_point); }

  inline void setPickedPointNull(void) { picked_point = -1; }

  inline int getPickedPoint(void) { return picked_point; }
  inline int getBrightnessId(void) { return num_colors; }
  inline GLfloat getCurAlpha(void) { return colors[getCurScalar()][3]; }		
  inline GLfloat getBrightness(void) { return brightness; }

  void setColorCode(int cc) { colorCode = cc; computeColors(); }
  void computeColors(void) { computeChromacity(); computeTF(); }

 private:
  int colorCode, picked_point, num_colors, tf_size;
  GLfloat brightness, layer_size;
  TColor *colors;

  void setColor(int c);
  void computeChromacity(void);
  void addColor(int index, TColor c) { colors[index] = c; }
  GLfloat gaussian(GLfloat x);

  inline GLfloat xOrtho(int x) { return (GLfloat)( ((x/(GLfloat)tfWinWidth) * TF_ORTHO_SIZE) + TF_ORTHO_MIN ); }
  inline GLfloat yOrtho(int y) { return (GLfloat)(((((GLfloat)tfWinHeight - y)/(GLfloat)tfWinHeight) * TF_ORTHO_SIZE) + TF_ORTHO_MIN); }
  inline int xScreen(GLfloat x) { return (int)(((x - TF_ORTHO_MIN)*(GLfloat)tfWinWidth)/TF_ORTHO_SIZE); }
  inline int yScreen(GLfloat y) { return tfWinHeight - (int)(((y - TF_ORTHO_MIN)*(GLfloat)tfWinHeight)/TF_ORTHO_SIZE); }

};
