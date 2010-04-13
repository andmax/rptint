/***********************************************************************************
 *                              Virtual Trackball
 *  Created by: Yalmar Ponce (yalmar@gmail.com|http://www.lcg.ufrj.br/Members/yalmar)
 *  Last Update: 06-05-2006
 *  --------------------------------------------------------------------------------
 *  This file is a part of work done by Yalmar. You is free to use this code in any 
 *  way you like. However, i expect you puts this note, how credits of the author
 *
 ***********************************************************************************/
 
#ifndef __F_TRACKBALL__
#define __F_TRACKBALL__

#include <GL/glut.h>
#include <math.h>

#include <iostream>

#define TRACKBALLSIZE  (0.8)
#define RENORMCOUNT 97

/***************************************************************************************************/

void trackball(double q[4], double p1x, double p1y, double p2x, double p2y);

/*
 * Given two quaternions, add them together to get a third quaternion.
 * Adding quaternions to get a compound rotation is analagous to adding
 * translations to get a compound translation.  When incrementally
 * adding rotations, the first argument here should be the new
 * rotation, the second and third the total rotation (which will be
 * over-written with the resulting new total rotation).
 */
void add_quats(double *q1, double *q2, double *dest);

/*
 * A useful function, builds a rotation matrix in Matrix based on
 * given quaternion.
 */
void build_rotmatrix(double m[4][4], double q[4]);

/*
 * This function computes a quaternion based on an axis (defined by
 * the given vector) and an angle about which to rotate.  The angle is
 * expressed in radians.  The result is put into the third argument.
 */
void axis_to_quat(double a[3], double phi, double q[4]);


//!This class generate a TrackBall
class FTrackBall
{
 private:
  double rot_matrix[4][4];
  double quat[4];
  double spin_quat[4];

  //! the zoom of the trackball
  float z;

  //! older coordinate of the mouse in x
  int oldx, poldx;
  //! older coordinate of the mouse in y
  int oldy, poldy;  
    
  //! The center of the trackball, where we look
  float center[3];

  //@{ \name Constructor
 public:
  /*! Default Constructor*/
  FTrackBall(float zoom = 29);
  //@}

  //@{ \name Functions
 public:
  /*! Set the zoom to the trackball*/
  void setZoom(float zoom);

  //sets the quaternion
  void setQuat(double q[4]);

  /*! Return the zoom of the trackball*/
  double getZoom(void) {return z;}

  /*! Set the old coordinates of the mouse position */
  void setOld(int  x, int  y);

  /** Calcul the Track-ball
   * You have to call this function when you move the mous
   * \param x : the x mouse position
   * \param y : the x mouse position
   * \param width : the width of the gl context
   * \param height : the height of the gl context*/
  void trackBallInAction(int x, int y, int width, int height);

  /** When you zoom
   * Most ofently use : with the right button of the mouse. 
   * \param y : the y mouse position*/		    
  void zoomInAction(int y);
		
  /** This function mult the matrix to the current opengl matrix
   * Most ofently use : before the rendering of the scene*/
  void multToOpenGL();

  void getRotMatrix(double **m);
    
  /*! Reset the trackball*/
  void reset(float zoom = 10);		
		
  /*! Return the older position of the pointer (x coord) */
  int getOldX();
		
  /*! Return the older position of the pointer (y coord) */
  int getOldY();   
    
  //! in order to achieve a smoth idle trackball  
  void prepareToIdle(int x, int y, int width, int height);
    
  /*! Idle mousetrackball */
  void idle();
        
  //@}

};

void trackball(double q[4], double p1x, double p1y, double p2x, double p2y);

/*
 * Given two quaternions, add them together to get a third quaternion.
 * Adding quaternions to get a compound rotation is analagous to adding
 * translations to get a compound translation.  When incrementally
 * adding rotations, the first argument here should be the new
 * rotation, the second and third the total rotation (which will be
 * over-written with the resulting new total rotation).
 */
void add_quats(double *q1, double *q2, double *dest);

/*
 * A useful function, builds a rotation matrix in Matrix based on
 * given quaternion.
 */
void build_rotmatrix(double m[4][4], double q[4]);

/*
 * This function computes a quaternion based on an axis (defined by
 * the given vector) and an angle about which to rotate.  The angle is
 * expressed in radians.  The result is put into the third argument.
 */
void axis_to_quat(double a[3], double phi, double q[4]);


/*


static double tb_project_to_sphere(double, double, double);
static void normalize_quat(double [4]);

void vzero(double *v);

void vset(double *v, double x, double y, double z);

void vsub(const double *src1, const double *src2, double *dst);
void vcopy(const double *v1, double *v2);
void vcross(const double *v1, const double *v2, double *cross);

double vlength(const double *v);
void vscale(double *v, double div);
void vnormal(double *v);

double vdot(const double *v1, const double *v2);
void vadd(const double *src1, const double *src2, double *dst);
void trackball(double q[4], double p1x, double p1y, double p2x, double p2y);
void axis_to_quat(double a[3], double phi, double q[4]);
static double tb_project_to_sphere(double r, double x, double y);

void add_quats(double q1[4], double q2[4], double dest[4]);
static void normalize_quat(double q[4]);
void build_rotmatrix(double m[4][4], double q[4]);

*/
#endif /* __F_TRACKBALL__ */

