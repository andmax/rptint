/***********************************************************************************
 *                              Virtual Trackball
 *  Created by: Yalmar Ponce (yalmar@gmail.com|http://www.lcg.ufrj.br/Members/yalmar)
 *  Last Update: 06-05-2006
 *  --------------------------------------------------------------------------------
 *  This file is a part of work done by Yalmar. You is free to use this code in any 
 *  way you like. However, i expect you puts this note, how credits of the author
 *
 ***********************************************************************************/
 
#include "ftrackball.h"

#define DELTA 0.0000001

/*
 * Local function prototypes (not defined in trackball.h)
 */
static double tb_project_to_sphere(double, double, double);
static void normalize_quat(double [4]);
static void vzero(double *v);

static void vset(double *v, double x, double y, double z);
static void vsub(const double *src1, const double *src2, double *dst);
static void vcopy(const double *v1, double *v2);
static void vcross(const double *v1, const double *v2, double *cross);
static double vlength(const double *v);
static void vscale(double *v, double div);
static void vnormal(double *v);
static double vdot(const double *v1, const double *v2);
static void vadd(const double *src1, const double *src2, double *dst);


/***********************************************************************************/
/*! Default Constructor*/
FTrackBall::FTrackBall(float zoom) { reset(zoom); }

/***********************************************************************************/

/*! set the zoom to the trackball*/
void FTrackBall::setZoom(float zoom) { z = zoom; }

/*! Set the old coordinate of the mouse position */
void FTrackBall::setOld(int x, int y) { poldx = oldx; oldx = x; poldy = oldy; oldy = y; }

/*! This function mult the matrix to the current opengl matrix*/
void FTrackBall::multToOpenGL()
{  
  //gluLookAt(center[0], center[1], center[2]+z, center[0], center[1], center[2], 0.0, 1.0, 0.0);
  build_rotmatrix(rot_matrix, quat);
  glMultMatrixd(&rot_matrix[0][0]);
}

void FTrackBall::getRotMatrix(double **m){
  build_rotmatrix(rot_matrix, quat);
  *m = &rot_matrix[0][0];
}

void FTrackBall::prepareToIdle(int x, int y, int width, int height){
  double vx = (double)(oldx - poldx);
  double vy = (double)(oldy - poldy);
  double len = sqrt(vx*vx + vy*vy);
  
  // it guarantees to not have a malfunction
  if (len < DELTA) len = DELTA;
    
  vx /= len;
  vy /= len;  
    
  trackball(spin_quat,
	    (2.0 * ((float)oldx - vx)  - (float)width) / (float)width,
	    ((float)height - 2.0 * ((float)oldy - vy)) / (float)height,
	    (2.0 * x      - (float)width) / (float)width,
	    ((float)height - 2.0 * y)     / (float)height);  
}

/*! Calcul the Track-ball*/
void FTrackBall::trackBallInAction(int x, int y, int width, int height)
{   
  trackball(spin_quat,
	    (2.0 * oldx  - (float)width) / (float)width,
	    ((float)height - 2.0 * oldy) / (float)height,
	    (2.0 * x      - (float)width) / (float)width,
	    ((float)height - 2.0 * y)     / (float)height);

  add_quats(spin_quat, quat, quat);
}

void FTrackBall::idle(){
  add_quats(spin_quat, quat, quat);
}

/*! When you zoom */
void FTrackBall::zoomInAction(int y) 
{ 
  z*= (float)exp(double((oldy - y) / 300.0));
}

/*! Reset the trackball*/
void FTrackBall::reset(float zoom)
{
  z = zoom;
  center[0] = 0.0f;
  center[1] = 0.5f;
  center[2] = 0.0f;
  quat[0] = quat[1] = quat[2] = 0;
  quat[3] = 1;
}

/*! Return the older position of the pointer (x coord) */
int FTrackBall::getOldX() { return oldx; }

/*! Return the older position of the pointer (y coord) */
int FTrackBall::getOldY() { return oldy; }

  //sets the quaternion
void FTrackBall::setQuat(double q[4])
{
  normalize_quat(q);
  quat[0] = q[0];
  quat[1] = q[1];
  quat[2] = q[2];
  quat[3] = q[3];
}

/*
 * Ok, simulate a track-ball.  Project the points onto the virtual
 * trackball, then figure out the axis of rotation, which is the cross
 * product of P1 P2 and O P1 (O is the center of the ball, 0,0,0)
 * Note:  This is a deformed trackball-- is a trackball in the center,
 * but is deformed into a hyperbolic sheet of rotation away from the
 * center.  This particular function was chosen after trying out
 * several variations.
 *
 * It is assumed that the arguments to this routine are in the range
 * (-1.0 ... 1.0)
 */
void trackball(double q[4], double p1x, double p1y, double p2x, double p2y)
{
  double a[3]; /* Axis of rotation */
  double phi;  /* how much to rotate about axis */
  double p1[3], p2[3], d[3];
  double t;

  if ( ((p1x < (p2x + DELTA)) && (p1x > (p2x - DELTA)) )
       && ((p1x < (p2x + DELTA)) && (p1x > (p2x - DELTA)) ) )
    {
      /* Zero rotation */
      vzero(q);
      q[3] = 1.0;
      return;
    }

  /*
   * First, figure out z-coordinates for projection of P1 and P2 to
   * deformed sphere
   */
  vset(p1,p1x,p1y,tb_project_to_sphere(TRACKBALLSIZE,p1x,p1y));
  vset(p2,p2x,p2y,tb_project_to_sphere(TRACKBALLSIZE,p2x,p2y));

  /*
   *  Now, we want the cross product of P1 and P2
   */
  vcross(p2,p1,a);

  /*
   *  Figure out how much to rotate around that axis.
   */
  vsub(p1,p2,d);
  t = vlength(d) / (2.0*TRACKBALLSIZE);

  /*
   * Avoid problems with out-of-control values...
   */
  if (t > 1.0) t = 1.0;
  if (t < -1.0) t = -1.0;
  phi = 2.0 * asin(t);

  axis_to_quat(a,phi,q);
}

/*
 *  Given an axis and angle, compute quaternion.
 */
void axis_to_quat(double a[3], double phi, double q[4])
{
  vnormal(a);
  vcopy(a,q);
  vscale(q,sin(phi/2.0));
  q[3] = cos(phi/2.0);
}

/*
 * Project an x,y pair onto a sphere of radius r OR a hyperbolic sheet
 * if we are away from the center of the sphere.
 */
static double tb_project_to_sphere(double r, double x, double y)
{
  double d, t, z;

  d = sqrt(x*x + y*y);
  if (d < r * 0.70710678118654752440) {    /* Inside sphere */
    z = sqrt(r*r - d*d);
  } else {           /* On hyperbola */
    t = r / 1.41421356237309504880;
    z = t*t / d;
  }
  return z;
}

/*
 * Given two rotations, e1 and e2, expressed as quaternion rotations,
 * figure out the equivalent single rotation and stuff it into dest.
 *
 * This routine also normalizes the result every RENORMCOUNT times it is
 * called, to keep error from creeping in.
 *
 * NOTE: This routine is written so that q1 or q2 may be the same
 * as dest (or each other).
 */

void add_quats(double q1[4], double q2[4], double dest[4])
{
  static int count=0;
  double t1[4], t2[4], t3[4];
  double tf[4];

  vcopy(q1,t1);
  vscale(t1,q2[3]);

  vcopy(q2,t2);
  vscale(t2,q1[3]);

  vcross(q2,q1,t3);
  vadd(t1,t2,tf);
  vadd(t3,tf,tf);
  tf[3] = q1[3] * q2[3] - vdot(q1,q2);

  dest[0] = tf[0];
  dest[1] = tf[1];
  dest[2] = tf[2];
  dest[3] = tf[3];

  if (++count > RENORMCOUNT) {
    count = 0;
    normalize_quat(dest);
  }
}

void vzero(double *v)
{
  v[0] = 0.0;
  v[1] = 0.0;
  v[2] = 0.0;
}

void vset(double *v, double x, double y, double z)
{
  v[0] = x;
  v[1] = y;
  v[2] = z;
}

void vsub(const double *src1, const double *src2, double *dst)
{
  dst[0] = src1[0] - src2[0];
  dst[1] = src1[1] - src2[1];
  dst[2] = src1[2] - src2[2];
}

void vcopy(const double *v1, double *v2)
{
  register int i;
  for (i = 0 ; i < 3 ; i++)
    v2[i] = v1[i];
}

void vcross(const double *v1, const double *v2, double *cross)
{
  double temp[3];

  temp[0] = (v1[1] * v2[2]) - (v1[2] * v2[1]);
  temp[1] = (v1[2] * v2[0]) - (v1[0] * v2[2]);
  temp[2] = (v1[0] * v2[1]) - (v1[1] * v2[0]);
  vcopy(temp, cross);
}

double vlength(const double *v)
{
  return (double)sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

void vscale(double *v, double div)
{
  v[0] *= div;
  v[1] *= div;
  v[2] *= div;
}

void vnormal(double *v)
{
  vscale(v,1.0/vlength(v));
}

double vdot(const double *v1, const double *v2)
{
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

void vadd(const double *src1, const double *src2, double *dst)
{
  dst[0] = src1[0] + src2[0];
  dst[1] = src1[1] + src2[1];
  dst[2] = src1[2] + src2[2];
}

/*
 * Quaternions always obey:  a^2 + b^2 + c^2 + d^2 = 1.0
 * If they don't add up to 1.0, dividing by their magnitued will
 * renormalize them.
 *
 * Note: See the following for more information on quaternions:
 *
 * - Shoemake, K., Animating rotation with quaternion curves, Computer
 *   Graphics 19, No 3 (Proc. SIGGRAPH'85), 245-254, 1985.
 * - Pletinckx, D., Quaternion calculus as a basic tool in computer
 *   graphics, The Visual Computer 5, 2-13, 1989.
 */
static void normalize_quat(double q[4])
{
  int i;
  double mag;

  mag = (q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
  for (i = 0; i < 4; i++) q[i] /= mag;
}

/*
 * Build a rotation matrix, given a quaternion rotation. *
 */
void build_rotmatrix(double m[4][4], double q[4])
{
  m[0][0] = 1.0 - 2.0 * (q[1] * q[1] + q[2] * q[2]);
  m[0][1] = 2.0 * (q[0] * q[1] - q[2] * q[3]);
  m[0][2] = 2.0 * (q[2] * q[0] + q[1] * q[3]);
  m[0][3] = 0.0;

  m[1][0] = 2.0 * (q[0] * q[1] + q[2] * q[3]);
  m[1][1]= 1.0 - 2.0 * (q[2] * q[2] + q[0] * q[0]);
  m[1][2] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
  m[1][3] = 0.0;

  m[2][0] = 2.0 * (q[2] * q[0] - q[1] * q[3]);
  m[2][1] = 2.0 * (q[1] * q[2] + q[0] * q[3]);
  m[2][2] = 1.0 - 2.0 * (q[1] * q[1] + q[0] * q[0]);
  m[2][3] = 0.0;

  m[3][0] = 0.0;
  m[3][1] = 0.0;
  m[3][2] = 0.0;
  m[3][3] = 1.0;


  //  for (int i = 0; i < 4; ++i)
//   std::cout << std::endl << "q : " << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << std::endl;
//   for (int i = 0; i < 3; ++i)
//     {
//       for (int j = 0; j < 3; ++j)
// 	std::cout << m[i][j] << " ";
//       std:: cout << std::endl;
//     }
}
