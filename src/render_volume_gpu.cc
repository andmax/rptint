/**
 *
 *    Render Volume GPU
 *
 *  File: render_volume_gpu.cc
 *
 *  Authors:
 *    Andre Maximo
 *    Ricardo Marroquim
 *
 *  Last Update: May 10, 2006
 *
 */

#include "volume.h"
   
#include <string>

#include <sstream>
#include <iostream>

#include <iomanip>
#include <math.h>

#include "ftrackball.h"

#define MIN_RENDER_TIME 1

using namespace std;

/// Global Variables

/// the volume
static volume* vol;

#include "main_utils.h" // temporally here...

// names

#ifndef NO_NVIDIA

const char titleModelWin[30] = "GPU Volume Rendering";
const char runingIn[4] = "GPU";

#else

const char titleModelWin[30] = "CPU Volume Rendering";
const char runingIn[4] = "CPU";

#endif

const char titleTFWin[30] = "Transfer Function";

static char *volume_name, *tfName, *info_file_name;
static fileType fType;
// windows
static int modelWindow, tfWindow;

// fps computation
static int time_frame = 0;

// rotation variables
static bool button_pressedModelWin[3] = {false, false, false};
static bool button_pressedTFWin = false;

// zoom box variables
GLint box_corners[2][2];

// trackball
static FTrackBall track;

// screen variable
GLint modelWinWidth=512, modelWinHeight=512,
  tfWinWidth=400, tfWinHeight=300;

GLint focusFace = 1;

// flag controls
bool integrating = false,
  debug_cout = true,
  debug_shaders = true,
  debug_setup = false,
  sorting = true;

static bool rotate_always = false,
  show_debug = true,
  show_volume = true,
  show_bbox = true,
  tf_visible = false,
  write_in_file = false;

// debug timers
static double update_time=0.0,
  sorting_time=0.0,
  setup_time=0.0,
  draw_time=0.0,
  total_time=0.0;

static unsigned long int time_elapsed = 0;
static bool write_once = false;

//static double gpu_fps = 0.0;

static uint num_tets = 0, num_verts = 0;

static uint dimensionX, dimensionY, dimensionZ;

static bool clear_white = true,
            text_white = false;

// information file
static ofstream info_file;

/// Static Local Functions

static void init(void);
static void reshapeModelWin(int w, int h);
static void reshapeTFWin(int w, int h);
static void idle(void);
static void glWrite(GLdouble x, GLdouble y, char *str);
static bool displayDebugModel(void);
static void displayDebugTF(void);
static void displayModelWin(void);
static void displayTFWin(void);
static void mouseModelWin(GLint button, GLint state, GLint x, GLint y);
static void motionModelWin(GLint x, GLint y);
static void keyboardModelWin(unsigned char key, int x, int y);
static void mouseTFWin(GLint button, GLint state, GLint x, GLint y);
static void motionTFWin(GLint x, GLint y);
static void keyboardTFWin(unsigned char key, int x, int y);
static GLuint processHits (GLuint, GLuint[]);
static void pickFace(int, int);
static GLint pickGrid(int, int);

/// Initialization (global)

void init(void) 
{
  vol->CreateTextures();
  vol->CreateShaders();
  vol->createArrays();

  vol->resetSelectionBox();
    
  num_verts = vol->getNumVerts();
  num_tets = vol->getNumTets();

  track.setZoom(1.0);

  if (write_in_file) {
    ostringstream oss_info;
    oss_info << "last_info_file.txt";

    info_file.open((char*)oss_info.str().c_str());
    info_file << "* " << runingIn << " Volume Rendering - Information file" << endl;
    info_file << "* Dimensions: " << dimensionX << " x " << dimensionY << " x " << dimensionZ << endl;
    info_file << "* Resolution: " << modelWinWidth << " x " << modelWinHeight << endl;
    info_file << "* Num Tets: " << num_tets << " ; Num Vertices: " << num_verts << endl;
    info_file << "* Psi-Gama Table: " << vol->getPsiGamaTableSize() << endl;
    info_file << "* Exponential Table: " << vol->getExpTexSize() << endl;
    info_file << "* Table (average times in seconds):" << endl;
    info_file << "---------------------------------------------------------------------------" << endl;
    if (debug_setup)
      info_file << "Time |   FPS   |  K Tets  |  Update   | Sorting   |     -     |   Setup   |" << endl;
    else
      info_file << "Time |   FPS   |  K Tets  |  Update   | Sorting   |     -     |   Draw    |" << endl;
    info_file << "---------------------------------------------------------------------------" << endl;
  }
}

/// Reshape (model)

void reshapeModelWin(int w, int h)
{
  glutSetWindow(modelWindow);
  glViewport(0, 0, modelWinWidth=w, modelWinHeight=h);
}

/// Reshape (tf)

void reshapeTFWin(int w, int h)
{
  glutSetWindow(tfWindow);
  glViewport(0, 0, tfWinWidth=w, tfWinHeight=h);
}

/// Idle (global)

void idle(void)
{
  if (rotate_always) {
    track.trackBallInAction(track.getOldX()+4, track.getOldY()+8,
			    modelWinWidth, modelWinHeight);
    glutSetWindow(modelWindow);
    glutPostRedisplay();
    return; // don't redisplay again
  }

  if (show_debug) {
    glutSetWindow(modelWindow);
    glutPostRedisplay(); // to update the debug information
  }
}

/// OpenGL Write (global)

void glWrite(GLdouble x, GLdouble y, char *str) {
  glRasterPos2d(x, y);
  for (char *s = str; *s; ++s)
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *s);
}

/// Display Debug Text (tf)

bool displayDebugModel(void)
{
  ostringstream oss_update, oss_sort, oss_setup, oss_draw;
  ostringstream oss_fps, oss_tps, oss_brightness;
  ostringstream oss_verts, oss_tets, oss_res;

  GLfloat actual_brightness = vol->tf.getBrightness();

  num_tets = vol->getNumTets();

  // trunc 3 decimals (and 2 for percentuals and rates)
  double update = (uint)(update_time*1000.0)/1000.0;
  double update_perc = (uint)((update_time/total_time)*10000.0)/100.0;
  double sort = (uint)(sorting_time*1000.0)/1000.0;
  double sort_perc = (uint)((sorting_time/total_time)*10000.0)/100.0;
  //  double setup = (uint)(setup_time*1000.0)/1000.0;
  //  double setup_perc = (uint)((setup_time/total_time)*10000.0)/100.0;
  double draw = (uint)(draw_time*1000.0)/1000.0;
  double draw_perc = (uint)((draw_time/total_time)*10000.0)/100.0;
  double fps = (uint)((1.0/total_time)*100.0)/100.0;
  //long double gpu_tets = gpu_fps*num_tets/1000.0; // K Tets / sec
  long double gpu_tets = ((long double)num_tets/total_time)/1000.0; // K Tets / sec
  long double tets = (unsigned long int)(gpu_tets*100.0)/100.0; // trunc 2 decimals



  oss_update << "CPU Update: " << update << " s (" << update_perc << "%)";

  oss_sort << "CPU Sort: " << sort << " s (" << sort_perc << " %)";

  //  oss_setup << "CPU Setup: " << setup << " s (" << setup_perc << " %)";

  if (debug_setup)
    oss_draw << "CPU Setup: " << draw << " s (" << draw_perc << " %)";
  else
    oss_draw << runingIn << " Draw: " << draw << " s (" << draw_perc << " %)";

  oss_fps << "FPS: " << fps;
  oss_tps << "Tet/s: " << tets << " K";
    
  oss_brightness << "Brightness: " << actual_brightness;

  oss_verts << "# Vertices: " << num_verts;	
  oss_tets << "# Tetrahedra: " << num_tets;
  oss_res << "Resolution: " << modelWinWidth << " x " << modelWinHeight;	

  if (text_white)
    glColor4f(1.0, 1.0, 1.0, 1.0);
  else
    glColor4f(0.0, 0.0, 0.0, 1.0);

  glWrite(-0.3, 1.0, volume_name);

  glWrite(-1.0, 0.8, (char*) oss_update.str().c_str());
  glWrite(-1.0, 0.7, (char*) oss_sort.str().c_str());
  //glWrite(-1.0, 0.7, (char*) oss_setup.str().c_str());
  glWrite(-1.0, 0.6, (char*) oss_draw.str().c_str());
  glWrite(-1.0, 0.5, (char*) oss_fps.str().c_str());
  glWrite(-1.0, 0.4, (char*) oss_tps.str().c_str());
  
  glWrite(-1.0, -0.7, (char*) oss_brightness.str().c_str());
  glWrite(-1.0, -0.8, (char*) oss_verts.str().c_str());
  glWrite(-1.0, -0.9, (char*) oss_tets.str().c_str());
  glWrite(-1.0, -1.0, (char*) oss_res.str().c_str());

  if (sorting)
    glWrite(0.7, -0.7, "* sorting *");
  if (rotate_always)
    glWrite(0.7, -0.8, "* rotating *");
  if (integrating)
    glWrite(0.7, -0.9, "* integrating *");

  if (write_in_file)
    {
      // total average values
      static long double total_fps = 0.0, total_tets = 0.0, total_update = 0.0,
	total_sorting = 0.0, total_setup = 0.0, total_draw = 0.0, num_sums = 0.0;
      // final average (fa) values
      static long double fa_total_fps = 0.0, fa_total_tets = 0.0, fa_total_update = 0.0,
	fa_total_sorting = 0.0, fa_total_setup = 0.0, fa_total_draw = 0.0, fa_num_sums = 0.0;

      total_fps += fps;
      total_tets += tets;
      total_update += update;
      total_sorting += sorting;
      //      total_setup += setup;
      total_draw += draw;
      num_sums += 1.0;

      if (write_once) { // write in file every 16 s
	
	info_file << setw(5) << time_elapsed << " "
		  << setw(9) << total_fps / num_sums << " "
		  << setw(9) << total_tets / num_sums << " "
		  << setw(11) << total_update / num_sums << " "
		  << setw(11) << total_sorting / num_sums << " "
		  << setw(11) << total_setup / num_sums << " "
		  << setw(11) << total_draw / num_sums << " "
		  << endl;

	fa_total_fps += total_fps;
	fa_total_tets += total_tets;
	fa_total_update += total_update;
	fa_total_sorting += total_sorting;
	fa_total_setup += total_setup;
	fa_total_draw += total_draw;
	fa_num_sums += num_sums;

	total_fps = 0.0; total_tets = 0.0; total_update = 0.0;
	total_sorting = 0.0; total_setup = 0.0; total_draw = 0.0;
	num_sums = 0;
	write_once = false;
      }

      if (time_elapsed >= 128) { // exit after 128 s
	// write the final average values
	
	info_file << "---------------------------------------------------------------------------" << endl;
	info_file << " avg "
		  << setw(9) << fa_total_fps / fa_num_sums << " "
		  << setw(9) << fa_total_tets / fa_num_sums << " "
		  << setw(11) << fa_total_update / fa_num_sums << " "
		  << setw(11) << fa_total_sorting / fa_num_sums << " "
		  << setw(11) << fa_total_setup / fa_num_sums << " "
		  << setw(11) << fa_total_draw / fa_num_sums << " "
		  << endl;
	info_file << "---------------------------------------------------------------------------" << endl;
	info_file << endl;

	glutDestroyWindow(modelWindow);
	glutDestroyWindow(tfWindow);
	return false;
      }
    }
  return true;
}

/// Display (model)
void displayModelWin(void)
{
  static int sta_sh=0, end_sh=0, now=0;
  static double total_dt=0;
  double update_dt=0, sorting_dt=0, setup_dt=0, draw_dt=0;

  glClear(GL_COLOR_BUFFER_BIT);   

  now = glutGet(GLUT_ELAPSED_TIME);
  total_dt = ( now - time_frame ) / 1000.0;
  if (total_dt >= MIN_RENDER_TIME) {
    time_frame = now;
    time_elapsed += MIN_RENDER_TIME;
    write_once = true;
  }

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  
  track.multToOpenGL();  

  double zoom = track.getZoom();
  glScaled(zoom, zoom, zoom);

  //--- Update ---
  //--- Update data ---
  if (show_debug)
    sta_sh = glutGet(GLUT_ELAPSED_TIME);
  
  vol->projectBasisHexahedra();
  
  if (show_debug) {
    end_sh = glutGet(GLUT_ELAPSED_TIME);
    update_dt = (end_sh - sta_sh) / 1000.0;	
  }
  
  //--- Draw ---

  if (show_debug)
    sta_sh = glutGet(GLUT_ELAPSED_TIME);

  glPushMatrix();

  if (show_volume)
    vol->draw();

  glPopMatrix();

  if (show_debug) {
    end_sh = glutGet(GLUT_ELAPSED_TIME);
    draw_dt = (end_sh - sta_sh) / 1000.0;
  }

  if (show_bbox)
    {
      glEnable(GL_BLEND);

      if (clear_white)
	glColor4d(0.0, 0.0, 0.0, 0.1);
      else
	glColor4d(1.0, 1.0, 1.0, 0.1);

      vol->drawBBox();

      glDisable(GL_BLEND);
    }

  //--- Text ---
  if (show_debug) {
    glPushMatrix();
    glLoadIdentity();
    if (total_dt >= MIN_RENDER_TIME)
    {
      update_time = update_dt;
      sorting_time = sorting_dt;
      setup_time = setup_dt;
      draw_time = draw_dt;
      total_time = update_time + sorting_time + setup_time + draw_time;
    }
    if (time_elapsed > 0) {
      if (!displayDebugModel())
	return;
    }
    glPopMatrix();
  }

  //draw selection box
  if (button_pressedModelWin[2] && (focusFace != -1))
    {

      GLfloat box[2][2];
      box[0][0] = (GLfloat)( ((box_corners[0][0]/(GLfloat)modelWinWidth) * (MAXORTHOSIZE - MINORTHOSIZE)) + MINORTHOSIZE );
      box[1][0] = (GLfloat)( ((box_corners[1][0]/(GLfloat)modelWinWidth) * (MAXORTHOSIZE - MINORTHOSIZE)) + MINORTHOSIZE );

      box[0][1] = (GLfloat)(((((GLfloat)modelWinHeight - box_corners[0][1])/(GLfloat)modelWinHeight) * (MAXORTHOSIZE - MINORTHOSIZE)) + MINORTHOSIZE );
      box[1][1] = (GLfloat)(((((GLfloat)modelWinHeight - box_corners[1][1])/(GLfloat)modelWinHeight) * (MAXORTHOSIZE - MINORTHOSIZE)) + MINORTHOSIZE );

      glMatrixMode(GL_PROJECTION);
      glPushMatrix();
      glLoadIdentity();
      glOrtho(MINORTHOSIZE, MAXORTHOSIZE, MINORTHOSIZE, MAXORTHOSIZE, MINORTHOSIZE, MAXORTHOSIZE);

      glMatrixMode(GL_MODELVIEW);
      glPushMatrix();
      glLoadIdentity();

      glEnable(GL_BLEND);
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      glColor4d(0.0, 0.0, 1.0, 0.2);

      glBegin(GL_QUADS);

      glVertex3f(box[0][0], box[0][1], 0.0); 
      glVertex3f(box[1][0], box[0][1], 0.0); 
      glVertex3f(box[1][0], box[1][1], 0.0);
      glVertex3f(box[0][0], box[1][1], 0.0);
 
      glEnd();
      glDisable(GL_BLEND);

      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      glLineWidth(2.0); 
      glEnable(GL_LINE_STIPPLE);
      glLineStipple(1, 3);
      glColor3d(0.0, 0.0, 1.0);

      glBegin(GL_QUADS);

      glVertex3f(box[0][0], box[0][1], 0.0); 
      glVertex3f(box[1][0], box[0][1], 0.0); 
      glVertex3f(box[1][0], box[1][1], 0.0);
      glVertex3f(box[0][0], box[1][1], 0.0);
 
      glEnd();

      glDisable(GL_LINE_STIPPLE);
      glLineWidth(1.0); 
      glPolygonMode(GL_FRONT, GL_FILL);

      glPopMatrix();
      glMatrixMode(GL_PROJECTION);
      glPopMatrix();

    }

  glutSwapBuffers();
}

/// Display Debug Text (tf)

void displayDebugTF(void)
{
  glColor3f(0.0, 0.0, 0.0);

  glWrite(-0.05, 1.12, "Alpha");
  glWrite(1.12, -0.01, "Scalar");
  glWrite(1.04, -0.16, "Brightness");

  int pPoint = vol->tf.getPickedPoint();

  if (pPoint != -1) {
    glColor3f(0.7, 0.3, 0.1);
    if (pPoint == vol->tf.getBrightnessId()) { // brightness
      GLfloat curBrightness = (int)(vol->tf.getBrightness()*100.0)/100.0;
      ostringstream oss_brightness;
      oss_brightness << curBrightness;
      glWrite(curBrightness/BRIGHTNESS_MAX - 0.03, -0.12,
	      (char*)oss_brightness.str().c_str());      
    }
    else { // control points
      GLfloat curAlpha = (int)(vol->tf.getCurAlpha()*100.0)/100.0;
      GLfloat curScalar = (int)(vol->tf.getCurScalar()*100.0)/100.0;
      ostringstream oss_alpha, oss_scalar;
      oss_alpha << curAlpha;
      oss_scalar << curScalar;
      glWrite(-0.12, curAlpha, (char*)oss_alpha.str().c_str());
      glWrite(curScalar/255.0 - 0.03, -0.05,
	      (char*)oss_scalar.str().c_str());
    }
  }
}

/// Display (tf)

void displayTFWin(void)
{
  glClearColor(1.0, 1.0, 1.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);
  
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  //--- Draw ---
  vol->tf.draw();

  //--- Text ---
  glColor3f(0.0, 0.0, 0.0);
  glWrite(0.25, 1.2, tfName);

  if (show_debug)
    displayDebugTF();

  glutSwapBuffers();
}

/// Mouse (model)

void mouseModelWin(GLint button, GLint state, GLint x, GLint y) {

  if ((button == GLUT_LEFT_BUTTON) && (state == GLUT_DOWN))
    {
      box_corners[0][0] = x;
      box_corners[0][1] = y;
      button_pressedModelWin[0] = true;     
    }
    
  if ((button == GLUT_MIDDLE_BUTTON) && (state == GLUT_DOWN))
    button_pressedModelWin[1] = true;

  if ((button == GLUT_LEFT_BUTTON) && (state == GLUT_UP))
    {
      if ((box_corners[0][0] == x) && (box_corners[0][1] == y))
	pickFace(x, y);
      button_pressedModelWin[0] = false;
      glutPostRedisplay();
    }
  if ((button == GLUT_MIDDLE_BUTTON) && (state == GLUT_UP))
    button_pressedModelWin[1] = false;

  if ((button == GLUT_RIGHT_BUTTON) && (state == GLUT_DOWN))
    {
      box_corners[0][0] = x;
      box_corners[0][1] = y;

      box_corners[1][0] = x;
      box_corners[1][1] = y;

      button_pressedModelWin[2] = true;
    }

  if ((button == GLUT_RIGHT_BUTTON) && (state == GLUT_UP))
    {
      if (focusFace != -1)
	{
	  box_corners[1][0] = x;
	  box_corners[1][1] = y;

	  GLfloat box[2][2];
	  box[0][0] = (GLfloat)( ((box_corners[0][0]/(GLfloat)modelWinWidth) * (MAXORTHOSIZE - MINORTHOSIZE)) + MINORTHOSIZE );
	  box[1][0] = (GLfloat)( ((box_corners[1][0]/(GLfloat)modelWinWidth) * (MAXORTHOSIZE - MINORTHOSIZE)) + MINORTHOSIZE );
	  
	  box[0][1] = (GLfloat)(((((GLfloat)modelWinHeight - box_corners[0][1])/(GLfloat)modelWinHeight) * (MAXORTHOSIZE - MINORTHOSIZE)) + MINORTHOSIZE );
	  box[1][1] = (GLfloat)(((((GLfloat)modelWinHeight - box_corners[1][1])/(GLfloat)modelWinHeight) * (MAXORTHOSIZE - MINORTHOSIZE)) + MINORTHOSIZE );

	  GLint stPoint = pickGrid(box_corners[0][0], box_corners[0][1]);
	  GLint endPoint = pickGrid(x, y);

// 	  cout << "points : " << box_corners[0][0] << " " << box_corners[0][1] << " " <<
// 	    x << " " << y << endl; 
// 	  cout << "pick grid : " << stPoint << " " << endPoint << endl;

	  //vol->setSelectionBox(box[0][0], box[0][1], box[1][0], box[1][1], focusFace);
	  
	  vol->setSelectionBox(stPoint, endPoint, focusFace);

	  num_tets = vol->getCurNumTets();
	}  
      button_pressedModelWin[2] = false;
    }
  track.setOld(x, y);
}

GLuint processHits (GLuint hits, GLuint buffer[])
{
  uint i;
  GLuint names, *ptr, minZ,*ptrNames;
  GLuint closest = 9;

  ptr = (GLuint *) buffer;
  minZ = 0xffffffff;
  for (i = 0; i < hits; i++) 
    {	
      names = *ptr;
      ptr++;
      if (*ptr < minZ) 
	{
	  closest = *(ptr+2);
	  minZ = *ptr;
	  ptrNames = ptr+2;
	} 
      ptr += names+2;
    }

  return closest;
}

/// Pick Face

void pickFace(int x, int y)
{
  GLuint selectBuf[512];
  GLint viewport[4];

  glSelectBuffer(512, selectBuf);
  glRenderMode(GL_SELECT);

  glMatrixMode(GL_PROJECTION);
  
  glPushMatrix();
  glLoadIdentity();


  glOrtho(MINORTHOSIZE, MAXORTHOSIZE, MINORTHOSIZE, MAXORTHOSIZE, MINORTHOSIZE, MAXORTHOSIZE);

  glGetIntegerv(GL_VIEWPORT, viewport);
  gluPickMatrix(x, viewport[3] - y, 3, 3, viewport);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  track.multToOpenGL();

  double zoom = track.getZoom();
  glScaled(zoom, zoom, zoom);
  glInitNames();

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  vol->drawBBox(true);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  uint hits;
	
  // restoring the original projection matrix
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glutSwapBuffers();
	
  // returning to normal rendering mode
  hits = glRenderMode(GL_RENDER);
	
  focusFace = -1;

  if (hits == 0)
    return; // if there aren't hits to be processed

  focusFace = processHits(hits, selectBuf);

  //  cout << "face : " << focusFace << endl;

  double axis[3] = {0.0, 0.0, 0.0};
  double q[4] = {0.0, 0.0, 0.0, 0.0};
  double angle = 0.0;
  switch (focusFace)
    {
    case 0:
      angle = M_PI;
      axis[1] = 1.0;
      break;
    case 1:
      angle = 0.0;
      axis[0] = 1.0;
      break;
    case 2:
      angle = -M_PI/2.0;
      axis[1] = 1.0;
      break;
    case 3:
      angle = M_PI/2.0;
      axis[1] = 1.0;
      break;
    case 4:
      angle = -M_PI/2.0;
      axis[0] = 1.0;
      break;
    case 5:
      angle = M_PI/2.0;
      axis[0] = 1.0;
      break;
    default:
      break;
    }
  
  axis_to_quat(axis, angle, q);
  
  track.reset(1.0);
  track.setQuat(q);

  
}

/// Pick Grid

GLint pickGrid(int x, int y)
{
  GLuint selectBuf[512];
  GLint viewport[4];

  glSelectBuffer(512, selectBuf);
  glRenderMode(GL_SELECT);

  glMatrixMode(GL_PROJECTION);
  
  glPushMatrix();
  glLoadIdentity();

  glOrtho(MINORTHOSIZE, MAXORTHOSIZE, MINORTHOSIZE, MAXORTHOSIZE, MINORTHOSIZE, MAXORTHOSIZE);
  
  glGetIntegerv(GL_VIEWPORT, viewport);
  gluPickMatrix(x, viewport[3] - y, 3, 3, viewport);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  double zoom = track.getZoom();
  glScaled(zoom, zoom, zoom);

  glInitNames();

  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  vol->drawGrid(focusFace);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  uint hits;
	
  // restoring the original projection matrix
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  glMatrixMode(GL_MODELVIEW);
  glutSwapBuffers();
	
  // returning to normal rendering mode
  hits = glRenderMode(GL_RENDER);

  if (hits == 0)
    return -1; // if there aren't hits to be processed

  return (GLint)processHits(hits, selectBuf);
}

/// Mouse (tf)

void mouseTFWin(GLint button, GLint state, GLint x, GLint y)
{
  if ((button == GLUT_LEFT_BUTTON) && (state == GLUT_DOWN))
    {
      vol->tf.pick(x, y);	  
      button_pressedTFWin = true;
      glutPostRedisplay();
    }

  if ((button == GLUT_LEFT_BUTTON) && (state == GLUT_UP))
    {
      vol->tf.updateTF(x, y);
      vol->tf.setPickedPointNull();
      vol->tf.computeTF();
      button_pressedTFWin = false;
      glutSetWindow(modelWindow);
      vol->createTransferFunctionTex();
      vol->updateHexaMask();
      glutPostRedisplay();
      glutSetWindow(tfWindow);
      glutPostRedisplay();
    }
}

/// Motion (model)

void motionModelWin(GLint x, GLint y)
{
  if (button_pressedModelWin[0]) { // rotate
    focusFace = -1;
    track.trackBallInAction(x, y, modelWinWidth, modelWinHeight);
    glutPostRedisplay();
  }
    
  if (button_pressedModelWin[1]) { // zoom
    track.zoomInAction(y);
    glutPostRedisplay();
  }

  if (button_pressedModelWin[2]) { //box selection

    box_corners[1][0] = x;
    box_corners[1][1] = y;
  }
  track.setOld(x, y);
}

/// Motion (tf)

void motionTFWin(GLint x, GLint y)
{
  if (button_pressedTFWin) 
    {
      vol->tf.updateTF(x, y);
      vol->tf.computeTF();
      glutSetWindow(modelWindow);
      vol->createTransferFunctionTex();
      glutPostRedisplay();
      glutSetWindow(tfWindow);
      glutPostRedisplay();
    }
}

/// keyboard (model)

void keyboardModelWin(unsigned char key, int x, int y)
{

  switch(key)
    {
    case '0':
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9':
      //char sk[1];
      //sk[0] = key;
      vol->tf.setColorCode(atoi((const char*)&key));
      glutSetWindow(tfWindow);
      glutPostRedisplay();
      glutSetWindow(modelWindow);
      vol->createTransferFunctionTex();
      glutPostRedisplay();
      break;

    case '+':
      vol->tf.updateBrightness(+0.2);
      glutSetWindow(tfWindow);
      glutPostRedisplay();
      glutSetWindow(modelWindow);
      glutPostRedisplay();
      break;
    case '-':
      vol->tf.updateBrightness(-0.2);
      glutSetWindow(tfWindow);
      glutPostRedisplay();
      glutSetWindow(modelWindow);
      glutPostRedisplay();
      break;
      //--- no ---
    case 'f':
    case 'F':
      glutFullScreen();
      break;
    case 'n':
    case 'N':
      sorting = !sorting;
      break;
      //-----------
    case 'o':
    case 'O':
      text_white = !text_white;
      clear_white = !clear_white;
      glutSetWindow(modelWindow);
      if (clear_white)
	glClearColor(1.0, 1.0, 1.0, 0.0);
      else
	glClearColor(0.0, 0.0, 0.0, 0.0);
      glutPostRedisplay();
      break;
    case 'v':
    case 'V':
      focusFace = 1;
      vol->resetSelectionBox();
      track.reset(1.0); 
      glutPostRedisplay();
      break;
    case 'r':
    case 'R':
      rotate_always = !rotate_always;
      glutPostRedisplay();
      break;
    case 'a':
    case 'A':
      show_volume = !show_volume;
      glutPostRedisplay();
      break;
    case 's':
    case 'S':
      show_debug = !show_debug;
      glutPostRedisplay();
      break;
    case 't':
    case 'T':
      tf_visible = !tf_visible;
      glutSetWindow(tfWindow);
      if (tf_visible) {
	glutShowWindow();
	glutPostRedisplay();
      }
      else
	glutHideWindow();
      glutSetWindow(modelWindow);
      break;
    case 'i':
    case 'I':
      integrating = !integrating;
      glutPostRedisplay();
      break;
    case 'b':
    case 'B':
      show_bbox = !show_bbox;
      glutPostRedisplay();
      break;
    case 'q':
    case 'Q':
    case 27:
      glutDestroyWindow(modelWindow);
      glutDestroyWindow(tfWindow);
      break;
    default:
      break;
    }
}

/// Keyboard (tf)

void keyboardTFWin(unsigned char key, int x, int y)
{
  switch(key)
    {
    case 'w':
    case 'W':
      vol->tf.writeTF(tfName);      
      break;
    case 'q':
    case 'Q':
    case 'h':
    case 'H':
    case 27:
      glutHideWindow();
      break;
    default:
      keyboardModelWin(key, x, y);
      break;
    }
}

/// Main

int main(int argc, char** argv)
{
  int sta_in=0, end_in=0;
  double init_time;

  vol = new volume;

  if (show_debug)
    sta_in = glutGet(GLUT_ELAPSED_TIME);

  // Extended options
  if ( (strcmp(argv[argc-1], "-t") == 0)) {
    write_in_file = true;
    rotate_always = true;
    show_debug = false;
    debug_setup = true;
    argc--;
  }

  if (argc == 1) { // for commodity...
    argv[1] = "fuel";
    argc = 3;
  }

  uint dim[3] = {0, 0, 0};

  fType = read_args(argc, argv, dim);

  volume_name = argv[1];
  tfName = argv[2];
  dimensionX = dim[0];
  dimensionY = dim[1];
  dimensionZ = dim[2];

  cout << "Volume name: " << volume_name << endl;

  vol->readMedBIN(volume_name, dim[0], dim[1], dim[2], fType);
  vol->tf.readTF(tfName);
  info_file_name = argv[1];

  // initialize opengl
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
  glutInitWindowSize(modelWinWidth, modelWinHeight);
  glutInitWindowPosition(50, 100);
  modelWindow = glutCreateWindow(titleModelWin);

  glDisable(GL_POINT_SMOOTH);
  glDisable(GL_LINE_SMOOTH);
  glDisable(GL_POLYGON_SMOOTH);

  glDisable(GL_CULL_FACE); 

  glDisable(GL_DEPTH_TEST);
  glDisable(GL_STENCIL_TEST);
  glDisable(GL_FOG);
  glDisable(GL_LIGHTING);
  glDisable(GL_LIGHT0);

  glDisable(GL_TEXTURE_1D);
  glDisable(GL_TEXTURE_2D);
  glDisable(GL_TEXTURE_3D);

  glPolygonMode(GL_FRONT, GL_FILL);
  if (clear_white)
    glClearColor(1.0, 1.0, 1.0, 0.0);
  else
    glClearColor(0.0, 0.0, 0.0, 0.0);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(MINORTHOSIZE, MAXORTHOSIZE, MINORTHOSIZE, MAXORTHOSIZE, MINORTHOSIZE, MAXORTHOSIZE);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  init();

  // user interaction
  glutMotionFunc(motionModelWin);
  glutMouseFunc(mouseModelWin);
  glutKeyboardFunc(keyboardModelWin);
  glutIdleFunc(idle);

    
  if (debug_cout) {
    end_in = glutGet(GLUT_ELAPSED_TIME);
    init_time = (end_in - sta_in) / 1000.0;
    cout << "*** Initialization time      : " << setw(8) << init_time << " s ***" << endl;
    cout << "*********************************************" << endl;
    cout << endl;
    time_frame = glutGet(GLUT_ELAPSED_TIME);
  }

  glutInitWindowSize(tfWinWidth, tfWinHeight);
  glutInitWindowPosition(600, 100);    
  glutReshapeFunc(reshapeModelWin);
  glutDisplayFunc(displayModelWin);

  tfWindow = glutCreateWindow(titleTFWin);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(TF_ORTHO_MIN, TF_ORTHO_MAX, TF_ORTHO_MIN, TF_ORTHO_MAX);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  if (!tf_visible)
    glutHideWindow();

  // user interaction
  glutMotionFunc(motionTFWin);
  glutMouseFunc(mouseTFWin);
  glutKeyboardFunc(keyboardTFWin);

  glutReshapeFunc(reshapeTFWin);
  glutDisplayFunc(displayTFWin);
    
  glutSetWindow(modelWindow);

  glutMainLoop();

  if (write_in_file)
    info_file.close();

  return 0;
}
