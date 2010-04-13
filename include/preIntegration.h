/**
 *
 *    Render Volume GPU
 *
 *  File: preIntegration.h
 *
 *  Authors:
 *    Andre Maximo
 *    Ricardo Marroquim
 *
 *  Last Update: Apr 29, 2006
 *
 */

/**
 * Note:
 * This header is used by volume.cc only in 
 * to compute the pre-integration table.
 */

/// Generate pre-integration 3D texture
/// Reference HARC 2002

GLfloat* volume::harcPreIntegration(void)
{
  GLfloat preIntBuffer[preIntTexSize][preIntTexSize][preIntTexSize][4];

  for (uint sfi = 0; sfi < preIntTexSize; ++sfi)
    for (uint sbi = 0; sbi < preIntTexSize; ++sbi)
      for (uint li = 0; li < preIntTexSize; ++li)
	{	    
	  GLfloat integratedColor[3] = {0.0, 0.0, 0.0};
	  GLfloat integratedAlpha = 0.0;

	  GLfloat sf = sfi / ((GLfloat)preIntTexSize - 1);
	  GLfloat sb = sbi / ((GLfloat)preIntTexSize - 1);
	  GLfloat l = li / ((GLfloat)preIntTexSize - 1);

	  GLfloat dw = 1.0 / ((GLfloat)preIntStep);
	  for(GLfloat w = 0.0; w < 1.0; w += dw)
	    {
	      GLfloat s = (1.0 - w)*sf + w * sb;
	      GLfloat color[3];
	      GLfloat tau;
	      GLfloat intExp = 0.0;

	      GLfloat dw_ = w / (GLfloat)preIntStep;
	      for(GLfloat w_ = 0.0; w_ < w; w_ += dw_)
		{
		  GLfloat s_ = (1 - w_)*sf + w_ * sb;
		  GLfloat tau_ = tf->GetAlpha(s_);
		  tau_ *= dw_ * l;
		  intExp += tau_;		   
		}

	      intExp = exp(-intExp);


	      for (int k = 0; k < 3; ++k)
		color[k] = (tf->Get(s))->GetColor(k);
	      tau = tf->GetAlpha(s);

	      for (int k = 0; k < 3; ++k)
		integratedColor[k] += tau * color[k] * l * intExp * dw;
	      integratedAlpha +=  tau * l * dw;
	    }

	  integratedAlpha = 1.0 - exp(-integratedAlpha);
	  //cout << sfi << " " << sbi << " " << li << " : ";
	  for (int k = 0; k < 3; ++k)
	    preIntBuffer[sfi][sbi][li][k] = integratedColor[k];
	  preIntBuffer[sfi][sbi][li][3] = integratedAlpha;
	  //cout << preIntBuffer[sfi][sbi][li][k] << " ";
	  //cout << endl;
	}


  cout << "pre integration buffer computed" << endl;

  return &preIntBuffer[0][0][0][0];
}

/// Generate pre-integration 3D texture
/// Reference A Fast High Accuracy Volume Renderer for Unstructured Data 2004

GLfloat* volume::morelandPreIntegration(void)
{
  GLfloat preIntBuffer[preIntTexSize][preIntTexSize][preIntTexSize][4];

  return &preIntBuffer[0][0][0][0];
}

/// Generate pre-integration 3D texture
/// Reference Harwared-Accelerated Volume and Isosurface Rendering based on Cell-Projection 2000

GLfloat* volume::rottgerPreIntegration(void)
{
  GLfloat preIntBuffer[preIntTexSize][preIntTexSize][preIntTexSize][4];

  for (uint sfi = 0; sfi < preIntTexSize; ++sfi)
    for (uint sbi = 0; sbi < preIntTexSize; ++sbi)
      for (uint li = 0; li < preIntTexSize; ++li)
	{	    
	  GLfloat integratedColor[3] = {0.0, 0.0, 0.0};
	  GLfloat integratedAlpha = 0.0;

	  GLfloat sf = sfi / ((GLfloat)preIntTexSize - 1);
	  GLfloat sb = sbi / ((GLfloat)preIntTexSize - 1);
	  GLfloat l = li / ((GLfloat)preIntTexSize - 1);

	  GLfloat dt = l / ((GLfloat)preIntStep);
	  for(GLfloat t = 0.0; t < l; t += dt)
	    {
	      GLfloat sl_t = sf + (t/l)*(sb - sf);
	      GLfloat color[3];
	      GLfloat rho_t;
	      GLfloat intExp = 0.0;

	      GLfloat du = t / (GLfloat)preIntStep;
	      for(GLfloat u = 0.0; u < t; u += du)
		{
		  GLfloat sl_u = sf + (u/l)*(sb - sf);
		  GLfloat rho_u = tf->GetAlpha(sl_u);
		  rho_u *= du;
		  intExp += rho_u; 
		}

	      intExp = exp(-intExp);


	      for (int k = 0; k < 3; ++k)
		color[k] = (tf->Get(sl_t))->GetColor(k);
	      rho_t = tf->GetAlpha(sl_t);

	      for (int k = 0; k < 3; ++k)
		integratedColor[k] += intExp * color[k] * rho_t * dt;
	      integratedAlpha +=  rho_t * dt;
	    }

	  integratedAlpha = 1.0 - exp(-integratedAlpha);
	  //cout << sfi << " " << sbi << " " << li << " : ";
	  for (int k = 0; k < 3; ++k)
	    preIntBuffer[sfi][sbi][li][k] = integratedColor[k];
	  preIntBuffer[sfi][sbi][li][3] = integratedAlpha;
	  //cout << preIntBuffer[sfi][sbi][li][k] << " ";
	  //cout << endl;
	}


  cout << "pre integration buffer computed" << endl;

  return &preIntBuffer[0][0][0][0];
}
