///GLSL

//-------- Integration --------------

const float M_SQRTPI = 1.77245385090551602792981;
const float M_SQRT1_2 = 0.70710678118654752440;
const float M_2_SQRTPI = 1.12837916709551257390;
const float M_1_SQRTPI = (0.5*M_2_SQRTPI);

const float e = 2.718281828;

float exp(float x)
{
  return pow(e, x);
}

float erf_fitting_function(in float u)
{
  return
    - 1.26551223 + u*(1.00002368 + u*(0.37409196 + u*(0.09678418 + 
        u*(-0.18628806 + u*(0.27886807 + u*(-1.13520398 + u*(1.48851587 +
        u*(-0.82215223 + u*0.17087277))))))));
}
/*
float erf(in float x)
{
  // Compute as described in Numerical Recipes in C++ by Press, et al.
  //   x = abs(x);        In this application, x should always be <= 0.
  float u = 1.0 / (1.0 + 0.5*x);
  float ans = u*exp(-x*x + erf_fitting_function(u));
  //   return (x >= 0 ? 1 - ans : ans - 1);    x should always be <= 0.
  return 1.0 - ans;
}

float erfi(in float x)
{
  return M_2_SQRTPI*exp(x*x)*dawson(x);
}
*/

// Compute Dawson's integral as described in Numerical Recipes in C++ by
// Press, et al.
float dawson(in float x)
{
  const float H = 0.4;
  const float NMAX = 6.0;
  float dawson_constant0 = 0.852144;
  float dawson_constant1 = 0.236928;
  float dawson_constant2 = 0.0183156;
  float dawson_constant3 = 0.000393669;
  float dawson_constant4 = 2.35258e-6;
  float dawson_constant5 = 3.90894e-9;

  float result;
  if (x > 0.2)
    {
/*  x = abs(x);       In this application, x should always be <= 0. */

    int n0 = 2 * floor( (0.5 / H)*x + 0.5);
    float xp = x - float(n0) * H;
    float e1 = exp((2.0 * H) * xp);
    float e2 = e1 * e1;
    float d1 = float(n0) + 1.0;
    float d2 = d1 - 2.0;
    float sum = 0.0;
    sum = dawson_constant0*(e1/d1 + 1.0/(d2*e1));
    d1 += 2.0;  d2 -= 2.0;  e1 *= e2;
    sum += dawson_constant1*(e1/d1 + 1.0/(d2*e1));
    d1 += 2.0;  d2 -= 2.0;  e1 *= e2;
    sum += dawson_constant2*(e1/d1 + 1.0/(d2*e1));
    d1 += 2.0;  d2 -= 2.0;  e1 *= e2;
    sum += dawson_constant3*(e1/d1 + 1.0/(d2*e1));
    d1 += 2.0;  d2 -= 2.0;  e1 *= e2;
    sum += dawson_constant4*(e1/d1 + 1.0/(d2*e1));
    d1 += 2.0;  d2 -= 2.0;  e1 *= e2;
    sum += dawson_constant5*(e1/d1 + 1.0/(d2*e1));
    result = M_1_SQRTPI*exp(-xp*xp)*sum;
    }
  else
    {
    float x2 = x * x;
    result = x * (1.0 - (2.0/3.0) * x2 * (1.0 - 0.4 * x2 * (1.0 - (2.0/7.0) * x2)));
    }
  return result;
}

float Psi(in float taub, in float tauf, in float l)
{
  float difftau = taub - tauf;
  bool useHomoTau = ((difftau > -0.0001) && (difftau < 0.0001));
  bool useErf = (difftau > 0.0);

  float Y;

  if (!useHomoTau)
    {
    float invsqrt2lengthdifftau = 1.0 / sqrt(2.0 * l * abs(difftau));
    float t = l*invsqrt2lengthdifftau;
    float frontterm = t*tauf;
    float backterm = t*taub;
    float expterm = exp(frontterm*frontterm-backterm*backterm);
    if (useErf)
      {
      /* Back more opaque. */
      float u = 1.0 / (1.0 + 0.5*frontterm);
      Y = u*exp(erf_fitting_function(u));
      u = 1.0 / (1.0 + 0.5*backterm);
      Y += -expterm*u*exp(erf_fitting_function(u));
      Y *= M_SQRTPI;
      }
    else
      {
      /* Front more opaque. */
      expterm = 1.0/expterm;
      Y = 2.0*(dawson(frontterm) - expterm*dawson(backterm));
      }
    Y *= invsqrt2lengthdifftau;
    }
  else
    {
      float tauD = taub*l;
      Y = (1.0 - exp(-tauD))/tauD;
    }

  return Y;
}

//-----------------------------------
