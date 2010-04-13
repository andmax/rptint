/// GLSL CODE

//--- LINEAR COLOR AVERAGE OPACITY ---

const float e = 2.718281828;

float exp(float x)
{
  return pow(e, x);
}

float loge(float x)
{
  return log2(x) / log2(e);
}

vec4 integrateRay(in vec4 cb, in vec4 cf, in float l)
{
  float dtau = -l * loge (1.0 - 0.5 * (cb.a + cf.a));
  float zeta = exp(-dtau);
  float alpha = 1.0 - zeta;
  float Psi = alpha / dtau;

  vec4 c;
  c.rgb = cf.rgb*(1.0 - Psi) + cb.rgb*(Psi - zeta);
  c.a = alpha;
  return c;
}

//------------------------------------
