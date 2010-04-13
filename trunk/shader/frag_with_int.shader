/// GLSL CODE

/// 2nd Fragment Shader * WITH INTEGRATION *

/// Role of Fragment Shader on the 2nd step:

uniform sampler1D tfTex;
uniform sampler1D expTex;
uniform float preIntTexSize;
uniform sampler2D psiGammaTableTex;
uniform float brightness;

void main(void)
{
    //----- Input data ( -, sf, sb, thickness ) -----
    float sf = gl_Color.r;
    float sb = gl_Color.g;
    float l = gl_Color.b;
    //float l = gl_Color.a; //thickness

    l *= brightness;

    if (l == 0.0) // no fragment color
	discard;

    vec4 colorFront = texture1D(tfTex, sf).rgba;
    vec4 colorBack = texture1D(tfTex, sb).rgba;

    vec4 color;

    vec2 tau = vec2(colorFront.a, colorBack.a) * l;

    vec2 halfVec = vec2(0.5);

    float zeta = texture1D(expTex, dot(tau, halfVec)).a; // using tex1d for exponential

    if (zeta == 1.0) // no fragment color
	discard;

    //zeta = exp(-1.0 * (dot(tau, halfVec))); // computing exponential

    vec2 gamma = tau / (1.0 + tau);

    float psi = texture2D(psiGammaTableTex, gamma + (halfVec / vec2(preIntTexSize))).a;

    //    color.rgb = colorFront.rgb*(1.0 - psi) + colorBack.rgb*(psi - zeta);
    color.rgb = colorFront.rgb*(1.0 - psi) + colorBack.rgb*(psi - zeta);
    color.a = 1.0 - zeta;

    //----- Output color -----
    gl_FragColor = color;
}
