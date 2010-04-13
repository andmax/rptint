/// GLSL CODE

/// 2nd Fragment Shader * NO INTEGRATION *

/// Role of Fragment Shader on the 2nd step:

uniform sampler1D tfTex;
uniform sampler1D expTex;
uniform float brightness;

void main(void)
{
    //----- Input data ( -, sf, sb, thickness ) -----
    float sf = gl_Color.r;
    float sb = gl_Color.g;
    float l = gl_Color.b; //thickness
    
    l *= brightness;

    if (l == 0.0) // no fragment color
	discard;
    
    float s_avg = (sf + sb) * 0.5; // average scalar
   


    vec4 c_avg = texture1D(tfTex, s_avg).rgba; // average color

    if (c_avg.a == 0.0) // no fragment color
	discard;

    c_avg.a = 1.0 - texture1D(expTex, l*c_avg.a).a; // exp alpha
    c_avg.rgb *= c_avg.a;

    //----- Output color -----
    gl_FragColor = c_avg;
}
