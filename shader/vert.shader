/// GLSL CODE

/// Vertex Shader

//uniform vec3 offset;
uniform vec3 scale;
uniform vec4 basisThickVerts[5];
uniform vec4 basisProjVerts[8];

void main(void)
{
  vec4 currPos = vec4(gl_Vertex.xyz, 1.0);
  //vec4 currPos = vec4(gl_Vertex.xyz - offset.xyz, 1.0);

  currPos = gl_ModelViewProjectionMatrix * currPos;

  currPos *= vec4(scale, 1.0);

  //pass id vars as normal
  int tetId = int(gl_Normal.x); 
  int vertId = int(gl_Normal.y);

  if (vertId == 8) //thick vertex
    currPos.xyz += basisThickVerts[tetId].xyz;
  else
    currPos.xyz += basisProjVerts[vertId].xyz;

  //gl_TexCoord[0] = gl_MultiTexCoord0;

  gl_Position = currPos;
  gl_FrontColor = gl_Color;
}
