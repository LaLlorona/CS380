#version 130

uniform mat4 uProjMatrix;
uniform mat4 uModelViewMatrix;
uniform mat4 uNormalMatrix;

//we are taking all the 3d information for the object 

in vec3 aPosition; 
in vec3 aNormal;

out vec3 vNormal;
out vec3 vPosition;

void main() {
  vNormal = vec3(uNormalMatrix * vec4(aNormal, 0.0)); 

  // send position (eye coordinates) to fragment shader
  //3d object of the eyeframe
  vec4 tPosition = uModelViewMatrix * vec4(aPosition, 1.0);
  vPosition = vec3(tPosition);
  //project the object to the screen 
  gl_Position = uProjMatrix * tPosition;
}