uniform float uVertexScale;
uniform float uCurrentWidth;
uniform float uCurrentHeight;

attribute vec2 aPosition;
attribute vec3 aColor;
attribute vec2 aTexCoord0, aTexCoord1;

varying vec3 vColor;
varying vec2 vTexCoord0, vTexCoord1;

void main() {
  

  if (uCurrentHeight >= uCurrentWidth) {
    gl_Position = vec4(aPosition.x , aPosition.y / uCurrentHeight * uCurrentWidth, 0,1);
  }
  else {
    gl_Position = vec4(aPosition.x  / uCurrentWidth * uCurrentHeight , aPosition.y, 0,1);
  }
 
  
  vColor = aColor;
  vTexCoord0 = aTexCoord0;
  vTexCoord1 = aTexCoord1;
}
