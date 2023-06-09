////////////////////////////////////////////////////////////////////////
//
//   KAIST, Spring 2023
//   CS380: Introduction to Computer Graphics
//   Instructor: Minhyuk Sung (mhsung@kaist.ac.kr)
//   Last Update: Juil Koo (63days@kaist.ac.kr)
//
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
// If your OS is LINUX, uncomment the line below.
//#include <tr1/memory>

#include <GL/glew.h>
#ifdef __MAC__
#   include <GLUT/glut.h>
#else
#   include <GL/glut.h>
#endif

#include "cvec.h"
#include "matrix4.h"
#include "rigtform.h"
#include "geometrymaker.h"
#include "ppm.h"
#include "glsupport.h"
#include "arcball.h"

using namespace std;      // for string, vector, iostream, and other standard C++ stuff
// If your OS is LINUX, uncomment the line below.
//using namespace tr1; // for shared_ptr

// G L O B A L S ///////////////////////////////////////////////////

// --------- IMPORTANT --------------------------------------------------------
// Before you start working on this assignment, set the following variable
// properly to indicate whether you want to use OpenGL 2.x with GLSL 1.0 or
// OpenGL 3.x+ with GLSL 1.3.
//
// Set g_Gl2Compatible = true to use GLSL 1.0 and g_Gl2Compatible = false to
// use GLSL 1.3. Make sure that your machine supports the version of GLSL you
// are using. In particular, on Mac OS X currently there is no way of using
// OpenGL 3.x with GLSL 1.3 when GLUT is used.
//
// If g_Gl2Compatible=true, shaders with -gl2 suffix will be loaded.
// If g_Gl2Compatible=false, shaders with -gl3 suffix will be loaded.
// To complete the assignment you only need to edit the shader files that get
// loaded
// ----------------------------------------------------------------------------
static const bool g_Gl2Compatible = true;


static const float g_frustMinFov = 60.0;  // A minimal of 60 degree field of view
static float g_frustFovY = g_frustMinFov; // FOV in y direction (updated by updateFrustFovY)

static const float g_frustNear = -0.1;    // near plane
static const float g_frustFar = -50.0;    // far plane
static const float g_groundY = -2.0;      // y coordinate of the ground
static const float g_groundSize = 10.0;   // half the ground length

static int g_windowWidth = 512;
static int g_windowHeight = 512;
static bool g_mouseClickDown = false;    // is the mouse button pressed
static bool g_mouseLClickButton, g_mouseRClickButton, g_mouseMClickButton;
static int g_mouseClickX, g_mouseClickY; // coordinates for mouse click event
static int g_activeShader = 0;

static double g_arcballScreenRadius;
static double g_currentArcballRadius;

struct ShaderState {
  GlProgram program;

  // Handles to uniform variables
  GLint h_uLight, h_uLight2;
  GLint h_uProjMatrix;
  GLint h_uModelViewMatrix;
  GLint h_uNormalMatrix;
  GLint h_uColor;

  // Handles to vertex attributes
  GLint h_aPosition;
  GLint h_aNormal;

  ShaderState(const char* vsfn, const char* fsfn) {
    readAndCompileShader(program, vsfn, fsfn);

    const GLuint h = program; // short hand

    // Retrieve handles to uniform variables
    h_uLight = safe_glGetUniformLocation(h, "uLight");
    h_uLight2 = safe_glGetUniformLocation(h, "uLight2");
    h_uProjMatrix = safe_glGetUniformLocation(h, "uProjMatrix");
    h_uModelViewMatrix = safe_glGetUniformLocation(h, "uModelViewMatrix");
    h_uNormalMatrix = safe_glGetUniformLocation(h, "uNormalMatrix");
    h_uColor = safe_glGetUniformLocation(h, "uColor");

    // Retrieve handles to vertex attributes
    h_aPosition = safe_glGetAttribLocation(h, "aPosition");
    h_aNormal = safe_glGetAttribLocation(h, "aNormal");

    if (!g_Gl2Compatible)
      glBindFragDataLocation(h, 0, "fragColor");
    checkGlErrors();
  }

};

static const int g_numShaders = 2;
static const char * const g_shaderFiles[g_numShaders][2] = {
  {"./shaders/basic-gl3.vshader", "./shaders/diffuse-gl3.fshader"},
  {"./shaders/basic-gl3.vshader", "./shaders/solid-gl3.fshader"}
};
static const char * const g_shaderFilesGl2[g_numShaders][2] = {
  {"./shaders/basic-gl2.vshader", "./shaders/diffuse-gl2.fshader"},
  {"./shaders/basic-gl2.vshader", "./shaders/solid-gl2.fshader"}
};
static vector<shared_ptr<ShaderState> > g_shaderStates; // our global shader states

// --------- Geometry

// Macro used to obtain relative offset of a field within a struct
#define FIELD_OFFSET(StructType, field) &(((StructType *)0)->field)

// A vertex with floating point position and normal
struct VertexPN {
  Cvec3f p, n;

  VertexPN() {}
  VertexPN(float x, float y, float z,
           float nx, float ny, float nz)
    : p(x,y,z), n(nx, ny, nz)
  {}

  // Define copy constructor and assignment operator from GenericVertex so we can
  // use make* functions from geometrymaker.h
  VertexPN(const GenericVertex& v) {
    *this = v;
  }

  VertexPN& operator = (const GenericVertex& v) {
    p = v.pos;
    n = v.normal;
    return *this;
  }
};

struct Geometry {
  GlBufferObject vbo, ibo;
  int vboLen, iboLen;

  Geometry(VertexPN *vtx, unsigned short *idx, int vboLen, int iboLen) {
    this->vboLen = vboLen;
    this->iboLen = iboLen;

    // Now create the VBO and IBO
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(VertexPN) * vboLen, vtx, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned short) * iboLen, idx, GL_STATIC_DRAW);
  }

  void draw(const ShaderState& curSS) {
    // Enable the attributes used by our shader
    safe_glEnableVertexAttribArray(curSS.h_aPosition);
    safe_glEnableVertexAttribArray(curSS.h_aNormal);

    // bind vbo
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    safe_glVertexAttribPointer(curSS.h_aPosition, 3, GL_FLOAT, GL_FALSE, sizeof(VertexPN), FIELD_OFFSET(VertexPN, p));
    safe_glVertexAttribPointer(curSS.h_aNormal, 3, GL_FLOAT, GL_FALSE, sizeof(VertexPN), FIELD_OFFSET(VertexPN, n));

    // bind ibo
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);

    //before drawing, we need to bind vertex buffer object + index buffer object 

    // draw!
    glDrawElements(GL_TRIANGLES, iboLen, GL_UNSIGNED_SHORT, 0);

    // Disable the attributes used by our shader
    safe_glDisableVertexAttribArray(curSS.h_aPosition);
    safe_glDisableVertexAttribArray(curSS.h_aNormal);
  }
};


// Vertex buffer and index buffer associated with the ground and cube geometry
static shared_ptr<Geometry> g_ground, g_cube, g_arcball;

// --------- Scene

static const Cvec3 g_light1(2.0, 3.0, 14.0), g_light2(-2, -3.0, -5.0);  // define two lights positions in world space
static RigTForm g_skyRbt = RigTForm().setTranslation(Cvec3(0.0, 0.25, 4.0));
// g_skyRbt.setTranslation(Cvec3(0.0, 0.25, 4.0));
// RigTForm.setTranslation(Cvec3(0.0, 0.25, 4.0));
static RigTForm g_objectRbt[1] = {RigTForm().setTranslation(Cvec3(-1,0,0))};  // currently only 1 obj is defined
static Cvec3f g_objectColors[1] = {Cvec3f(1, 0, 0)};

//this part is implemented by kkm. fuck!
static RigTForm g_object2Rbt[1] = {RigTForm().setTranslation(Cvec3(1,0,0))};  // defining other object
static Cvec3f g_object2Colors[1] = {Cvec3f(0, 1, 0)};

static RigTForm g_arcballRbt = RigTForm();
static Cvec3f g_arcballColor = Cvec3f(1, 1, 1);

static RigTForm g_myRbt = RigTForm().setTranslation(Cvec3(0.0, 0.25, 4.0));

static int cameraStatus = 0;
static int manipulationStatus = 0;
static int worldSkyFrameStatus = 0; // 0: sky sky, 1: world sky

enum CAMERA_STATUS {
  CAMERA_SKY,
  CAMERA_CUBE1,
  CAMERA_CUBE2,
};

enum MANIPULATION_STATUS {
  MANIPULATION_CUBE1,
  MANIPULATION_CUBE2,
  MANIPULATION_SKY,
};

enum WORLDSKYFRAME_STATUS {
  SKYSKY,
  WORLDSKY
};


///////////////// END OF G L O B A L S //////////////////////////////////////////////////




static void initGround() {
  // A x-z plane at y = g_groundY of dimension [-g_groundSize, g_groundSize]^2
  VertexPN vtx[4] = {
    VertexPN(-g_groundSize, g_groundY, -g_groundSize, 0, 1, 0), //first three are position information, next three are normal information
    VertexPN(-g_groundSize, g_groundY,  g_groundSize, 0, 1, 0),
    VertexPN( g_groundSize, g_groundY,  g_groundSize, 0, 1, 0),
    VertexPN( g_groundSize, g_groundY, -g_groundSize, 0, 1, 0),
  };
  unsigned short idx[] = {0, 1, 2, 0, 2, 3}; //index for drawing a triangle. 
  g_ground.reset(new Geometry(&vtx[0], &idx[0], 4, 6)); //vertex information + index information
}

static void initCubes() {
  int ibLen, vbLen;
  getCubeVbIbLen(vbLen, ibLen);

  // Temporary storage for cube geometry
  vector<VertexPN> vtx(vbLen);
  vector<unsigned short> idx(ibLen);

  makeCube(1, vtx.begin(), idx.begin());
  g_cube.reset(new Geometry(&vtx[0], &idx[0], vbLen, ibLen));
}

static void initArcball() {
  int ibLen, vbLen;
  int slices = 10;
  int stacks = 10;

  //radius 
  // this part should be changed;
  // float radius = 0.25 * min(g_windowHeight, g_windowWidth) * 0.001; 
  // arcballInitialZ = g_arcballRbt.getTranslation()[2];
  // arcballInitialZ = g_skyRbt.getTranslation()[2] - g_arcballRbt.getTranslation()[2];
  float radius = 1;
  // currentArcballSize = radius;
  getSphereVbIbLen(slices, stacks, vbLen, ibLen);

  vector<VertexPN> vtx(vbLen);
  vector<unsigned short> idx(ibLen);

  makeSphere(radius, slices, stacks, vtx.begin(), idx.begin());
  g_arcball.reset(new Geometry(&vtx[0], &idx[0], vbLen, ibLen));




}

// takes a projection matrix and send to the the shaders
static void sendProjectionMatrix(const ShaderState& curSS, const Matrix4& projMatrix) {
  GLfloat glmatrix[16];
  projMatrix.writeToColumnMajorMatrix(glmatrix); // send projection matrix
  safe_glUniformMatrix4fv(curSS.h_uProjMatrix, glmatrix);
}

// takes MVM and its normal matrix to the shaders
static void sendModelViewNormalMatrix(const ShaderState& curSS, const Matrix4& MVM, const Matrix4& NMVM) {
  GLfloat glmatrix[16];
  MVM.writeToColumnMajorMatrix(glmatrix); // send MVM
  safe_glUniformMatrix4fv(curSS.h_uModelViewMatrix, glmatrix);

  NMVM.writeToColumnMajorMatrix(glmatrix); // send NMVM
  safe_glUniformMatrix4fv(curSS.h_uNormalMatrix, glmatrix);
}

// update g_frustFovY from g_frustMinFov, g_windowWidth, and g_windowHeight
static void updateFrustFovY() {
  if (g_windowWidth >= g_windowHeight)
    g_frustFovY = g_frustMinFov;
  else {
    const double RAD_PER_DEG = 0.5 * CS380_PI/180;
    g_frustFovY = atan2(sin(g_frustMinFov * RAD_PER_DEG) * g_windowHeight / g_windowWidth, cos(g_frustMinFov * RAD_PER_DEG)) / RAD_PER_DEG;
  }
}

static Matrix4 makeProjectionMatrix() {
  return Matrix4::makeProjection(
           g_frustFovY, g_windowWidth / static_cast <double> (g_windowHeight),
           g_frustNear, g_frustFar);
}

static void drawStuff() {
  // short hand for current shader state
  const ShaderState& curSS = *g_shaderStates[g_activeShader];

  // build & send proj. matrix to vshader
  const Matrix4 projmat = makeProjectionMatrix();
  sendProjectionMatrix(curSS, projmat);

  // use the skyRbt as the eyeRbt
  if (cameraStatus == 0) {
    g_myRbt = g_skyRbt;
  }
  else if (cameraStatus == 1) {
    g_myRbt = g_objectRbt[0];
  }
  else {
    g_myRbt = g_object2Rbt[0];
  }

  const RigTForm eyeRbt = g_myRbt; //eyeframe. rbt means rigit body transformation. in this assignment, we need to change this eyerbt.
  const RigTForm invEyeRbt = inv(eyeRbt);

  //for testing
  // const Matrix4 eyeRbt = g_myRbt; //eyeframe. rbt means rigit body transformation. in this assignment, we need to change this eyerbt.
  // const Matrix4 invEyeRbt = inv(eyeRbt);
 

  const Cvec3 eyeLight1 = Cvec3(invEyeRbt * Cvec4(g_light1, 1)); // g_light1 position in eye coordinates
  const Cvec3 eyeLight2 = Cvec3(invEyeRbt * Cvec4(g_light2, 1)); // g_light2 position in eye coordinates
  safe_glUniform3f(curSS.h_uLight, eyeLight1[0], eyeLight1[1], eyeLight1[2]);
  safe_glUniform3f(curSS.h_uLight2, eyeLight2[0], eyeLight2[1], eyeLight2[2]);

  // draw ground
  // ===========
  //
  const RigTForm groundRbt = RigTForm();  // identity
  Matrix4 MVM = rigTFormToMatrix(invEyeRbt * groundRbt); //E-1O
  Matrix4 NMVM = normalMatrix(MVM);
  sendModelViewNormalMatrix(curSS, MVM, NMVM);
  safe_glUniform3f(curSS.h_uColor, 0.1, 0.95, 0.1); // set color
  g_ground->draw(curSS);

  // draw cubes
  // ==========
  MVM = rigTFormToMatrix(invEyeRbt * g_objectRbt[0]); 
  NMVM = normalMatrix(MVM);
  sendModelViewNormalMatrix(curSS, MVM, NMVM);
  safe_glUniform3f(curSS.h_uColor, g_objectColors[0][0], g_objectColors[0][1], g_objectColors[0][2]);
  g_cube->draw(curSS);

  //draw another cube
  MVM = rigTFormToMatrix(invEyeRbt * g_object2Rbt[0]); 
  NMVM = normalMatrix(MVM);
  sendModelViewNormalMatrix(curSS, MVM, NMVM);
  safe_glUniform3f(curSS.h_uColor, g_object2Colors[0][0], g_object2Colors[0][1], g_object2Colors[0][2]);
  g_cube->draw(curSS);

  //draw arcball

  bool canDrawArcball = false;
  if (cameraStatus == CAMERA_SKY) {
    if (manipulationStatus != MANIPULATION_SKY) {
      canDrawArcball = true;
    }
    else {
      if (worldSkyFrameStatus == WORLDSKY) {
        canDrawArcball = true;
      }
    }
  }
  if (cameraStatus == CAMERA_CUBE1) {
    if (manipulationStatus != MANIPULATION_CUBE1) {
      canDrawArcball = true;
    }
  }
  if (cameraStatus == CAMERA_CUBE2) {
    if (manipulationStatus != MANIPULATION_CUBE2) {
      canDrawArcball = true;
    }
  }
  
  if (canDrawArcball) {
    double scaleFactor = getScreenToEyeScale((invEyeRbt * g_arcballRbt).getTranslation()(2), g_frustFovY, g_windowHeight);
    // double scaleFactor = invEyeRbt.getTranslation()[2];
    g_currentArcballRadius = scaleFactor * g_arcballScreenRadius;
    cout << "currentArcballSize is " << g_currentArcballRadius << "\n";
    MVM = rigTFormToMatrix(invEyeRbt * g_arcballRbt) * Matrix4().makeScale(Cvec3(g_currentArcballRadius,g_currentArcballRadius,g_currentArcballRadius));
      
    NMVM = normalMatrix(MVM);
    sendModelViewNormalMatrix(curSS, MVM, NMVM);

      
    safe_glUniform3f(curSS.h_uColor, g_arcballColor[0], g_arcballColor[1], g_arcballColor[2]);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    g_arcball->draw(curSS);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    
      
      
    
    
  }
  
}

static void display() {
  glUseProgram(g_shaderStates[g_activeShader]->program);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);                   // clear framebuffer color&depth

  drawStuff();

  glutSwapBuffers();                                    // show the back buffer (where we rendered stuff)

  checkGlErrors();
}

static void reshape(const int w, const int h) {
  g_windowWidth = w;
  g_windowHeight = h;
  g_arcballScreenRadius = 0.25 * min(g_windowHeight, g_windowWidth);
  glViewport(0, 0, w, h);
  cerr << "Size of window is now " << w << "x" << h << endl;
  updateFrustFovY();
  glutPostRedisplay();
}

Cvec3 isIntersectWithArcball(double clickX, double clickY) {
  Cvec2 clickInput(clickX, clickY);

  Cvec3 ray = getModelViewRay(clickInput, g_frustFovY, g_windowWidth, g_windowHeight);

  Cvec3 u = ray;
  Cvec3 o = g_myRbt.getTranslation();
  Cvec3 center = g_arcballRbt.getTranslation();
  double r = g_currentArcballRadius;

  double a = norm2(ray);
  double b = 2 * dot(ray, (o - center));
  double c = norm2(o - center) - r * r;

  bool isIntersect = b * b - 4 * a * c > 0;

  if (isIntersect) {
    double d = (-b - sqrt(b * b - 4 * a * c)) / (2 * a);
    Cvec3 result = o + u * d;
    return result;
  }
  else {
    Cvec3();
  }

  
}

bool IsZeroVector(Cvec3& a) {
  return ((a(0) < CS380_EPS) && (a(1) < CS380_EPS) && (a(2) < CS380_EPS));
}

static void motion(const int x, const int y) {

  const double dx = x - g_mouseClickX;
  const double dy = g_windowHeight - y - 1 - g_mouseClickY;

  //this part is checking clicked position information

  Cvec2 clickInput(g_mouseClickX, g_mouseClickY);

 

  Cvec3 p1 = isIntersectWithArcball(g_mouseClickX, g_mouseClickY);
  Cvec3 p2 = isIntersectWithArcball(g_mouseClickX + dx, g_mouseClickY + dy);
  

  Cvec3 v1 = p1 - g_arcballRbt.getTranslation();
  Cvec3 v2 = p2 - g_arcballRbt.getTranslation();

  bool needToFollowArcball = false;
  Quat targetQuaternion = Quat();

  

  if (!IsZeroVector(p1)) {
    needToFollowArcball = true;
      
    Quat quatV2 = Quat(0, v2(0), v2(1), v2(2));
    Quat quatV1 = Quat(0, v1(0), v1(1), v1(2));

    targetQuaternion =  quatV2 * (quatV1 * -1);
    
    
    
    
  }
  else {
    cout << "do not intersect with arcball\n";
  }

  



  RigTForm A;
  RigTForm AMAI;

  if (cameraStatus == CAMERA_SKY) { // when current eyeframe is sky frame
    if (manipulationStatus == MANIPULATION_CUBE1) { // when the current manipulation is cube1

      A.setTranslation(g_objectRbt[0].getTranslation());
      A.setRotation(g_skyRbt.getRotation());
    }
    else if (manipulationStatus == MANIPULATION_CUBE2) { // when the current manipulation is the cube 2

      A.setTranslation(g_object2Rbt[0].getTranslation());
      A.setRotation(g_skyRbt.getRotation());

      
    }
    else { //when the current manipulation is the skyframe => we need to handle input m in this condition

      if (worldSkyFrameStatus == SKYSKY) { // sky sky frame


        A.setTranslation(g_skyRbt.getTranslation());
        A.setRotation(g_skyRbt.getRotation());
      }
      else { // world sky frame

        A.setTranslation(RigTForm().getTranslation());
        A.setRotation(g_skyRbt.getRotation());

        //this code 
      }
    }
  }

  else if (cameraStatus == CAMERA_CUBE1) { // when the current eyeframe is cube 1
    if (manipulationStatus == MANIPULATION_CUBE1) { //manipulation is cub1

      A.setTranslation(g_objectRbt[0].getTranslation());
      A.setRotation(g_objectRbt[0].getRotation());
    }
    else if (manipulationStatus == MANIPULATION_CUBE2) { //manipulation is cube2

      A.setTranslation(g_object2Rbt[0].getTranslation());
      A.setRotation(g_objectRbt[0].getRotation());
    }
    else { //manipulation is skycamera - NOT ALLOWED!
      return;

    }
  }

  else { //when the current eyeframe is cube2
    if (manipulationStatus == MANIPULATION_CUBE1) {

      A.setTranslation(g_objectRbt[0].getTranslation());
      A.setRotation(g_object2Rbt[0].getRotation());
    }
    else if (manipulationStatus == MANIPULATION_CUBE2) {

      A.setTranslation(g_object2Rbt[0].getTranslation());
      A.setRotation(g_object2Rbt[0].getRotation());
    }
    else { //manipulation is sky camera - NOT ALLOWED!
      return;

    }
  }

  RigTForm m;
  if (g_mouseLClickButton && !g_mouseRClickButton) { // left button down?
    if ((cameraStatus == CAMERA_SKY && manipulationStatus == MANIPULATION_SKY) || (cameraStatus == CAMERA_CUBE1 && manipulationStatus == MANIPULATION_CUBE1
    || (cameraStatus == CAMERA_CUBE2 && manipulationStatus == MANIPULATION_CUBE2)
    )) {
      RigTForm xRotation = RigTForm().setRotation(RigTForm().getRotation().makeXRotation(dy));
      RigTForm yRotation = RigTForm().setRotation(RigTForm().getRotation().makeYRotation(-dx));
      m = xRotation * yRotation;

      if (needToFollowArcball) {
        m = RigTForm().setRotation(targetQuaternion);
      }

    }
    else {
      RigTForm xRotation = RigTForm().setRotation(RigTForm().getRotation().makeXRotation(-dy));
      RigTForm yRotation = RigTForm().setRotation(RigTForm().getRotation().makeYRotation(dx));
      m = xRotation * yRotation;

      if (needToFollowArcball) {
        m = RigTForm().setRotation(targetQuaternion);
      }
    }
    
  }
  else if (g_mouseRClickButton && !g_mouseLClickButton) { // right button down?
    if (cameraStatus == CAMERA_SKY && manipulationStatus == MANIPULATION_SKY && worldSkyFrameStatus == WORLDSKY) {
      m = RigTForm().setTranslation(Cvec3(-dx, -dy, 0) * 0.01);
    }
    else {
      m = RigTForm().setTranslation(Cvec3(dx, dy, 0) * 0.01);
    }
  }
  else if (g_mouseMClickButton || (g_mouseLClickButton && g_mouseRClickButton)) {  // middle or (left and right) button down?
    m = RigTForm().setTranslation(Cvec3(0, 0, -dy) * 0.01);
  }

  AMAI = A * m * inv(A);

  if (g_mouseClickDown) {
    if (needToFollowArcball) {
      g_arcballRbt = AMAI * g_arcballRbt;
    }

    if (manipulationStatus == MANIPULATION_CUBE1) {
   

      g_objectRbt[0] = AMAI * g_objectRbt[0]; // Simply right-multiply is WRONG

      if (cameraStatus != CAMERA_CUBE1) {

        g_arcballRbt.setTranslation(g_objectRbt[0].getTranslation());
      }
    }
    else if (manipulationStatus == MANIPULATION_CUBE2) {

      g_object2Rbt[0] = AMAI * g_object2Rbt[0]; // Simply right-multiply is WRONG

      if (cameraStatus != CAMERA_CUBE2) {
        g_arcballRbt.setTranslation(g_object2Rbt[0].getTranslation());
      }
     
    }
    else {

      g_skyRbt = AMAI * g_skyRbt;

      if (worldSkyFrameStatus == WORLDSKY) {
        g_arcballRbt.setTranslation(RigTForm().getTranslation());
      }

    }
    glutPostRedisplay(); // we always redraw if we changed the scene
  }

  g_mouseClickX = x;
  g_mouseClickY = g_windowHeight - y - 1;
}


static void mouse(const int button, const int state, const int x, const int y) {
  g_mouseClickX = x;
  g_mouseClickY = g_windowHeight - y - 1;  // conversion from GLUT window-coordinate-system to OpenGL window-coordinate-system

  g_mouseLClickButton |= (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN);
  g_mouseRClickButton |= (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN);
  g_mouseMClickButton |= (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN);

  g_mouseLClickButton &= !(button == GLUT_LEFT_BUTTON && state == GLUT_UP);
  g_mouseRClickButton &= !(button == GLUT_RIGHT_BUTTON && state == GLUT_UP);
  g_mouseMClickButton &= !(button == GLUT_MIDDLE_BUTTON && state == GLUT_UP);

  g_mouseClickDown = g_mouseLClickButton || g_mouseRClickButton || g_mouseMClickButton;
}


static void keyboard(const unsigned char key, const int x, const int y) {
  switch (key) {
  case 27:
    exit(0);                                  // ESC
  case 'h':
    cout << " ============== H E L P ==============\n\n"
    << "h\t\thelp menu\n"
    << "s\t\tsave screenshot\n"
    << "f\t\tToggle flat shading on/off.\n"
    << "o\t\tCycle object to edit\n"
    << "v\t\tCycle view\n"
    << "drag left mouse to rotate\n" << endl;
    break;
  case 's':
    glFlush();
    writePpmScreenshot(g_windowWidth, g_windowHeight, "out.ppm");
    break;
  case 'f':
    g_activeShader ^= 1;
    break;
  case 'v':
    cout << "v key pressed";
    if (cameraStatus == 0) {
   
      cameraStatus = 1;
    } 
    else if (cameraStatus == 1) {

      cameraStatus = 2;
    }
    else {
 
      cameraStatus = 0;
    }
    
    break;
  case 'o':
    cout << "o key was pressed!";
    if (manipulationStatus == 0) {
   
      manipulationStatus = 1;
    } 
    else if (manipulationStatus == 1) {

      manipulationStatus = 2;
    }
    else {
 
      manipulationStatus = 0;
    }
    //o 키가 눌렸을 때 현재 manipulated 되고 있는 object 를 바꿔야 함
    //
    break;
  case 'm':
    

    if (worldSkyFrameStatus == 0) {
      worldSkyFrameStatus = 1;
    }
    else {
      worldSkyFrameStatus = 0;
    }

    Matrix4 baseMatrix;
    for (int y = 0 ; y < 4; y ++) {
      for (int x = 0; x < 4; x++) {
        baseMatrix(y, x) = y * 4.0 + double(x) + 1.0;
      }
    }
    cout << "given matrix is \n";
    PrintMatrix(baseMatrix);
    cout << "transfact matrix is \n";
    PrintMatrix(transFact(baseMatrix));
    cout << "linfact matrix is \n";
    PrintMatrix(linFact(baseMatrix));

    cout <<"-1 product of transfact matrix is\n";
    PrintMatrix(transFact(baseMatrix * -1));

    cout << "-1 product of linfact matrix is \n";
    PrintMatrix(linFact(baseMatrix * -1));
    
  }
  glutPostRedisplay();
}

static void initGlutState(int argc, char * argv[]) {
  glutInit(&argc, argv);                                  // initialize Glut based on cmd-line args
  glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH);  //  RGBA pixel channels and double buffering
  glutInitWindowSize(g_windowWidth, g_windowHeight);      // create a window
  glutCreateWindow("CS380: Assignment 3");                       // title the window

  glutDisplayFunc(display);                               // display rendering callback
  glutReshapeFunc(reshape);                               // window reshape callback
  glutMotionFunc(motion);                                 // mouse movement callback
  glutMouseFunc(mouse);                                   // mouse click callback
  glutKeyboardFunc(keyboard);
}

static void initGLState() {
  glClearColor(128./255., 200./255., 255./255., 0.);
  glClearDepth(0.);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glCullFace(GL_BACK);
  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_GREATER);
  glReadBuffer(GL_BACK);
  if (!g_Gl2Compatible)
    glEnable(GL_FRAMEBUFFER_SRGB);
}

static void initShaders() {
  g_shaderStates.resize(g_numShaders);
  for (int i = 0; i < g_numShaders; ++i) {
    if (g_Gl2Compatible)
      g_shaderStates[i].reset(new ShaderState(g_shaderFilesGl2[i][0], g_shaderFilesGl2[i][1]));
    else
      g_shaderStates[i].reset(new ShaderState(g_shaderFiles[i][0], g_shaderFiles[i][1]));
  }
}

static void initGeometry() {
  initGround();
  initCubes(); //we can initialize other objects in here 
  initArcball();
}

int main(int argc, char * argv[]) {
  try {
    initGlutState(argc,argv);

    glewInit(); // load the OpenGL extensions

    cout << (g_Gl2Compatible ? "Will use OpenGL 2.x / GLSL 1.0" : "Will use OpenGL 3.x / GLSL 1.3") << endl;
    if ((!g_Gl2Compatible) && !GLEW_VERSION_3_0)
      throw runtime_error("Error: card/driver does not support OpenGL Shading Language v1.3");
    else if (g_Gl2Compatible && !GLEW_VERSION_2_0)
      throw runtime_error("Error: card/driver does not support OpenGL Shading Language v1.0");

    initGLState();
    initShaders();
    initGeometry();

    glutMainLoop();
    return 0;
  }
  catch (const runtime_error& e) {
    cout << "Exception caught: " << e.what() << endl;
    return -1;
  }
}
