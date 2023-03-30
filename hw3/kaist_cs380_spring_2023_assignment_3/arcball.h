#ifndef ARCBALL_H
#define ARCBALL_H

#include <iostream>
#include "cvec.h"
#include "matrix4.h"

// Return the screen space projection in terms of pixels of a 3d point
// given in eye-frame coordinates. 
//
// Ideally you should never call this for a point behind the Z=0 plane,
// sinch such a point wouldn't be visible.  
//
// But if you do pass in a point behind Z=0 plane, we'll just
// print a warning, and return the center of the screen.
inline Cvec2 getScreenSpaceCoord(const Cvec3& p,
                                 const Matrix4& projection,
                                 double frustNear, double frustFovY,
                                 int screenWidth, int screenHeight) {
  if (p[2] > -CS380_EPS) {
    std::cerr << "WARNING: getScreenSpaceCoord of a point near or behind Z=0 plane. Returning screen-center instead." << std::endl;
    return Cvec2((screenWidth-1)/2.0, (screenHeight-1)/2.0);
  }
  Cvec4 q = projection * Cvec4(p, 1);
  Cvec3 clipCoord = Cvec3(q) / q[3];
  return Cvec2(clipCoord[0] * screenWidth / 2.0 + (screenWidth - 1)/2.0,
               clipCoord[1] * screenHeight / 2.0 + (screenHeight - 1)/2.0);
}

// Return the scale between 1 unit in screen pixels and 1 unit in the eye-frame
// (or world-frame, since we always use rigid transformations to represent one
// frame with resepec to another frame)
//
// Ideally you should never call this using a z behind the Z=0 plane,
// sinch such a point wouldn't be visible.  
//
// But if you do pass in a point behind Z=0 plane, we'll just
// print a warning, and return 1
inline double getScreenToEyeScale(double z, double frustFovY, int screenHeight) {
  if (z > -CS380_EPS) {
    std::cerr << "WARNING: getScreenToEyeScale on z near or behind Z=0 plane. Returning 1 instead." << std::endl;
    return 1;
  }
  return -(z * tan(frustFovY * CS380_PI/360.0)) * 2 / screenHeight;
}

inline Cvec3 getModelViewRay(const Cvec2& p, double frustFovY,
                             int screenWidth, int screenHeight) {
  const double aspectRatio = screenWidth / static_cast <double> (screenHeight);
  const double ang = frustFovY * 0.5 * CS380_PI/180;
  const double f = std::abs(std::sin(ang)) < CS380_EPS ? 0 : 1/std::tan(ang);

  const double clipcoord_x = (p[0] - (screenWidth - 1)/2.0) * 2.0 / screenWidth;
  const double clipcoord_y = (p[1] - (screenHeight - 1)/2.0) * 2.0 / screenHeight;

  Cvec3 d(-clipcoord_x * aspectRatio / f, -clipcoord_y / f, 1.0);
  d.normalize();
  return d;
}

#endif

