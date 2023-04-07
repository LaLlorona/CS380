#ifndef RIGTFORM_H
#define RIGTFORM_H

#include <iostream>
#include <cassert>

#include "matrix4.h"
#include "quat.h"

class RigTForm {
  Cvec3 t_; // translation component
  Quat r_;  // rotation component represented as a quaternion

public:
  RigTForm() : t_(0) {
    assert(norm2(Quat(1,0,0,0) - r_) < CS380_EPS2);
  }

  RigTForm(const Cvec3& t, const Quat& r) {
    //TODO
    t_ = t;
    r_ = r;
  }

  explicit RigTForm(const Cvec3& t) {
    // TODO
    t_ = t;
    r_ = Quat();
  }

  explicit RigTForm(const Quat& r) {
    // TODO
    t_ = Cvec3();
    r_ = r;
  }

  Cvec3 getTranslation() const {
    return t_;
  }

  Quat getRotation() const {
    return r_;
  }

  RigTForm& setTranslation(const Cvec3& t) {
    t_ = t;
    return *this;
  }

  RigTForm& setRotation(const Quat& r) {
    r_ = r;
    return *this;
  }

  

  Cvec4 operator * (const Cvec4& a) const {
    // TODO
    //full afine matrix 에다가 Cvec4를 곱한다. 
    //왜 필요한가? - Matrix4 에 Cvec4를 곱하는 함수가 존재하기 때문.

    // return r_ * (a + Cvec4(t_, 0));
    return r_ * a + Cvec4(t_, 0);
  }

  RigTForm operator * (const RigTForm& a) const {
    // TODO
    RigTForm result;

    Cvec3 t2 = a.getTranslation();
    Quat r2 = a.getRotation();

    Cvec3 t1 = t_;
    Quat r1 = r_;

    //this part should be changed

    // Cvec4 t2Cvec4(0, t2(0), t2(1), t2(2));

    Quat cvec3ToCvec4Quat(0, t2(0), t2(1), t2(2));

    Quat r1t2Quat = r1 * (cvec3ToCvec4Quat * inv(r1));
    Cvec4 r1t2(r1t2Quat(0), r1t2Quat(1), r1t2Quat(2), r1t2Quat(3)); 


    Cvec4 t1_r1t2 = Cvec4(0, t1(0), t1(1), t1(2)) + r1t2;

    Cvec3 resultTrans(t1_r1t2(1), t1_r1t2(2), t1_r1t2(3));

    Quat resultQuaterion = r1 * r2;

    result.setTranslation(resultTrans);
    result.setRotation(resultQuaterion);

    return result;
  }
};



inline RigTForm inv(const RigTForm& tform) {
  // TODO
  RigTForm result;

  Quat r = tform.getRotation();
  Cvec3 t = tform.getTranslation();


  //this part should be changed 
  // Cvec4 tCvec4(0, t(0), t(1), t(2));

  Quat tCvec4Quat(0, t(0), t(1),t(2));

  Quat minusRInvT = inv(r) * (tCvec4Quat * inv(inv(r))) * -1;

  Cvec4 minusRInvTCvec4(minusRInvT(0), minusRInvT(1), minusRInvT(2), minusRInvT(3));

  Cvec3 transResultCvec3(minusRInvTCvec4(1), minusRInvTCvec4(2), minusRInvTCvec4(3));

  result.setTranslation(transResultCvec3);
  result.setRotation(inv(r));

  return result;

}

inline void PrintTForm(const RigTForm& tform) {
  std::cout << "transpart is :";
  for (int i = 0; i < 3; i++) {
    std::cout << tform.getTranslation()(i) << " ";
  }
  std::cout <<"\n";

  std::cout <<"rotation part is:";

  for (int i =0 ; i < 4; i++) {
    std::cout << tform.getRotation()[i] << " ";
  }
  std::cout <<"\n";
}

inline RigTForm transFact(const RigTForm& tform) {
  return RigTForm(tform.getTranslation());
}

inline RigTForm linFact(const RigTForm& tform) {
  return RigTForm(tform.getRotation());
}


inline Matrix4 rigTFormToMatrix(const RigTForm& tform) {
  // TODO
  

  
  //translation part
  Cvec3 translation = tform.getTranslation();
  Matrix4 transMatrix = cvec3ToMatrix4(translation);

  //rotation part
  Quat rotation = tform.getRotation();
  Matrix4 rotationMatrix = quatToMatrix(rotation);
  
  return transMatrix * rotationMatrix;
}





#endif
