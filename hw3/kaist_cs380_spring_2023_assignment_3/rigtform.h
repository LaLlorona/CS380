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
  }

  RigTForm operator * (const RigTForm& a) const {
    // TODO
    Matrix4 thisToMatrix4 = rigTFormToMatrix((RigTForm)*this);
    Matrix4 aToMatrix4 = rigTFormToMatrix(a);

    Matrix4 t1 = transFact(thisToMatrix4);
    Matrix4 r1 = linFact(thisToMatrix4);

    Matrix4 t2 = transFact(aToMatrix4);
    Matrix4 r2 = linFact(aToMatrix4);

    Matrix4 resultTrans = transFact(t1 + r1 * t2);
    Matrix4 resultRot = linFact(r1 * r2);

    //question: how can I convert matrix4 to quaternion?????????
    // 굳이 Matrix4 로 변환한 필요가 없다. 

    //r1: quaternion 으로 정의됨. 



    

  }
};

inline RigTForm inv(const RigTForm& tform) {
  // TODO
  RigTForm result;

  Quat tformRot = tform.getRotation();
  Cvec3 tformTrans = tform.getTranslation();
  Cvec4 tformTransCvec4(tformTrans(0), tformTrans(1), tformTrans(2), 1);

  // Cvec4 transResult = quatToMatrix(inv(tformRot)) * cvec3ToMatrix4(tformTrans) * -1;

  Matrix4 rInv = quatToMatrix(inv(tformRot));

  Cvec4 transResult = rInv * tformTransCvec4 * -1;

  Cvec3 transResultCvec3(transResult(0), transResult(1), transResult(2));
  result.setTranslation(transResultCvec3);
  result.setRotation(inv(tformRot));

  return result;




  

  // Cvec3 transResultCvec3(transResult(0), transResult(1), transResult(2));



  

  // Matrix4 tformToMatrix4 = rigTFormToMatrix(tform);

  // Matrix4 t = transFact(tformToMatrix4);
  // Matrix4 r = linFact(tformToMatrix4);

  // Matrix4 resultTrans = transFact(inv(r) * t * -1);
  // Matrix4 resultRot = linFact(inv(r));

  //question again: Matrix4 에서 쿼터니온으로 어떻게 바꾸????????????????????


}

inline RigTForm transFact(const RigTForm& tform) {
  return RigTForm(tform.getTranslation());
}

inline RigTForm linFact(const RigTForm& tform) {
  return RigTForm(tform.getRotation());
}


inline Matrix4 rigTFormToMatrix(const RigTForm& tform) {
  // TODO
  
  Matrix4 rotMatrix = Matrix4();
  
  //translation part
  Cvec3 translation = tform.getTranslation();
  Matrix4 transMatrix = cvec3ToMatrix4(translation);

  //rotation part
  Quat rotation = tform.getRotation();
  Matrix4 rotationMatrix = quatToMatrix(rotation);
  
  return transMatrix * rotMatrix;
}

inline Matrix4 cvec3ToMatrix4(const Cvec3& a) {
  Matrix4 result = Matrix4();
  for (int y = 0; y < 3; y++) {
    result(y, 3) = a[y];
  }
  return result;
}

#endif
