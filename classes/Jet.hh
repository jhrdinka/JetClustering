#ifndef JET_H
#define JET_H

#include "TLorentzVector.h"
#include "TVector3.h"

class Jet {
 private:
  TLorentzVector _mom;
  float _tau1;
  float _tau2;
  float _tau3;
  float _tau21;
  float _tau32;
  float _sdmass;
  float _vertexZ;

  Jet* _ref = NULL;

 public:
  // constructors
  Jet();
  Jet(const TLorentzVector& p4, float vertexZ = 0);
  Jet(const Jet& j);

  void setTau1(float tau1);
  void setTau2(float tau2);
  void setTau3(float tau3);
  void setTau21(float tau21);
  void setTau32(float tau32);
  void setSDmass(float sdmass);

  void setRef(Jet*);

  float eta() const;
  float phi() const;
  float pt() const;
  float energy() const;
  float px() const;
  float py() const;
  float pz() const;
  float mass() const;

  float tau1() const;
  float tau2() const;
  float tau3() const;
  float tau21() const;
  float tau32() const;

  float massSD() const;

  TLorentzVector p4() const;

  float vertexZ() const;

  Jet* ref();

  void print() const;
};

#endif
