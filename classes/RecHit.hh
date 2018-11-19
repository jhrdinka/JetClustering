#ifndef RECHIT_H
#define RECHIT_H

#include "RecHitCalibration.hh"
#include "TLorentzVector.h"

class RecHit {
 private:
  float _thickness;
  float _layer;
  float _puOffset;
  TLorentzVector _mom;
  TLorentzVector _pos;
  int _bits;
  float _time;

 public:
  // constructors
  RecHit(RecHit&);
  RecHit(TLorentzVector p4, TLorentzVector pos, int layer, int bits,
         float time);
  RecHit(TLorentzVector p4, TLorentzVector pos, int layer, float thickness,
         int bits, float time);

  void setP4(TLorentzVector p4);

  void setThickness(float thickness);
  void setPuOffset(float);

  bool isAboveThreshold(RecHitCalibration rc, float ecut) const;
  bool isAbovePuNoise();

  float eta() const;
  float phi() const;
  float pt() const;
  float energy() const;
  float x() const;
  float y() const;
  float z() const;
  float px() const;
  float py() const;
  float pz() const;
  float thickness() const;
  int layer() const;
  TLorentzVector p4() const;
  TLorentzVector pos() const;
  int bits() const;
  float time() const;
};

#endif
