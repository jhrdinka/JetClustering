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

 public:
  // constructors
  RecHit(RecHit&);
  RecHit(TLorentzVector p4, TLorentzVector pos, int layer, int bits);
  RecHit(TLorentzVector p4, TLorentzVector pos, int layer, float thickness,
         int bits);

  void setP4(TLorentzVector p4);

  void setThickness(float thickness);
  void setPuOffset(float);

  bool isAboveThreshold(RecHitCalibration rc, float ecut);
  bool isAbovePuNoise();

  float eta();
  float phi();
  float pt();
  float energy();
  float x();
  float y();
  float z();
  float px();
  float py();
  float pz();
  float thickness();
  int layer();
  TLorentzVector p4();
  TLorentzVector pos();
  int bits();
};

#endif
