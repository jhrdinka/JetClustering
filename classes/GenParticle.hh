#ifndef GENPARTICLE_H
#define GENPARTICLE_H

#include "TLorentzVector.h"

class GenParticle {
 private:
  int _pdgid;
  int _status;
  TLorentzVector _mom;
  TVector3 _vertex;
  float _charge;

 public:
  // constructor
  GenParticle(const TLorentzVector& p4, int pdgid, unsigned status,
              const TVector3& vertex, float charge);
  GenParticle(const GenParticle& g);

  bool isClusterable() const;

  float eta() const;
  float phi() const;
  float pt() const;
  float energy() const;
  float px() const;
  float py() const;
  float pz() const;
  float mass() const;
  float charge() const;
  TLorentzVector p4() const;
  int pdgid() const;
  int status() const;
  TVector3 vertex() const;
};

#endif
