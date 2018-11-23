#include "GenParticle.hh"
#include <iostream>

GenParticle::GenParticle(const TLorentzVector& p4, int pdgid, unsigned status,
                         const TVector3& vertex, float charge) {
  _mom = p4;
  _pdgid = pdgid;
  _status = status;
  _vertex = vertex;

  _mom.SetZ(p4.Z() + vertex.Z());

  _charge = charge;
}

GenParticle::GenParticle(const GenParticle& g) {
  _pdgid = g._pdgid;
  _status = g._status;
  _mom = g._mom;
  _vertex = g._vertex;
  _charge = g._charge;
}

bool GenParticle::isClusterable() const {
  bool pass = true;
  // primary
  if (_status != 1) pass = false;
  // do not use neutrinos
  if (fabs(_pdgid) == 12) pass = false;
  if (fabs(_pdgid) == 14) pass = false;
  if (fabs(_pdgid) == 16) pass = false;
  return pass;
}

float GenParticle::eta() const { return _mom.Eta(); }
float GenParticle::phi() const { return _mom.Phi(); }
float GenParticle::pt() const { return _mom.Pt(); }
float GenParticle::energy() const { return _mom.Energy(); }
float GenParticle::px() const { return _mom.Px(); }
float GenParticle::py() const { return _mom.Py(); }
float GenParticle::pz() const { return _mom.Pz(); }
float GenParticle::mass() const { return _mom.M(); }
float GenParticle::charge() const { return _charge; }
TLorentzVector GenParticle::p4() const { return _mom; }
int GenParticle::status() const { return _status; }
int GenParticle::pdgid() const { return _pdgid; }
TVector3 GenParticle::vertex() const { return _vertex; }
