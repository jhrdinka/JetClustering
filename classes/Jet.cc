#include "Jet.hh"
#include <string.h>
#include <iomanip>
#include <iostream>
#include <sstream>

Jet::Jet(){};

Jet::Jet(const TLorentzVector& p4, float vertexZ)
    : _mom(p4), _vertexZ(vertexZ) {
  // shift caused by vertex
  _mom.SetZ(p4.Z() + vertexZ);
}

Jet::Jet(const Jet& j) {
  _mom = j._mom;
  _tau1 = j._tau1;
  _tau2 = j._tau2;
  _tau3 = j._tau3;
  _tau21 = j._tau21;
  _tau32 = j._tau32;
  _sdmass = j._sdmass;
  _ref = j._ref;
  _vertexZ = j._vertexZ;
}

void Jet::setTau1(float tau1) { _tau1 = tau1; }
void Jet::setTau2(float tau2) { _tau2 = tau2; }
void Jet::setTau3(float tau3) { _tau3 = tau3; }
void Jet::setTau21(float tau21) { _tau21 = tau21; }
void Jet::setTau32(float tau32) { _tau32 = tau32; }

void Jet::setRef(Jet* jet) { _ref = jet; }

void Jet::setSDmass(float sdmass) { _sdmass = sdmass; }

float Jet::eta() const { return _mom.Eta(); }
float Jet::phi() const { return _mom.Phi(); }
float Jet::pt() const { return _mom.Pt(); }
float Jet::energy() const { return _mom.Energy(); }
float Jet::px() const { return _mom.Px(); }
float Jet::py() const { return _mom.Py(); }
float Jet::pz() const { return _mom.Pz(); }
float Jet::mass() const { return _mom.M(); }
TLorentzVector Jet::p4() const { return _mom; }
float Jet::tau1() const { return _tau1; }
float Jet::tau2() const { return _tau2; }
float Jet::tau3() const { return _tau3; }
float Jet::tau21() const { return _tau21; }
float Jet::tau32() const { return _tau32; }
float Jet::massSD() const { return _sdmass; }
float Jet::vertexZ() const { return _vertexZ; }

Jet* Jet::ref() { return _ref; }

void Jet::print() const {
  std::string commentStr = "";
  std::cout << commentStr << std::setprecision(4) << std::setw(5) << _mom.Pt()
            << std::setprecision(4) << std::setw(10) << _mom.Eta()
            << std::setprecision(4) << std::setw(11) << _mom.Phi()
            << std::setprecision(4) << std::setw(11) << _mom.M()
            << std::setprecision(4) << std::setw(11) << _mom.E();
  std::cout << std::setprecision(4) << std::setw(11) << _tau21;
  std::cout << std::setprecision(4) << std::setw(11) << _tau32;
  std::cout << std::setprecision(4) << std::setw(11) << _sdmass;
  std::cout << std::endl;
}
