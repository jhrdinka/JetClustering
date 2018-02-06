#include <iostream>
#include "GenParticle.hh"

GenParticle::GenParticle(TLorentzVector p4, int pdgid, int status){
   _mom = p4;
   _pdgid = pdgid;
   _status = status;
}

GenParticle::GenParticle(GenParticle &g){
   _pdgid = g._pdgid;
   _status = g._status;
   _mom = g._mom;
}

bool GenParticle::isClusterable(){
   bool pass = true;
   if(_status != 1) pass = false;
   if(fabs(_pdgid) == 12) pass = false;
   if(fabs(_pdgid) == 14) pass = false;
   if(fabs(_pdgid) == 16) pass = false;
   return pass; 
}


float GenParticle::eta() {return _mom.Eta();}
float GenParticle::phi() {return _mom.Phi();}
float GenParticle::pt() {return _mom.Pt();}
float GenParticle::energy() {return _mom.Energy();}
float GenParticle::px() {return _mom.Px();}
float GenParticle::py() {return _mom.Py();}
float GenParticle::pz() {return _mom.Pz();}
float GenParticle::mass() {return _mom.M();}
TLorentzVector GenParticle::p4() {return _mom;}
int GenParticle::status(){return _status;}
int GenParticle::pdgid(){return _pdgid;}
