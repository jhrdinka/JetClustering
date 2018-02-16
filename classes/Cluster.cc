#include <iostream>
#include "Cluster.hh"

Cluster::Cluster(Cluster& r){
    _puOffset = r._puOffset;
    _mom = r._mom;
    _pos = r._pos;

};

Cluster::Cluster(TLorentzVector p4, TLorentzVector pos){
   _mom = p4;
   _pos = pos;
}

void Cluster::setP4(TLorentzVector p4){
   _mom = p4;
}

float Cluster::eta() {return _mom.Eta();}
float Cluster::phi() {return _mom.Phi();}
float Cluster::pt() {return _mom.Pt();}
float Cluster::energy() {return _mom.E();}
float Cluster::x() {return _pos.X();}
float Cluster::y() {return _pos.Y();}
float Cluster::z() {return _pos.Z();}
float Cluster::px() {return _mom.Px();}
float Cluster::py() {return _mom.Py();}
float Cluster::pz() {return _mom.Pz();}
TLorentzVector Cluster::p4() {return _mom;}
TLorentzVector Cluster::pos() {return _pos;}

