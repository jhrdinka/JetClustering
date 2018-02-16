#ifndef CLUSTER_H
#define CLUSTER_H

#include "TLorentzVector.h"

class Cluster{

    private:
        float _puOffset;
        TLorentzVector _mom;
        TLorentzVector _pos;

    public:
        
        // constructors
        Cluster(Cluster&);
        Cluster(TLorentzVector p4, TLorentzVector pos);

        void setP4(TLorentzVector p4);

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
        TLorentzVector p4();
        TLorentzVector pos();
};

#endif 
