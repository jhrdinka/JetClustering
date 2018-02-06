#ifndef GENPARTICLE_H
#define GENPARTICLE_H

#include "TLorentzVector.h"

class GenParticle{

    private:
        int   _pdgid;
        int   _status;
        TLorentzVector _mom;


    public:
        
        // constructor
        GenParticle(TLorentzVector p4, int pdgid, int status);
        GenParticle(GenParticle& g);

        bool isClusterable();
        
        float eta();
        float phi();
        float pt();
        float energy();
        float px();
        float py();
        float pz();
        float mass();
        TLorentzVector p4();
        int   pdgid();
        int   status();


};

#endif 
