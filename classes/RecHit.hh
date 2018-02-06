#ifndef RECHIT_H
#define RECHIT_H

#include "TLorentzVector.h"
#include "RecHitCalibration.hh"

class RecHit{

    private:
        float _thickness;
        float _layer;
        float _puOffset;
        TLorentzVector _mom;
        TLorentzVector _pos;


    public:
        
        // constructors
        RecHit(RecHit&);
        RecHit(TLorentzVector p4, TLorentzVector pos, int layer);
        RecHit(TLorentzVector p4, TLorentzVector pos, int layer, float thickness);

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
        int   layer();
        TLorentzVector p4();
        TLorentzVector pos();
};

#endif 
