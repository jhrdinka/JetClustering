#ifndef JET_H
#define JET_H

#include "TLorentzVector.h"

class Jet{

    private:
        TLorentzVector _mom;
        float _tau1;
        float _tau2;
        float _tau3;
        float _tau21;
        float _tau32;
        float _sdmass;

        Jet* _ref = NULL;

    public:
        
        // constructors
        Jet();
        Jet(TLorentzVector p4);
        Jet(Jet& j);

        void setTau1(float tau1);
        void setTau2(float tau2);
        void setTau3(float tau3);
        void setTau21(float tau21);
        void setTau32(float tau32);
        void setSDmass(float sdmass);

        void setRef(Jet*);

        float eta();
        float phi();
        float pt();
        float energy();
        float px();
        float py();
        float pz();
        float mass();

        float tau1();
        float tau2();
        float tau3();
        float tau21();
        float tau32();

        float massSD();

        TLorentzVector p4();

        Jet* ref();

        void print();

};

#endif 
