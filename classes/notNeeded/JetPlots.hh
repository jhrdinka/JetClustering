#ifndef JETPLOTS_H
#define JETPLOTS_H

#include "Jet.hh"
#include "JetCollection.hh"
#include <vector>
#include <TH1F.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TDirectory.h>
#include <iostream>

using namespace std;

class JetPlots{

    private:

        TString _name;

        TH1F _pt;
        TH1F _eta;
        TH1F _phi;
        TH1F _mass;
        TH1F _energy;
        TH1F _tau1;
        TH1F _tau2;
        TH1F _tau3;
        TH1F _tau21;
        TH1F _tau32;
        TH1F _massSD;

        vector<pair<pair<float,float>, TH1F>> _ptbins;

    public:
        
        // constructors
        JetPlots(const TString, const vector<float> ptvals);
        void fill(JetCollection&);
        void write();

};

#endif 
