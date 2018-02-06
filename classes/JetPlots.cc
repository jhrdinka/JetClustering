#include "JetPlots.hh"
#include <iostream>

JetPlots::JetPlots(const TString name){

   _pt      = TH1F(name+"_pt",name+"_pt",100,0.,1000.);
   _eta     = TH1F(name+"_eta",name+"_eta",100,-5.0,5.0);
   _phi     = TH1F(name+"_phi",name+"_phi",100,-TMath::Pi(),-TMath::Pi());
   _mass    = TH1F(name+"_mass",name+"_mass",100,0.,500.);
   _energy  = TH1F(name+"_energy",name+"_energy",100,0.,10000.);
   _tau1    = TH1F(name+"_tau1",name+"_tau1",100,0.,1.0);
   _tau2    = TH1F(name+"_tau2",name+"_tau2",100,0.,1.0);
   _tau3    = TH1F(name+"_tau3",name+"_tau3",100,0.,1.0);
   _tau21   = TH1F(name+"_tau21",name+"_tau21",100,0.,1.0);
   _tau32   = TH1F(name+"_tau32",name+"_tau32",100,0.,1.0);
   _massSD  = TH1F(name+"_massSD",name+"_massSD",100,0.,500.);

}

void JetPlots::fill(JetCollection& coll){
   for (unsigned i = 0; i < coll.size(); i++) {
      
      // fill first two only for now
      
      Jet jet = *(coll.at(i));
      _pt      .Fill(jet.pt());
      _eta     .Fill(jet.eta());
      _phi     .Fill(jet.phi());
      _mass    .Fill(jet.mass());
      _energy  .Fill(jet.energy());
      _tau1    .Fill(jet.tau1());
      _tau2    .Fill(jet.tau2());
      _tau3    .Fill(jet.tau3());
      _tau21   .Fill(jet.tau21());
      _tau32   .Fill(jet.tau32());
      _massSD  .Fill(jet.massSD());
   }
}

void JetPlots::write(){
   _pt      .Write();
   _eta     .Write();
   _phi     .Write();
   _mass    .Write();
   _energy  .Write();
   _tau1    .Write();
   _tau2    .Write();
   _tau3    .Write();
   _tau21   .Write();
   _tau32   .Write();
   _massSD  .Write();

}
