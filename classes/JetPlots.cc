#include "JetPlots.hh"

JetPlots::JetPlots(const TString name, const vector<float> ptvals){
   
   _name    = name;   
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


  //create dictionary for resolution plots
  for(vector<float>::const_iterator it = ptvals.begin(); it != ptvals.end()-1; it++) {
      float ptmin = *it;
      float ptmax = *(it+1);
      TString ptbin_str;
      ptbin_str.Form("_%.0f_%.0f",ptmin,ptmax);
      ptbin_str = name + ptbin_str;
      _ptbins.push_back(make_pair(make_pair(ptmin, ptmax), TH1F(ptbin_str,ptbin_str,100,0.,4.0)));
  }

}

void JetPlots::fill(JetCollection& coll){
   for (unsigned i = 0; i < coll.size(); i++) {

       // fill first two only for now
       if(i > 1) break;

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

       // fill resolution plots
       for(vector<pair<pair<float,float>, TH1F>>::iterator it = _ptbins.begin(); it != _ptbins.end()-1; it++) {
           float ptmin=(it->first).first;
           float ptmax=(it->first).second;
           TH1F histo = it->second;
           cout<<ptmin<<","<<ptmax<<endl;
           
           Jet *genjet = jet.ref();
           if(genjet){
               if(genjet->pt() > ptmin && genjet->pt() < ptmax){
                   (it->second).Fill(jet.pt()/genjet->pt());
               }
               
           }
           
       }
   
   }
}

void JetPlots::write(){
   
   gDirectory->mkdir(_name);
   gDirectory->cd(_name);
   
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

   // now fill resolution plots
   gDirectory->mkdir("reso");
   gDirectory->cd("reso");
   
   for(vector<pair<pair<float,float>, TH1F>>::const_iterator it = _ptbins.begin(); it != _ptbins.end()-1; it++) {
       (it->second).Write();
   }

   gDirectory->cd("../..");
}
