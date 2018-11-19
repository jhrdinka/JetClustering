#include "JetAnalysis.hh"
#include <iostream>

JetAnalysis::JetAnalysis(const TString treeName, float ptMin, float ptMax)
    : _treeName(treeName), _pT(0.), _eta(0), _phi(0.), _energy() {
  _outputTree = new TTree(_treeName, "JetAnalysis");
  if (!_outputTree) throw std::bad_alloc();

  // Initial parameters
  _outputTree->Branch("pT", &_pT);
  _outputTree->Branch("eta", &_eta);
  _outputTree->Branch("phi", &_phi);
  _outputTree->Branch("energy", &_energy);
  _outputTree->SetDirectory(0);

  _eta_h = new TH1F("jet_eta_h", "jets in #eta", 120., -6., 6.);
  _pt_h = new TH1F("jet_pt_h", "jets p_{T}", 120., ptMin - 100., ptMax + 100.);
  _ptOverEta_h =
      new TProfile("jet_ptOverEta_h", "gjet p_{T} over #eta", 120., -6., 6.);
  _ptOverPhi_h =
      new TProfile("jet_ptOverPhi_h", "jet p_{T} over #phi", 120., -M_PI, M_PI);

  _eta_h->SetDirectory(0);
  _pt_h->SetDirectory(0);
  _ptOverEta_h->SetDirectory(0);
  _ptOverPhi_h->SetDirectory(0);
}

void JetAnalysis::fill(const Jet& jet) {
  _pT = jet.pt();
  _eta = jet.eta();
  _phi = jet.phi();
  _energy = jet.energy();
  _outputTree->Fill();

  _eta_h->Fill(jet.eta());
  _pt_h->Fill(jet.pt());
  _ptOverEta_h->Fill(jet.eta(), jet.pt());
  _ptOverPhi_h->Fill(jet.phi(), jet.pt());
}

void JetAnalysis::write() {
  //  gDirectory->mkdir(_treeName);
  //  gDirectory->cd(_treeName);

  _outputTree->Write();

  _eta_h->Write();
  _pt_h->Write();
  _ptOverEta_h->Write();
  _ptOverPhi_h->Write();
  // _deltaR_h->Write();
}

void JetAnalysis::normalize(float norm) {
  float scalor = 1. / norm;
  if (norm == 0.) {
    throw std::runtime_error(
        "TrackClusterAnalysis::normalizeMergedClusters:Trying to divide by "
        "zero!");
  }
  _eta_h->Scale(scalor);
  _pt_h->Scale(scalor);
}
