#include "GenParticleAnalysis.hh"
#include <iostream>

GenParticleAnalysis::GenParticleAnalysis(const std::string& analysisName,
                                         float ptMin, float ptMax,
                                         bool particlesInJets)
    : _particlesInJets(particlesInJets) {
  _eta_h = new TH1F((analysisName + "_eta_h").c_str(), "particles in #eta",
                    120., -6., 6.);
  _pt_h = new TH1F((analysisName + "_pt_h").c_str(), "particle p_{T}", 600.,
                   ptMin, ptMax);
  _ptOverEta_h = new TProfile((analysisName + "_ptOverEta_h").c_str(),
                              "particle p_{T} over #eta", 120., -6., 6.);
  _ptOverPhi_h = new TProfile((analysisName + "_ptOverPhi_h").c_str(),
                              "particle p_{T} over #phi", 120., -M_PI, M_PI);

  _eta_h->SetDirectory(0);
  _pt_h->SetDirectory(0);
  _ptOverEta_h->SetDirectory(0);
  _ptOverPhi_h->SetDirectory(0);

  if (_particlesInJets) {
    _tracks_deltaR = new TH1F((analysisName + "_tracks_deltaR").c_str(),
                              "Tracks within jet", 120., 0., 0.4);
  }
}

void GenParticleAnalysis::fill(const GenParticle& genParticle) {
  // only write charge particles
  if (genParticle.charge()) {
    _eta_h->Fill(genParticle.eta());
    _pt_h->Fill(genParticle.pt());
    _ptOverEta_h->Fill(genParticle.eta(), genParticle.pt());
    _ptOverPhi_h->Fill(genParticle.phi(), genParticle.pt());
  }
}

void GenParticleAnalysis::fillParticleAndJets(const GenParticle& genParticle,
                                              JetCollection& jets,
                                              float jetCone) {
  // only write charge particles
  if (!_particlesInJets) {
    throw std::runtime_error(
        "You need to turn on 'particlesInJets'-flag if you want to use "
        "this function!");
  }
  if (genParticle.charge()) {
    _eta_h->Fill(genParticle.eta());
    _pt_h->Fill(genParticle.pt());
    _ptOverEta_h->Fill(genParticle.eta(), genParticle.pt());
    _ptOverPhi_h->Fill(genParticle.phi(), genParticle.pt());
    // go through jets and check if particle is within one of them
    float deltaR = jetCone;
    for (size_t i = 0; i < jets.size(); i++) {
      // first access jet
      auto pJet = jets.at(i)->p4();
      // delta R
      float dR = pJet.DeltaR(genParticle.p4());
      // delta R of vertex (check if particle is actually within the cone)
      // correct z vertex
      TVector3 vertex(genParticle.vertex().X(), genParticle.vertex().Y(),
                      (genParticle.vertex().Z() - jets.at(i)->vertexZ()));

      float dR_vertex = pJet.Vect().DeltaR(vertex);

      // assign to closer jet
      if (dR < deltaR) {
        if ((vertex.Perp() > 20.) && (dR_vertex < jetCone)) {
          deltaR = dR;
        } else {
          if (std::fabs(jets.at(i)->vertexZ() - vertex.z()) < 1.) deltaR = dR;
        }
      }
    }
    if (deltaR < jetCone) {
      _tracks_deltaR->Fill(deltaR);
    }
  }
}

void GenParticleAnalysis::write() {
  _eta_h->Write();
  _pt_h->Write();
  _ptOverEta_h->Write();
  _ptOverPhi_h->Write();
  if (_particlesInJets) {
    _tracks_deltaR->Write();
  }
}

void GenParticleAnalysis::normalize(float norm) {
  float scalor = 1. / norm;
  if (norm == 0.) {
    throw std::runtime_error(
        "TrackClusterAnalysis::normalizeMergedClusters:Trying to divide by "
        "zero!");
  }
  _eta_h->Scale(scalor);
  _pt_h->Scale(scalor);
  if (_particlesInJets) {
    _tracks_deltaR->Scale(scalor);
    ;
  }
}
