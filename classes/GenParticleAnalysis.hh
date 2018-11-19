#ifndef GENPARTICLEANALYSIS_H
#define GENPARTICLEANALYSIS_H

#include <TDirectory.h>
#include <TString.h>
#include <TTree.h>
#include <ios>
#include <stdexcept>
#include "GenParticle.hh"
#include "JetCollection.hh"
#include "TH1F.h"
#include "TProfile.h"

class GenParticleAnalysis {
 private:
  /// flag, if tracks within jets should be written out
  bool _particlesInJets = false;
  /// eta distribution
  TH1F* _eta_h;
  /// pt distribution
  TH1F* _pt_h;
  /// pt distribution over eta
  TProfile* _ptOverEta_h;
  /// pt distribution over phi
  TProfile* _ptOverPhi_h;
  /// number of tracks within jet cone (only filled if 'particlesInJets'
  /// switched on)
  TH1F* _tracks_deltaR;

 public:
  // constructors
  GenParticleAnalysis(const std::string& analysisName, float ptMin, float ptMax,
                      bool particlesInJets = false);
  void fill(const GenParticle& genParticle);
  void fillParticleAndJets(const GenParticle& genParticle, JetCollection& jets,
                           float jetCone);
  void write();
  void normalize(float norm);
};

#endif  // GENPARTICLEANALYSIS_H
