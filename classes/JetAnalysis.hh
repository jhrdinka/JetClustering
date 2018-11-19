#ifndef JETANALYSIS_H
#define JETANALYSIS_H

#include <TDirectory.h>
#include <TString.h>
#include <TTree.h>
#include <ios>
#include <stdexcept>
#include "Jet.hh"
#include "TH1F.h"
#include "TProfile.h"

using namespace std;

class JetAnalysis {
 private:
  /// Name of the outputtree
  TString _treeName;
  /// Jet pT
  float _pT;
  /// Jet eta
  float _eta;
  /// Jet pz
  float _phi;
  /// Jet energy
  float _energy;
  /// The output tree
  TTree* _outputTree;
  /// eta distribution
  TH1F* _eta_h;
  /// pt distribution
  TH1F* _pt_h;
  /// pt distribution over eta
  TProfile* _ptOverEta_h;
  /// pt distribution over phi
  TProfile* _ptOverPhi_h;

 public:
  // constructors
  JetAnalysis(const TString treeName, float ptMin, float ptMax);
  void fill(const Jet& jet);
  void write();
  void normalize(float norm = 1);
};

#endif  // JETANALYSIS_H
