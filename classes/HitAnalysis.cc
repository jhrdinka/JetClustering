#include "HitAnalysis.hh"

HitAnalysis::HitAnalysis(const TString treeName)
    : _treeName(treeName),
      _layerID(0),
      _r(0.),
      _z(0.),
      _deltaR({}),
      _mean(0.),
      _min(0.),
      _max(0.) {
  _outputTree = new TTree(_treeName, "HitAnalysis");
  if (!_outputTree) throw std::bad_alloc();

  // Initial parameters
  _outputTree->Branch("layerID", &_layerID);
  _outputTree->Branch("r", &_r);
  _outputTree->Branch("z", &_z);
  _outputTree->Branch("deltaR", &_deltaR);
  _outputTree->Branch("mean", &_mean);
  _outputTree->Branch("min", &_min);
  _outputTree->Branch("max", &_max);
}

void HitAnalysis::fill(const int& layerID, const float& r, const float& z,
                       const std::vector<float>& deltaR, float mean, float min,
                       float max) {
  _layerID = layerID;
  _r = r;
  _z = z;
  _deltaR = deltaR;
  _mean = mean;
  _min = min;
  _max = max;

  _outputTree->Fill();
}

void HitAnalysis::write() {
  //  gDirectory->mkdir(_treeName);
  //  gDirectory->cd(_treeName);

  _outputTree->Write();

  // gDirectory->cd("../..");
}
