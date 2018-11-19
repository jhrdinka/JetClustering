#include "HitAnalysis.hh"
#include <iostream>

HitAnalysis::HitAnalysis(const TString treeName)
    : _treeName(treeName),
      _layerRZ(0.),
      _layerID(0),
      _pt(0.),
      _distancesRZ({}),
      _meanDrz(0.),
      _minDrz(0.),
      _maxDrz(0.),
      _distancesS({}),
      _meanDs(0.),
      _minDs(0.),
      _maxDs(0.),
      _deltaR({}),
      _hitsWithinCluster({}),
      _hitsWithinPitch({}) {
  _outputTree = new TTree(_treeName, "HitAnalysis");
  if (!_outputTree) throw std::bad_alloc();

  // Initial parameters
  _outputTree->Branch("layerRZ", &_layerRZ);
  _outputTree->Branch("layerID", &_layerID);
  _outputTree->Branch("pt", &_pt);
  _outputTree->Branch("distancesRZ", &_distancesRZ);
  _outputTree->Branch("meanDrz", &_meanDrz);
  _outputTree->Branch("minDrz", &_minDrz);
  _outputTree->Branch("maxDrz", &_maxDrz);
  _outputTree->Branch("distancesS", &_distancesS);
  _outputTree->Branch("meanDs", &_meanDs);
  _outputTree->Branch("minDs", &_minDs);
  _outputTree->Branch("maxDs", &_maxDs);
  _outputTree->Branch("deltaR", &_deltaR);
  _outputTree->Branch("hitsWithinCluster", &_hitsWithinCluster);
  _outputTree->Branch("hitsWithinPitch", &_hitsWithinPitch);
}

void HitAnalysis::fill(const float& layerRZ, const int& layerID,
                       const float& pt, const std::vector<float>& deltaRZ,
                       const float& meanRZ, const float& minRZ,
                       const float& maxRZ, const std::vector<float>& deltaS,
                       const float& meanS, const float& minS, const float& maxS,
                       const std::vector<float>& deltaR,
                       const std::vector<int>& hitsWithinCluster,
                       const std::vector<int>& hitsWithinPitch) {
  _layerRZ = layerRZ;
  _layerID = layerID;
  _pt = pt;
  _distancesRZ = deltaRZ;
  _meanDrz = meanRZ;
  _minDrz = minRZ;
  _maxDrz = maxRZ;
  _distancesS = deltaS;
  _meanDs = meanS;
  _minDs = minS;
  _maxDs = maxS;
  _deltaR = deltaR;
  _hitsWithinCluster = hitsWithinCluster;
  _hitsWithinPitch = hitsWithinPitch;

  _outputTree->Fill();
}

void HitAnalysis::write() {
  //  gDirectory->mkdir(_treeName);
  //  gDirectory->cd(_treeName);

  _outputTree->Write();

  // gDirectory->cd("../..");
}
