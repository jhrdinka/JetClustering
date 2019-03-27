#include "TrackClusterAnalysis.hh"
#include <iostream>

TrackClusterAnalysis::TrackClusterAnalysis(const std::string& analysisName)
    : _name(analysisName) {
  _outputTree_cluster =
      new TTree((_name + "_TrackClusterAnalysis_cluster").c_str(),
                "TrackClusterAnalysis_cluster");
  if (!_outputTree_cluster) throw std::bad_alloc();

  // tree cluster - to possibly redo histograms
  _outputTree_cluster->Branch("deltaR_cluster", &_deltaR_cluster);
  _outputTree_cluster->Branch("eta", &_eta);
  _outputTree_cluster->Branch("eta_jet", &_eta_jet);
  _outputTree_cluster->Branch("phi", &_phi);
  _outputTree_cluster->Branch("phi_jet", &_phi_jet);
  _outputTree_cluster->Branch("nTracksPerCluster", &_nTracksPerCluster);
  _outputTree_cluster->Branch("nTracksPerClusterPTCut",
                              &_nTracksPerClusterPTCut);
  _outputTree_cluster->Branch("mergedClusters", &_mergedCluster);
  _outputTree_cluster->Branch("mergedClustersPTCut", &_mergedClusterPTCut);
  _outputTree_cluster->Branch("clusterSizeSingle", &_clusterSizeSingle);
  _outputTree_cluster->Branch("clusterSizeMulti", &_clusterSizeMulti);
  _outputTree_cluster->Branch("clusterSizeSinglePTCut",
                              &_clusterSizeSinglePTCut);
  _outputTree_cluster->Branch("clusterSizeMultiPTCut", &_clusterSizeMultiPTCut);
  _outputTree_cluster->Branch("isPU", &_isPU);
  _outputTree_cluster->SetDirectory(0);

  // histograms and profiles
  _deltaR_h = new TH1F((_name + "_deltaR_cluster").c_str(), "deltaR_cluster",
                       120., 0., 0.4);

  _nTracksPerCluster_dR =
      new TProfile((_name + "_nTracksPerCluster_dR").c_str(),
                   "_nTracksPerCluster_dR", 120, 0., 0.4);

  _nTracksPerCluster_PU_dR =
      new TProfile((_name + "_nTracksPerCluster_PU_dR").c_str(),
                   "_nTracksPerCluster_PU_dR", 120, 0., 0.4);

  _nTracksPerClusterPTCut_dR =
      new TProfile((_name + "_nTracksPerClusterPTCut_dR").c_str(),
                   "_nTracksPerCluster_dR", 120, 0., 0.4);

  _nTracksPerClusterPTCut_PU_dR =
      new TProfile((_name + "_nTracksPerClusterPTCut_PU_dR").c_str(),
                   "_nTracksPerCluster_PU_dR", 120, 0., 0.4);

  _mergedClusters_dR = new TH1F((_name + "_mergedClusters_dR").c_str(),
                                "_mergedClusters_dR", 120, 0., 0.4);
  _mergedClustersPTCut_dR =
      new TH1F((_name + "_mergedClustersPTCut_dR").c_str(),
               "_mergedClusters_dR", 120, 0., 0.4);

  _mergedClusterRate_dR =
      new TProfile((_name + "_mergedClusterRate_dR").c_str(),
                   "_mergedClusterrate_dR", 120, 0., 0.4);
  _mergedClusterRatePTCut_dR =
      new TProfile((_name + "_mergedClusterRatePTCut_dR").c_str(),
                   "_mergedClusterRate_dR", 120, 0., 0.4);

  _mergedClusterRate_PU_dR =
      new TProfile((_name + "_mergedClusterRate_PU_dR").c_str(),
                   "_mergedClusterrate_PU_dR", 120, 0., 0.4);
  _mergedClusterRatePTCut_PU_dR =
      new TProfile((_name + "_mergedClusterRatePTCut_PU_dR").c_str(),
                   "_mergedClusterRate_PU_dR", 120, 0., 0.4);

  _hits_dR = new TH1F((_name + "_hits_dR").c_str(), "_hits_dR", 120, 0., 0.4);

  _clusterSize_single_dR =
      new TProfile((_name + "_clusterSize_single_dR").c_str(),
                   "_clusterSize_single_dR", 120, 0, 0.4);

  _clusterSize_multi_dR =
      new TProfile((_name + "_clusterSize_multi_dR").c_str(),
                   "_clusterSize_multi_dR", 120, 0, 0.4);

  _clusterSizePTCut_single_dR =
      new TProfile((_name + "_clusterSizeTCut_single_dR").c_str(),
                   "_clusterSize_single_dR", 120, 0, 0.4);

  _clusterSizePTCut_multi_dR =
      new TProfile((_name + "_clusterSizeTCut_multi_dR").c_str(),
                   "_clusterSize_multi_dR", 120, 0, 0.4);

  _clusterSize_single = new TH1F((_name + "_clusterSize_single").c_str(),
                                 "_clusterSize_single", 120, 0, 120);

  _clusterSize_multi = new TH1F((_name + "_clusterSize_multi").c_str(),
                                "_clusterSize_multi", 120, 0, 120);

  _clusterSizePTCut_single =
      new TH1F((_name + "_clusterSizeTCut_single").c_str(),
               "_clusterSize_single", 120, 0, 120);

  _clusterSizePTCut_multi = new TH1F((_name + "_clusterSizeTCut_multi").c_str(),
                                     "_clusterSize_multi", 120, 0, 120);

  _hits_PU_dR = new TH1F((_name + "_hits_PU_Jet_h").c_str(), "_hits_PU_Jet_h",
                         120, 0., 0.4);

  _deltaR_PU_h =
      new TH1F((_name + "_deltaR_PU_h").c_str(), "_deltaR_PU_h", 120, 0., 0.4);
}

TrackClusterAnalysis::~TrackClusterAnalysis() {
  delete _deltaR_h;
  delete _outputTree_cluster;
  delete _nTracksPerCluster_dR;
  delete _nTracksPerClusterPTCut_dR;
  delete _nTracksPerCluster_PU_dR;
  delete _nTracksPerClusterPTCut_PU_dR;
  delete _mergedClusters_dR;
  delete _mergedClustersPTCut_dR;
  delete _mergedClusterRate_dR;
  delete _mergedClusterRatePTCut_dR;
  delete _mergedClusterRate_PU_dR;
  delete _mergedClusterRatePTCut_PU_dR;
  delete _hits_dR;
  delete _clusterSize_single_dR;
  delete _clusterSize_multi_dR;
  delete _clusterSizePTCut_single_dR;
  delete _clusterSizePTCut_multi_dR;
  delete _clusterSize_single;
  delete _clusterSize_multi;
  delete _clusterSizePTCut_single;
  delete _clusterSizePTCut_multi;
  delete _deltaR_PU_h;
  delete _hits_PU_dR;
}

void TrackClusterAnalysis::fill_cluster(float deltaR_cluster, float eta,
                                        float eta_jet, float phi, float phi_jet,
                                        float nTracksPerCluster,
                                        float nTracksPerClusterPTCut,
                                        unsigned short nCellsPerCluster,
                                        bool isPU) {
  // fill cluster tree
  _deltaR_cluster = deltaR_cluster;
  _eta = eta;
  _eta_jet = eta_jet;
  _phi = phi;
  _phi_jet = phi_jet;
  _nTracksPerCluster = nTracksPerCluster;
  _nTracksPerClusterPTCut = nTracksPerClusterPTCut;
  _mergedCluster = false;
  _mergedClusterPTCut = false;
  _isPU = isPU;
  // fill histograms
  if (!isPU) {
    _nTracksPerCluster_dR->Fill(deltaR_cluster, nTracksPerCluster);
    _nTracksPerClusterPTCut_dR->Fill(deltaR_cluster, nTracksPerClusterPTCut);
  }
  _nTracksPerCluster_PU_dR->Fill(deltaR_cluster, nTracksPerCluster);
  _nTracksPerClusterPTCut_PU_dR->Fill(deltaR_cluster, nTracksPerClusterPTCut);

  if (nTracksPerCluster > 1) {
    _mergedClusters_dR->Fill(deltaR_cluster);
    if (!isPU) _mergedClusterRate_dR->Fill(deltaR_cluster, 1.);
    _mergedClusterRate_PU_dR->Fill(deltaR_cluster, 1.);
    _clusterSize_multi_dR->Fill(deltaR_cluster, nCellsPerCluster);
    _clusterSize_multi->Fill(nCellsPerCluster);
    _mergedCluster = true;
    _clusterSizeMulti = nCellsPerCluster;
  } else {
    if (!isPU) _mergedClusterRate_dR->Fill(deltaR_cluster, 0.);
    _mergedClusterRate_PU_dR->Fill(deltaR_cluster, 0.);
    _clusterSize_single_dR->Fill(deltaR_cluster, nCellsPerCluster);
    _clusterSize_single->Fill(nCellsPerCluster);
    _clusterSizeSingle = nCellsPerCluster;
  }
  if (nTracksPerClusterPTCut > 1) {
    _mergedClustersPTCut_dR->Fill(deltaR_cluster);
    if (!isPU) _mergedClusterRatePTCut_dR->Fill(deltaR_cluster, 1.);
    _mergedClusterRatePTCut_PU_dR->Fill(deltaR_cluster, 1.);
    _clusterSizePTCut_multi_dR->Fill(deltaR_cluster, nCellsPerCluster);
    _clusterSizePTCut_multi->Fill(nCellsPerCluster);
    _mergedClusterPTCut = true;
    _clusterSizeMultiPTCut = nCellsPerCluster;
  } else {
    if (!isPU) _mergedClusterRatePTCut_dR->Fill(deltaR_cluster, 0.);
    _mergedClusterRatePTCut_PU_dR->Fill(deltaR_cluster, 0.);
    _clusterSizePTCut_single_dR->Fill(deltaR_cluster, nCellsPerCluster);
    _clusterSizePTCut_single->Fill(nCellsPerCluster);
    _clusterSizeSinglePTCut = nCellsPerCluster;
  }
  if (isPU) {
    _deltaR_PU_h->Fill(deltaR_cluster);
  }

  _deltaR_h->Fill(deltaR_cluster);
  for (unsigned short i = 0; i < nCellsPerCluster; i++) {
    _hits_dR->Fill(deltaR_cluster);
    if (isPU) {
      _hits_PU_dR->Fill(deltaR_cluster);
    }
  }
  _outputTree_cluster->Fill();
}
/*
void TrackClusterAnalysis::fill_module(float deltaR_module, float nClusters,
                                       float occupancy,
                                       unsigned long long moduleID) {
  _deltaR_module = deltaR_module;
  _nClusters = nClusters;
  _occupancy = occupancy;
  _moduleID = moduleID;
  _occ_dR->Fill(deltaR_module, occupancy);
  _occ_dR_h->Fill(deltaR_module, occupancy);

  _outputTree_module->Fill();
}

void TrackClusterAnalysis::fill_allClusters(float eta, float phi,
                                            unsigned long long moduleID) {
  _etaAll = eta;
  _phiAll = phi;
  _moduleIDAll = moduleID;
  _outputTree_clusterAll->Fill();
}
*/
void TrackClusterAnalysis::write() {
  _nTracksPerCluster_dR->Write();
  _nTracksPerCluster_PU_dR->Write();
  _nTracksPerClusterPTCut_dR->Write();
  _nTracksPerClusterPTCut_PU_dR->Write();

  _mergedClusters_dR->Write();
  _mergedClustersPTCut_dR->Write();

  _mergedClusterRate_dR->Write();
  _mergedClusterRatePTCut_dR->Write();
  _mergedClusterRate_PU_dR->Write();
  _mergedClusterRatePTCut_PU_dR->Write();

  // _outputTree_module->Write();
  _outputTree_cluster->Write();
  // _outputTree_clusterAll->Write();
  _deltaR_h->Write();
  _hits_dR->Write();

  _deltaR_PU_h->Write();
  _hits_PU_dR->Write();

  _clusterSize_single_dR->Write();
  _clusterSize_multi_dR->Write();
  _clusterSizePTCut_single_dR->Write();
  _clusterSizePTCut_multi_dR->Write();
  _clusterSize_single->Write();
  _clusterSize_multi->Write();
  _clusterSizePTCut_single->Write();
  _clusterSizePTCut_multi->Write();
}

void TrackClusterAnalysis::normalize(float norm) {
  float scalor = 1. / norm;
  if (norm == 0.) {
    throw std::runtime_error(
        "TrackClusterAnalysis::normalizeMergedClusters:Trying to divide by "
        "zero!");
  }
  _mergedClusters_dR->Scale(scalor);
  _mergedClustersPTCut_dR->Scale(scalor);
  _deltaR_h->Scale(scalor);
  _hits_dR->Scale(scalor);
  _clusterSize_single->Scale(scalor);
  _clusterSize_multi->Scale(scalor);
  _clusterSizePTCut_single->Scale(scalor);
  _clusterSizePTCut_multi->Scale(scalor);
  _deltaR_PU_h->Scale(scalor);
  _hits_PU_dR->Scale(scalor);
}
