#include "TrackClusterAnalysis.hh"
#include <iostream>

TrackClusterAnalysis::TrackClusterAnalysis(const std::string& analysisName)
    : _name(analysisName) {
  _outputTree_cluster =
      new TTree((_name + "_TrackClusterAnalysis_cluster").c_str(),
                "TrackClusterAnalysis_cluster");
  _outputTree_module =
      new TTree((_name + "_TrackClusterAnalysis_module").c_str(),
                "TrackClusterAnalysis_module");
  _outputTree_clusterAll =
      new TTree((_name + "_TrackClusterAnalysis_clusterAll").c_str(),
                "TrackClusterAnalysis_clusterAll");
  if (!_outputTree_clusterAll) throw std::bad_alloc();
  if (!_outputTree_module) throw std::bad_alloc();
  if (!_outputTree_clusterAll) throw std::bad_alloc();

  // tree module - to possibly redo histograms
  _outputTree_module->Branch("deltaR_module", &_deltaR_module);
  _outputTree_module->Branch("nClusters", &_nClusters);
  _outputTree_module->Branch("occupancy", &_occupancy);
  _outputTree_module->Branch("moduleID", &_moduleID);
  _outputTree_module->SetDirectory(0);
  // tree cluster - to possibly redo histograms
  _outputTree_cluster->Branch("deltaR_cluster", &_deltaR_cluster);
  _outputTree_cluster->Branch("eta", &_eta);
  _outputTree_cluster->Branch("eta_jet", &_eta_jet);
  _outputTree_cluster->Branch("phi", &_phi);
  _outputTree_cluster->Branch("phi_jet", &_phi_jet);
  _outputTree_cluster->Branch("nTracksPerCluster", &_nTracksPerCluster);
  _outputTree_cluster->Branch("nTracksPerClusterPTCut",
                              &_nTracksPerClusterPTCut);
  _outputTree_cluster->Branch("nTracksPerClusterLayerCut",
                              &_nTracksPerClusterLayerCut);
  _outputTree_cluster->Branch("mergedClusters", &_mergedCluster);
  _outputTree_cluster->Branch("mergedClustersPTCut", &_mergedClusterPTCut);
  _outputTree_cluster->Branch("mergedClustersLayerCut",
                              &_mergedClusterLayerCut);
  _outputTree_cluster->Branch("moduleID", &_moduleID_cluster);
  _outputTree_cluster->SetDirectory(0);

  _outputTree_clusterAll->Branch("etaAll", &_etaAll);
  _outputTree_clusterAll->Branch("phiAll", &_phiAll);
  _outputTree_clusterAll->Branch("moduleID", &_moduleIDAll);
  _outputTree_clusterAll->SetDirectory(0);

  // histograms and profiles
  _deltaR_h = new TH1F((_name + "_deltaR_cluster").c_str(), "deltaR_cluster",
                       120., 0., 0.4);
  _occ_dR = new TProfile((_name + "_occ_dR").c_str(), "occupancies over dR",
                         120, 0., 0.4);
  _occ_dR_h = new TH2F((_name + "_occ_dR_h").c_str(), "occupancies over dR",
                       120, 0., 0.4, 120., 0., 0.1);
  _nTracksPerCluster_dR =
      new TProfile((_name + "_nTracksPerCluster_dR").c_str(),
                   "_nTracksPerCluster_dR", 120, 0., 0.4);
  _nTracksPerClusterPTCut_dR =
      new TProfile((_name + "_nTracksPerClusterPTCut_dR").c_str(),
                   "_nTracksPerCluster_dR", 120, 0., 0.4);
  _nTracksPerClusterLayerCut_dR =
      new TProfile((_name + "_nTracksPerClusterLayerCut_dR").c_str(),
                   "_nTracksPerCluster_dR", 120, 0., 0.4);

  _mergedClusters_dR = new TH1F((_name + "_mergedClusters_dR").c_str(),
                                "_mergedClusters_dR", 120, 0., 0.4);
  _mergedClustersPTCut_dR =
      new TH1F((_name + "_mergedClustersPTCut_dR").c_str(),
               "_mergedClusters_dR", 120, 0., 0.4);
  _mergedClustersLayerCut_dR =
      new TH1F((_name + "_mergedClustersLayerCut_dR").c_str(),
               "_mergedClusters_dR", 120, 0., 0.4);

  _mergedClusterRate_dR =
      new TProfile((_name + "_mergedClusterRate_dR").c_str(),
                   "_mergedClusterrate_dR", 120, 0., 0.4);
  _mergedClusterRatePTCut_dR =
      new TProfile((_name + "_mergedClusterRatePTCut_dR").c_str(),
                   "_mergedClusterRate_dR", 120, 0., 0.4);
  _mergedClusterRateLayerCut_dR =
      new TProfile((_name + "_mergedClusterRateLayerCut_dR").c_str(),
                   "_mergedClusterRate_dR", 120, 0., 0.4);

  _hits_dR = new TH1F((_name + "_hits_dR").c_str(), "_hits_dR", 120, 0., 0.4);
}

TrackClusterAnalysis::~TrackClusterAnalysis() {
  delete _deltaR_h;
  delete _occ_dR;
  delete _occ_dR_h;
  delete _outputTree_cluster;
  delete _outputTree_module;
  delete _outputTree_clusterAll;
  delete _nTracksPerCluster_dR;
  delete _nTracksPerClusterPTCut_dR;
  delete _nTracksPerClusterLayerCut_dR;
  delete _mergedClusters_dR;
  delete _mergedClustersPTCut_dR;
  delete _mergedClustersLayerCut_dR;
  delete _mergedClusterRate_dR;
  delete _mergedClusterRatePTCut_dR;
  delete _mergedClusterRateLayerCut_dR;
  delete _hits_dR;
}

void TrackClusterAnalysis::fill_cluster(float deltaR_cluster, float eta,
                                        float eta_jet, float phi, float phi_jet,
                                        float nTracksPerCluster,
                                        float nTracksPerClusterPTCut,
                                        float nTracksPerClusterLayerCut,
                                        unsigned long long moduleID,
                                        unsigned short nCellsPerCluster) {
  // fill cluster tree
  _deltaR_cluster = deltaR_cluster;
  _eta = eta;
  _eta_jet = eta_jet;
  _phi = phi;
  _phi_jet = phi_jet;
  _nTracksPerCluster = nTracksPerCluster;
  _nTracksPerClusterPTCut = nTracksPerClusterPTCut;
  _nTracksPerClusterLayerCut = nTracksPerClusterLayerCut;
  _moduleID_cluster = moduleID;
  _mergedCluster = false;
  _mergedClusterPTCut = false;
  _mergedClusterLayerCut = false;

  // fill histograms
  _nTracksPerCluster_dR->Fill(deltaR_cluster, nTracksPerCluster);
  _nTracksPerClusterPTCut_dR->Fill(deltaR_cluster, nTracksPerClusterPTCut);
  _nTracksPerClusterLayerCut_dR->Fill(deltaR_cluster,
                                      nTracksPerClusterLayerCut);

  if (nTracksPerCluster > 1) {
    _mergedClusters_dR->Fill(deltaR_cluster);
    _mergedClusterRate_dR->Fill(deltaR_cluster, 1.);
    _mergedCluster = true;
  } else {
    _mergedClusterRate_dR->Fill(deltaR_cluster, 0.);
  }
  if (nTracksPerClusterPTCut > 1) {
    _mergedClustersPTCut_dR->Fill(deltaR_cluster);
    _mergedClusterRatePTCut_dR->Fill(deltaR_cluster, 1.);
    _mergedClusterPTCut = true;
  } else {
    _mergedClusterRatePTCut_dR->Fill(deltaR_cluster, 0.);
  }
  if (nTracksPerClusterLayerCut > 1) {
    _mergedClustersLayerCut_dR->Fill(deltaR_cluster);
    _mergedClusterRateLayerCut_dR->Fill(deltaR_cluster, 1.);
    _mergedClusterLayerCut = true;
  } else {
    _mergedClusterRateLayerCut_dR->Fill(deltaR_cluster, 0.);
  }

  _deltaR_h->Fill(deltaR_cluster);

  for (unsigned short i = 0; i < nCellsPerCluster; i++) {
    _hits_dR->Fill(deltaR_cluster);
  }
  _outputTree_cluster->Fill();
}

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

void TrackClusterAnalysis::write() {
  _occ_dR->Write();
  _nTracksPerCluster_dR->Write();
  _nTracksPerClusterPTCut_dR->Write();
  _nTracksPerClusterLayerCut_dR->Write();

  _mergedClusters_dR->Write();
  _mergedClustersPTCut_dR->Write();
  _mergedClustersLayerCut_dR->Write();

  _mergedClusterRate_dR->Write();
  _mergedClusterRatePTCut_dR->Write();
  _mergedClusterRateLayerCut_dR->Write();

  _occ_dR_h->Write();

  _outputTree_module->Write();
  _outputTree_cluster->Write();
  _outputTree_clusterAll->Write();
  _deltaR_h->Write();
  _hits_dR->Write();
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
  _mergedClustersLayerCut_dR->Scale(scalor);
  _deltaR_h->Scale(scalor);
  _hits_dR->Scale(scalor);
}
