#ifndef TRACKCLUSTERANALYSIS_H
#define TRACKCLUSTERANALYSIS_H

#include <TDirectory.h>
#include <TString.h>
#include <TTree.h>
#include <ios>
#include <stdexcept>
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

using namespace std;

class TrackClusterAnalysis {
 private:
  std::string _name;
  /// difference in r (angular region)
  float _deltaR_module;
  /// difference in r (angular region)
  float _deltaR_cluster;
  /// number of clusters
  float _nClusters;
  /// channel occupancy
  float _occupancy;
  /// eta distribution
  float _eta;
  /// eta distribution
  float _eta_jet;
  /// eta distribution
  float _phi;
  /// eta distribution
  float _phi_jet;
  /// moduleIDs for the trees
  ULong64_t _moduleID;
  ULong64_t _moduleID_cluster;
  ULong64_t _moduleIDAll;
  /// number of tracks per cluster for differnt cuts
  float _nTracksPerCluster;          // no cut
  float _nTracksPerClusterPTCut;     // cut on pt
  float _nTracksPerClusterLayerCut;  // cut if it found on a second layer
  /// flag indicating if it is a merged cluster for differn cuts
  bool _mergedCluster;          // no cut
  bool _mergedClusterPTCut;     // cut on pt
  bool _mergedClusterLayerCut;  // cut if it found on a second layer
  /// The output trees
  TTree* _outputTree_cluster;
  TTree* _outputTree_module;
  TTree* _outputTree_clusterAll;
  // all clusters - also the ones outside the jets
  float _etaAll;
  float _phiAll;
  // number of clusters over dR within a jet
  TH1F* _deltaR_h;
  // occupancy per module over dR
  TProfile* _occ_dR;
  // 2D histogram of occupancy per module over dR
  TH2F* _occ_dR_h;
  // number of tracks over dR
  TProfile* _nTracksPerCluster_dR;
  // number of tracks over dR
  TProfile* _nTracksPerClusterPTCut_dR;
  // number of tracks over dR
  TProfile* _nTracksPerClusterLayerCut_dR;
  // mergedClusters - number of merged clusters within this jet
  TH1F* _mergedClusters_dR;
  // mergedClusters
  TH1F* _mergedClustersPTCut_dR;
  // mergedClusters
  TH1F* _mergedClustersLayerCut_dR;
  // mergedClusterRate - number of merged clusters within this jet divided by
  // the total number of clusters of this jet
  TProfile* _mergedClusterRate_dR;
  // mergedClusterRate
  TProfile* _mergedClusterRatePTCut_dR;
  // mergedClusterRate
  TProfile* _mergedClusterRateLayerCut_dR;
  // number of hits (activated cells) over dR within the jet cone
  TH1F* _hits_dR;

 public:
  // constructors
  TrackClusterAnalysis(const std::string& analysisName);
  ~TrackClusterAnalysis();
  void fill_cluster(float deltaR_cluster, float eta, float eta_jet, float phi,
                    float phi_jet, float nTracksPerCluster,
                    float nTracksPerClusterPTCut,
                    float nTracksPerClusterLayerCut,
                    unsigned long long moduleID,
                    unsigned short nCellsPerCluster);
  void fill_module(float deltaR_module, float nClusters, float occupancy,
                   unsigned long long moduleID);

  void fill_allClusters(float eta, float phi, unsigned long long moduleID);

  void write();

  void normalize(float norm);
};

#endif  // TRACKCLUSTERANALYSIS_H
