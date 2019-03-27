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
  /// number of tracks per cluster for differnt cuts
  float _nTracksPerCluster;       // no cut
  float _nTracksPerClusterPTCut;  // cut on pt
  /// flag indicating if it is a merged cluster for differn cuts
  bool _mergedCluster;       // no cut
  bool _mergedClusterPTCut;  // cut on pt
  /// cluster size single particle cluster
  int _clusterSizeSingle;
  /// cluster size multi particle cluster
  int _clusterSizeMulti;
  /// cluster size single particle cluster
  int _clusterSizeSinglePTCut;
  /// cluster size multi particle cluster
  int _clusterSizeMultiPTCut;
  /// if cluster is from PU
  bool _isPU;

  /// The output trees
  TTree* _outputTree_cluster = nullptr;
  // number of clusters over dR within a jet
  TH1F* _deltaR_h = nullptr;
  // number of tracks over dR
  TProfile* _nTracksPerCluster_dR = nullptr;
  // number of tracks over dR
  TProfile* _nTracksPerCluster_PU_dR = nullptr;
  // number of tracks over dR
  TProfile* _nTracksPerClusterPTCut_dR = nullptr;
  // number of tracks over dR
  TProfile* _nTracksPerClusterPTCut_PU_dR = nullptr;
  // mergedClusters - number of merged clusters within this jet
  TH1F* _mergedClusters_dR = nullptr;
  // mergedClusters
  TH1F* _mergedClustersPTCut_dR = nullptr;
  // mergedClusterRate - number of merged clusters within this jet divided by
  // the total number of clusters of this jet
  TProfile* _mergedClusterRate_dR = nullptr;
  // mergedClusterRate - number of merged clusters within this jet divided by
  // the total number of clusters of this jet
  TProfile* _mergedClusterRate_PU_dR = nullptr;
  // mergedClusterRate
  TProfile* _mergedClusterRatePTCut_dR = nullptr;
  // mergedClusterRate
  TProfile* _mergedClusterRatePTCut_PU_dR = nullptr;
  // number of hits (activated cells) over dR within the jet cone
  TH1F* _hits_dR = nullptr;
  // cluster size of single particle clusters over dR
  TProfile* _clusterSize_single_dR = nullptr;
  // cluster size of multi-particle clusters over dR
  TProfile* _clusterSize_multi_dR = nullptr;
  // cluster size of single particle clusters over dR
  TProfile* _clusterSizePTCut_single_dR = nullptr;
  // cluster size of multi-particle clusters over dR
  TProfile* _clusterSizePTCut_multi_dR = nullptr;
  // cluster size of single particle clusters per jet
  TH1F* _clusterSize_single = nullptr;
  // cluster size of multi-particle clusters per jet
  TH1F* _clusterSize_multi = nullptr;
  // cluster size of single particle clusters per jet
  TH1F* _clusterSizePTCut_single = nullptr;
  // cluster size of multi-particle clusters per jet
  TH1F* _clusterSizePTCut_multi = nullptr;
  // number of clusters of pile up only over dR within a jet
  TH1F* _deltaR_PU_h = nullptr;
  // number of hits (activated cells) of pile up only over dR within the jet
  // cone
  TH1F* _hits_PU_dR = nullptr;

 public:
  // constructors
  TrackClusterAnalysis(const std::string& analysisName);
  ~TrackClusterAnalysis();
  void fill_cluster(float deltaR_cluster, float eta, float eta_jet, float phi,
                    float phi_jet, float nTracksPerCluster,
                    float nTracksPerClusterPTCut,
                    unsigned short nCellsPerCluster, bool isPU);
  /* void fill_module(float deltaR_module, float nClusters, float occupancy,
                    unsigned long long moduleID);

   void fill_allClusters(float eta, float phi, unsigned long long moduleID);
 */
  void write();

  void normalize(float norm);
};

#endif  // TRACKCLUSTERANALYSIS_H
