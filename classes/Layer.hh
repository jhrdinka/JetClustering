#ifndef LAYER_H
#define LAYER_H

#include "TH1F.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TrackCluster.hh"

class Layer {
 public:
  /// default constructor
  Layer();
  /// constructor
  Layer(const std::vector<float>& clusterPosX,
        const std::vector<float>& clusterPosY,
        const std::vector<float>& clusterPosZ);
  /// copy constructor
  Layer(const Layer& Layer);

  /// desctructor
  ~Layer() = default;

  // initialize
  void initialize(size_t layerNum);

  /// add clusters to this layer
  void addClusters(
      std::vector<float> clusterPosX, std::vector<float> clusterPosY,
      std::vector<float> clusterPosZ, float deltaR, float occupancy,
      const std::vector<std::vector<unsigned>>& trackIDsPerCluster);

  /// access cluster position in x
  const std::vector<float>& clusterPositionX() const;

  /// access cluster position in y
  const std::vector<float>& clusterPositionY() const;

  /// access cluster position in z
  const std::vector<float>& clusterPositionZ() const;

  /// number of clusters on this Layer
  size_t nClusters() const;

  void print() const;

  // normalize
  void normalize(float norm);

 private:
  /// the cluster position in x
  std::vector<float> _clusterPosX;
  /// the cluster position in y
  std::vector<float> _clusterPosY;
  /// the cluster position in z
  std::vector<float> _clusterPosZ;
  // number of clusters over dR
  TH1F* _deltaR_h;
  // occupancy over dR
  TProfile* _occ_dR;
  // number of tracks over dR
  TProfile* _nTracksPerCluster_dR;
};

inline const std::vector<float>& Layer::clusterPositionX() const {
  return (_clusterPosX);
}

inline const std::vector<float>& Layer::clusterPositionY() const {
  return (_clusterPosY);
}

inline const std::vector<float>& Layer::clusterPositionZ() const {
  return (_clusterPosZ);
}

inline size_t Layer::nClusters() const { return _clusterPosX.size(); }

#endif  // LAYER_H
