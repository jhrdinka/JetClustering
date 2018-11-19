#ifndef MODULE_H
#define MODULE_H

#include <iostream>
#include "TVector3.h"
#include "TrackCluster.hh"

class Module {
 public:
  /// constructor
  Module(long long int moduleID, const TVector3& position, int nChannelsOn,
         int nChannels, const std::vector<float>& clusterPosX,
         const std::vector<float>& clusterPosY,
         const std::vector<float>& clusterPosZ,
         const std::vector<std::vector<unsigned>>& trackIDsPerCluster,
         const std::vector<short int>& nCellsPerCluster);
  /// copy constructor
  Module(const Module& module);

  /// desctructor
  ~Module() = default;

  /// access module ID
  long long int moduleID() const;

  /// access position
  TVector3 position() const;

  /// access number of activated channels
  int nChannelsOn() const;

  /// access number of channels
  int nChannels() const;

  /// get occupancy for this module
  float occupancy() const;

  /// get the merged cluster rate for this module
  //  float mergedClusterRate() const;

  /// access cluster position in x
  const std::vector<float>& clusterPositionX() const;

  /// access cluster position in y
  const std::vector<float>& clusterPositionY() const;

  /// access cluster position in z
  const std::vector<float>& clusterPositionZ() const;

  /// the total number of clusters(also with no track assigned, i.e. only
  /// produced by secondaries not leaving the detector material)
  size_t nClusters() const;

  /// the trackIDs per cluster (row: cluster,column: trackIDS of that cluster)
  const std::vector<std::vector<unsigned>>& trackIDsPerCluster() const;

  /// number of cells per cluster
  const std::vector<short int>& nCellsPerCluster() const;

 private:
  /// the module ID
  long long int _moduleID;
  /// the module position
  TVector3 _position;
  /// the number of channels turned on
  int _nChannelsOn;
  /// the number of channels of this module
  int _nChannels;
  /// information on  cluster level
  /// the cluster position in x
  std::vector<float> _clusterPosX;
  /// the cluster position in y
  std::vector<float> _clusterPosY;
  /// the cluster position in z
  std::vector<float> _clusterPosZ;
  /// number of tracks per cluster
  std::vector<std::vector<unsigned>> _trackIDsPerCluster;
  /// number of cells per cluster
  std::vector<short int> _nCellsPerCluster;
};

inline long long int Module::moduleID() const { return _moduleID; }

inline TVector3 Module::position() const { return _position; }

inline int Module::nChannelsOn() const { return _nChannelsOn; }

inline int Module::nChannels() const { return _nChannels; }

inline float Module::occupancy() const {
  return (float(_nChannelsOn) / float(_nChannels));
}

inline const std::vector<float>& Module::clusterPositionX() const {
  return (_clusterPosX);
}

inline const std::vector<float>& Module::clusterPositionY() const {
  return (_clusterPosY);
}

inline const std::vector<float>& Module::clusterPositionZ() const {
  return (_clusterPosZ);
}

inline size_t Module::nClusters() const { return _clusterPosX.size(); }

inline const std::vector<std::vector<unsigned>>& Module::trackIDsPerCluster()
    const {
  return _trackIDsPerCluster;
}

inline const std::vector<short int>& Module::nCellsPerCluster() const {
  return (_nCellsPerCluster);
}

#endif  //  MODULE_H