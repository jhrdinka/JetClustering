#include "Module.hh"
#include <iostream>
#include <stdexcept>

Module::Module(unsigned long long moduleID, const TVector3& position,
               unsigned nChannelsOn, unsigned nChannels,
               const std::vector<float>& clusterPosX,
               const std::vector<float>& clusterPosY,
               const std::vector<float>& clusterPosZ,
               const std::vector<std::vector<unsigned>>& trackIDsPerCluster,
               const std::vector<unsigned short>& nCellsPerCluster)
    : _moduleID(moduleID),
      _position(position),
      _nChannelsOn(nChannelsOn),
      _nChannels(nChannels),
      _clusterPosX(clusterPosX),
      _clusterPosY(clusterPosY),
      _clusterPosZ(clusterPosZ),
      _trackIDsPerCluster(trackIDsPerCluster),
      _nCellsPerCluster(nCellsPerCluster) {
  // todo check if sizes of all vectors are the same
  size_t nClusters = clusterPosX.size();
  if ((nClusters != clusterPosY.size()) || (nClusters != clusterPosZ.size()) ||
      (nClusters != trackIDsPerCluster.size()) ||
      (nClusters != nCellsPerCluster.size())) {
    throw std::runtime_error(
        "In Module constructor: Number of clusters not consistent!");
  }
}

Module::Module(const Module& module) {
  _moduleID = module._moduleID;
  _position = module._position;
  _nChannelsOn = module._nChannelsOn;
  _nChannels = module._nChannels;
  _clusterPosX = module._clusterPosX;
  _clusterPosY = module._clusterPosY;
  _clusterPosZ = module._clusterPosZ;
  _trackIDsPerCluster = module._trackIDsPerCluster;
  _nCellsPerCluster = module._nCellsPerCluster;
}
