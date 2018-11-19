#include "Layer.hh"
#include <iostream>

Layer::Layer() {}

void Layer::initialize(size_t layerNum) {
  _deltaR_h = new TH1F(("deltaR_cluster" + std::to_string(layerNum)).c_str(),
                       "deltaR_cluster", 100., 0., 0.4);
  _occ_dR = new TProfile(("_occ_dR" + std::to_string(layerNum)).c_str(),
                         "occupancies over dR", 100, 0., 0.4);
  _nTracksPerCluster_dR =
      new TProfile(("_nTracksPerCluster_dR" + std::to_string(layerNum)).c_str(),
                   "_nTracksPerCluster_dR", 100, 0., 0.4);
}

Layer::Layer(const std::vector<float>& clusterPosX,
             const std::vector<float>& clusterPosY,
             const std::vector<float>& clusterPosZ)
    : _clusterPosX(clusterPosX),
      _clusterPosY(clusterPosY),
      _clusterPosZ(clusterPosZ) {
  /* _clusters_eta_h =
       new TH1F("_clusters_eta_h", "clusters in #eta", 100., -6., 6.);
   _clusters_phi_h =
       new TH1F("_clusters_phi_h", "clusters in #phi", 100., -M_PI, M_PI);*/
}

Layer::Layer(const Layer& Layer) {
  _clusterPosX = Layer._clusterPosX;
  _clusterPosY = Layer._clusterPosY;
  _clusterPosZ = Layer._clusterPosZ;
}

void Layer::print() const {
  /* for (size_t i = 0; i < _clusterPosX.size(); i++) {
     TVector3 clusterPos(_clusterPosX.at(i), _clusterPosY.at(i),
                         _clusterPosZ.at(i));

     _clusters_eta_h->Fill(clusterPos.Eta());
     _clusters_phi_h->Fill(clusterPos.Phi());
   }
   _clusters_eta_h->Write();
   _clusters_phi_h->Write();
 */
  _deltaR_h->Write();
  _occ_dR->Write();
  _nTracksPerCluster_dR->Write();
}

void Layer::addClusters(
    std::vector<float> clusterPosX, std::vector<float> clusterPosY,
    std::vector<float> clusterPosZ, float deltaR, float occupancy,
    const std::vector<std::vector<unsigned>>& trackIDsPerCluster) {
  _clusterPosX.insert(_clusterPosX.end(), clusterPosX.begin(),
                      clusterPosX.end());
  _clusterPosY.insert(_clusterPosY.end(), clusterPosY.begin(),
                      clusterPosY.end());
  _clusterPosZ.insert(_clusterPosZ.end(), clusterPosZ.begin(),
                      clusterPosZ.end());

  _occ_dR->Fill(deltaR, occupancy);
  _deltaR_h->Fill(deltaR);
  for (auto& trackIDs : trackIDsPerCluster) {
    _nTracksPerCluster_dR->Fill(deltaR, trackIDs.size());
  }
}

void Layer::normalize(float norm) { _deltaR_h->Scale(norm); }
/*
void Layer::fillDistances(TProfile* dz_dR_p, TProfile* dphi_dR_p) const {
  for (size_t i = 0; i < (_clusterPosX.size() - 1); i++) {
    TVector3 clusterPos1(_clusterPosX.at(i), _clusterPosY.at(i),
                         _clusterPosZ.at(i));
    for (size_t j = i + 1; j < _clusterPosX.size(); j++) {
      TVector3 clusterPos2(_clusterPosX.at(j), _clusterPosY.at(j),
                           _clusterPosZ.at(j));
      float dz = std::fabs(clusterPos1.Z() - clusterPos2.Z());
      float dphi = std::fabs(clusterPos1.Phi() - clusterPos2.Phi());

      dz_dR_p->Fill(deltaR, dz);
      dphi_dR_p->Fill(deltaR, dphi);
    }
  }
}*/
