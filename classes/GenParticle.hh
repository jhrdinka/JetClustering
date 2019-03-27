#ifndef GENPARTICLE_H
#define GENPARTICLE_H

#include <set>
#include "TLorentzVector.h"

class GenParticle {
 private:
  TLorentzVector _mom;
  int _pdgid;
  int _status;
  TVector3 _vertex;
  float _charge;
  float _deltaR;
  bool _isPU = false;

 public:
  // constructor
  GenParticle(const TLorentzVector& p4, int pdgid, unsigned status,
              const TVector3& vertex, float charge, float deltaR, bool isPU);
  GenParticle(const GenParticle& g);

  bool isClusterable() const;

  float eta() const;
  float phi() const;
  float pt() const;
  float energy() const;
  float px() const;
  float py() const;
  float pz() const;
  float mass() const;
  float charge() const;
  TLorentzVector p4() const;
  int pdgid() const;
  int status() const;
  TVector3 vertex() const;
  // deltaR to closest jet
  float deltaR() const;
  bool isPU() const;

  // all layers hit along the track
  std::set<unsigned long long> layerIDs;
  // total number of clusters
  size_t nClusters = 0;
  // the other trackIDs belonging to each cluster
  std::vector<std::set<unsigned>> shared_trackIDs;
};

#endif
