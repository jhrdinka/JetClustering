#ifndef TRACKANALYSIS_H
#define TRACKANALYSIS_H

#include <TDirectory.h>
#include <TString.h>
#include <TTree.h>
#include <ios>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include "GenParticle.hh"
#include "Jet.hh"
#include "TH1F.h"
#include "TProfile.h"

using namespace std;

// Analyse tracks: e.g. number of shared hits along track

class TrackAnalysis {
 private:
  /// Name of the output tree
  TString _treeName;
  /// minimum number if hits on different layers to form a track
  size_t _hitRequirement = 8;  // see flavour tagging note Estel
  /// The output tree
  TTree* _outputTree = nullptr;
  /// number of shared clusters per track (track requirement)
  TH1F* _sharedClusters_h = nullptr;
  /// number of clusters per track (track requirement)
  TH1F* _nclusters_h = nullptr;
  /// number of clusters on different layers (track requirement)
  TH1F* _nLayersHit_h = nullptr;
  /// all generated particles
  TH1F* _genPart_dR_h = nullptr;
  /// tracks (at least _hitRequirement number of hits)
  TH1F* _tracks_dR_h = nullptr;
  /// number of tracks (at least _hitRequirement number of hits) with 0 shared
  /// cluster
  TH1F* _tracks_0shared_dR_h = nullptr;
  /// number of tracks (at least _hitRequirement number of hits) with 1 shared
  /// cluster
  TH1F* _tracks_1shared_dR_h = nullptr;
  /// number of tracks (at least _hitRequirement number of hits) with 2 shared
  /// cluster
  TH1F* _tracks_2shared_dR_h = nullptr;
  /// number of tracks (at least _hitRequirement number of hits) with 3 shared
  /// cluster
  TH1F* _tracks_3shared_dR_h = nullptr;
  /// number of tracks (at least _hitRequirement number of hits) with 4 shared
  /// cluster
  TH1F* _tracks_4shared_dR_h = nullptr;
  /// number of tracks (at least _hitRequirement number of hits) with more than
  /// shared cluster
  TH1F* _tracks_manyshared_dR_h = nullptr;

  // With PU
  /// number of shared clusters per track (track requirement)
  TH1F* _sharedClusters_PU_h = nullptr;
  /// number of clusters per track (track requirement)
  TH1F* _nclusters_PU_h = nullptr;
  /// number of clusters on different layers (track requirement)
  TH1F* _nLayersHit_PU_h = nullptr;
  /// all generated particles
  TH1F* _genPart_dR_PU_h = nullptr;
  /// tracks (at least _PU_hitRequirement number of hits)
  TH1F* _tracks_dR_PU_h = nullptr;
  /// number of tracks (at least _PU_hitRequirement number of hits) with 0
  /// shared
  /// cluster
  TH1F* _tracks_0shared_dR_PU_h = nullptr;
  /// number of tracks (at least _PU_hitRequirement number of hits) with 1
  /// shared
  /// cluster
  TH1F* _tracks_1shared_dR_PU_h = nullptr;
  /// number of tracks (at least _PU_hitRequirement number of hits) with 2
  /// shared
  /// cluster
  TH1F* _tracks_2shared_dR_PU_h = nullptr;
  /// number of tracks (at least _PU_hitRequirement number of hits) with 3
  /// shared
  /// cluster
  TH1F* _tracks_3shared_dR_PU_h = nullptr;
  /// number of tracks (at least _PU_hitRequirement number of hits) with 4
  /// shared
  /// cluster
  TH1F* _tracks_4shared_dR_PU_h = nullptr;

  /// number of tracks (at least _PU_hitRequirement number of hits) with more
  /// than 4
  /// shared
  /// cluster
  TH1F* _tracks_manyshared_dR_PU_h = nullptr;

 public:
  // constructors
  TrackAnalysis(const TString treeName = "TrackAnalysis",
                size_t hitRequirement = 8);
  void fill(std::unordered_map<unsigned, GenParticle>& tracks, float coneR);
  void write();
  void normalize(float norm = 1);
};

#endif  // TRACKANALYSIS_H
