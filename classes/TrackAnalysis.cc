#include "TrackAnalysis.hh"
#include <iostream>

TrackAnalysis::TrackAnalysis(const TString treeName, size_t hitRequirement)
    : _treeName(treeName), _hitRequirement(hitRequirement) {
  _sharedClusters_h =
      new TH1F("track_sharedClusters_h", "number of shared clusters per track",
               50, 0., 50.);
  _nclusters_h = new TH1F("track_nClusters_h", "number of clusters per track",
                          50, 0., 50.);
  _nLayersHit_h =
      new TH1F("track_nLayersHit_h", "number of different layers hit per track",
               50, 0., 50.);

  _genPart_dR_h = new TH1F("track_genPart_dR_h",
                           "total number of generated particles", 600, 0., 0.4);
  _tracks_dR_h = new TH1F(
      "track_tracks_dR_h",
      "total number of tracks fulfilling the track requirement", 600, 0., 0.4);
  _tracks_0shared_dR_h =
      new TH1F("track_0shared_tracks_dR_h",
               "total number of tracks fulfilling the track requirement with 0 "
               "shared clusters along the track",
               600, 0., 0.4);
  _tracks_1shared_dR_h =
      new TH1F("track_1shared_tracks_dR_h",
               "total number of tracks fulfilling the track requirement with 1 "
               "shared clusters along the track",
               600, 0., 0.4);
  _tracks_2shared_dR_h =
      new TH1F("track_2shared_tracks_dR_h",
               "total number of tracks fulfilling the track requirement with 2 "
               "shared clusters along the track",
               600, 0., 0.4);
  _tracks_3shared_dR_h =
      new TH1F("track_3shared_tracks_dR_h",
               "total number of tracks fulfilling the track requirement with 3 "
               "shared clusters along the track",
               600, 0., 0.4);
  _tracks_4shared_dR_h =
      new TH1F("track_4shared_tracks_dR_h",
               "total number of tracks fulfilling the track requirement with 4 "
               "shared clusters along the track",
               600, 0., 0.4);

  _tracks_manyshared_dR_h =
      new TH1F("track_manyshared_tracks_dR_h",
               "total number of tracks fulfilling the track requirement with 4 "
               "shared clusters along the track",
               600, 0., 0.4);

  // with PU
  _sharedClusters_PU_h =
      new TH1F("track_sharedClusters_PU_h",
               "number of shared clusters per track", 50, 0., 50.);
  _nclusters_PU_h = new TH1F("track_nClusters_PU_h",
                             "number of clusters per track", 50, 0., 50.);
  _nLayersHit_PU_h =
      new TH1F("track_nLayersHit_PU_h",
               "number of different layers hit per track", 50, 0., 50.);
  _genPart_dR_PU_h =
      new TH1F("track_genPart_dR_PU_h", "total number of generated particles",
               600, 0., 0.4);
  _tracks_dR_PU_h = new TH1F(
      "track_tracks_dR_PU_h",
      "total number of tracks fulfilling the track requirement", 600, 0., 0.4);
  _tracks_0shared_dR_PU_h =
      new TH1F("track_0shared_tracks_dR_PU_h",
               "total number of tracks fulfilling the track requirement with 0 "
               "shared clusters along the track",
               600, 0., 0.4);
  _tracks_1shared_dR_PU_h =
      new TH1F("track_1shared_tracks_dR_PU_h",
               "total number of tracks fulfilling the track requirement with 1 "
               "shared clusters along the track",
               600, 0., 0.4);
  _tracks_2shared_dR_PU_h =
      new TH1F("track_2shared_tracks_dR_PU_h",
               "total number of tracks fulfilling the track requirement with 2 "
               "shared clusters along the track",
               600, 0., 0.4);
  _tracks_3shared_dR_PU_h =
      new TH1F("track_3shared_tracks_dR_PU_h",
               "total number of tracks fulfilling the track requirement with 3 "
               "shared clusters along the track",
               600, 0., 0.4);
  _tracks_4shared_dR_PU_h =
      new TH1F("track_4shared_tracks_dR_PU_h",
               "total number of tracks fulfilling the track requirement with 4 "
               "shared clusters along the track",
               600, 0., 0.4);
  _tracks_manyshared_dR_PU_h =
      new TH1F("track_manyshared_tracks_dR_PU_h",
               "total number of tracks fulfilling the track requirement with 4 "
               "shared clusters along the track",
               600, 0., 0.4);
}

void TrackAnalysis::fill(std::unordered_map<unsigned, GenParticle>& tracks,
                         float coneR) {
  std::cout << "filling trackAnalysis with " << tracks.size()
            << " number of tracks." << std::endl;

  // first choose all tracks fulfilling track requirement
  size_t nChosenTracks = 0;
  size_t signalTracks = 0;
  for (auto& trackInfo : tracks) {
    float deltaR = trackInfo.second.deltaR();

    // is it within jet cone
    if (deltaR < coneR) {
      // all gen particles
      _genPart_dR_PU_h->Fill(deltaR);
      if (!trackInfo.second.isPU()) {
        _genPart_dR_h->Fill(deltaR);
        //    std::cout << "within cone" << std::endl;
      }
    }

    // is it a track?
    if (trackInfo.second.layerIDs.size() >= _hitRequirement) {
      // check how many shared clusters are along the track
      size_t nClusters_shared_signal = 0;
      size_t nClusters_shared = 0;
      nChosenTracks++;
      // go through clusters per track
      for (const auto& sharedTrackIDs : trackInfo.second.shared_trackIDs) {
        size_t nTracks_shared_signal = 0;
        size_t nTracks_shared = 0;
        // go through tracks per cluster
        for (const auto& sharedTrackID : sharedTrackIDs) {
          // current trackID
          if (sharedTrackID == trackInfo.first) {
            throw std::runtime_error(
                "TrackID of this found in shared track IDs");
          }

          auto search = tracks.find(sharedTrackID);
          if (search != tracks.end() &&
              search->second.layerIDs.size() >= _hitRequirement) {
            if (!search->second.isPU()) {
              nTracks_shared_signal++;
            }
            nTracks_shared++;
          }
        }
        if (nTracks_shared_signal) nClusters_shared_signal++;
        if (nTracks_shared) nClusters_shared++;
      }  // clusters per track
      //
      // fill histograms filled for each track
      if (!trackInfo.second.isPU()) {
        //  std::cout << "fill nClusters" << std::endl;
        _nclusters_h->Fill(trackInfo.second.nClusters);
        _nLayersHit_h->Fill(trackInfo.second.layerIDs.size());
        _sharedClusters_h->Fill(nClusters_shared_signal);
      }
      _nclusters_PU_h->Fill(trackInfo.second.nClusters);
      _nLayersHit_PU_h->Fill(trackInfo.second.layerIDs.size());
      _sharedClusters_PU_h->Fill(nClusters_shared);

      // is it within jet cone
      if (deltaR < coneR) {
        if (!trackInfo.second.isPU()) {
          //    std::cout << "signal track" << std::endl;
          //     std::cout << "  shared clusters: " << nClusters_shared_signal
          //                << ", with PU: " << nClusters_shared << std::endl;
          signalTracks++;
          _tracks_dR_h->Fill(deltaR);
          switch (nClusters_shared_signal) {
            case 0:
              _tracks_0shared_dR_h->Fill(deltaR);
              //   std::cout << "    0" << std::endl;
              break;
            case 1:
              _tracks_1shared_dR_h->Fill(deltaR);
              //   std::cout << "    1" << std::endl;
              break;
            case 2:
              _tracks_2shared_dR_h->Fill(deltaR);
              //    std::cout << "    2" << std::endl;
              break;
            case 3:
              _tracks_3shared_dR_h->Fill(deltaR);
              // std::cout << "    3" << std::endl;
              break;
            case 4:
              _tracks_4shared_dR_h->Fill(deltaR);
              //   std::cout << "    4" << std::endl;
              break;
            default:
              _tracks_manyshared_dR_h->Fill(deltaR);
              //   std::cout << "    >4" << std::endl;
          }  // number of shared clusters

        }  // track itself is signal
        // total PU
        _tracks_dR_PU_h->Fill(deltaR);
        switch (nClusters_shared) {
          case 0:
            _tracks_0shared_dR_PU_h->Fill(deltaR);
            //   std::cout << "    0" << std::endl;
            break;
          case 1:
            _tracks_1shared_dR_PU_h->Fill(deltaR);
            //   std::cout << "    1" << std::endl;
            break;
          case 2:
            _tracks_2shared_dR_PU_h->Fill(deltaR);
            //    std::cout << "    2" << std::endl;
            break;
          case 3:
            _tracks_3shared_dR_PU_h->Fill(deltaR);
            // std::cout << "    3" << std::endl;
            break;
          case 4:
            _tracks_4shared_dR_PU_h->Fill(deltaR);
            //   std::cout << "    4" << std::endl;
            break;
          default:
            _tracks_manyshared_dR_PU_h->Fill(deltaR);
            //   std::cout << "    >4" << std::endl;
        }  // number of shared clusters
      }    // inside cone
    }      // hit requirement
  }
  std::cout << "Number of tracks fulfilling the track requirement: "
            << nChosenTracks << std::endl;
  std::cout << "SignalTracks chosen within deltaR: " << signalTracks
            << std::endl;
}
/*for (auto& trackInfo : chosenTracks) {
  // check if tracks shared within the cluster fulfill technical
  // track requirement
  size_t nClusters_shared_signal = 0;
  size_t nClusters_shared = 0;
  // go through clusters
  for (const auto& sharedTrackIDs : trackInfo.second.shared_trackIDs) {
    size_t nTracks_shared_signal = 0;
    size_t nTracks_shared = 0;
    // go through tracks per cluster
    for (const auto& sharedTrackID : sharedTrackIDs) {
      // current trackID
      if (sharedTrackID == trackInfo.first) {
        continue;
      }
      auto search = tracks.find(sharedTrackID);
      if (search != tracks.end() &&
          search->second.layerIDs.size() >= _hitRequirement) {
        if (!search->second.isPU()) {
          nTracks_shared_signal++;
        }
        nTracks_shared++;
      }
    }
         if (!trackInfo.second.isPU) {
           std::cout << "Number of tracks in cluster: " <<
       sharedTrackIDs.size()
                     << std::endl;
           std::cout << "  fullfilling hitrequirement: " << nTracks_shared
                     << std::endl;
         }
if (nTracks_shared) nClusters_shared++;
if (nTracks_shared_signal) nClusters_shared_signal++;
}
}


  for (auto& trackInfo : tracks) {
// check if tracks shared within the cluster fulfill technical
// track requirement
size_t nClusters_shared_signal = 0;
size_t nClusters_shared = 0;
// go through clusters
for (const auto& sharedTrackIDs : trackInfo.second.shared_trackIDs) {
  size_t nTracks_shared_signal = 0;
  size_t nTracks_shared = 0;
  // go through tracks per cluster
  for (const auto& sharedTrackID : sharedTrackIDs) {
    // current trackID
    if (sharedTrackID == trackInfo.first) {
      continue;
    }
    auto search = tracks.find(sharedTrackID);
    if (search != tracks.end() &&
        search->second.layerIDs.size() >= _hitRequirement) {
      if (!search->second.isPU()) {
        nTracks_shared_signal++;
      }
      nTracks_shared++;
    }
  }
       if (!trackInfo.second.isPU) {
         std::cout << "Number of tracks in cluster: " <<
     sharedTrackIDs.size()
                   << std::endl;
         std::cout << "  fullfilling hitrequirement: " << nTracks_shared
                   << std::endl;
       }
  if (nTracks_shared) nClusters_shared++;
  if (nTracks_shared_signal) nClusters_shared_signal++;
}

// current deltaR
float deltaR = trackInfo.second.deltaR();
// number of clusters
if (!trackInfo.second.isPU()) {
  std::cout << "_______signalTrack" << std::endl;
  std::cout << "number of shared tracks: " << nClusters_shared
            << ", purely signal: " << nClusters_shared_signal << std::endl;
  // signal only
  _nclusters_h->Fill(trackInfo.second.nClusters);
  _sharedClusters_h->Fill(nClusters_shared_signal);
  _nLayersHit_h->Fill(trackInfo.second.layerIDs.size());
  // inside cone
  if (deltaR < coneR) {
    _genPart_dR_h->Fill(deltaR);
    // track only if it fulfills the minimum hit requirement (on different
    // layers)
    if (trackInfo.second.layerIDs.size() >= _hitRequirement) {
      _tracks_dR_h->Fill(deltaR);

      //   std::cout << "----Number of shared clusters for this track: "
      //             << trackInfo.second.nClusters_shared << std::endl;

      switch (nClusters_shared_signal) {
        case 0:
          _tracks_0shared_dR_h->Fill(deltaR);
          //   std::cout << "    0" << std::endl;
          break;
        case 1:
          _tracks_1shared_dR_h->Fill(deltaR);
          //   std::cout << "    1" << std::endl;
          break;
        case 2:
          _tracks_2shared_dR_h->Fill(deltaR);
          //    std::cout << "    2" << std::endl;
          break;
        case 3:
          _tracks_3shared_dR_h->Fill(deltaR);
          // std::cout << "    3" << std::endl;
          break;
        case 4:
          _tracks_4shared_dR_h->Fill(deltaR);
          //   std::cout << "    4" << std::endl;
          break;
        default:
          _tracks_manyshared_dR_PU_h->Fill(deltaR);
          //   std::cout << "    >4" << std::endl;
      }  // number of shared clusters
    }    // track
  }
}  // !PU
_nclusters_PU_h->Fill(trackInfo.second.nClusters);
_sharedClusters_PU_h->Fill(nClusters_shared);
_nLayersHit_PU_h->Fill(trackInfo.second.layerIDs.size());
// inside cone
if (deltaR < coneR) {
  _genPart_dR_PU_h->Fill(deltaR);
  // track only if it fulfills the minimum hit requirement (on different
  // layers)
  if (trackInfo.second.layerIDs.size() >= _hitRequirement) {
    _tracks_dR_PU_h->Fill(deltaR);
    switch (nClusters_shared) {
      case 0:
        _tracks_0shared_dR_PU_h->Fill(deltaR);
        break;
      case 1:
        _tracks_1shared_dR_PU_h->Fill(deltaR);
        break;
      case 2:
        _tracks_2shared_dR_PU_h->Fill(deltaR);
        break;
      case 3:
        _tracks_3shared_dR_PU_h->Fill(deltaR);
        break;
      case 4:
        _tracks_4shared_dR_PU_h->Fill(deltaR);
        break;
      default:
        _tracks_manyshared_dR_PU_h->Fill(deltaR);
    }  // number of shared clusters
  }    // track
}
}  // going through tracks
}*/

void TrackAnalysis::write() {
  std::cout << "Writing trackanalysis" << std::endl;
  _sharedClusters_h->Write();
  _nclusters_h->Write();
  _nLayersHit_h->Write();
  _genPart_dR_h->Write();
  _tracks_dR_h->Write();
  _tracks_0shared_dR_h->Write();
  _tracks_1shared_dR_h->Write();
  _tracks_2shared_dR_h->Write();
  _tracks_3shared_dR_h->Write();
  _tracks_4shared_dR_h->Write();
  _tracks_manyshared_dR_h->Write();

  _sharedClusters_PU_h->Write();
  _nclusters_PU_h->Write();
  _nLayersHit_PU_h->Write();
  _genPart_dR_PU_h->Write();
  _tracks_dR_PU_h->Write();
  _tracks_0shared_dR_PU_h->Write();
  _tracks_1shared_dR_PU_h->Write();
  _tracks_2shared_dR_PU_h->Write();
  _tracks_3shared_dR_PU_h->Write();
  _tracks_4shared_dR_PU_h->Write();
  _tracks_manyshared_dR_PU_h->Write();
}

void TrackAnalysis::normalize(float norm) {
  if (norm == 0.) {
    throw std::runtime_error(
        "TrackClusterAnalysis::normalizeMergedClusters:Trying to divide by "
        "zero!");
  }
  float scale = 1. / norm;
  std::cout << "Normalizing track analysis: " << norm << std::endl;
  _sharedClusters_h->Scale(scale);
  _nclusters_h->Scale(scale);
  _nLayersHit_h->Scale(scale);
  _genPart_dR_h->Scale(scale);
  _tracks_dR_h->Scale(scale);
  _tracks_0shared_dR_h->Scale(scale);
  _tracks_1shared_dR_h->Scale(scale);
  _tracks_2shared_dR_h->Scale(scale);
  _tracks_3shared_dR_h->Scale(scale);
  _tracks_4shared_dR_h->Scale(scale);
  _tracks_manyshared_dR_h->Scale(scale);

  _sharedClusters_PU_h->Scale(scale);
  _nclusters_PU_h->Scale(scale);
  _nLayersHit_PU_h->Scale(scale);
  _genPart_dR_PU_h->Scale(scale);
  _tracks_dR_PU_h->Scale(scale);
  _tracks_0shared_dR_PU_h->Scale(scale);
  _tracks_1shared_dR_PU_h->Scale(scale);
  _tracks_2shared_dR_PU_h->Scale(scale);
  _tracks_3shared_dR_PU_h->Scale(scale);
  _tracks_4shared_dR_PU_h->Scale(scale);
  _tracks_manyshared_dR_PU_h->Scale(scale);
}
