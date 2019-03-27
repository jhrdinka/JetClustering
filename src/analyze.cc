#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <istream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// ROOT includes
#include <TFile.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TTree.h>
#include "TCanvas.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TSystemDirectory.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

// fastjet includes
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/SoftKiller.hh"  // In external code, this should be fastjet/contrib/SoftKiller.hh
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"

// my classes
#include "Cluster.hh"
#include "ClusterCollection.hh"
#include "GenParticle.hh"
#include "GenParticleAnalysis.hh"
#include "GenParticleCollection.hh"
#include "Jet.hh"
#include "JetAnalysis.hh"
#include "JetCollection.hh"
#include "Layer.hh"
#include "Module.hh"
#include "RecHit.hh"
#include "RecHitCalibration.hh"
#include "RecHitCollection.hh"
#include "TrackAnalysis.hh"
#include "TrackClusterAnalysis.hh"

#include <chrono>
#include <ctime>

bool debug = false;
// bool debug = true;

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

//----------------------------------------------------------------------
struct selection {
  double ptmin;
  double ptmax;
  double absetamin;
  double absetamax;
};

//----------------------------------------------------------------------

void produceJets(GenParticleCollection &constituents, JetCollection &jets,
                 float r, const selection &cuts, float ptMinCut, float ptMaxCut,
                 size_t &nJets, bool doSubstructure = false);

template <class Sequence>
void convertJets(Sequence seq, const vector<PseudoJet> &pseudojets, float r,
                 GenParticleCollection &constituents, JetCollection &jets,
                 float ptMinCut, float ptMaxCut, size_t &nJets,
                 bool doSubstructure = false);

std::vector<std::vector<RecHit *>> hitsInJets(JetCollection &recojets,
                                              RecHitCollection &recHits,
                                              float dr, float zMin = -1650.,
                                              float zMax = 1650.,
                                              float rMin = 0.,
                                              float rMax = 160.);

bool createJetsPerEvent(std::vector<JetCollection> &jetsPerEvent,
                        unsigned &eventCount, unsigned &nEvents, double R,
                        const selection &cuts, double ptMinCut, double ptMaxCut,
                        size_t &nJets, bool doSubstructure, TFile *inFile,
                        const std::string &genInfoTreeName = "genInfo");

bool createSimParticlesMap(
    std::vector<std::unordered_map<unsigned, GenParticle>>
        &simParticlesPerEvent,
    unsigned &eventCount, unsigned &nEvents, unsigned pileUpOffset,
    TFile *inFile, const std::vector<JetCollection> &jetsPerEvent, float maxR,
    const std::string &simInfoTreeName = "simInfo");

bool isClusterable(int pdgid, int status);

unsigned long long globalLayerID(const ULong64_t &moduleID);

bool selectPT(float pT, float pTMin, float pTMax);

//----------------------------------------------------------------------
int main(int argc, char *argv[]) {
  // Check the number of parameters
  if (argc < 8) {
    // Tell the user how to run the program
    std::cerr << "Usage: " << argv[0] << " [inputdir/] "
              << " [output.root] "
              << " [Nevts] "
              << " [pTmin] "
              << " [pTmax] "
              << " [file offset] "
              << " [nFiles] " << std::endl;
    return 1;
  }

  //-----------------------------------------------------------------------
  //----------------------------Definitions--------------------------------
  //-----------------------------------------------------------------------

  // Global parameters and definitions
  std::string dirname = argv[1];
  unsigned nEvents = std::atoi(argv[3]);
  float ptMinCut = std::atof(argv[4]);
  float ptMaxCut = std::atof(argv[5]);
  unsigned fileOffset = std::atof(argv[6]);
  unsigned nFiles = std::atof(argv[7]);
  std::cout << "Begin analysis of '" << nEvents << "' Events "
            << "of directory '" << dirname << "' , with pTCuts: '[" << ptMinCut
            << "," << ptMaxCut << "]'." << std::endl;
  const char *ext = ".root";
  std::string treeName = "events";
  // Analysis defintions
  // JetAnalysis
  /// JetAnalysis jetAnalysis("JetAnalysis", ptMinCut, ptMaxCut);
  // GenParticleAnalysis
  /// GenParticleAnalysis genParticleAnalysis("GenPart", ptMinCut, ptMaxCut,
  /// true);
  // SimParticleAnalysis
  /// GenParticleAnalysis simParticleAnalysis("SimPart", ptMinCut, ptMaxCut,
  /// true);
  // cluster analysis
  TrackClusterAnalysis analysis("all");
  TrackClusterAnalysis pixelAnalysis("pixel");
  TrackClusterAnalysis macroPixelAnalysis("macroPixel");
  TrackClusterAnalysis stripAnalysis("Strip");
  TrackAnalysis trackAnalysis;
  TH1F *_hits_PU_h = new TH1F("all_hits_PU_h", "_hits_PU_h", 1200, 0., 100.);
  TH1F *_hits_PU_pixel_h =
      new TH1F("pixel_hits_PU_h", "pixel_hits_PU_h", 1200, 0., 0.1);
  TH1F *_hits_PU_macroPixel_h =
      new TH1F("macroPixel_hits_PU_h", "macroPixel_hits_PU_h", 1200, 0., 100.);
  TH1F *_hits_PU_strip_h =
      new TH1F("strip_hits_PU_h", "strip_hits_PU_h", 1200, 0., 100.);

  // total number of clusters
  size_t nClusters = 0;
  // total number of jets (for normalization)
  size_t nJets = 0;
  // Global parameters
  // jet parameters
  constexpr float R = 0.4;
  constexpr bool doSubstructure = false;
  constexpr unsigned pileUpOffset = 2.5e6;
  constexpr double pTCut =
      0.015;  // 15MeV pT cut: minimum requirement to pass at least 2 layers
  constexpr float siliconTolerance = 0.1;  // tolerance that particle is within
                                           // silicon (O(silicon thickness = 100
                                           // mum))
  const long long unsigned mask = 0xf;
  // jet selection cuts
  selection cuts;
  cuts.ptmin = 2.5;
  cuts.ptmax = 20000.;
  cuts.absetamin = -2.5;
  cuts.absetamax = 2.5;

  unsigned eventCount = 0;
  unsigned offsetFileCount = 0;
  unsigned fileCount = 0;
  // when jet is taken
  unsigned validEventCount = 0;

  //-----------------------------------------------------------------------
  //------------------------Read in gen particles---------------------------
  //-----------------------------------------------------------------------

  // go through files, events
  TSystemDirectory dir(dirname.data(), dirname.data());
  TList *files = dir.GetListOfFiles();
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    // file offset
    while ((file = (TSystemFile *)next()) && (offsetFileCount <= fileOffset)) {
      offsetFileCount++;
    }

    while ((file = (TSystemFile *)next()) && (fileCount < nFiles) &&
           (eventCount < nEvents)) {
      fileCount++;
      fname = dirname + file->GetName();
      if (!file->IsDirectory() && fname.EndsWith(ext)) {
        // open file
        TFile *inFile = new TFile(fname.Data());
        if (!inFile || inFile->IsZombie() || !inFile->IsOpen() ||
            (inFile->GetSize() < 1000)) {
          std::cerr << "Could not open file: " << fname << std::endl;
          delete inFile;
          continue;
        }

        // jets per event
        std::vector<JetCollection> jetsPerEvent;
        if (!createJetsPerEvent(jetsPerEvent, eventCount, nEvents, R, cuts,
                                ptMinCut, ptMaxCut, nJets, doSubstructure,
                                inFile)) {
          inFile->Close();
          delete inFile;
          continue;
        }
        std::cout << "Number of valid jets for this event: "
                  << jetsPerEvent.size() << std::endl;
        // map indicating if the particle created a track above threshold
        std::vector<std::unordered_map<unsigned, GenParticle>>
            simParticlesPerEvent;
        if (!createSimParticlesMap(simParticlesPerEvent, eventCount, nEvents,
                                   pileUpOffset, inFile, jetsPerEvent, R)) {
          inFile->Close();
          delete inFile;
          continue;
        }

        //-----------------------------------------------------------------------
        //--------------------------read-in-clusters-----------------------------
        //-----------------------------------------------------------------------
        std::cout << "Reading in clusters......" << std::endl;

        // number of total hits per Event
        std::vector<unsigned> nHitsPerEvent(simParticlesPerEvent.size());
        // number of pile-up hits per event
        std::vector<unsigned> nPUHitsPerEvent(simParticlesPerEvent.size());
        // number of total hits per Event
        std::vector<unsigned> nHitsPerEvent_pixel(simParticlesPerEvent.size());
        // number of pile-up hits per event
        std::vector<unsigned> nPUHitsPerEvent_pixel(
            simParticlesPerEvent.size());
        // number of total hits per Event
        std::vector<unsigned> nHitsPerEvent_macroPixel(
            simParticlesPerEvent.size());
        // number of pile-up hits per event
        std::vector<unsigned> nPUHitsPerEvent_macroPixel(
            simParticlesPerEvent.size());
        // number of total hits per Event
        std::vector<unsigned> nHitsPerEvent_strip(simParticlesPerEvent.size());
        // number of pile-up hits per event
        std::vector<unsigned> nPUHitsPerEvent_strip(
            simParticlesPerEvent.size());

        // first read in cluster and module information
        TTreeReader clusterReader(treeName.c_str(), inFile);
        TTreeReaderValue<UInt_t> eventNr(clusterReader, "event_nr");
        TTreeReaderValue<ULong64_t> moduleID(clusterReader, "moduleID");
        TTreeReaderValue<std::vector<float>> clusterX(clusterReader, "g_x");
        TTreeReaderValue<std::vector<float>> clusterY(clusterReader, "g_y");
        TTreeReaderValue<std::vector<float>> clusterZ(clusterReader, "g_z");
        TTreeReaderValue<std::vector<UShort_t>> nCells(clusterReader, "nCells");
        TTreeReaderValue<std::vector<UShort_t>> numTracksPerCluster(
            clusterReader, "nTracksPerCluster");
        TTreeReaderValue<std::vector<UInt_t>> tracksPerCluster(
            clusterReader, "trackIDsPerCluster");

        // collect modules per event
        if (clusterReader.GetEntries(true) < 1) {
          eventCount--;
          continue;
        }
        // go through all modules
        while (clusterReader.Next()) {
          if ((*eventNr) < nEvents) {
            // trackInfo for this event
            // access jets of corresponding event
            auto &genjets = jetsPerEvent.at(*eventNr);
            unsigned trackIDIndex = 0;

            // access sim particles for current event
            auto &simParticles = simParticlesPerEvent.at(*eventNr);
            // calculate global layer ID
            auto layerID = globalLayerID(*moduleID);
            //------------now fill information per cluster------------
            for (size_t j = 0; j < clusterX->size(); j++) {
              // fill cluster Position
              TVector3 clusterPosition(clusterX->at(j), clusterY->at(j),
                                       clusterZ->at(j));

              // set parameters
              float deltaR_cluster = R;
              float eta_jet = 0.;
              float phi_jet = 0.;
              // number of tracks belonging to this cluster
              unsigned short nTracks = numTracksPerCluster->at(j);
              size_t nValidTracks = 0;
              bool isPU = 0;
              size_t puTracks = 0;
              // all shared trackIDs for this cluster
              std::set<unsigned> shared_trackIDs;
              shared_trackIDs.insert(
                  tracksPerCluster->begin() + trackIDIndex,
                  tracksPerCluster->begin() + trackIDIndex + nTracks);

              for (size_t iTracks = 0; iTracks < nTracks; iTracks++) {
                unsigned trackID = tracksPerCluster->at(trackIDIndex);

                auto search = simParticles.find(trackID);
                if (search != simParticles.end()) {
                  // to be counted particles need to have certain pT and should
                  // not be created in the sensitive silicon
                  bool outsideSilicon =
                      (std::fabs(search->second.vertex().Perp() -
                                 clusterPosition.Perp()) <= siliconTolerance)
                          ? false
                          : true;
                  if (search->second.pt() >= pTCut && outsideSilicon) {
                    nValidTracks++;
                    // store track information if it is a valid track access
                    // track information
                    // if there is at least one track from PU
                    if (search->second.isPU()) puTracks++;
                    // fill track info
                    search->second.nClusters++;
                    search->second.shared_trackIDs.push_back(shared_trackIDs);
                    // erase trackID from current track
                    search->second.shared_trackIDs.back().erase(trackID);
                    search->second.layerIDs.insert(layerID);
                  }
                } else {
                  std::cerr << "fileName. " << fname << std::endl;
                  std::cerr
                      << "trackID: " << trackID
                      << ", number of sim particles: " << simParticles.size()
                      << std::endl;
                  std::cerr << "nTracksPerCluster::TrackID not found in sim "
                               "particles!"
                            << std::endl;
                  return 0;
                }

                trackIDIndex++;
              }
              if (puTracks) isPU = true;

              if (nValidTracks == 0) {
                nValidTracks = 1;
              }
              // count hits per event with signal and PU only
              nHitsPerEvent.at(*eventNr) += nCells->at(j);
              if (isPU) {
                nPUHitsPerEvent.at(*eventNr) += nCells->at(j);
              }

              // go through jets
              for (size_t i = 0; i < genjets.size(); i++) {
                // first access jet
                auto jet = genjets.at(i);
                // vertex correction
                TVector3 clusterPos(clusterPosition.X(), clusterPosition.Y(),
                                    clusterPosition.Z() - jet->vertexZ());
                // access LorentzVector
                auto pJet = jet->p4().Vect();
                // calculate deltaR
                float dR = pJet.DeltaR(clusterPos);
                if (dR < deltaR_cluster) {
                  deltaR_cluster = dR;
                  eta_jet = jet->p4().Eta();
                  phi_jet = jet->p4().Phi();
                }
              }
              // only write out, when inside jet
              if (deltaR_cluster < R) {
                analysis.fill_cluster(deltaR_cluster, clusterPosition.Eta(),
                                      eta_jet, clusterPosition.Phi(), phi_jet,
                                      nTracks, nValidTracks, nCells->at(j),
                                      isPU);
                nClusters++;
              }

              // fill region info per cluster and barrel
              const unsigned long long &systemID = ((*moduleID) & mask);

              if ((systemID == 0) || (systemID == 1)) {
                //-------------------barrel-------------------
                // fill region information
                if (layerID < 4) {
                  nHitsPerEvent_pixel.at(*eventNr) += nCells->at(j);
                  if (isPU) {
                    nPUHitsPerEvent_pixel.at(*eventNr) += nCells->at(j);
                  }
                  if (deltaR_cluster < R) {
                    pixelAnalysis.fill_cluster(
                        deltaR_cluster, clusterPosition.Eta(), eta_jet,
                        clusterPosition.Phi(), phi_jet, nTracks, nValidTracks,
                        nCells->at(j), isPU);
                  }
                } else if (layerID < 8) {
                  nHitsPerEvent_macroPixel.at(*eventNr) += nCells->at(j);
                  if (isPU) {
                    nPUHitsPerEvent_macroPixel.at(*eventNr) += nCells->at(j);
                  }
                  // macroPixel
                  if (deltaR_cluster < R) {
                    macroPixelAnalysis.fill_cluster(
                        deltaR_cluster, clusterPosition.Eta(), eta_jet,
                        clusterPosition.Phi(), phi_jet, nTracks, nValidTracks,
                        nCells->at(j), isPU);
                  }
                } else {
                  nHitsPerEvent_strip.at(*eventNr) += nCells->at(j);
                  if (isPU) {
                    nPUHitsPerEvent_strip.at(*eventNr) += nCells->at(j);
                  }
                  // strip
                  if (deltaR_cluster < R) {
                    stripAnalysis.fill_cluster(
                        deltaR_cluster, clusterPosition.Eta(), eta_jet,
                        clusterPosition.Phi(), phi_jet, nTracks, nValidTracks,
                        nCells->at(j), isPU);
                  }
                }  // different regions
              }    // barrel
            }      // go through clusters
          }        // check number of events
        }          // go through reader entries (modules)

        std::cout << "number of hits for this event: " << nHitsPerEvent.size()
                  << std::endl;
        for (size_t nEvent = 0; nEvent < nHitsPerEvent.size(); nEvent++) {
          // calculate only when we have a valid event - when jets passed the
          // cut
          if (jetsPerEvent.at(nEvent).size()) {
            validEventCount++;
            // the relative hit difference per event between signal + PU and
            // PU
            // only
            double relHitDiff =
                nHitsPerEvent.at(nEvent) - nPUHitsPerEvent.at(nEvent);
            if (relHitDiff < 0. || nHitsPerEvent.at(nEvent) == 0. ||
                nPUHitsPerEvent.at(nEvent) == 0.) {
              std::cout << "Number of hits for this event: "
                        << nHitsPerEvent.at(nEvent)
                        << ", number of PU hits: " << nPUHitsPerEvent.at(nEvent)
                        << std::endl;
              throw std::runtime_error(
                  "Relative hit difference < 0! Or number of hits == 0.!");
            }
            relHitDiff /= nHitsPerEvent.at(nEvent);
            _hits_PU_h->Fill(relHitDiff);
            // pixel
            double relHitDiff_pixel = nHitsPerEvent_pixel.at(nEvent) -
                                      nPUHitsPerEvent_pixel.at(nEvent);
            if (relHitDiff_pixel < 0. || nHitsPerEvent_pixel.at(nEvent) == 0. ||
                nPUHitsPerEvent_pixel.at(nEvent) == 0.) {
              throw std::runtime_error(
                  "Relative hit difference pixel < 0! Or number of hits == "
                  "0.!");
            }
            relHitDiff_pixel /= nHitsPerEvent_pixel.at(nEvent);
            _hits_PU_pixel_h->Fill(relHitDiff_pixel);
            // macroPixel
            double relHitDiff_macroPixel =
                nHitsPerEvent_macroPixel.at(nEvent) -
                nPUHitsPerEvent_macroPixel.at(nEvent);
            if (relHitDiff_macroPixel < 0. ||
                nHitsPerEvent_macroPixel.at(nEvent) == 0. ||
                nPUHitsPerEvent_macroPixel.at(nEvent) == 0.) {
              throw std::runtime_error(
                  "Relative hit difference macroPixel < 0! Or number of hits "
                  "== 0.!");
            }
            relHitDiff_macroPixel /= nHitsPerEvent_macroPixel.at(nEvent);
            _hits_PU_macroPixel_h->Fill(relHitDiff_macroPixel);
            // strip
            double relHitDiff_strip = nHitsPerEvent_strip.at(nEvent) -
                                      nPUHitsPerEvent_strip.at(nEvent);
            if (relHitDiff_strip < 0. || nHitsPerEvent_strip.at(nEvent) == 0. ||
                nPUHitsPerEvent_strip.at(nEvent) == 0.) {
              throw std::runtime_error(
                  "Relative hit difference strip < 0! Or number of hits == "
                  "0.!");
            }
            relHitDiff_strip /= nHitsPerEvent_strip.at(nEvent);
            _hits_PU_strip_h->Fill(relHitDiff_strip);

            std::cout << "Fill track analysis..." << std::endl;
            //-------------------trackAnalysis------------------
            trackAnalysis.fill(simParticlesPerEvent.at(nEvent), R);
          }  // check if event is valid (jets passed cuts)

        }  // go through events

        inFile->Close();
        delete inFile;
      }
    }
  }

  std::cout << "Before writing: " << std::endl;
  // store plots in output file
  TFile outfile(argv[2], "RECREATE");

  // normalize
  std::cout << "number of files: " << fileCount << std::endl;
  std::cout << "eventCount: " << eventCount << std::endl;
  std::cout << "number of jets: " << nJets << std::endl;
  /// jetAnalysis.normalize(nJets);
  /// genParticleAnalysis.normalize(nJets);
  /// simParticleAnalysis.normalize(nJets);
  analysis.normalize(nJets);
  // @todo can this vary? the number jets per region?
  pixelAnalysis.normalize(nJets);
  macroPixelAnalysis.normalize(nJets);
  stripAnalysis.normalize(nJets);
  trackAnalysis.normalize(nJets);
  std::cout << "Total number of clusters: " << nClusters << std::endl;

  std::cout << "valid Events: " << validEventCount << ", nJets: " << nJets
            << std::endl;

  _hits_PU_h->Scale(1. / validEventCount);
  _hits_PU_h->Write();

  _hits_PU_pixel_h->Scale(1. / validEventCount);
  _hits_PU_pixel_h->Write();

  _hits_PU_macroPixel_h->Scale(1. / validEventCount);
  _hits_PU_macroPixel_h->Write();

  _hits_PU_strip_h->Scale(1. / validEventCount);
  _hits_PU_strip_h->Write();

  analysis.write();
  trackAnalysis.write();
  /// jetAnalysis.write();
  /// genParticleAnalysis.write();
  /// simParticleAnalysis.write();

  std::cout << "After writing" << std::endl;
  // per region
  TDirectory *dir_pixel = outfile.mkdir("pixel");
  dir_pixel->cd();
  pixelAnalysis.write();
  outfile.cd();

  TDirectory *dir_macroPixel = outfile.mkdir("macroPixel");
  dir_macroPixel->cd();
  macroPixelAnalysis.write();
  outfile.cd();

  TDirectory *dir_strip = outfile.mkdir("strip");
  dir_strip->cd();
  stripAnalysis.write();
  outfile.cd();
  outfile.Close();

  return 0;
}

//-----------------------------------------------------------------------
//-----------------------Read in gen particles-------------------------
//-----------------------------------------------------------------------
bool createJetsPerEvent(std::vector<JetCollection> &jetsPerEvent,
                        unsigned &eventCount, unsigned &nEvents, double R,
                        const selection &cuts, double ptMinCut, double ptMaxCut,
                        size_t &nJets, bool doSubstructure, TFile *inFile,
                        const std::string &genInfoTreeName) {
  std::cout << "before gen" << std::endl;
  TTreeReader reader(genInfoTreeName.c_str(), inFile);
  TTreeReaderValue<std::vector<float>> genpart_eta(reader, "gen_eta");
  TTreeReaderValue<std::vector<float>> genpart_phi(reader, "gen_phi");
  TTreeReaderValue<std::vector<float>> genpart_pt(reader, "gen_pt");
  TTreeReaderValue<std::vector<float>> genpart_energy(reader, "gen_energy");
  TTreeReaderValue<std::vector<int>> genpart_charge(reader, "gen_charge");
  TTreeReaderValue<std::vector<unsigned>> genpart_status(reader, "gen_status");
  TTreeReaderValue<std::vector<int>> genpart_pdgid(reader, "gen_pdgid");
  TTreeReaderValue<std::vector<float>> genpart_vertexX(reader, "gen_vertexX");
  TTreeReaderValue<std::vector<float>> genpart_vertexY(reader, "gen_vertexY");
  TTreeReaderValue<std::vector<float>> genpart_vertexZ(reader, "gen_vertexZ");
  std::cout << "Processing event: " << eventCount << " of " << nEvents
            << std::endl;
  std::cout << "reading file: " << inFile->GetName() << std::endl;
  if (reader.GetEntries(true) < 1) {
    return false;
  }
  // go through events
  while (reader.Next() && (eventCount < nEvents)) {
    eventCount++;
    // create gen particle collection for this event
    GenParticleCollection genparts;
    // read in genparticles
    for (size_t i = 0; i < (*genpart_pdgid).size(); i++) {
      // only use clusterable particles (exclude neutrinos, secondaries
      // (should not happen anyway))
      if (isClusterable((*genpart_pdgid)[i], (*genpart_status)[i])) {
        // set the lorentz vector of the gen particle
        TLorentzVector genpart_p4;
        genpart_p4.SetPtEtaPhiE((*genpart_pt)[i], (*genpart_eta)[i],
                                (*genpart_phi)[i], (*genpart_energy)[i]);
        TVector3 vertex((*genpart_vertexX)[i], (*genpart_vertexY)[i],
                        (*genpart_vertexZ)[i]);
        // add gen particle to collection
        const GenParticle *genpart = new GenParticle(
            genpart_p4, (*genpart_pdgid)[i], (*genpart_status)[i], vertex,
            (*genpart_charge)[i], R, false);

        genparts.Add(genpart);
      }
    }  //** end go through gen particles

    //-----------------------------------------------------------------------
    //----------------------------Produce-jets-------------------------------
    //-----------------------------------------------------------------------
    JetCollection genjets;
    // produce jet collections (anti-kT R = 0.4)
    produceJets(genparts, genjets, R, cuts, ptMinCut, ptMaxCut, nJets,
                doSubstructure);
    // add to genjets
    jetsPerEvent.push_back(genjets);

    //-----------------------------------------------------------------------
    //--------------------------genjets-analysis-----------------------------
    //-----------------------------------------------------------------------
    /*  for (size_t i = 0; i < genjets.size(); i++) {
        // first access jet
        auto jet = genjets.at(i);
        jetAnalysis.fill(*jet);
      }*/
    //-----------------------------------------------------------------------
    //--------------------------genpart-analysis-----------------------------
    //-----------------------------------------------------------------------
    // go through each gen particle
    /// for (size_t ipart = 0; ipart < genparts.size(); ipart++) {
    ///   genParticleAnalysis.fillParticleAndJets(*(genparts.at(ipart)),
    ///                                           genjets, R);
    ///  }

  }  //*end go through events of genparticle tree of current file
  return true;
}

//-----------------------------------------------------------------------
//-----------------------Read-in-sim-particle-map------------------------
//-----------------------------------------------------------------------
bool createSimParticlesMap(
    std::vector<std::unordered_map<unsigned, GenParticle>>
        &simParticlesPerEvent,
    unsigned &eventCount, unsigned &nEvents, unsigned pileUpOffset,
    TFile *inFile, const std::vector<JetCollection> &jetsPerEvent, float maxR,
    const std::string &simInfoTreeName) {
  std::cout << "before sim" << std::endl;
  std::vector<float> *sim_eta = new std::vector<float>;
  std::vector<float> *sim_phi = new std::vector<float>;
  std::vector<float> *sim_pt = new std::vector<float>;
  std::vector<float> *sim_energy = new std::vector<float>;
  std::vector<int> *sim_charge = new std::vector<int>;
  std::vector<unsigned> *sim_bits = new std::vector<unsigned>;
  std::vector<unsigned> *sim_status = new std::vector<unsigned>;
  std::vector<int> *sim_pdgid = new std::vector<int>;
  std::vector<float> *sim_vertexX = new std::vector<float>;
  std::vector<float> *sim_vertexY = new std::vector<float>;
  std::vector<float> *sim_vertexZ = new std::vector<float>;

  TTree *simInfo = (TTree *)inFile->Get(simInfoTreeName.c_str());

  if (!simInfo) {
    eventCount--;
    return false;
  }

  if (simInfo->GetEntries() < 1) {
    eventCount--;
    return false;
  }

  simInfo->SetBranchAddress("sim_eta", &sim_eta);
  simInfo->SetBranchAddress("sim_phi", &sim_phi);
  simInfo->SetBranchAddress("sim_pt", &sim_pt);
  simInfo->SetBranchAddress("sim_energy", &sim_energy);
  simInfo->SetBranchAddress("sim_status", &sim_status);
  simInfo->SetBranchAddress("sim_pdgid", &sim_pdgid);
  simInfo->SetBranchAddress("sim_charge", &sim_charge);
  simInfo->SetBranchAddress("sim_vertexX", &sim_vertexX);
  simInfo->SetBranchAddress("sim_vertexY", &sim_vertexY);
  simInfo->SetBranchAddress("sim_vertexZ", &sim_vertexZ);
  simInfo->SetBranchAddress("sim_bits", &sim_bits);

  size_t puCounter = 0;
  for (unsigned iEvent = 0;
       (iEvent < simInfo->GetEntries()) && (iEvent < nEvents); iEvent++) {
    simInfo->GetEntry(iEvent);
    // access jets of corresponding event
    const auto &genjets = jetsPerEvent[iEvent];
    //    auto currentJets = jetsPerEvent.at(iEvent);
    std::unordered_map<unsigned, GenParticle> simParticleMap;
    for (size_t i_simPart = 0; i_simPart < sim_bits->size(); i_simPart++) {
      //	if ((*sim_pt)[i_simPart]>0.){
      // the four momentum
      TLorentzVector simpart_p4;
      simpart_p4.SetPtEtaPhiE((*sim_pt)[i_simPart], (*sim_eta)[i_simPart],
                              (*sim_phi)[i_simPart], (*sim_energy)[i_simPart]);
      // the vertex
      TVector3 vertex((*sim_vertexX)[i_simPart], (*sim_vertexY)[i_simPart],
                      (*sim_vertexZ)[i_simPart]);

      bool isPU = false;
      if (sim_bits->at(i_simPart) >= pileUpOffset) {
        puCounter++;
        isPU = true;
      }
      // deltaR to closest jet
      // go through jets (maximum 2) for each track and take
      // closest one
      float deltaR_track = maxR;
      // go throuh jets for this event (maximal two)
      for (size_t i = 0; i < genjets.size(); i++) {
        // calculate deltaR
        float dR = genjets.at(i)->p4().DeltaR(simpart_p4);
        if (dR < deltaR_track) {
          deltaR_track = dR;
        }
      }

      // store the sim particle in the map
      GenParticle simPart(simpart_p4, (*sim_pdgid)[i_simPart],
                          (*sim_status)[i_simPart], vertex,
                          (*sim_charge)[i_simPart], deltaR_track, isPU);

      simParticleMap.insert(
          std::pair<unsigned, GenParticle>(sim_bits->at(i_simPart), simPart));

      ///     simParticleAnalysis.fillParticleAndJets(simPart,
      ///     currentJets, R);
      //	}
    }
    std::cout << "#simparticles for this event: " << simParticleMap.size()
              << ", sim_bits->size(): " << sim_bits->size() << ", with "
              << puCounter << "particles coming from PU" << std::endl;
    simParticlesPerEvent.push_back(simParticleMap);
  }

  delete sim_eta;
  delete sim_phi;
  delete sim_pt;
  delete sim_energy;
  delete sim_charge;
  delete sim_bits;
  delete sim_status;
  delete sim_pdgid;
  delete sim_vertexX;
  delete sim_vertexY;
  delete sim_vertexZ;

  return true;
}

//------------------------------------------------------------------------------------------------------

void produceJets(GenParticleCollection &constituents, JetCollection &jets,
                 float r, const selection &cuts, float ptMinCut, float ptMaxCut,
                 size_t &nJets, bool doSubstructure) {
  // first convert constituents into fastjet pseudo-jets
  vector<PseudoJet> input_particles;
  for (unsigned i = 0; i < constituents.size(); i++) {
    // Constituent pj = constituents.at(i);
    auto pseudojet =
        PseudoJet(constituents.at(i)->px(), constituents.at(i)->py(),
                  constituents.at(i)->pz(), constituents.at(i)->energy());
    pseudojet.set_user_index(i);
    input_particles.push_back(pseudojet);
  }
  // Initial clustering with anti-kt algorithm
  JetDefinition jetDef = JetDefinition(antikt_algorithm, r, E_scheme, Best);

  Selector select_eta = SelectorAbsEtaRange(cuts.absetamin, cuts.absetamax);
  Selector select_pt = SelectorPtRange(cuts.ptmin, cuts.ptmax);
  Selector select_jets = select_eta && select_pt;

  ClusterSequence clust_seq(input_particles, jetDef);

  vector<PseudoJet> akjets = sorted_by_pt(clust_seq.inclusive_jets());
  // apply cuts
  akjets = select_jets(akjets);

  // eventually apply substructure and store in custom dataformat
  convertJets(clust_seq, akjets, r, constituents, jets, ptMinCut, ptMaxCut,
              nJets, doSubstructure);
}

//------------------------------------------------------------------------------------------------------
template <class Sequence>
void convertJets(Sequence seq, const vector<PseudoJet> &pseudojets, float r,
                 GenParticleCollection &constituents, JetCollection &jets,
                 float ptMinCut, float ptMaxCut, size_t &nJets,
                 bool doSubstructure) {
  // select only the two jets with the highest pT per event
  float jetPT0 = 0;
  float jetPT1 = 1;
  Jet jet0;
  Jet jet1;
  // momentum of current jet
  TLorentzVector p4;
  for (unsigned j = 0; j < pseudojets.size(); j++) {
    // get the jet for analysis
    PseudoJet this_jet = pseudojets[j];
    // get corresponding particle
    p4.SetPtEtaPhiM(this_jet.pt(), this_jet.eta(), this_jet.phi(),
                    std::max(this_jet.m(), 0.0));

    // map of all jet vertex z position to their pTs of this jet
    std::unordered_map<float, float> jetVertices;
    // go through all subjets and find the vertex with most weight (highest pT)
    // size_t counter = 0;
    for (auto &subjet : this_jet.constituents()) {
      // the z vertex of the constituent
      float vertexZ = constituents.at(subjet.user_index())->vertex().Z();
      //  std::cout << "vertex: " << vertexZ << std::endl;
      // find the corresponding pT
      auto search = jetVertices.find(vertexZ);
      if (search != jetVertices.end()) {
        search->second += subjet.pt();
        //  counter++;
      } else {
        jetVertices[vertexZ] = subjet.pt();
      }
    }
    // find the vertex with highest pt and assign it to this jet
    auto finalVertexZ = std::max_element(
        jetVertices.begin(), jetVertices.end(),
        [](const pair<float, float> &a, const pair<float, float> &b) {
          return a.second < b.second;
        });
    /*  std::cout << "doubled: " << counter
                << " number of vertices per jet: " << jetVertices.size()
                << " final vertex: " << finalVertexZ->first
                << ", pT: " << finalVertexZ->second << std::endl;*/
    // current het
    Jet jet(p4, finalVertexZ->first);

    if (doSubstructure) {
      // N-subjettiness
      double beta = 1.0;

      Nsubjettiness nSub1_beta1(1, OnePass_KT_Axes(),
                                NormalizedMeasure(beta, r));
      Nsubjettiness nSub2_beta1(2, OnePass_KT_Axes(),
                                NormalizedMeasure(beta, r));
      Nsubjettiness nSub3_beta1(3, OnePass_KT_Axes(),
                                NormalizedMeasure(beta, r));
      NsubjettinessRatio nSub21_beta1(2, 1, OnePass_KT_Axes(),
                                      NormalizedMeasure(beta, r));
      NsubjettinessRatio nSub32_beta1(3, 2, OnePass_KT_Axes(),
                                      NormalizedMeasure(beta, r));

      jet.setTau1(nSub1_beta1(this_jet));
      jet.setTau2(nSub2_beta1(this_jet));
      jet.setTau3(nSub3_beta1(this_jet));
      jet.setTau21(nSub21_beta1(this_jet));
      jet.setTau32(nSub32_beta1(this_jet));

      // soft drop
      beta = 0.0;
      double zcut = 0.1;
      SoftDrop softDrop(beta, zcut);
      PseudoJet softdrop_jet = softDrop(this_jet);
      p4.SetPtEtaPhiM(softdrop_jet.pt(), softdrop_jet.eta(), softdrop_jet.phi(),
                      std::max(softdrop_jet.m(), 0.0));
      jet.setSDmass(p4.M());

      // store jet in collection
    }

    // jet with highest pT
    if (jet.pt() > jetPT0) {
      jetPT0 = jet.pt();
      jet0 = jet;
    } else if (jet.pt() > jetPT1) {
      jetPT1 = jet.pt();
      jet1 = jet;
    }
  }
  std::cout << "Final two jets of event: " << std::endl;
  std::cout << " pT: " << jetPT0 << "," << jetPT1 << std::endl;
  if (selectPT(jetPT0, ptMinCut, ptMaxCut)) {
    nJets++;
    jets.Add(new Jet(jet0));
  }
  if (selectPT(jetPT1, ptMinCut, ptMaxCut)) {
    nJets++;
    jets.Add(new Jet(jet1));
  }
}

bool selectPT(float pT, float pTMin, float pTMax) {
  return ((pT >= pTMin) && (pT <= pTMax));
}

// return global layerID given the cellID or moduleID
unsigned long long globalLayerID(const ULong64_t &moduleID) {
  //
  // masks to decode on which layer we are
  // @todo could be done nicer with decoder or at least defining the
  // masks somewhere else
  const long long unsigned mask = 0xf;
  const long long unsigned layerMaskBarrel = 0x1f0;
  const long long unsigned posNegMask = 0x10;
  const long long unsigned layerMaskEC = 0x3e0;
  const unsigned nLayers_barrel = 12;
  const unsigned nLayers_ec = 20;
  long long unsigned systemID = (moduleID & mask);
  long long unsigned globalLayerID = 0;
  if ((systemID == 0) || (systemID == 1)) {
    //-------------------barrel-------------------
    // calculate global barrel layer ID
    globalLayerID = ((moduleID & layerMaskBarrel) >> 4) + 6 * systemID;
    if (globalLayerID > 11) {
      throw std::runtime_error("error wrong globalLayerID in barrel");
    }

  } else {
    //-------------------endcaps-------------------
    // layer offset for endcaps
    unsigned layerOffset = 0;
    if (systemID == 3) {
      // outer
      layerOffset = 5;
    } else if (systemID == 4) {
      // fwd
      layerOffset = 11;
    }
    // distinguish positive and negative side
    // calculate global ec pos layer ID
    globalLayerID = ((moduleID & layerMaskEC) >> 5) + layerOffset;
    // add barrel layers to offset
    if (!((moduleID & posNegMask) >> 4)) {
      globalLayerID += nLayers_ec;
      globalLayerID += nLayers_barrel;
      //-------------------pos endcap-------------------
      if (globalLayerID > 52 || globalLayerID < 32) {
        throw std::runtime_error("error wrong globalLayerID in positive EC");
      }

    } else {
      globalLayerID = (19 - globalLayerID);
      globalLayerID += nLayers_barrel;
      //-------------------neg endcap-------------------
      if (globalLayerID > 32 || globalLayerID < 12) {
        throw std::runtime_error("error wrong globalLayerID in negative EC");
      }
    }
  }
  return globalLayerID;
}

//------------------------------------------------------------------------------------------------------

bool isClusterable(int pdgid, int status) {
  bool pass = true;
  // primary
  if (status != 1) pass = false;
  // do not use neutrinos
  if (fabs(pdgid) == 12) pass = false;
  if (fabs(pdgid) == 14) pass = false;
  if (fabs(pdgid) == 16) pass = false;
  return pass;
}
