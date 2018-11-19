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
#include <TTree.h>
#include "TCanvas.h"
#include "TGraph.h"
#include "TMatrixD.h"
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
#include "SimParticle.hh"
#include "TrackClusterAnalysis.hh"

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

bool isClusterable(int pdgid, int status);

std::vector<std::vector<unsigned>> translateMatrix(
    const TMatrixD &trackIDsMatrix, const long long int &globalLayerID,
    std::vector<std::unordered_set<int>> &trackIDsPerLayer,
    std::unordered_set<int> &validTrackIDs);

// hands back as a pair the valid
std::pair<size_t, size_t> nTrackAndMergedClusters(
    const std::vector<std::vector<unsigned>> &trackIDsPerCluster,
    const std::unordered_set<int> &validTrackIDs);

std::pair<size_t, size_t> nTracksPerCluster(
    const std::vector<unsigned> &clusterTrackIDs,
    const std::unordered_map<int, GenParticle> &simParticles, float pTCut,
    const std::unordered_set<int> &validTrackIDs_layer,
    size_t &nValidTracksLayer);

bool selectPT(float pT, float pTMin, float pTMax);

//----------------------------------------------------------------------
int main(int argc, char *argv[]) {
  // Check the number of parameters
  if (argc < 6) {
    // Tell the user how to run the program
    std::cerr << "Usage: " << argv[0] << " [inputdir/] "
              << " [output.root] "
              << " [Nevts] "
              << " [pTmin] "
              << " [pTmax] " << std::endl;
    return 1;
  }

  //-----------------------------------------------------------------------
  //----------------------------Definitions--------------------------------
  //-----------------------------------------------------------------------

  // Global parameters and definitions
  std::string dirname = argv[1];
  int nEvents = std::atoi(argv[3]);
  float ptMinCut = std::atof(argv[4]);
  float ptMaxCut = std::atof(argv[5]);
  std::cout << "Begin analysis of '" << nEvents << "' Events "
            << "of directory '" << dirname << "' , with pTCuts: '[" << ptMinCut
            << "," << ptMaxCut << "]'." << std::endl;
  const char *ext = ".root";
  std::string treeName = "events";
  // Analysis defintions
  // JetAnalysis
  JetAnalysis jetAnalysis("JetAnalysis", ptMinCut, ptMaxCut);
  // GenParticleAnalysis
  GenParticleAnalysis genParticleAnalysis("GenPart", ptMinCut, ptMaxCut, true);
  // SimParticleAnalysis
  GenParticleAnalysis simParticleAnalysis("SimPart", ptMinCut, ptMaxCut, true);
  // cluster analysis
  TrackClusterAnalysis analysis("all");
  TrackClusterAnalysis pixelAnalysis("pixel");
  TrackClusterAnalysis macroPixelAnalysis("macroPixel");
  TrackClusterAnalysis stripAnalysis("Strip");
  // total number of clusters
  size_t nClusters = 0;
  // total number of jets (for normalization)
  size_t nJets = 0;
  // layer analysis
  std::vector<Layer> layers(12);
  for (size_t i = 0; i < 12; i++) {
    layers.at(i).initialize(i);
  }

  // Global parameters
  // jet parameters
  constexpr float R = 0.4;
  constexpr bool doSubstructure = false;  // @todo what is that exactly?
  // jet selection cuts
  selection cuts;
  cuts.ptmin = 2.5;
  cuts.ptmax = 20000.;
  cuts.absetamin = -2.5;  // @todo maybe change to 2.5
  cuts.absetamax = 2.5;

  int eventCount = 0;

  // masks for layers
  long long int mask = 0xf;
  long long int layerMaskBarrel = 0x1f0;
  long long int posNegMask = 0x10;
  long long int layerMaskEC = 0x3e0;
  // number of layers
  const int nLayers_barrel = 12;
  const int nLayers_ec = 20;

  //-----------------------------------------------------------------------
  //------------------------Read in gen paricles---------------------------
  //-----------------------------------------------------------------------

  // go through files, events
  TSystemDirectory dir(dirname.data(), dirname.data());
  TList *files = dir.GetListOfFiles();
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file = (TSystemFile *)next()) && (eventCount < nEvents)) {
      fname = dirname + file->GetName();
      if (!file->IsDirectory() && fname.EndsWith(ext)) {
        // open file
        TFile *inFile = new TFile(fname.Data());
        if (!inFile) {
          std::cerr << "Could not open file: " << fname << std::endl;
        }

        //-----------------------------------------------------------------------
        //-----------------------Read in gen particles-------------------------
        //-----------------------------------------------------------------------
        std::cout << "before gen" << std::endl;
        TTreeReader reader("genInfo", inFile);
        TTreeReaderValue<std::vector<float>> genpart_eta(reader, "gen_eta");
        TTreeReaderValue<std::vector<float>> genpart_phi(reader, "gen_phi");
        TTreeReaderValue<std::vector<float>> genpart_pt(reader, "gen_pt");
        TTreeReaderValue<std::vector<float>> genpart_energy(reader,
                                                            "gen_energy");
        TTreeReaderValue<std::vector<int>> genpart_charge(reader, "gen_charge");
        TTreeReaderValue<std::vector<int>> genpart_status(reader, "gen_status");
        TTreeReaderValue<std::vector<int>> genpart_pdgid(reader, "gen_pdgid");
        TTreeReaderValue<std::vector<float>> genpart_vertexX(reader,
                                                             "gen_vertexX");
        TTreeReaderValue<std::vector<float>> genpart_vertexY(reader,
                                                             "gen_vertexY");
        TTreeReaderValue<std::vector<float>> genpart_vertexZ(reader,
                                                             "gen_vertexZ");

        // jets per event
        std::vector<JetCollection> jetsPerEvent;

        std::cout << "Processing event: " << eventCount << " of " << nEvents
                  << std::endl;
        std::cout << "reading file: " << fname << std::endl;
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
                  (*genpart_charge)[i]);

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
          for (size_t i = 0; i < genjets.size(); i++) {
            // first access jet
            auto jet = genjets.at(i);
            jetAnalysis.fill(*jet);
          }
          //-----------------------------------------------------------------------
          //--------------------------genpart-analysis-----------------------------
          //-----------------------------------------------------------------------
          // go through each gen particle
          for (size_t ipart = 0; ipart < genparts.size(); ipart++) {
            genParticleAnalysis.fillParticleAndJets(*(genparts.at(ipart)),
                                                    genjets, R);
          }

        }  //*end go through events of genparticle tree of current file

        //-----------------------------------------------------------------------
        //-----------------------Read-in-sim-particle-map------------------------
        //-----------------------------------------------------------------------
        std::cout << "before sim" << std::endl;
        // map indicating if the particle created a track above threshold
        std::vector<std::unordered_map<int, GenParticle>> simParticlesPerEvent;
        std::vector<float> *sim_eta = new std::vector<float>;
        std::vector<float> *sim_phi = new std::vector<float>;
        std::vector<float> *sim_pt = new std::vector<float>;
        std::vector<float> *sim_energy = new std::vector<float>;
        std::vector<int> *sim_charge = new std::vector<int>;
        std::vector<int> *sim_bits = new std::vector<int>;
        std::vector<int> *sim_status = new std::vector<int>;
        std::vector<int> *sim_pdgid = new std::vector<int>;
        std::vector<float> *sim_vertexX = new std::vector<float>;
        std::vector<float> *sim_vertexY = new std::vector<float>;
        std::vector<float> *sim_vertexZ = new std::vector<float>;

        TTree *simInfo = (TTree *)inFile->Get("simInfo");
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

        for (int iEvent = 0;
             (iEvent < simInfo->GetEntries()) && (iEvent < nEvents); iEvent++) {
          simInfo->GetEntry(iEvent);
          auto currentJets = jetsPerEvent.at(iEvent);
          std::unordered_map<int, GenParticle> simParticleMap;
          for (size_t i_simPart = 0; i_simPart < sim_bits->size();
               i_simPart++) {
            //	if ((*sim_pt)[i_simPart]>0.){
            // the four momentum
            TLorentzVector simpart_p4;
            simpart_p4.SetPtEtaPhiE((*sim_pt)[i_simPart], (*sim_eta)[i_simPart],
                                    (*sim_phi)[i_simPart],
                                    (*sim_energy)[i_simPart]);
            // the vertex
            TVector3 vertex((*sim_vertexX)[i_simPart],
                            (*sim_vertexY)[i_simPart],
                            (*sim_vertexZ)[i_simPart]);
            // store the simparticle in the map
            GenParticle simPart(simpart_p4, (*sim_pdgid)[i_simPart],
                                (*sim_status)[i_simPart], vertex,
                                (*sim_charge)[i_simPart]);
            simParticleMap.insert(
                std::pair<int, GenParticle>(sim_bits->at(i_simPart), simPart));

            simParticleAnalysis.fillParticleAndJets(simPart, currentJets, R);
            //	}
          }
          std::cout << "#simparticles for this event: " << simParticleMap.size()
                    << ", sim_bits->size(): " << sim_bits->size() << std::endl;
          simParticlesPerEvent.push_back(simParticleMap);
        }

        //-----------------------------------------------------------------------
        //--------------------------read-in-clusters-----------------------------
        //-----------------------------------------------------------------------
        // first read in cluster and module information
        TTreeReader clusterReader(treeName.c_str(), inFile);
        TTreeReaderValue<int> eventNr(clusterReader, "event_nr");
        TTreeReaderValue<long long int> moduleID(clusterReader, "moduleID");
        TTreeReaderValue<int> nChannels(clusterReader, "nChannels");
        TTreeReaderValue<int> nChannelsOn(clusterReader, "nChannelsOn");
        TTreeReaderValue<float> moduleX(clusterReader, "s_x");
        TTreeReaderValue<float> moduleY(clusterReader, "s_y");
        TTreeReaderValue<float> moduleZ(clusterReader, "s_z");
        TTreeReaderValue<std::vector<float>> clusterX(clusterReader, "g_x");
        TTreeReaderValue<std::vector<float>> clusterY(clusterReader, "g_y");
        TTreeReaderValue<std::vector<float>> clusterZ(clusterReader, "g_z");
        TTreeReaderValue<std::vector<short int>> nCells(clusterReader,
                                                        "nCells");

        // todo current workaround because of
        // https://root-forum.cern.ch/t/problem-reading-in-tmatrixd-with-ttreereader-while-working-using-simple-tree/31113
        //     TTreeReaderValue<TMatrixD>
        //     trackIDsPerCluster(clusterReader,"trackIDsPerCluster");
        TMatrixD *trackIDsMatrix = new TMatrixD();
        TTree *clusterTree = (TTree *)inFile->Get(treeName.c_str());
        clusterTree->SetBranchAddress("trackIDsPerCluster", &trackIDsMatrix);
        // collect modules per event
        std::vector<std::vector<Module>> modulesPerEvent(jetsPerEvent.size());

        // the global map for valid trackIDs (trackIDs which have at least
        // occured on one other layer
        //***-----LayerCut-begin
        std::vector<std::unordered_set<int>> validTrackIDsPerEvent(
            jetsPerEvent.size());
        // the map for each layer and event (0-11 barrel,12-32 nEC, 32-52 pEC)
        // of valid trackIDs (tracks that are seen in at least one other layer,
        // to substract secondaries only created within the material)
        std::vector<std::vector<std::unordered_set<int>>>
        trackIDsPerLayerAndEvent(jetsPerEvent.size(),
                                 std::vector<std::unordered_set<int>>(
                                     nLayers_barrel + 2 * nLayers_ec));
        //***-----LayerCut-end
        // go through all modules
        while (clusterReader.Next()) {
          if ((*eventNr) < nEvents) {
            //***-----LayerCut-begin
            // first get the global layerID of the current module to check if
            // track has occured on another layer
            long long int systemID = ((*moduleID) & mask);
            long long int globalLayerID = 0;
            if ((systemID == 0) || (systemID == 1)) {
              //-------------------barrel-------------------
              // calculate global barrel layer ID
              globalLayerID =
                  (((*moduleID) & layerMaskBarrel) >> 4) + 6 * systemID;
              if (globalLayerID > 11 || globalLayerID < 0) {
                throw std::runtime_error("error wrong globalLayerID in barrel");
              }
            } else {
              //-------------------endcaps-------------------
              // layer offset for endcaps
              int layerOffset = 0;
              if (systemID == 3) {
                // outer
                layerOffset = 5;
              } else if (systemID == 4) {
                // fwd
                layerOffset = 11;
              }
              // distinguish positive and negative side
              // calculate global ec pos layer ID
              globalLayerID = (((*moduleID) & layerMaskEC) >> 5) + layerOffset;
              // add barrel layers to offset
              if (!(((*moduleID) & posNegMask) >> 4)) {
                globalLayerID += nLayers_ec;
                globalLayerID += nLayers_barrel;
                //-------------------pos endcap-------------------
                if (globalLayerID > 52 || globalLayerID < 32) {
                  throw std::runtime_error(
                      "error wrong globalLayerID in positive EC");
                }

              } else {
                globalLayerID = (19 - globalLayerID);
                globalLayerID += nLayers_barrel;
                //-------------------neg endcap-------------------
                if (globalLayerID > 32 || globalLayerID < 12) {
                  throw std::runtime_error(
                      "error wrong globalLayerID in negative EC");
                }
              }
            }
            //***-----LayerCut-end
            auto &currentLayers = trackIDsPerLayerAndEvent.at(*eventNr);
            auto &currentValidTrackIDSet = validTrackIDsPerEvent.at(*eventNr);
            // access the current entry of the trackIDs matrix (work-around)
            clusterTree->GetEntry(clusterReader.GetCurrentEntry());
            auto trackIDsPerCluster = std::move(
                translateMatrix(*trackIDsMatrix, globalLayerID, currentLayers,
                                currentValidTrackIDSet));
            // save module information
            TVector3 modulePosition((*moduleX), (*moduleY), (*moduleZ));
            modulesPerEvent[(*eventNr)].push_back(
                Module((*moduleID), modulePosition, (*nChannelsOn),
                       (*nChannels), (*clusterX), (*clusterY), (*clusterZ),
                       std::move(trackIDsPerCluster), (*nCells)));

          }  // check number of events
        }
        //-----------------------------------------------------------------------
        //--------------------------cluster-analysis-----------------------------
        //-----------------------------------------------------------------------
        // go through events
        for (size_t iEvent = 0; iEvent < jetsPerEvent.size(); iEvent++) {
          // access jets of corresponding event
          auto &genjets = jetsPerEvent.at(iEvent);
          // access modules of corresponding event
          auto &modules = modulesPerEvent.at(iEvent);

          size_t nTracksPerEvent = 0;
          size_t nValidTracksPerEvent = 0;
          size_t nValidTracksPerEvent_layer = 0;
          // go through modules and calculate deltaR to each jet and take the
          // smallest one
          for (auto &module : modules) {
            // the systemID
            long long int systemID = (module.moduleID() & mask);
            // calculate global barrel layer ID
            long long int barrelLayerID =
                ((module.moduleID() & layerMaskBarrel) >> 4) + 6 * systemID;
            // calculate deltaR
            float defaultCone = R;
            float deltaR_module = defaultCone;

            // access the modulePosition
            auto modulePosition = module.position();
            // access the trackIDs of this module
            auto trackIDsPerCluster = module.trackIDsPerCluster();

            // calculate deltaR to each jet and take the smallest one
            // go through jets (should be two only)
            for (size_t i = 0; i < genjets.size(); i++) {
              // first access jet
              auto jet = genjets.at(i);
              // access LorentzVector
              auto pJet = jet->p4().Vect();
              // module position corrected to the vertex of the current jet
              TVector3 modulePos(modulePosition.X(), modulePosition.Y(),
                                 (modulePosition.Z() - jet->vertexZ()));
              // calculate deltaR
              float dR = pJet.DeltaR(modulePos);
              if (dR < deltaR_module) {
                deltaR_module = dR;
              }
            }  // go through jets

            const std::vector<float> &clusterPositionX =
                module.clusterPositionX();
            const std::vector<float> &clusterPositionY =
                module.clusterPositionY();
            const std::vector<float> &clusterPositionZ =
                module.clusterPositionZ();

            const std::vector<short int> &nCellsPerCluster =
                module.nCellsPerCluster();

            //------------now fill information per cluster------------
            for (size_t j = 0; j < clusterPositionX.size(); j++) {
              // fill cluster analysis
              TVector3 clusterPosition(clusterPositionX.at(j),
                                       clusterPositionY.at(j),
                                       clusterPositionZ.at(j));
              // all clusters also outside jets
              analysis.fill_allClusters(clusterPosition.Eta(),
                                        clusterPosition.Phi(),
                                        module.moduleID());
              // set parameters
              float deltaR_cluster = defaultCone;
              float eta_jet = 0.;
              float phi_jet = 0.;
              // number of tracks in this cluster
              size_t nTracksLayer = 0;
              auto nTracks = nTracksPerCluster(
                  trackIDsPerCluster.at(j), simParticlesPerEvent.at(iEvent),
                  0.015, validTrackIDsPerEvent.at(iEvent), nTracksLayer);
              nTracksPerEvent += nTracks.first;
              nValidTracksPerEvent += nTracks.second;
              nValidTracksPerEvent_layer += nTracksLayer;
              // go through jets
              for (size_t i = 0; i < genjets.size(); i++) {
                // first access jet
                auto jet = genjets.at(i);

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
              }  // go through jets
              // only write out, when inside jet
              if (deltaR_cluster < defaultCone) {
                analysis.fill_cluster(deltaR_cluster, clusterPosition.Eta(),
                                      eta_jet, clusterPosition.Phi(), phi_jet,
                                      nTracks.first, nTracks.second,
                                      nTracksLayer, module.moduleID(),
                                      nCellsPerCluster.at(j));
                nClusters++;
              }

              // fill region info per cluster
              if ((systemID == 0) || (systemID == 1)) {
                //-------------------barrel-------------------
                // fill region information
                if (barrelLayerID < 4) {
                  // pixel
                  pixelAnalysis.fill_allClusters(clusterPosition.Eta(),
                                                 clusterPosition.Phi(),
                                                 module.moduleID());
                  if (deltaR_cluster < defaultCone) {
                    pixelAnalysis.fill_cluster(
                        deltaR_cluster, clusterPosition.Eta(), eta_jet,
                        clusterPosition.Phi(), phi_jet, nTracks.first,
                        nTracks.second, nTracksLayer, module.moduleID(),
                        nCellsPerCluster.at(j));
                  }
                } else if (barrelLayerID < 8) {
                  // macroPixel
                  macroPixelAnalysis.fill_allClusters(clusterPosition.Eta(),
                                                      clusterPosition.Phi(),
                                                      module.moduleID());
                  if (deltaR_cluster < defaultCone) {
                    macroPixelAnalysis.fill_cluster(
                        deltaR_cluster, clusterPosition.Eta(), eta_jet,
                        clusterPosition.Phi(), phi_jet, nTracks.first,
                        nTracks.second, nTracksLayer, module.moduleID(),
                        nCellsPerCluster.at(j));
                  }
                } else {
                  // strip
                  stripAnalysis.fill_allClusters(clusterPosition.Eta(),
                                                 clusterPosition.Phi(),
                                                 module.moduleID());
                  if (deltaR_cluster < defaultCone) {
                    stripAnalysis.fill_cluster(
                        deltaR_cluster, clusterPosition.Eta(), eta_jet,
                        clusterPosition.Phi(), phi_jet, nTracks.first,
                        nTracks.second, nTracksLayer, module.moduleID(),
                        nCellsPerCluster.at(j));
                  }
                }
              }

            }  // go through clusters of this module

            // fill module information
            if (deltaR_module < defaultCone) {
              analysis.fill_module(deltaR_module, module.nClusters(),
                                   module.occupancy(), module.moduleID());
              // fill layer info
              if ((systemID == 0) || (systemID == 1)) {
                //-------------------barrel-------------------
                // create layer information
                layers.at(barrelLayerID)
                    .addClusters(module.clusterPositionX(),
                                 module.clusterPositionY(),
                                 module.clusterPositionZ(), deltaR_module,
                                 module.occupancy(), trackIDsPerCluster);
                // fill region information
                if (barrelLayerID < 4) {
                  // pixel
                  pixelAnalysis.fill_module(deltaR_module, module.nClusters(),
                                            module.occupancy(),
                                            module.moduleID());
                } else if (barrelLayerID < 8) {
                  // macroPixel
                  macroPixelAnalysis.fill_module(
                      deltaR_module, module.nClusters(), module.occupancy(),
                      module.moduleID());
                } else {
                  // strip
                  stripAnalysis.fill_module(deltaR_module, module.nClusters(),
                                            module.occupancy(),
                                            module.moduleID());
                }
              }
            }
          }  // go through  modules
          std::cout << "nTracksPerEvent: " << nTracksPerEvent
                    << ", nValidTracksPerEvent: " << nValidTracksPerEvent
                    << ", validLayer: " << nValidTracksPerEvent_layer
                    << std::endl;

        }  // go through events
      }
    }
  }

  std::cout << "Before writing: " << std::endl;
  // store plots in output file
  TFile outfile(argv[2], "RECREATE");

  // normalize
  std::cout << "eventCount: " << eventCount << std::endl;
  std::cout << "number of jets: " << nJets << std::endl;
  jetAnalysis.normalize(nJets);
  genParticleAnalysis.normalize(nJets);
  simParticleAnalysis.normalize(nJets);
  analysis.normalize(nJets);
  // @todo can this vary? the number jets per region?
  pixelAnalysis.normalize(nJets);
  macroPixelAnalysis.normalize(nJets);
  stripAnalysis.normalize(nJets);
  std::cout << "Total number of clusters: " << nClusters << std::endl;

  analysis.write();
  jetAnalysis.write();
  genParticleAnalysis.write();
  simParticleAnalysis.write();

  float nClusters_layer[12];
  float nlayers[12] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

  TDirectory *dir_layers = outfile.mkdir("layers");
  dir_layers->cd();
  size_t layerCounter = 0;
  for (auto &layer : layers) {
    nClusters_layer[layerCounter] = layer.nClusters();
    layerCounter++;
    //  layer.fillDistances(_dz_dR, _dphi_dR);
    layer.normalize(eventCount);
    layer.print();
  }
  outfile.cd();

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

  TGraph *nClusters_layers = new TGraph(12, nlayers, nClusters_layer);
  nClusters_layers->SetTitle("#clusters per barrel in jet");
  nClusters_layers->GetXaxis()->SetTitle("layer number");
  nClusters_layers->GetYaxis()->SetTitle("#clusters in jet");
  nClusters_layers->Write();

  outfile.Close();

  return 0;
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
  //  std::cout << "Final two jets of event: " << std::endl;
  //  std::cout << " pT: " << jetPT0 << "," << jetPT1 << std::endl;
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

//------------------------------------------------------------------------------------------------------
std::vector<std::vector<unsigned>> translateMatrix(
    const TMatrixD &trackIDsMatrix, const long long int &globalLayerID,
    std::vector<std::unordered_set<int>> &trackIDsPerLayer,
    std::unordered_set<int> &validTrackIDs) {
  std::vector<std::vector<unsigned>> trackIDs;
  for (int i = 0; i < trackIDsMatrix.GetNrows(); i++) {
    std::vector<unsigned> trackIDsPerCluster;
    for (size_t j = 0; j < trackIDsMatrix.GetNcols(); j++) {
      unsigned trackID = trackIDsMatrix[i][j];
      if (trackID) {
        trackIDsPerCluster.push_back(trackID);
        // check if it is already marked as valid
        auto searchValid = validTrackIDs.find(trackID);
        if (searchValid == validTrackIDs.end()) {
          // not marked as valid yet, check if it appears in another layer
          size_t layerCounter = 0;
          for (auto &layer : trackIDsPerLayer) {
            // do not search in current layer
            if (layerCounter != globalLayerID) {
              auto searchLayer = layer.find(trackID);
              if (searchLayer != layer.end()) {
                // found on another layer - mark as valid
                validTrackIDs.insert(trackID);
              }
            }
            layerCounter++;
          }  // go through layers
          // insert to current layer
          trackIDsPerLayer.at(globalLayerID).insert(trackID);
        }
      }
    }  // go through columns
    trackIDs.push_back(trackIDsPerCluster);
  }  // go through rows
  return trackIDs;
}

//------------------------------------------------------------------------------------------------------

std::pair<size_t, size_t> nTrackAndMergedClusters(
    const std::vector<std::vector<unsigned>> &trackIDsPerCluster,
    const std::unordered_set<int> &validTrackIDs) {
  // The number of clusters which have a track assigned
  size_t nValidClusters = 0;
  // the number of merged clusters
  size_t nMergedClusters = 0;
  // go through the clusters of this module and count
  for (auto &clusterTrackIDs : trackIDsPerCluster) {
    // the valid track IDs of this cluster
    size_t nValidTracks = 0;
    // go through the trackIDs of the current cluster
    for (auto trackID : clusterTrackIDs) {
      // check if it was assigned to be valid
      auto search = validTrackIDs.find(trackID);
      if (search != validTrackIDs.end()) {
        nValidTracks++;
      }
    }  // go through trackIDs of the cluster
    // if there is at least one valid track in the cluster, the cluster is
    // valid
    if (nValidTracks) {
      if (nValidTracks > 1) {
        // if there is more than one track in the cluster, the cluster is
        // merged
        nMergedClusters++;
      }
      nValidClusters++;
    }
  }  // go through clusters
  return std::make_pair(nValidClusters, nMergedClusters);
}
//------------------------------------------------------------------------------------------------------

std::pair<size_t, size_t> nTracksPerCluster(
    const std::vector<unsigned> &clusterTrackIDs,
    const std::unordered_map<int, GenParticle> &simParticles, float pTCut,
    const std::unordered_set<int> &validTrackIDs_layer,
    size_t &nValidTracksLayer) {
  // the valid track IDs of this cluster
  size_t nTracks = 0;
  size_t nValidTracks = 0;
  // go through the trackIDs of the current cluster
  if (!clusterTrackIDs.size()) std::cout << "no size" << std::endl;
  for (auto trackID : clusterTrackIDs) {
    // check if it was assigned to be valid
    nTracks++;
    auto search = simParticles.find(trackID);
    if (search != simParticles.end()) {
      if (search->second.pt() >= pTCut) {
        nValidTracks++;
      }
    } else {
      throw std::runtime_error(
          "nTracksPerCluster::TrackID not found in simparticles!");
    }
    // check if it was found on more than one layer
    auto searchLayer = validTrackIDs_layer.find(trackID);
    if (searchLayer != validTrackIDs_layer.end()) {
      nValidTracksLayer++;
    }
  }  // go through trackIDs of the cluster
  // number of tracks can not be zero, set to one if no track of this cluster
  // passes the cut
  if (nValidTracksLayer == 0) {
    nValidTracksLayer = 1;
  }
  if (nValidTracks == 0) {
    nValidTracks = 1;
  }
  if (nTracks == 0) {
    throw std::runtime_error(
        "nTracksPerCluster::Number of tracks participating to this cluster "
        "zero!");
  }
  return std::make_pair(nTracks, nValidTracks);
}
