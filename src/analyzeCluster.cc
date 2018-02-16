#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <glob.h>

// ROOT includes
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TLorentzVector.h>

// fastjet includes
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/SoftKiller.hh" // In external code, this should be fastjet/contrib/SoftKiller.hh
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh" 

// my classes
#include "Cluster.hh"
#include "ClusterCollection.hh"
#include "GenParticle.hh"
#include "GenParticleCollection.hh"
#include "Jet.hh"
#include "JetCollection.hh"
#include "JetPlots.hh"

bool Debug = false;
//bool debug = true;

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

template <class Collection>
void produceJets(Collection& input_particles, JetCollection & jets, const float& r, const selection cuts, const bool doPuSubtraction = false, const bool doSubstructure = false);

template <class Sequence>
void ConvertJets(Sequence seq, vector<PseudoJet> pseudojets, const float& r, JetCollection& jets, const bool doSubstructure = false);
void MatchJets(JetCollection& genjets, JetCollection& recojets, float dr);
void PrintJets(JetCollection & jets);

vector<std::string> GlobVector(const std::string& pattern){
  glob_t glob_result;
  glob(pattern.c_str(),GLOB_TILDE,NULL,&glob_result);
  vector<std::string> files;
  for(unsigned int i=0;i<glob_result.gl_pathc;++i){
    files.push_back(string(glob_result.gl_pathv[i]));
  }
  globfree(&glob_result);
  return files;
}

//----------------------------------------------------------------------
int Main(int argc, char* argv[]){


 // Check the number of parameters
  if (argc < 5) {
      // Tell the user how to run the program
      std::cerr << "Usage: " << argv[0] << " [input.root] " << " [output.root] " <<" [Nevts] [CMS/FCC]" <<std::endl;
      return 1;
  }

  TString runType = argv[4];
  if(runType != "CMS" && runType != "FCC"){
     cerr << "Should specifiy CMS or FCC as last argument"<<endl;
     return 1;
  }

  // declare histograms                                                                                              
  vector<float> ptvals;
  ptvals = {10., 20., 30.,50., 75., 100., 150., 200., 300., 500., 750., 1000., 1500., 2000., 5000.};

  JetPlots gen_plots  = JetPlots("gen", ptvals);
  JetPlots reco_plots = JetPlots("reco", ptvals);

  std::stringstream inputPath;
  inputPath << argv[1] << "/*ntuple.root";
  std::cout << inputPath.str() << std::endl;
  vector<std::string> files = GlobVector(inputPath.str());

  // -- Loop over files
  for (const auto& ifile : files){
    std::cout << ifile << std::endl;
    // ---   Tree stuff declarations
    TFile *f = new TFile(ifile.c_str());
    TTree *t = (TTree*)f->Get("ana/hgc");

    vector<Float_t> *cluster_pt        = 0;
    vector<Float_t> *cluster_eta       = 0;
    vector<Float_t> *cluster_phi       = 0;
    vector<Float_t> *cluster_energy    = 0;
    vector<Float_t> *cluster_x         = 0;
    vector<Float_t> *cluster_y         = 0;
    vector<Float_t> *cluster_z         = 0;
    //vector<Float_t> *cluster_thickness = 0;
    //vector<Float_t> *cluster_layer     = 0;

    vector<Float_t> *genpart_pt       = 0;
    vector<Float_t> *genpart_eta      = 0;
    vector<Float_t> *genpart_phi      = 0;
    vector<Float_t> *genpart_energy   = 0;
    vector<Float_t> *genpart_status   = 0;
    vector<Float_t> *genpart_pdgid    = 0;

    t->SetBranchAddress("cluster_eta",    &cluster_eta);
    t->SetBranchAddress("cluster_phi",    &cluster_phi);
    t->SetBranchAddress("cluster_pt",     &cluster_pt);
    t->SetBranchAddress("cluster_energy", &cluster_energy);
    t->SetBranchAddress("cluster_x",      &cluster_x);
    t->SetBranchAddress("cluster_y",      &cluster_y);
    t->SetBranchAddress("cluster_z",      &cluster_z);
    //if(runType == "CMS")t->SetBranchAddress("rechit_thickness", &rechit_thickness);
    //t->SetBranchAddress("rechit_layer", &cluster_layer);

    t->SetBranchAddress("gen_eta", &genpart_eta);
    t->SetBranchAddress("gen_phi", &genpart_phi);
    t->SetBranchAddress("gen_pt", &genpart_pt);
    t->SetBranchAddress("gen_energy", &genpart_energy);
    t->SetBranchAddress("gen_status", &genpart_status);
    t->SetBranchAddress("gen_pdgid", &genpart_pdgid);

    // generic declarations
    TLorentzVector cluster_p4, cluster_pos;
    TLorentzVector genpart_p4;

    //read all entries and fill the histograms
    Long64_t nentries = t->GetEntries();
    Long64_t nmax = stoi(argv[3]);
    Int_t nrun = TMath::Min(nentries, nmax);

    for (Long64_t i=0;i<nrun;i++) {
      t->GetEntry(i);

      cout<<" ---- processing event : "<<i<<endl;

      // ---  prepare genparts
      GenParticleCollection genparts;
      unsigned genpart_size = genpart_pt->size();

      for (unsigned i = 0; i < genpart_size; i++) {
	// initialize genpart
	genpart_p4.SetPtEtaPhiE(genpart_pt->at(i), genpart_eta->at(i), genpart_phi->at(i), genpart_energy->at(i));
	GenParticle *genpart = genparts.AddGenParticle(genpart_p4, genpart_pdgid->at(i), genpart_status->at(i));
      }

      GenParticleCollection clean_genparts;
      for (unsigned i = 0; i < genparts.size(); i++) {
        GenParticle *g = genparts.at(i);
        if (!g->isClusterable()) continue;
        clean_genparts.Add(new GenParticle(*g));
      }

      // ---  prepare clusters ----------------------------------------------
    
      ClusterCollection clusters;
      unsigned cluster_size = cluster_pt->size();
      for (unsigned i = 0; i < cluster_size; i++) {
	//for (unsigned i = 0; i < 1000; i++) {
	// initialize cluster
	cluster_p4.SetPtEtaPhiE(cluster_pt->at(i), cluster_eta->at(i), cluster_phi->at(i), cluster_energy->at(i));
	cluster_pos.SetXYZT(cluster_x->at(i), cluster_y->at(i), cluster_z->at(i), 0.0);
	Cluster *cluster = clusters.AddCluster(cluster_p4, cluster_pos);
      }
      
      if(Debug) cout<<"cluster size: "<<clusters.size()<<endl;
      ClusterCollection clean_clusters;
      for (unsigned i = 0; i < clusters.size(); i++) {
	Cluster *r = clusters.at(i);
	clean_clusters.Add(new Cluster(*r));
      }
      if(Debug) cout<<"clean cluster size: "<<clean_clusters.size()<<endl;
    
    
      // ---------- Produce jets ------------------------------------------------
    
      // declare jet collections
      JetCollection genjets;
      JetCollection recojets;
    
      // produce jet collections (anti-kT R = 0.4)
      bool doSubstructure  = false;
      bool doPuSubtraction = false;
    
      selection cuts;
      cuts.ptmin  = 2.5;
      cuts.ptmax  = 5000.;
      cuts.absetamin = 0.0;
      cuts.absetamax = 1.3;
    
      produceJets(clean_genparts, genjets, 0.4, cuts, false, doSubstructure);
      produceJets(clean_clusters, recojets, 0.4, cuts, doPuSubtraction, doSubstructure);

      // match reco to gen (need this in order to make resolution plots)
      MatchJets(genjets, recojets, 0.4);

      if (Debug) { 
        cout<<" ------  gen jets ------"<<endl;
        PrintJets(genjets);
        cout<<" ------  reco jets ------ "<<endl;
        PrintJets(recojets);
      }

      // fill plots
    
      gen_plots.fill(genjets);
      reco_plots.fill(recojets);
    
    } // end event loop
  }  
  //store plots in output file
  TFile outfile(argv[2],"RECREATE");
  
  // store plots in output root tree
  gen_plots.write();
  reco_plots.write();
  
  outfile.Close();

  return 0;
}


//------------------------------------------------------------------------------------------------------
template <class Collection>
void produceJets(Collection& constituents, JetCollection& jets, const float& r, const selection cuts, const bool doPuSubtraction = false, const bool doSubstructure = false) {
      
   // first convert constituents into fastjet pseudo-jets
   vector <PseudoJet> input_particles;
   for (unsigned i = 0; i < constituents.size(); i++) {
      //Constituent pj = constituents.at(i);
      input_particles.push_back( PseudoJet( constituents.at(i)->px(), constituents.at(i)->py(), constituents.at(i)->pz(), constituents.at(i)->energy()));
   }
   
   // Initial clustering with anti-kt algorithm
   JetAlgorithm algorithm = antikt_algorithm;
   double jet_rad = r; // jet radius for anti-kt algorithm
   JetDefinition jetDef = JetDefinition(algorithm,jet_rad,E_scheme,Best);
   
   vector<PseudoJet> akjets;
   
   Selector select_eta   = SelectorAbsEtaRange(cuts.absetamin, cuts.absetamax);
   Selector select_pt    = SelectorPtRange(cuts.ptmin, cuts.ptmax);
   Selector select_jets  = select_eta && select_pt;
   
   if(doPuSubtraction) {
       
       // jet area correction parameters, are used only if doPU = true
       float etamin = 0;
       float etamax = 2.5;
       float spacing = 0.50;
       Selector selector = SelectorAbsRapRange(etamin, etamax);

       AreaDefinition areaDef(active_area, GhostedAreaSpec(selector));
       ClusterSequenceArea clust_seq(input_particles,jetDef, areaDef);
       
       vector<PseudoJet> antikt_jets  = clust_seq.inclusive_jets();

       // now apply PU subtraction
       RectangularGrid grid(-etamax, etamax, spacing, spacing, selector);
       GridMedianBackgroundEstimator gmbge(grid);
       gmbge.set_particles(input_particles);
       Subtractor subtractor(&gmbge);
       akjets = subtractor(antikt_jets);
       
       // apply cuts
       akjets = select_jets(sorted_by_pt(akjets));
       
       // eventually apply substructure and store in custom dataformat
       ConvertJets(clust_seq, akjets, r, jets, doSubstructure);

       // soft killer PU (could be used later)

       /*double grid_size = 0.4;
       SoftKiller soft_killer(-3.0, 3.0, grid_size, grid_size) ;
       double pt_threshold;
       vector<PseudoJet> soft_killed_event;
       soft_killer.apply(input_particles, soft_killed_event, pt_threshold);

       ClusterSequence clust_seq_kill(soft_killed_event, jetDef);  
       vector<PseudoJet> kill_jets = clust_seq_kill.inclusive_jets();

       kill_jets = sel_jets(kill_jets);

       if(debug) cout << setprecision(4);
       if(debug) cout << "Soft Killer applied a pt threshold of " << pt_threshold << endl;

       // run things and print the result
       //----------------------------------------------------------
       if(debug) cout << "# original hard jets" << endl;
       for (unsigned int i=0; i<antikt_jets.size(); i++){
         const PseudoJet &jet = antikt_jets[i];
         if(debug) cout << "pt = " << jet.pt()
              << ", rap = " << jet.rap()
              << ", mass = " << jet.m() << endl;
       }
       if(debug) cout << endl;

       if(debug) cout << "# jets after applying the soft killer" << endl;
       for (unsigned int i=0; i<kill_jets.size(); i++){
         const PseudoJet &jet = kill_jets[i];
         if(debug) cout << "pt = " << jet.pt()
              << ", rap = " << jet.rap()
              << ", mass = " << jet.m() << endl;
       }
       if(debug) cout << endl;
       }
       */
   }
   else {
       
       ClusterSequence clust_seq(input_particles,jetDef);
       akjets  = sorted_by_pt(clust_seq.inclusive_jets());

       // apply cuts
       akjets = select_jets(akjets);

       // eventually apply substructure and store in custom dataformat
       ConvertJets(clust_seq, akjets, r, jets, doSubstructure);
   }
   
}



//------------------------------------------------------------------------------------------------------
template <class Sequence>
void ConvertJets(Sequence seq, vector<PseudoJet> pseudojets, const float& r, JetCollection& jets, const bool doSubstructure = false){

   TLorentzVector p4;
   for (unsigned j = 0; j < pseudojets.size() ; j++) { 

       // get the jet for analysis
       PseudoJet this_jet = pseudojets[j];

       p4.SetPtEtaPhiM(this_jet.pt(), this_jet.eta(), this_jet.phi(), std::max(this_jet.m(),0.0));
       Jet jet(p4);

       if (doSubstructure) {
       
           // N-subjettiness
           double beta = 1.0;

           Nsubjettiness         nSub1_beta1(1,   OnePass_KT_Axes(), NormalizedMeasure(beta, r));
           Nsubjettiness         nSub2_beta1(2,   OnePass_KT_Axes(), NormalizedMeasure(beta, r));
           Nsubjettiness         nSub3_beta1(3,   OnePass_KT_Axes(), NormalizedMeasure(beta, r));
           NsubjettinessRatio    nSub21_beta1(2,1, OnePass_KT_Axes(), NormalizedMeasure(beta, r));
           NsubjettinessRatio    nSub32_beta1(3,2, OnePass_KT_Axes(), NormalizedMeasure(beta, r));

           jet.setTau1(nSub1_beta1(this_jet));
           jet.setTau2(nSub2_beta1(this_jet));
           jet.setTau3(nSub3_beta1(this_jet));
           jet.setTau21(nSub21_beta1(this_jet));
           jet.setTau32(nSub32_beta1(this_jet));

           // soft drop
           beta = 0.0;
           double zcut = 0.1;
           SoftDrop softDrop(beta,zcut);
           PseudoJet softdrop_jet = softDrop(this_jet);
           p4.SetPtEtaPhiM(softdrop_jet.pt(), softdrop_jet.eta(), softdrop_jet.phi(), std::max(softdrop_jet.m(),0.0));
           jet.setSDmass(p4.M());

           // store jet in collection
       }
       
       jets.Add(new Jet(jet));
   }   
}


//------------------------------------------------------------------------------------------------------
void MatchJets(JetCollection& genjets, JetCollection& recojets, float dr){

    for (unsigned i = 0; i < genjets.size() ; i++) { 
       float dr0 = 999.;
       Jet *gj = genjets.at(i);
       // will be best matching recojet
       Jet *rj0; 
       for (unsigned j = 0; j < recojets.size() ; j++) { 
          Jet *rj = recojets.at(j);
          float dr_gr = rj->p4().DeltaR(gj->p4());
          if( dr_gr < dr0 ){
             rj0 = rj;
             dr0 = dr_gr;
          } 
       }
       // assign genjet ref. to best matching recojet (and vice versa)
       if (dr0 < dr){ 
         rj0->setRef(gj);
         gj->setRef(rj0);
       }
    }
}

//------------------------------------------------------------------------------------------------------
void PrintJets(JetCollection& jets){
    cout<<" -- Print Jet collection -- "<<endl; 
    for (unsigned j = 0; j < jets.size() ; j++){
       jets.at(j)->print();
       if(jets.at(j)->ref()) {
           (jets.at(j)->ref())->print();
       }
    }
}
