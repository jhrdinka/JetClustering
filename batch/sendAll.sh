python submitJetClustering.py -i /eos/experiment/fcc/hh/simulation/samples/tracker_calobarrel/diJets/20GeV/NTUP -n 100 --njobs 100 -o output/FccDiJet20GeV -q 1nh
python submitJetClustering.py -i /eos/experiment/fcc/hh/simulation/samples/tracker_calobarrel/diJets/50GeV/NTUP -n 100 --njobs 100 -o output/FccDiJet50GeV -q 1nh
python submitJetClustering.py -i /eos/experiment/fcc/hh/simulation/samples/tracker_calobarrel/diJets/100GeV/NTUP -n 100 --njobs 100 -o output/FccDiJet100GeV  -q 1nh
python submitJetClustering.py -i /eos/experiment/fcc/hh/simulation/samples/tracker_calobarrel/diJets/200GeV/NTUP -n 100 --njobs 100 -o output/FccDiJet200GeV  -q 1nh
python submitJetClustering.py -i /eos/experiment/fcc/hh/simulation/samples/tracker_calobarrel/diJets/500GeV/NTUP -n 100 --njobs 100 -o output/FccDiJet500GeV  -q 1nh
python submitJetClustering.py -i /eos/experiment/fcc/hh/simulation/samples/tracker_calobarrel/diJets/1000GeV/NTUP -n 100 --njobs 100 -o output/FccDiJet1000GeV  -q 1nh
python submitJetClustering.py -i /eos/experiment/fcc/hh/simulation/samples/tracker_calobarrel/diJets/2000GeV/NTUP -n 100 --njobs 100 -o output/FccDiJet2000GeV  -q 1nh


python submitJetClustering.py -o output/FccDiJet20GeV --collect
python submitJetClustering.py -o output/FccDiJet50GeV --collect
python submitJetClustering.py -o output/FccDiJet100GeV --collect
python submitJetClustering.py -o output/FccDiJet200GeV --collect
python submitJetClustering.py -o output/FccDiJet500GeV --collect
python submitJetClustering.py -o output/FccDiJet1000GeV --collect
python submitJetClustering.py -o output/FccDiJet2000GeV --collect

hadd -f output/merged_dijet.root output/FccDiJet*GeV/out/FccDiJet*GeV.root
