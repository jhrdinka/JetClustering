[]() Package installation
--------------------------
First load the required environment on lxplus:
```
source init.sh
```
Install fastjet (to be done only once):
```
./install-fastjet.sh
```

[]() Run instructions
----------------------


First load the required environment on lxplus:
```
source init.sh
```
Compile analysis code (can be found and modified in ```src/analyze.cc```):

```
make -j 4
```
Run analysis:
```
./analyze [input_file] [output_file] [number_of_events] [CMS/FCC]
```
e.g:
```
./analyze /eos/experiment/fcc/hh/simulation/samples/tracker_calobarrel/ditop/500GeV/NTUP/output_helsens_20171011151211690.root test.root 10 FCC
```


[]() LSF submission
--------------------

The script ```batch/submitJetClustering.py``` allows to run this script on LSF queues.
Choose a new output directory, in which one out/ and std/ will be build.
Specify algorithm: "antiKt", "antiKt_cluster" or "simpleCone" and if input are cluster add: -t
```
python batch/submitJetClustering.py -i [NTUP_dir] -n [nevts_per_job] -o [output_dir] --njobs [number_of_jobs] -q [queue]
```
example:
``
 python submitJetClustering.py -i /eos/experiment/fcc/hh/simulation/samples/v02_pre/physics/ljets/bFieldOn/etaTo1.5/20GeV/reco/topoClusters/ -n 100 --algorithm antiKt_cluster -o 20GeVljets --njobs 166 -q 8nh -t
``
Jobs will be collected for certain jet alorithm and pt, in directory /eos/experiment/fcc/users/c/cneubuse/JetClustering/:

```
python submitJetClustering.py --collect antiKt --collectPt 20
```

Already merged files per pt, can be collected via:

```
python submitJetClustering.py --collect antiKt --allPts
```

for default pts = ( 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000 ) GeV

e.g:
```
python submitJetClustering.py -i /eos/experiment/fcc/hh/simulation/samples/tracker_calobarrel/diW/500GeV/NTUP -n 1 -o output/test --njobs 10 -q 1nh
```
and:

``` 
python submitJetClustering.py -o output/test --collect
```
