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

Compile analysis code:
```
make -j 4
```
Run analysis:
```
./analyze [input_file] [output_file] [number_of_events] [CMS/FCC]
```
e.g:
```
./analyze /eos/experiment/fcc/hh/simulation/samples/tracker_calobarrel/ditop/500GeV/NTUP/output_helsens_20171011151211690.root test.root 10 FCC```
```


[]() LSF submission
--------------------

The script ```batch/submitJetClustering.py``` allows to run this script on LSF queues:

```
python batch/submitJetClustering.py -i [NTUP_dir] -n [nevts_per_job] -o [output_dir] --njobs [number_of_jobs] -q [queue]
```

Jobs can be collected via:


```
python batch/submitJetClustering.py -o [output_dir] --collect
```

e.g:
```
python batch/submitJetClustering.py -i /eos/experiment/fcc/hh/simulation/samples/tracker_calobarrel/diW/500GeV/NTUP -n 1 -o output/test --njobs 10 -q 1nh
```
and:

``` 
python batch/submitJetClustering.py -o output/test --collect
```
