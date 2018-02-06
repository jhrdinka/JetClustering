First load the required environment on lxplus:
```
source init.sh
```
Install fastjet (to be done only once):
```
./install-fastjet.sh
```
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
