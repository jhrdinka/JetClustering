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
./analyze /eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/RunLHE_selvaggi_TTbarGunPt500_0PU_20171014/NTUP/RunLHE_x100_NTUP_20.root test.root 20
```
