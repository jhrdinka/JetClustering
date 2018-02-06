First load the environment
```
source init.sh
```
Install fastjet (to be done only once)
```
./install-fastjet.sh
```
Compile analysis code:
```
make
```
Run analysis:
```
./analyze [input_file] [output_file] []number_of_events] [CMS/FCC]
```

Ex:
```
./analyze /eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/RunLHE_selvaggi_TTbarGunPt500_0PU_20171014/NTUP/RunLHE_x100_NTUP_20.root test.root 20
```
