unset LD_LIBRARY_PATH
unset PYTHONHOME
unset PYTHONPATH

export PATH="/afs/cern.ch/work/s/selvaggi/private/fastjet-3.3.0/fastjet-install-cgal/bin:$PATH"
source /cvmfs/sft.cern.ch/lcg/releases/LCG_88/gcc/4.9.1/x86_64-slc6/setup.sh 
source /cvmfs/sft.cern.ch/lcg/releases/LCG_88/Boost/1.62.0/x86_64-slc6-gcc49-opt/Boost-env.sh 
source /cvmfs/sft.cern.ch/lcg/releases/LCG_88/Python/2.7.13/x86_64-slc6-gcc49-opt/Python-env.sh 
source /cvmfs/sft.cern.ch/lcg/releases/LCG_88/ROOT/6.08.06/x86_64-slc6-gcc49-opt/ROOT-env.sh; 

INPUT_FILE=${1}
OUTPUT_FILE=${2}
NEVENTS=${3}
RUNTYPE=${4}

cp /afs/cern.ch/user/c/cneubuse/FCC/JetClustering/analyze .

./analyze ${INPUT_FILE} ${OUTPUT_FILE} ${NEVENTS} ${RUNTYPE}
 
