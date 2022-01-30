# Drift Tubes 2021 Test Beam offline analysis
## Setup

On lxplus/Bari ReCAS:

```
source /cvmfs/sft.cern.ch/lcg/views/LCG_98python3/x86_64-centos7-gcc10-opt/setup.sh
source /cvmfs/cms.cern.ch/slc7_amd64_gcc700/lcg/root/6.12.07-gnimlf7/bin/thisroot.sh 
source setDCDataReaderEnv.sh
```

## Instructions

These macros find voltage amplitude peaks without any filter algorithm is applied to the waveform.
For each sample and each channel it is able to count how many events with an actual signal we have.
Config files and executables are created to run on more than one ROOT file (not available here, too much large in size).

```
bash submit_root_to_histos_root.sh
```

It will produce:
- executable files `submit_executable_conversion*.sh` per each `run_*.root` file that has to be converted in a root file containing histos 
with the most important physical variables
- config files that can be run as job in recas for accelerating the process of the histos ROOT file generation
```
bash submit_recas.sh 
```
It will produce:
- executable file `submit_executable.sh`
- config file that can be run as job in Recas
```
bash submit_executable.sh
```

It will produce:
- by using the plots.txt, per each .root file in the first column, per each channel in the second column,
per each event in the third column some physical quantities which are related to 
1) number of peaks distributions
2) maxima
3) integral of the wavefunction 
4) wave functions w & w/o peak arrows
5) minima
6) number of events per each channel which passes a voltage amplitude requirement for the waveform maximum of 5 mV




