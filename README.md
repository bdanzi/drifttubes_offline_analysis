<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#Drift-Tubes-2021-Test-Beam-offline-analysis">Drift Tubes 2021 Test Beam offline analysis</a>
      <ul>
        <li><a href="#setup">Setup</a></li>
         <li> <a href="#instructions">Instructions</a></li>
        <li><a href="#channels-correspondance">Channels correspondance</a></li>
       
      </ul>
    </li>
  </ol>
</details>

# Drift Tubes 2021 Test Beam offline analysis

The `data_testbeam.xlsx` file contains the details of the different configurations associated to the ROOT files used for this data analysis.
The focus is on the 165 GeV/c momentum muon beam runs from 11st of November. 

## Setup

On Bari ReCAS and in the `testbeam_analysis\` directory of this repository:

```bash
$ source /cvmfs/sft.cern.ch/lcg/views/LCG_98python3/x86_64-centos7-gcc10-opt/setup.sh
$ source setDCDataReaderEnv.sh
$ bash compile.sh
$ ./read_data . 4 0 -10 1
```
On lxplus and in the `testbeam_analysis\` directory of this repository:

```bash
$ source setDCDataReaderEnv.sh
$ bash compile.sh
$ ./read_data . 4 0 -10 1
```
where in the last code line:

- 4 is the run number

- 0 -10 is the number of events to be processed

- 1 is a kind of option that can be fixed.

## Instructions

These macros find voltage amplitude peaks without any filter algorithm is applied to the waveform.
For each sample and each channel it is able to count how many events with an actual signal we have.
Config files and executables are created to run on more than one ROOT file (not available here, too much large in size).

```bash
$ bash submit_root_to_histos_root.sh
```

It will produce in `executables\`:
- executable files `submit_executable_conversion*.sh` per each `run_*.root` file that has to be converted in a root file containing histos 
with the most important physical variables
- config files that can be run as job in recas for accelerating the process of the histos ROOT file generation
```bash
$ bash submit_recas.sh 
```
It will produce in `executables\`::
- executable file `submit_executable.sh`
- config file that can be run as job in Recas
```bash
$ bash submit_executable.sh
```

It will produce:
- by using the plots.txt, per each .root file in the first column, per each channel in the second column,
per each event in the third column some physical quantities which are related to 
1) number of peaks, First Time Peak and Last Time Peak distributions
2) maxima
3) integral of the wavefunction 
4) wave functions w & w/o peak arrows
5) minima
6) number of events per each channel which passes a voltage amplitude requirement for the waveform maximum of 5 mV


## Channels correspondance

- 0,1,2,3 : Trigger Counters
- 4,5,6,7,8,9 : 6 Drift Tubes of 1 cm cell size respectively:
  - Channel 4 Wire diameter of 10 micrometer 
  - Channel 5 Wire diameter of 15 micrometer 
  - Channel 6 and 7 Wire diameter of 20 micrometer 
  - Channel 8 and 9 Wire diameter of 25 micrometer 
- 10,11,12 : 3 Drift Tubes of 2 cm cell size respectively:
  - Channel 10 Wire diameter of 20 micrometer 
  - Channel 11 Wire diameter of 25 micrometer 
  - Channel 12 Wire diameter of 40 micrometer 
- 13,14 : 2 Drift Tubes of 3 cm cell size respectively:
  - Channel 13 Wire diameter of 25 micrometer 
  - Channel 14 Wire diameter of 40 micrometer 

<img width="964" alt="Channel Schematics" src="https://github.com/bdanzi/drifttubes_offline_analysis/blob/master/Schermata%202022-01-30%20alle%2015.42.17.png">
Credits for the Channel Schematics: Franco Grancagnolo (https://agenda.infn.it/event/28676/)


 




