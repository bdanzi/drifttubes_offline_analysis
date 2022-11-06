<!-- PROJECT SHIELDS -->
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![MIT][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#Drift-Tubes-2021-Test-Beam-offline-analysis">Drift Tubes 2021 Test Beam offline analysis</a>
      <ul>
        <li><a href="#authors">Authors</a></li>
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

## Authors

- [Federica Cuna](https://github.com/federicacuna) (University and INFN Lecce)
- [Brunella D'Anzi](https://github.com/bdanzi) (University and INFN Bari)
- [Nicola De Filippis](https://github.com/ndefilip) (Politecnico and INFN Bari)
- [Walaa Elmentenawee](https://github.com/elmetenawee) (University and INFN Bari)
- [Matteo Greco](https://github.com/matteogreco2) (University and INFN Lecce)


## Setup

On Bari ReCAS and in the `testbeam_analysis\` directory of this repository:

```bash
$ source /cvmfs/sft.cern.ch/lcg/views/LCG_98python3/x86_64-centos7-gcc10-opt/setup.sh
$ source setDCDataReaderEnv.sh
$ bash compile.sh
$ ./read_data . 91 0 -1 1.2 4 0 1.0 0.5 25.0 1024
```
On lxplus and in the `testbeam_analysis\` directory of this repository:

```bash
$ source setDCDataReaderEnv.sh
$ bash compile.sh
$ ./read_data . 91 0 -1 1.2 4 0 1.0 0.5 25.0 1024
```
where in the last code line:

- 91 is the run number

- 0 -10 is the number of events to be processed (-1 to process the whole data set)

- 1.2 is the sampling rate

- 4, 0, 1.0, 0.5 are the cuts on the signal amplitude, first and second derivatives before starting the search for electron peaks (see FindPeak-algo.C for further details)

- 25.0 is the time interval (ns) starting from the beginning of the waveform in which the baseline and the rms is computed

- 1024 are the waveform bins available

## Instructions

These macros find voltage amplitude peaks without any filter algorithm is applied to the waveform.
For each sample and each channel it is able to count how many events with an actual signal we have.
Config files and executables are created to run on more than one ROOT file (not available here, too much large in size).

```bash
$ submit_root_to_histos_root_22test_variableCuts.sh
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
$ submit_executable_oldTestBeam.sh
```

It will produce:
- by using the plots_oldTestBeam.txt, per each .root file in the first column, per each channel in the second column,
per each event in the third column some physical quantities which are related to 
1) Number of Electron peaks, First Peak Time of Arrival and Last Peak Time of arrival distributions
2) Maxima
3) Charge Integral of the wavefunction (in pC)
4) Wave functions w & w/o peak arrows
5) Number of events per each channel which passes a voltage amplitude requirement for the waveform maximum of 5 mV
6) Number of Cluster distributions
7) Time Difference between two consecutive clusters
8) Time Difference between two consecutive Electrons
9) Cluster population (Number of Electron Peaks per Primary Ionization Cluster)

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

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/bdanzi/drifttubes_offline_analysis.svg?style=for-the-badge
[contributors-url]: https://github.com/bdanzi/drifttubes_offline_analysis/contributors

[forks-shield]: https://img.shields.io/github/forks/bdanzi/ML_COURSEBARI.svg?style=for-the-badge
[forks-url]: https://github.com/bdanzi/drifttubes_offline_analysis/network/members

[stars-shield]: https://img.shields.io/github/stars/bdanzi/ML_COURSEBARI.svg?style=for-the-badge
[stars-url]: https://github.com/bdanzi/drifttubes_offline_analysis/stargazers

[issues-shield]: https://img.shields.io/github/issues/bdanzi/drifttubes_offline_analysis.svg?style=for-the-badge
[issues-url]: https://github.com/bdanzi/drifttubes_offline_analysis/issues

[license-shield]: https://img.shields.io/github/license/bdanzi/drifttubes_offline_analysis.svg?style=for-the-badge
[license-url]: https://github.com/bdanzi/drifttubes_offline_analysis/blob/main/LICENSE.txt

[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/brunella-d-anzi


 




