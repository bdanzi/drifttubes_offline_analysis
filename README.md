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
The focus is on the 165 GeV/c momentum muon beam runs from 11st of November. A flag isNov2021TestBeam is used to set the condition of this beam test. IsJuly2022TestBeam is used to analyze data of July 2022 test beam summarized in `run_table.xlsx`.

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
$ bash submit_first.sh
```
On lxplus and in the `testbeam_analysis\` directory of this repository:

```bash
$ source setDCDataReaderEnv.sh
$ bash compile.sh
$ bash submit_first.sh
```
where the last line of code runs in parallel the `./read_data $path_to_your_run_files $i 0 -1 $sampling $N_1 $N_2 $N_3 $N_4 $binTimeInterval $dim $cut_scale` command:

- i is the run number

- 0 -10 are the number of events to be processed (-1 to process the whole data set)

- $sampling is the sampling rate to be used in analysis

- $N_1 $N_2 $N_3 $N_4 are the cuts respectively on the signal amplitude, first and second derivatives before starting the search for electron peaks (see FindPeak-algo.C for further details)

- $binTimeInterval is the time interval (ns) starting from the beginning of the waveform in which the baseline and the rms is computed

- $dim are the waveform bins available from the data acquisition window

- $cut_scale is the threshold for electron peaks to be in the same primary ionization cluster

## Instructions

These macros find voltage amplitude peaks without any filter algorithm is applied to the waveform.
For each sample and each channel it is able to count how many events with an actual signal we have.
Config files and executables are created to run on more than one ROOT file (not available here, too much large in size).

```bash
$ submit_second.sh
```

It will produce in `executables\`:
- executable files `submit_executable_conversion*.sh` per each `run_*.root` file that has to be converted in a root file containing histos 
with the most important physical variables
- config files that can be run as job in recas for accelerating the process of the histos ROOT file generation

Moreover, it will produce:
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

<img width="964" alt="Channel Schematics of November 2021 Beam Test" src="https://github.com/bdanzi/drifttubes_offline_analysis/blob/master/Schermata%202022-01-30%20alle%2015.42.17.png">
Credits for the Channel Schematics: Franco Grancagnolo (https://agenda.infn.it/event/28676/)

<img width="964" alt="Channel Schematics of July 2022 Beam Test" src="https://github.com/bdanzi/TestBeam2022/blob/main/Schermata%202022-09-10%20alle%2020.18.11.png">
Credits for the Channel Schematics: Federica Cuna (https://indico.ihep.ac.cn/event/16837/)

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


 




