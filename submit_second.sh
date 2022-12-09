#!/bin/bash
source /lustrehome/bdanzi/cmsset_default.sh
cd CMSSW_10_6_22/src/
eval `scramv1 runtime -sh`
cd ../..
source setDCDataReaderEnv.sh
root -l -b ReadSglChannel_test_NovTestBeam.C -q
#root -l -b hist.c -q
#!/bin/bash
