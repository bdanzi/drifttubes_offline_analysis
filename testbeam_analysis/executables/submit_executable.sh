source /lustrehome/bdanzi/cmsset_default.sh
cd /lustrehome/bdanzi/latest/CMSSW_10_2_22/src/HiggsAnalysis/HiggsToZZ4Leptons/test/
eval `scramv1 runtime -sh`
cd /lustrehome/bdanzi/offline_analysis/testbeam_analysis/
source setDCDataReaderEnv.sh
root -l -b ReadSglChannel.C -q
root -l -b hist.c -q
