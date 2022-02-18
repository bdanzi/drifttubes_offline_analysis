#create sh script for training
path="/lustrehome/bdanzi/offline_analysis/testbeam_analysis/executables"
pathf="/lustrehome/bdanzi/offline_analysis/testbeam_analysis"
rm ${path}/test.cfg
rm ${path}/submit_executable.sh
rm ${path}/submitted_recas.log
rm ${path}/submitted_recas.error
rm ${path}/submitted_recas.out
echo '#!/bin/bash' >> /lustrehome/bdanzi/offline_analysis/testbeam_analysis/submit_executable.sh
echo "source /lustrehome/bdanzi/cmsset_default.sh" >> ${path}/submit_executable.sh
echo "cd /lustrehome/bdanzi/latest/CMSSW_10_2_22/src/HiggsAnalysis/HiggsToZZ4Leptons/test/" >> ${path}/submit_executable.sh
echo 'eval `scramv1 runtime -sh`' >> ${path}/submit_executable.sh
echo "cd ${pathf}/" >> ${path}/submit_executable.sh
echo "source setDCDataReaderEnv.sh" >> ${path}/submit_executable.sh
echo "root -l -b ReadSglChannel.C -q" >> ${path}/submit_executable.sh
echo "root -l -b hist.c -q" >> ${path}/submit_executable.sh
chmod +x ${path}/submit_executable.sh


#creates the cfg file to call the sh
echo "universe = vanilla" > ${path}/test.cfg
echo "Executable = ${path}/submit_executable.sh" >> ${path}/test.cfg                      
echo "RequestMemory = 3000" >> ${path}/test.cfg
echo "WhenToTransferOutput = ON_EXIT" >> ${path}/test.cfg
echo 'Requirements = TARGET.OpSys == "LINUX" && (TARGET.Arch != "DUMMY" )' >> ${path}/test.cfg
echo "Log = ${path}/submitted_recas.log" >> ${path}/test.cfg
echo "Output = ${path}/submitted_recas.out" >> ${path}/test.cfg
echo "Error = ${path}/submitted_recas.error" >>${path}/test.cfg
echo "Queue " >> ${path}/test.cfg
###then submit the job
echo "condor_submit -name ettore test.cfg"
cd ${path}
#condor_submit -name ettore test.cfg


