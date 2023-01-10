#create sh script for histos root file creation
pathf="/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis"
path="/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/executables"
runinputsevents_1_gsample=("89" "90" "91" "92") # New TestBeam
runinputsevents_2_gsample=("10" "39" "41" "63") # New TestBeam
runinputsevents_1p5_gsample=("86" "87" "88" "93" "94" "95" "96" "97" "98" "99" "127" "117") # New TestBeam
runinputsevents_1p2_gsample=("91" "94" "99" "98" "96") # Old TestBeam
dim=1024
binTimeInterval=25.0

for N_1 in 4
do
    for N_2 in 0
    do
        for N_3 in 1.0
        do
            for N_4 in 0.5
            do
            for gsample in 3
                do
                    runinputs=()
                    if [ $gsample -eq 0 ]
                    then
                	sampling=1.5
                        runinputs=("${runinputsevents_1p5_gsample[@]}")
                    fi
                    if [ $gsample -eq 1 ]
                    then
                	sampling=1
                        runinputs=("${runinputsevents_1_gsample[@]}")
                    fi

                    if [ $gsample -eq 2 ]
                    then
                	sampling=2
                        runinputs=("${runinputsevents_2_gsample[@]}")
                    fi

                    if [ $gsample -eq 3 ]
                    then
                	sampling=1.2
                        runinputs=("${runinputsevents_1p2_gsample[@]}")
                    fi

                    for i in ${runinputs[@]}
                    do
                	    rm ${path}/log$i
                        rm ${path}/test$i.cfg 
                        rm ${path}/test_conversion$i.cfg
                        rm ${path}/submit_executable_conversion$i.sh
                        rm ${path}/submitted_recas_conversion$i.log
                        rm ${path}/submitted_recas_conversion$i.error
                        rm ${path}/submitted_recas_conversion$i.out
                        echo '#!/bin/bash' >> ${path}/submit_executable_conversion$i.sh
                        echo "source /lustrehome/bdanzi/cmsset_default.sh" >> ${path}/submit_executable_conversion$i.sh
                        echo "cd ${pathf}/CMSSW_10_6_22/src/" >> ${path}/submit_executable_conversion$i.sh
                        echo 'eval `scramv1 runtime -sh`' >> ${path}/submit_executable_conversion$i.sh
                        echo "cd ${pathf}/" >> ${path}/submit_executable_conversion$i.sh
                        echo "source setDCDataReaderEnv.sh" >> ${path}/submit_executable_conversion$i.sh
                	    echo "./read_data . $i 0 -1 $sampling $N_1 $N_2 $N_3 $N_4 $binTimeInterval $dim >&  ${path}/log$i &" >> ${path}/submit_executable_conversion$i.sh
                        chmod +x ${path}/submit_executable_conversion$i.sh


                        #creates the cfg file to call the sh
                        echo "universe = vanilla" > ${path}/test_conversion$i.cfg
                        echo "Executable = ${path}/submit_executable_conversion$i.sh" >> ${path}/test_conversion$i.cfg   
                        echo "RequestMemory = 500000" >> ${path}/test_conversion$i.cfg
                	    echo "RequestDisk = 500000" >> ${path}/test_conversion$i.cfg
                        echo "WhenToTransferOutput = ON_EXIT" >> ${path}/test_conversion$i.cfg  
                        echo 'Requirements = TARGET.OpSys == "LINUX" && (TARGET.Arch != "DUMMY" )' >> ${path}/test_conversion$i.cfg  
                        echo "Log = ${path}/submitted_recas_conversion$i.log" >> ${path}/test_conversion$i.cfg  
                        echo "Output = ${path}/submitted_recas_conversion$i.out" >> ${path}/test_conversion$i.cfg  
                        echo "Error = ${path}/submitted_recas_conversion$i.error" >>${path}/test_conversion$i.cfg  
                        echo "Queue " >> ${path}/test_conversion$i.cfg  
                        ###then submit the job
                        echo "condor_submit -name ettore test_conversion$i.cfg"
                        cd ${path}
                       # condor_submit -name ettore test_conversion$i.cfg
                       	bash submit_executable_conversion$i.sh 
                    done
                done
            done
        done
    done
done
echo "DONE!"
