#create sh script for histos root file creation
pathf="/lustrehome/bdanzi/offline_analysis/testbeam_analysis"
path="/lustrehome/bdanzi/offline_analysis/testbeam_analysis/executables"
runinputs10000events=("89" "90" "91" "92")
runinputs5000events=("86" "87" "88" "93" "94" "95" "96" "97" "98" "99" "127" "117")
#runinputs10000events=("1" "2" "3" "4" "5" "6" "7" "8" "21" "64" "89" "90" "91" "92" "100" "101")
#runinputs5000events=("23" "24" "26" "27" "28" "29" "30" "31" "32" "33" "34" "35" "36" "37" "38" "42" "43" "45" "46" "47" "48" "49" "50" "51" "52" "53" "54" "55" "56" "57" "58" "59" "60" "61" "62" "63" "65" "66" "67" "68" "69" "70" "71" "72" "73" "74" "75" "79" "80" "81" "82" "83" "84" "85" "86" "87" "88" "93" "94" "95" "96" "97" "98" "99" "102" "103" "104" "105" "106" "107" "108" "109" "110" "111" "112" "115" "116" "117" "118" "119" "120" "121" "122" "123" "124" "125" "126" "127" "128" "129" "130" "131" "132" "133" "142")
for nevents in 10000 5000
do
    runinputs=()
    if [ $nevents -eq 10000 ]
    then
        runinputs=("${runinputs10000events[@]}")
    fi
    if [ $nevents -eq 5000 ]
    then
        runinputs=("${runinputs5000events[@]}")
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
        echo "cd /lustrehome/bdanzi/latest/CMSSW_10_2_22/src/" >> ${path}/submit_executable_conversion$i.sh
        echo 'eval `scramv1 runtime -sh`' >> ${path}/submit_executable_conversion$i.sh
        echo "cd ${pathf}/" >> ${path}/submit_executable_conversion$i.sh
        echo "source setDCDataReaderEnv.sh" >> ${path}/submit_executable_conversion$i.sh
        echo "./read_data . $i 0 -1 1>&  ${path}/log$i &" >> ${path}/submit_executable_conversion$i.sh
        chmod +x ${path}/submit_executable_conversion$i.sh
	
	
        #creates the cfg file to call the sh
        echo "universe = vanilla" > ${path}/test_conversion$i.cfg
        echo "Executable = ${path}/submit_executable_conversion$i.sh" >> ${path}/test_conversion$i.cfg   
        echo "RequestMemory = 10000" >> ${path}/test_conversion$i.cfg  
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
echo "DONE!"
