#create sh script for histos root file creation
pathf="/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis"
path="/lustrehome/bdanzi/TestBeam2022/analysistestbeam2022/analysis_2021/drifttubes_offline_analysis/testbeam_analysis/executables"
#runinputsevents_1_gsample=() # Old TestBeam 
# runinputsevents_2_gsample=("10" "39" "41" "63") # 2.0 Gsample July 2022 Test Beam 45° run_10 90/10 80/20 merged from 39 and 41, run_63 85/15
#runinputsevents_2_gsample=("40" "41" "42" "43") # Study July 2022 Test Beam Runs with different HV : Nom -10V, Nom, Nom +10V, Nom +20V
runinputsevents_2_gsample=("10") #Study July 2022 Test Beam Runs with different Sampling rate: 1.5 2.0#
runinputsevents_1p5_gsample=("9") #Study July 2022 Test Beam Runs with different Sampling rate: 1.5 2.0
#runinputsevents_2_gsample=("44" "45" "41" "46") # Study July 2022 Test Beam Scan in Angle 0° 30° 45° 60°
#runinputsevents_2_gsample=("41" "63" "10") # Study July 2022 Test Beam Different Gas Mixtures 80/20 85/15 90/10
#runinputsevents_1p5_gsample=("9" ) # 1.5 Gsample July 2022 Test Beam
#runinputsevents_1p2_gsample=("99" "98")
runinputsevents_1p2_gsample=("99" "98" "96" "94" "91") # Nov 2021 TestBeam 1.2 Gsample 90/10 165 GeV 0° 15° 30° 45° 60°
#runinputsevents_1p2_gsample=("117" "127") # Nov 2021 TestBeam 1.2 Gsample 80/20 165 GeV 0° 15° 30° 45° 60°
#runinputsevents_1p2_gsample=("99") # Old TestBeam
dim=1024
binTimeInterval=25.0
rm plots_oldTestBeam.txt
# Original cuts for Nov 2021 cut_scale 0.26 N1 10.0 N_2 1.0 N_3 0.2 N_4 0.2 gsample 3
# Last cuts for July 2022 cut scale  0.2 N_1 10.0 N_2 1.2 N_3 0.0 N_4 0.0 gsample 0 2
for cut_scale in 0.26 #0.20 # 0.26 0.22 0.24 0.20 0.26 0.28 #0.1 0.4 0.5 #con N_2 1.0 # second derivative in the peak candidate
do
    for N_1 in 8.0 #6 # Amplitude cut
    do
        for N_2 in 1.1 #1.5 0.8 0.5 #1.0 amplitude on the average before and after the peak candidate
        do
            for N_3 in 0.1 #0.1 0.4 0.5 #con N_2 1.0 # first derivative in the peak candidate
            do
                for N_4 in 0.1 #0.1 0.4 0.5 #con N_2 1.0 # second derivative in the peak candidate
                do
                for gsample in 0 2
                    do
                        runinputs=()
                        if [ $gsample -eq 0 ]
                        then
                    	sampling=1.5
                            runinputs=("${runinputsevents_1p5_gsample[@]}")
                        fi
                        if [ $gsample -eq 1 ]
                        then
                    	sampling=1.0
                            runinputs=("${runinputsevents_1_gsample[@]}")
                        fi

                        if [ $gsample -eq 2 ]
                        then
                    	sampling=2.0
                            runinputs=("${runinputsevents_2_gsample[@]}")
                        fi

                        if [ $gsample -eq 3 ]
                        then
                    	sampling=1.2
                            runinputs=("${runinputsevents_1p2_gsample[@]}")
                        fi

                        for i in ${runinputs[@]}
                        do

                    	rm ${path}/log_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i
                        rm ${path}/test_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.cfg 
                        rm ${path}/test_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.cfg
                        rm ${path}/submit_executable_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.sh
                        rm ${path}/submitted_recas_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.log
                        rm ${path}/submitted_recas_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.error
                        rm ${path}/submitted_recas_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.out
                        echo '#!/bin/bash' >> ${path}/submit_executable_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.sh
                        echo "source /lustrehome/bdanzi/cmsset_default.sh" >> ${path}/submit_executable_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.sh
                        echo "cd ${pathf}/CMSSW_10_6_22/src/" >> ${path}/submit_executable_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.sh
                        echo 'eval `scramv1 runtime -sh`' >> ${path}/submit_executable_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.sh
                        echo "cd ${pathf}/" >> ${path}/submit_executable_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.sh
                        echo "source setDCDataReaderEnv.sh" >> ${path}/submit_executable_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.sh
                    	#echo "./read_data . $i 0 -1 $sampling $N_1 $N_2 $N_3 $N_4 $binTimeInterval $dim $cut_scale >&  ${path}/log_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i &" >> ${path}/submit_executable_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.sh
                        echo "./read_data /lustre/cms/store/user/bdanzi/TestBeam20212022Analysis $i 0 -1 $sampling $N_1 $N_2 $N_3 $N_4 $binTimeInterval $dim $cut_scale >&  ${path}/log_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i &" >> ${path}/submit_executable_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.sh
                        chmod +x ${path}/submit_executable_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.sh
                        echo "histosTB_run_${i}_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}.root ${N_1} ${N_2} ${N_3} ${N_4} ${cut_scale} ${sampling} histosTB_run_${i}.root 16 10000" >> ${pathf}/plots_oldTestBeam.txt

                        #creates the cfg file to call the sh
                        echo "universe = vanilla" > ${path}/test_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.cfg
                        echo "Executable = ${path}/submit_executable_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.sh" >> ${path}/test_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.cfg   
                        echo "RequestMemory = 500000" >> ${path}/test_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.cfg
                    	echo "RequestDisk = 500000" >> ${path}/test_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.cfg
                        echo "WhenToTransferOutput = ON_EXIT" >> ${path}/test_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.cfg  
                        echo 'Requirements = TARGET.OpSys == "LINUX" && (TARGET.Arch != "DUMMY" )' >> ${path}/test_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.cfg  
                        echo "Log = ${path}/submitted_recas_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_$i.log" >> ${path}/test_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.cfg  
                        echo "Output = ${path}/submitted_recas_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_$i.out" >> ${path}/test_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.cfg  
                        echo "Error = ${path}/submitted_recas_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_$i.error" >> ${path}/test_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.cfg  
                        echo "Queue " >> ${path}/test_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.cfg  
                        ###then submit the job
                        echo "condor_submit -name ettore test_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.cfg"
                        cd ${path}
                        condor_submit -name ettore test_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.cfg
                        bash submit_executable_conversion_N1_${N_1}_N2_${N_2}_N3_${N_3}_N4_${N_4}_cut_scale_${cut_scale}_sampling_${sampling}_$i.sh 
                        done
                    done
                done
            done
        done
    done
done
echo "DONE!"
