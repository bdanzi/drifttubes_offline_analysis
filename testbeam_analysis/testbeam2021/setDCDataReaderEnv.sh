#source /home/federica/Utils/root_v6.22/root/bin/thisroot.sh
export WAVEDATA_DIR=/eos/user/g/grecom/TestBeamData_Nov2021/testbeam2021/WaveData
export WAVEDATA_HEADERS=${WAVEDATA_DIR}/inc
export WAVEDATA_LINKDEF=${WAVEDATA_DIR}/linkdef_inc
export WAVEDATA_LIB=${WAVEDATA_DIR}/lib/libWaveData.so
export WAVEDATA_PCM=${WAVEDATA_DIR}/WaveDataDict_rdict.pcm
export LD_LIBRARY_PATH=$WAVEDATA_DIR:${WAVEDATA_DIR}/lib:$LD_LIBRARY_PATH

#export LD_LIBRARY_PATH=$WAVEDATA_DIR:${WAVEDATA_DIR}:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=${ROOTSYS}/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/opt/root-6.14.02/lib:$LD_LIBRARY_PATH
