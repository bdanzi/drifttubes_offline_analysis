#!/bin/bash

#export WAVEDATA_DIR=WaveData
#export LD_LIBRARY_PATH=$WAVEDATA_DIR:${WAVEDATA_DIR}/lib:$LD_LIBRARY_PATH
APPOUT="read_data"
APPOUT1=''
APPOUT2=''
COPT=''
if [ $1 ]; then
  COPT=" -D_OSC"
  APPOUT="read_osc_data"
  APPOUT1="calbElec"
  APPOUT2="simulatedSign"
  APPOUT3="fcTranf"
fi
  echo "--------/---/-/- $APPOUT "
g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${WAVEDATA_DIR}/.. $COPT compile_read_data.C read_data.C  `root-config --glibs` `root-config --libs` `root-config --cflags` -l Spectrum -l Minuit -L ${WAVEDATA_DIR}/lib -l WaveData -o $APPOUT

if [ $APPOUT1 ]; then
  echo "--------/----- $APPOUT1 "
  g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${WAVEDATA_DIR}/.. $COPT compile_read_data.C read_data-EL.C  `root-config --glibs` `root-config --libs` `root-config --cflags` -l Spectrum -l Minuit -L ${WAVEDATA_DIR}/lib -l WaveData -o $APPOUT1
fi

if [ $APPOUT3 ]; then
  echo "---///////----- $APPOUT3 "
  g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${WAVEDATA_DIR}/.. $COPT  compile_read_data.C read_data-FT.C  `root-config --glibs` `root-config --libs` `root-config --cflags` -l Spectrum -l Minuit -L ${WAVEDATA_DIR}/lib -l WaveData -o $APPOUT3
fi

if [ $APPOUT2 ]; then
  echo "-----/////----- $APPOUT2 "
  g++ -I $ROOTSYS/include -I $ROOFITSYS/include $COPT compile_simSign.C  `root-config --glibs` `root-config --libs` `root-config --cflags` -l Spectrum -l Minuit -o $APPOUT2
fi


