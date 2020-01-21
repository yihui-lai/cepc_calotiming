#!/bin/bash

cd /data/users/eno/dualReadout/cepc_calotiming/
export UNIQUE_ID=sce
export CONDOR_PROCESS=sce
export RUN_DIR=`pwd`

#
# header 
#


START_TIME=`/bin/date`
echo "started at $START_TIME"

echo ""
echo "parameter set:"
echo "UNIQUE_ID: $UNIQUE_ID"
echo "CONDOR_PROCESS: $CONDOR_PROCESS"
echo "RUN_DIR: $RUN_DIR"
echo "PARAMETER_SET: $PARAMETER_SET"


#
# setup software environment at UMD
#
source /cvmfs/sft.cern.ch/lcg/contrib/gcc/6.3/x86_64-slc6/setup.sh
source /cvmfs/geant4.cern.ch/geant4/10.5/x86_64-slc6-gcc63-opt-MT/CMake-setup.sh
export CXX=`which g++`
export CC=`which gcc`
export PATH=$PATH:/cvmfs/sft.cern.ch/lcg/contrib/CMake/3.11.1/Linux-x86_64/bin
source /cvmfs/sft.cern.ch/lcg/views/LCG_95/x86_64-slc6-gcc62-opt/setup.sh
source $ROOTSYS/bin/thisroot.sh



#
# run 
#


 ./CEPC_CaloTiming -c template.cfg -m run.mac -o testtest > run.log

exitcode=$?


#
# end run
#

echo ""
END_TIME=`/bin/date`
echo "finished at $END_TIME"
exit $exitcode
