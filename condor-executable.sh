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
source ./g4env.sh


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
