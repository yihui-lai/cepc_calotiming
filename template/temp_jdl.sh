#!/bin/sh

#define parameters which are passed in.
path=$1
particle=$2
energy=$3

#define the template.
cat  << EOF
universe = vanilla
Executable = condor-executable_${energy}GeV_${particle}.sh
should_transfer_files = NO
Requirements = TARGET.FileSystemDomain == "privnet"
Output = ${path}/condor_output/sce_\$(cluster)_\$(process).stdout
Error =  ${path}/condor_output/sce_\$(cluster)_\$(process).stderr
Log =    ${path}/condor_output/sce_\$(cluster)_\$(process).condor
Arguments = SCE
Queue 1

EOF
