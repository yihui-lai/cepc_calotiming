path="/data/users/yihuilai/dualReadout/cepc_calotiming"
number="1000"
Particle="electron" #pion"
pa="e-" #pi-"

echo "" > submit_${Particle}.sh
echo "the particle is $Particle ($pa), the number is $number /n"
for energy in 1 2 5 10 20 50 100
do
  source temp_mac.sh $pa $energy $number > run_${energy}GeV_${Particle}.mac
  source temp_jdl.sh $path ${Particle} $energy > condor-jobs_${energy}GeV_${Particle}.jdl
  source temp_condor_sh.sh $path ${Particle} $energy > condor-executable_${energy}GeV_${Particle}.sh
  chmod 777 condor-executable_${energy}GeV_${Particle}.sh
  echo "condor_submit condor-jobs_${energy}GeV_${Particle}.jdl" >> submit_${Particle}.sh
done


