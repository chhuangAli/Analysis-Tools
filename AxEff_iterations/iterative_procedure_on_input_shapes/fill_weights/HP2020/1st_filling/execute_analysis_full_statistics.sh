#RUNNUMBERS=$(cat /home/chunlu/local/runList_LHC17p/runList.txt)
RUNNUMBERS=$(cat /home/chunlu/local/runList_LHC15o/runList.txt)
#RUNNUMBERS=$(cat /home/chunlu/local/runList_LHC18q_r/muon_calo_pass3/LHC18q_r_runList.txt)

now=$(date +"%T")
echo "Starting time : $now"
for RUN in $RUNNUMBERS; do
     echo $RUN
     ./GOrun $RUN
done
echo "Current time : $now"
