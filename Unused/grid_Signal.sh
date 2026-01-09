#!/usr/bin/env bash
# Prompt the user to enter a value for date
echo -n "Enter the Date: "
read date
echo -n "Attempt: "
read Attempt
read -p "Signal (s), Download Signal (ds), Merge Signal (ms): " op

echo -e "\nDate Choosen: $date"
echo -e "Attempt: $Attempt"
echo -e "Option Choosen: $op\n"

if [ "$op" = "s" ]; then

gbasf2 /home/belle2/amubarak/Ds2D0enue_Analysis/KEK_Server/Steering_Scripts/Signal_Grid/Dsp-Reconstruction.py -p Dsp_${date}24_${Attempt}_Signal -s light-2406-ragdoll -i /belle/MC/release-06-00-08/DB00002100/MC15ri_b/prod00046485/s00/e1003/4S/r00000/2653230000/mdst --force
echo -e "\nD_s+ Project Submitted\n"
gbasf2 /home/belle2/amubarak/Ds2D0enue_Analysis/KEK_Server/Steering_Scripts/Signal_Grid/Dsm-Reconstruction.py -p Dsm_${date}24_${Attempt}_Signal -s light-2406-ragdoll -i /belle/MC/release-06-00-08/DB00002100/MC15ri_b/prod00046486/s00/e1003/4S/r00000/2653230001/mdst --force
echo -e "\nD_s- Project Submitted\n"

elif [ "$op" = "ds" ]; then

echo -e "Grid Signal Download"

(nohup gb2_ds_get --force Dsp_${date}24_${Attempt}_Signal) > C03-Grid/Terminal_Output/Dsp_${date}24_${Attempt}_Signal.out & #If things slow down add "--new"
(nohup gb2_ds_get --force Dsm_${date}24_${Attempt}_Signal) > C03-Grid/Terminal_Output/Dsm_${date}24_${Attempt}_Signal.out & #If things slow down add "--new"
echo -e "Download in Progress"
# sleep 2m # Waits 6 minutes.
jobs

elif [ "$op" = "ms" ]; then

echo -e "Merge Signal File Start"
rm C01-Simulated_Events/Ds2D0enu-Signal.root -f
hadd -f C01-Simulated_Events/Ds2D0enu-Signal.root Dsp_${date}24_${Attempt}_Signal/sub00/*.root Dsm_${date}24_${Attempt}_Signal/sub00/*.root
rm -r Dsp_${date}24_${Attempt}_Signal/ Dsm_${date}24_${Attempt}_Signal/ -f
rm C03-Grid/Terminal_Output/Dsp_${date}24_${Attempt}_Signal.out C03-Grid/Terminal_Output/Dsm_${date}24_${Attempt}_Signal.out -f 

fi

rm /home/belle2/amubarak/Ds2D0enue_Analysis/KEK_Server/Sources/grid_Signal.sh -f
cp -f grid_Signal.sh Ds2D0enue_Analysis/KEK_Server/Sources
echo -e "\nDone\n"

# Grid:
# gb2_load
# Password: .BgPPHs4dKc22_G
# source /cvmfs/belle.kek.jp/grid/gbasf2/pro/bashrc