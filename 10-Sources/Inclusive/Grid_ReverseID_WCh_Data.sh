#!/usr/bin/env bash
# Prompt the user to enter a value for date
echo -e "Reversal electronID Cut for Data"
echo -n "Enter the Date: "
read date
echo -n "Attempt: "
read Attempt
read -p "Data (p), Download Data (dd), Merge Data (md) : " op

echo -e "\nDate Choosen: $date"
echo -e "Attempt: $Attempt"
echo -e "Option Choosen: $op\n"

if [ "$op" = "p" ]; then

gbasf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Data/Ds2D0e-Reconstruction_ReverseID_WCh.py -p Ds_${date}25_${Attempt}_Data_ReverseID -s light-2503-ceres -i /belle/collection/Data/RSC_proc13prompt_4S_skim_17230100_v3 --force
echo -e "\nData Project Submitted\n"

p

elif [ "$op" = "dd" ]; then

echo -e "Download Data"

(nohup gb2_ds_get --force --new Ds_${date}25_${Attempt}_Data_ReverseID) > /home/belle2/amubarak/C02-Grid/Terminal_Output/Ds_${date}25_${Attempt}_Data_ReverseID.out & #If things slow down add "--new"
echo -e "Download in Progress"

elif [ "$op" = "md" ]; then

echo -e "Merge Data Start"
hadd -f /group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh/Ds2D0e-Generic_Ds_${date}25_${Attempt}_Data_ReverseID_WCh.root Ds_${date}25_${Attempt}_Data_ReverseID/sub00/*.root
# rm -r Ds_${date}25_${Attempt}_ccbar/ -f
# rm /home/belle2/amubarak/C03-Grid/Terminal_Output/Ds_${date}25_${Attempt}_ccbar.out -f

fi

# rm /home/belle2/amubarak/Ds2D0enue_Analysis/09-Sources/Inclusive/Grid_ReverseID_Data.sh -f
# cp -f /home/belle2/amubarak/Grid_ReverseID_Data.sh /home/belle2/amubarak/Ds2D0enue_Analysis/09-Sources/Inclusive
echo -e "\nDone\n"

# Grid:
# gb2
# Password: AzerDSRT4116!
# source /cvmfs/belle.kek.jp/grid/gbasf2/pro/bashrc