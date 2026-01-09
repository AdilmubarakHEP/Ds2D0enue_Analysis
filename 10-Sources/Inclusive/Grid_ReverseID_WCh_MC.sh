#!/usr/bin/env bash
# Prompt the user to enter a value for date
echo -e "Reversal electronID Cut for MC"
echo -n "Enter the Date: "
read date
echo -n "Attempt: "
read Attempt
read -p "Generic Charm Background (pcb), All Generic Background (pab), Download Generic Charm Background (dcb), Download All Generic Background (dab), Merge Generic Charm Background (mcb), Merge All Generic Background (mab) : " op

echo -e "\nDate Choosen: $date"
echo -e "Attempt: $Attempt"
echo -e "Option Choosen: $op\n"

if [ "$op" = "pcb" ]; then

# gbasf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction_ReverseID_WCh.py -p Ds_${date}25_${Attempt}_ccbar_ReverseID -s light-2503-ceres --input_dslist /home/belle2/amubarak/Ds2D0enue_Analysis/03-Grid/Inputs/grid_ccbar.txt --force
# echo -e "\nCharm Background Project Submitted\n"
gbasf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction_ReverseID_WCh.py -p Ds_${date}25_${Attempt}_ccbar_ReverseID -s light-2503-ceres -i /belle/collection/MC/RSC_MC15rd_4S_ccbar_skim_17230100_v3 --force
echo -e "\nCharm Background Project Submitted\n"

p

elif [ "$op" = "pab" ]; then

gbasf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction_ReverseID_WCh.py -p Ds_${date}25_${Attempt}_ccbar_ReverseID -s light-2503-ceres -i /belle/collection/MC/RSC_MC15rd_4S_ccbar_skim_17230100_v3 --force
echo -e "1). ccbar submitted\n"
gbasf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction_ReverseID_WCh.py -p Ds_${date}25_${Attempt}_BB_ReverseID -s light-2503-ceres -i /belle/collection/MC/RSC_MC15rd_4S_BB_skim_17230100_v3 --force
echo -e "2). BB submitted\n"
gbasf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction_ReverseID_WCh.py -p Ds_${date}25_${Attempt}_ddbar_ReverseID -s light-2503-ceres -i /belle/collection/MC/RSC_MC15rd_4S_ddbar_skim_17230100_v3 --force
echo -e "3). ddbar submitted\n"
# gbasf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction_ReverseID_WCh.py -p Ds_${date}25_${Attempt}_mixed_ReverseID -s light-2503-ceres -i /belle/collection/MC/MC15rd_mixed_exp20-26_4S_v2 --force
# echo -e "4). mixed submitted\n"
gbasf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction_ReverseID_WCh.py -p Ds_${date}25_${Attempt}_ssbar_ReverseID -s light-2503-ceres -i /belle/collection/MC/RSC_MC15rd_4S_ssbar_skim_17230100_v3 --force
echo -e "4). ssbar submitted\n"
gbasf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction_ReverseID_WCh.py -p Ds_${date}25_${Attempt}_taupair_ReverseID -s light-2503-ceres -i /belle/collection/MC/RSC_MC15rd_4S_taupair_skim_17230100_v3 --force
echo -e "5). taupair submitted\n"
gbasf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction_ReverseID_WCh.py -p Ds_${date}25_${Attempt}_uubar_ReverseID -s light-2503-ceres -i /belle/collection/MC/RSC_MC15rd_4S_uubar_skim_17230100_v3 --force
echo -e "6). uubar submitted\n"

elif [ "$op" = "dcb" ]; then

echo -e "Download Generic Charm Background"

(nohup gb2_ds_get --force --new Ds_${date}25_${Attempt}_ccbar_ReverseID) > /home/belle2/amubarak/C02-Grid/Terminal_Output/Ds_${date}25_${Attempt}_ccbar_ReverseID.out & #If things slow down add "--new"
echo -e "Download in Progress"

elif [ "$op" = "dab" ]; then

echo -e "Download All Generic Background"

(nohup gb2_ds_get --force --new Ds_${date}25_${Attempt}_ccbar_ReverseID) > /home/belle2/amubarak/C02-Grid/Terminal_Output/Ds_${date}25_${Attempt}_ccbar_ReverseID.out & #If things slow down add "--new"
echo -e "1). ccbar submitted\n"
(nohup gb2_ds_get --force --new Ds_${date}25_${Attempt}_BB_ReverseID) > /home/belle2/amubarak/C02-Grid/Terminal_Output/Ds_${date}25_${Attempt}_BB_ReverseID.out & #If things slow down add "--new"
echo -e "2). BB submitted\n"
(nohup gb2_ds_get --force --new Ds_${date}25_${Attempt}_ddbar_ReverseID) > /home/belle2/amubarak/C02-Grid/Terminal_Output/Ds_${date}25_${Attempt}_ddbar_ReverseID.out & #If things slow down add "--new"
echo -e "3). ddbar submitted\n"
(nohup gb2_ds_get --force --new Ds_${date}25_${Attempt}_ssbar_ReverseID) > /home/belle2/amubarak/C02-Grid/Terminal_Output/Ds_${date}25_${Attempt}_ssbar_ReverseID.out & #If things slow down add "--new"
echo -e "4). ssbar submitted\n"
(nohup gb2_ds_get --force --new Ds_${date}25_${Attempt}_taupair_ReverseID) > /home/belle2/amubarak/C02-Grid/Terminal_Output/Ds_${date}25_${Attempt}_taupair_ReverseID.out & #If things slow down add "--new"
echo -e "5). taupair submitted\n"
(nohup gb2_ds_get --force --new Ds_${date}25_${Attempt}_uubar_ReverseID) > /home/belle2/amubarak/C02-Grid/Terminal_Output/Ds_${date}25_${Attempt}_uubar_ReverseID.out & #If things slow down add "--new"
echo -e "6). uubar submitted\n"

echo -e "Download in Progress"

elif [ "$op" = "mcb" ]; then

echo -e "Merge Generic Charm Background Start"
hadd -f /group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh/Ds2D0e-Generic_Ds_${date}25_${Attempt}_ccbar_ReverseID.root Ds_${date}25_${Attempt}_ccbar_ReverseID/sub00/*.root Ds_${date}25_${Attempt}_ccbar_ReverseID/sub01/*.root
# hadd -f /group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh/Ds2D0e-Generic_Ds_${date}25_${Attempt}_ccbar.root Ds_${date}25_${Attempt}_ccbar/sub00/*.root Ds_${date}25_${Attempt}_ccbar/sub01/*.root Ds_${date}25_${Attempt}_ccbar/sub02/*.root Ds_${date}25_${Attempt}_ccbar/sub03/*.root Ds_${date}25_${Attempt}_ccbar/sub04/*.root Ds_${date}25_${Attempt}_ccbar/sub05/*.root
# rm -r Ds_${date}25_${Attempt}_ccbar/ -f
# rm /home/belle2/amubarak/C02-Grid/Terminal_Output/Ds_${date}25_${Attempt}_ccbar.out -f

elif [ "$op" = "mab" ]; then

# echo -e "Merge All Generic Background Start"
hadd -f /group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh/Ds2D0e-Generic_Ds_${date}25_${Attempt}_ccbar_ReverseID_WCh.root Ds_${date}25_${Attempt}_ccbar_ReverseID/sub00/*.root Ds_${date}25_${Attempt}_ccbar_ReverseID/sub01/*.root
# hadd -f /group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh/Ds2D0e-Generic_Ds_${date}25_${Attempt}_BB_ReverseID_WCh.root Ds_${date}25_${Attempt}_BB_ReverseID/sub00/*.root
# hadd -f /group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh/Ds2D0e-Generic_Ds_${date}25_${Attempt}_ddbar_ReverseID_WCh.root Ds_${date}25_${Attempt}_ddbar_ReverseID/sub00/*.root
# hadd -f /group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh/Ds2D0e-Generic_Ds_${date}25_${Attempt}_ssbar_ReverseID_WCh.root Ds_${date}25_${Attempt}_ssbar_ReverseID/sub00/*.root
# hadd -f /group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh/Ds2D0e-Generic_Ds_${date}25_${Attempt}_taupair_ReverseID_WCh.root Ds_${date}25_${Attempt}_taupair_ReverseID/sub00/*.root
# hadd -f /group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID_WCh/Ds2D0e-Generic_Ds_${date}25_${Attempt}_uubar_ReverseID_WCh.root Ds_${date}25_${Attempt}_uubar_ReverseID/sub00/*.root Ds_${date}25_${Attempt}_uubar_ReverseID/sub01/*.root
# rm /home/belle2/amubarak/C02-Grid/Terminal_Output/*.out -f

fi

# rm /home/belle2/amubarak/Ds2D0enue_Analysis/09-Sources/Inclusive/Grid_ReverseID_MC.sh -f
# cp -f /home/belle2/amubarak/Grid_ReverseID_MC.sh /home/belle2/amubarak/Ds2D0enue_Analysis/09-Sources/Inclusive
echo -e "\nDone\n"

# Grid:
# gb2
# Password: AzerDSRT4116!
# source /cvmfs/belle.kek.jp/grid/gbasf2/pro/bashrc