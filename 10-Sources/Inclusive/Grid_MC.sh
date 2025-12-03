#!/usr/bin/env bash
# Prompt the user to enter a value for date
echo -n "Enter the Date: "
read date
echo -n "Attempt: "
read Attempt
read -p "Generic Charm Background (pcb), All Generic Background (pab), Download Generic Charm Background (dcb), Download All Generic Background (dab), Merge Generic Charm Background (mcb), Merge All Generic Background (mab) : " op

echo -e "\nDate Choosen: $date"
echo -e "Attempt: $Attempt"
echo -e "Option Choosen: $op\n"

if [ "$op" = "pcb" ]; then

# gbasf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -p Ds_${date}25_${Attempt}_ccbar -s light-2409-toyger --input_dslist /home/belle2/amubarak/Ds2D0enue_Analysis/03-Grid/Inputs/grid_ccbar.txt --force
# echo -e "\nCharm Background Project Submitted\n"
gbasf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -p Ds_${date}25_${Attempt}_ccbar -s light-2409-toyger -i /belle/collection/MC/RSC_MC15rd_4S_ccbar_skim_17230100_v3 --force
echo -e "\nCharm Background Project Submitted\n"

p

elif [ "$op" = "pab" ]; then

gbasf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -p Ds_${date}25_${Attempt}_ccbar -s light-2409-toyger -i /belle/collection/MC/RSC_MC15rd_4S_ccbar_skim_17230100_v3 --force
echo -e "1). ccbar submitted\n"
gbasf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -p Ds_${date}25_${Attempt}_BB -s light-2409-toyger -i /belle/collection/MC/RSC_MC15rd_4S_BB_skim_17230100_v3 --force
echo -e "2). BB submitted\n"
gbasf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -p Ds_${date}25_${Attempt}_ddbar -s light-2409-toyger -i /belle/collection/MC/RSC_MC15rd_4S_ddbar_skim_17230100_v3 --force
echo -e "3). ddbar submitted\n"
# gbasf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -p Ds_${date}25_${Attempt}_mixed -s light-2409-toyger -i /belle/collection/MC/MC15rd_mixed_exp20-26_4S_v2 --force
# echo -e "4). mixed submitted\n"
gbasf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -p Ds_${date}25_${Attempt}_ssbar -s light-2409-toyger -i /belle/collection/MC/RSC_MC15rd_4S_ssbar_skim_17230100_v3 --force
echo -e "4). ssbar submitted\n"
gbasf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -p Ds_${date}25_${Attempt}_taupair -s light-2409-toyger -i /belle/collection/MC/RSC_MC15rd_4S_taupair_skim_17230100_v3 --force
echo -e "5). taupair submitted\n"
gbasf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -p Ds_${date}25_${Attempt}_uubar -s light-2409-toyger -i /belle/collection/MC/RSC_MC15rd_4S_uubar_skim_17230100_v3 --force
echo -e "6). uubar submitted\n"

elif [ "$op" = "dcb" ]; then

echo -e "Download Generic Charm Background"

(nohup gb2_ds_get --force --new Ds_${date}25_${Attempt}_ccbar) > /home/belle2/amubarak/C02-Grid/Terminal_Output/Ds_${date}25_${Attempt}_ccbar.out & #If things slow down add "--new"
echo -e "Download in Progress"

elif [ "$op" = "dab" ]; then

echo -e "Download All Generic Background"

(nohup gb2_ds_get --force --new Ds_${date}25_${Attempt}_ccbar) > /home/belle2/amubarak/C02-Grid/Terminal_Output/Ds_${date}25_${Attempt}_ccbar.out & #If things slow down add "--new"
echo -e "1). ccbar submitted\n"
(nohup gb2_ds_get --force --new Ds_${date}25_${Attempt}_BB) > /home/belle2/amubarak/C02-Grid/Terminal_Output/Ds_${date}25_${Attempt}_BB.out & #If things slow down add "--new"
echo -e "2). BB submitted\n"
(nohup gb2_ds_get --force --new Ds_${date}25_${Attempt}_ddbar) > /home/belle2/amubarak/C02-Grid/Terminal_Output/Ds_${date}25_${Attempt}_ddbar.out & #If things slow down add "--new"
echo -e "3). ddbar submitted\n"
(nohup gb2_ds_get --force --new Ds_${date}25_${Attempt}_ssbar) > /home/belle2/amubarak/C02-Grid/Terminal_Output/Ds_${date}25_${Attempt}_ssbar.out & #If things slow down add "--new"
echo -e "4). ssbar submitted\n"
(nohup gb2_ds_get --force --new Ds_${date}25_${Attempt}_taupair) > /home/belle2/amubarak/C02-Grid/Terminal_Output/Ds_${date}25_${Attempt}_taupair.out & #If things slow down add "--new"
echo -e "5). taupair submitted\n"
(nohup gb2_ds_get --force --new Ds_${date}25_${Attempt}_uubar) > /home/belle2/amubarak/C02-Grid/Terminal_Output/Ds_${date}25_${Attempt}_uubar.out & #If things slow down add "--new"
echo -e "6). uubar submitted\n"

echo -e "Download in Progress"

elif [ "$op" = "mcb" ]; then

echo -e "Merge Generic Charm Background Start"
hadd -f /group/belle/users/amubarak/02-Grid/Sample_Grid/Ds2D0e-Generic_Ds_${date}25_${Attempt}_ccbar.root Ds_${date}25_${Attempt}_ccbar/sub00/*.root Ds_${date}25_${Attempt}_ccbar/sub01/*.root
# hadd -f /group/belle/users/amubarak/02-Grid/Sample_Grid/Ds2D0e-Generic_Ds_${date}25_${Attempt}_ccbar.root Ds_${date}25_${Attempt}_ccbar/sub00/*.root Ds_${date}25_${Attempt}_ccbar/sub01/*.root Ds_${date}25_${Attempt}_ccbar/sub02/*.root Ds_${date}25_${Attempt}_ccbar/sub03/*.root Ds_${date}25_${Attempt}_ccbar/sub04/*.root Ds_${date}25_${Attempt}_ccbar/sub05/*.root
# rm -r Ds_${date}25_${Attempt}_ccbar/ -f
# rm /home/belle2/amubarak/C02-Grid/Terminal_Output/Ds_${date}25_${Attempt}_ccbar.out -f

elif [ "$op" = "mab" ]; then

echo -e "Merge All Generic Background Start"
# hadd -f /group/belle/users/amubarak/02-Grid/Sample_Grid/Ds2D0e-Generic_Ds_${date}25_${Attempt}_ccbar.root Ds_${date}25_${Attempt}_ccbar/sub00/*.root Ds_${date}25_${Attempt}_ccbar/sub01/*.root
# hadd -f /group/belle/users/amubarak/02-Grid/Sample_Grid/Ds2D0e-Generic_Ds_${date}25_${Attempt}_BB.root Ds_${date}25_${Attempt}_BB/sub00/*.root
# hadd -f /group/belle/users/amubarak/02-Grid/Sample_Grid/Ds2D0e-Generic_Ds_${date}25_${Attempt}_ddbar.root Ds_${date}25_${Attempt}_ddbar/sub00/*.root
# hadd -f /group/belle/users/amubarak/02-Grid/Sample_Grid/Ds2D0e-Generic_Ds_${date}25_${Attempt}_ssbar.root Ds_${date}25_${Attempt}_ssbar/sub00/*.root
# hadd -f /group/belle/users/amubarak/02-Grid/Sample_Grid/Ds2D0e-Generic_Ds_${date}25_${Attempt}_taupair.root Ds_${date}25_${Attempt}_taupair/sub00/*.root
# hadd -f /group/belle/users/amubarak/02-Grid/Sample_Grid/Ds2D0e-Generic_Ds_${date}25_${Attempt}_uubar.root Ds_${date}25_${Attempt}_uubar/sub00/*.root Ds_${date}25_${Attempt}_uubar/sub01/*.root
rm -r Ds_${date}25_${Attempt}_ccbar/ Ds_${date}25_${Attempt}_BB/ Ds_${date}25_${Attempt}_ddbar/ Ds_${date}25_${Attempt}_ssbar/ Ds_${date}25_${Attempt}_taupair/ Ds_${date}25_${Attempt}_uubar/ -f
# rm /home/belle2/amubarak/C02-Grid/Terminal_Output/*.out -f

fi

# rm /home/belle2/amubarak/Ds2D0enue_Analysis/09-Sources/Inclusive/Grid_MC.sh -f
# cp -f Grid_MC.sh /home/belle2/amubarak/Ds2D0enue_Analysis/09-Sources/Inclusive
echo -e "\nDone\n"

# Grid:
# gb2
# Password: AzerDSRT4116!
# source /cvmfs/belle.kek.jp/grid/gbasf2/pro/bashrc