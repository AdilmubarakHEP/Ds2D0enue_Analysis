#!/usr/bin/env bash
# Prompt the user to enter a value for date
echo -n "Enter the Date: "
read date
echo -n "Attempt: "
read Attempt
read -p "Generic Background (s or m) : " op

echo -e "\nDate Choosen: $date"
echo -e "Attempt: $Attempt"
echo -e "Option Choosen: $op\n"

if [ "$op" = "s" ]; then

mkdir Ds_${date}24_${Attempt}_uubar
mkdir Ds_${date}24_${Attempt}_uubar/sub00 Ds_${date}24_${Attempt}_uubar/sub01

fnames_1=(/group/belle2/dataprod/MC/MC15ri/uubar/sub00/*.root)
for ((i=0;i<=1001;i++))
do
        bsub basf2 Ds2D0enue_Analysis/02-Steering_Scripts/Ds2D0e-Reconstruction.py -i ${fnames_1[i]} -o Ds_${date}24_${Attempt}_uubar/sub00/ntuple_$i.root
done
fnames_2=(/group/belle2/dataprod/MC/MC15ri/uubar/sub01/*.root)
for ((i=0;i<=275;i++))
do
        bsub basf2 Ds2D0enue_Analysis/02-Steering_Scripts/Ds2D0e-Reconstruction.py -i ${fnames_2[i]} -o Ds_${date}24_${Attempt}_uubar/sub01/ntuple_$i.root
done

elif [ "$op" = "m" ]; then

echo -e "Merge Generic Charm Background Start"
hadd -f C03-Grid/Completed/Ds2D0e-Generic_Ds_${date}24_${Attempt}_uubar.root Ds_${date}24_${Attempt}_uubar/sub00/*.root Ds_${date}24_${Attempt}_uubar/sub01/*.root
rm -r Ds_${date}24_${Attempt}_uubar/ -f

fi

echo -e "\nDone\n"