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

mkdir Dstar0_${date}24_${Attempt}_ccbar
mkdir Dstar0_${date}24_${Attempt}_ccbar/sub00 Dstar0_${date}24_${Attempt}_ccbar/sub01

fnames_1=(/group/belle2/dataprod/MC/MC15ri/ccbar/sub00/*.root)
for ((i=0;i<=1001;i++))
do
        bsub basf2 Ds2D0enue_Analysis/02-Steering_Scripts/Background/Dstar0-Reconstruction.py -i ${fnames_1[i]} -o Dstar0_${date}24_${Attempt}_ccbar/sub00/ntuple_$i.root
done
fnames_2=(/group/belle2/dataprod/MC/MC15ri/ccbar/sub01/*.root)
for ((i=0;i<=157;i++))
do
        bsub basf2 Ds2D0enue_Analysis/02-Steering_Scripts/Background/Dstar0-Reconstruction.py -i ${fnames_2[i]} -o Dstar0_${date}24_${Attempt}_ccbar/sub01/ntuple_$i.root
done

elif [ "$op" = "m" ]; then

echo -e "Merge Generic Charm Background Start"
hadd -f Dstar02D0e-Generic_Dstar0_${date}24_${Attempt}_ccbar.root Dstar0_${date}24_${Attempt}_ccbar/sub00/*.root Dstar0_${date}24_${Attempt}_ccbar/sub01/*.root
rm -r Dstar0_${date}24_${Attempt}_ccbar/ -f

fi

echo -e "\nDone\n"