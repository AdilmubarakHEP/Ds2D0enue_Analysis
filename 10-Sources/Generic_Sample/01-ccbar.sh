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

mkdir /group/belle2/users2022/amubarak/Ds_${date}25_${Attempt}_ccbar
mkdir /group/belle2/users2022/amubarak/Ds_${date}25_${Attempt}_ccbar/sub00 /group/belle2/users2022/amubarak/Ds_${date}25_${Attempt}_ccbar/sub01

fnames_1=(/group/belle2/dataprod/MC/MC15ri/ccbar/sub00/*.root)
for ((i=0;i<=1001;i++))
do
        bsub basf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -i ${fnames_1[i]} -o /group/belle2/users2022/amubarak/Ds_${date}25_${Attempt}_ccbar/sub00/ntuple_${Attempt}_$i.root
done
fnames_2=(/group/belle2/dataprod/MC/MC15ri/ccbar/sub01/*.root)
for ((i=0;i<=157;i++))
do
        bsub basf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -i ${fnames_2[i]} -o /group/belle2/users2022/amubarak/Ds_${date}25_${Attempt}_ccbar/sub01/ntuple_${Attempt}_$i.root
done

elif [ "$op" = "m" ]; then

echo -e "Merge Generic Charm Background Start"
hadd -f /group/belle2/users2022/amubarak/02-Grid/Sample_KEKCC/Ds2D0e-Generic_Ds_${date}25_${Attempt}_ccbar.root /group/belle2/users2022/amubarak/Ds_${date}25_${Attempt}_ccbar/sub00/*.root /group/belle2/users2022/amubarak/Ds_${date}25_${Attempt}_ccbar/sub01/*.root
rm -r /group/belle2/users2022/amubarak/Ds_${date}25_${Attempt}_ccbar/ -f

fi

echo -e "\nDone\n"