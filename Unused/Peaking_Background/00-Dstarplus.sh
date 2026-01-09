#!/usr/bin/env bash
# Prompt the user to enter a value for date
echo -n "Total Files: "
read Files
read -p "Submit or Merge (s or m) : " op

echo -e "Total Files: $Files"
echo -e "Option Choosen: $op\n"

if [ "$op" = "s" ]; then

mkdir /group/belle2/users2022/amubarak/Dstarplus

fnames_1=(/group/belle2/users2022/amubarak/00-Generation/Dstar+/*.root)
for ((i=1;i<=${Files};i++))
do
        bsub basf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Background/Dstarplus-Reconstruction.py -i ${fnames_1[i]} -o /group/belle2/users2022/amubarak/Dstarplus/ntuple_$i.root
done

elif [ "$op" = "m" ]; then

echo -e "Merge Generic Charm Background Start"
hadd -f /home/belle2/amubarak/C01-Simulated_Events/Dstarplus-Background.root /group/belle2/users2022/amubarak/Dstarplus/*.root
rm -r /group/belle2/users2022/amubarak/Dstarplus/ -f

fi

echo -e "\nDone\n"