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

mkdir Ds_${date}25_${Attempt}_ccbar
mkdir Ds_${date}25_${Attempt}_ccbar/sub00 Ds_${date}25_${Attempt}_ccbar/sub01

fnames_1=(/group/belle2/dataprod/MC/MC15ri/ccbar/sub00/*.root)
for ((i=0;i<=1001;i++))
do
        bsub basf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -i ${fnames_1[i]} -o Ds_${date}25_${Attempt}_ccbar/sub00/ntuple_$i.root
done
fnames_2=(/group/belle2/dataprod/MC/MC15ri/ccbar/sub01/*.root)
for ((i=0;i<=157;i++))
do
        bsub basf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -i ${fnames_2[i]} -o Ds_${date}25_${Attempt}_ccbar/sub01/ntuple_$i.root
done

mkdir Ds_${date}25_${Attempt}_charged
mkdir Ds_${date}25_${Attempt}_charged/sub00

fnames_1=(/group/belle2/dataprod/MC/MC15ri/charged/sub00/*.root)
for ((i=0;i<=601;i++))
do
        bsub basf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -i ${fnames_1[i]} -o Ds_${date}25_${Attempt}_charged/sub00/ntuple_$i.root
done

mkdir Ds_${date}25_${Attempt}_ddbar
mkdir Ds_${date}25_${Attempt}_ddbar/sub00

fnames_1=(/group/belle2/dataprod/MC/MC15ri/ddbar/sub00/*.root)
for ((i=0;i<=320;i++))
do
        bsub basf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -i ${fnames_1[i]} -o Ds_${date}25_${Attempt}_ddbar/sub00/ntuple_$i.root
done

mkdir Ds_${date}25_${Attempt}_mixed
mkdir Ds_${date}25_${Attempt}_mixed/sub00

fnames_1=(/group/belle2/dataprod/MC/MC15ri/mixed/sub00/*.root)
for ((i=0;i<=560;i++))
do
        bsub basf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -i ${fnames_1[i]} -o Ds_${date}25_${Attempt}_mixed/sub00/ntuple_$i.root
done

mkdir Ds_${date}25_${Attempt}_ssbar
mkdir Ds_${date}25_${Attempt}_ssbar/sub00

fnames_1=(/group/belle2/dataprod/MC/MC15ri/ssbar/sub00/*.root)
for ((i=0;i<=305;i++))
do
        bsub basf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -i ${fnames_1[i]} -o Ds_${date}25_${Attempt}_ssbar/sub00/ntuple_$i.root
done

mkdir Ds_${date}25_${Attempt}_uubar
mkdir Ds_${date}25_${Attempt}_uubar/sub00 Ds_${date}25_${Attempt}_uubar/sub01

fnames_1=(/group/belle2/dataprod/MC/MC15ri/uubar/sub00/*.root)
for ((i=0;i<=1001;i++))
do
        bsub basf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -i ${fnames_1[i]} -o Ds_${date}25_${Attempt}_uubar/sub00/ntuple_$i.root
done
fnames_2=(/group/belle2/dataprod/MC/MC15ri/uubar/sub01/*.root)
for ((i=0;i<=275;i++))
do
        bsub basf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -i ${fnames_2[i]} -o Ds_${date}25_${Attempt}_uubar/sub01/ntuple_$i.root
done

elif [ "$op" = "m" ]; then

echo -e "Merge Generic Charm Background Start"
hadd -f /group/belle2/users2022/amubarak/02-Grid/Sample_KEKCC/Ds2D0e-Generic_Ds_${date}25_${Attempt}_ccbar.root Ds_${date}25_${Attempt}_ccbar/sub00/*.root Ds_${date}25_${Attempt}_ccbar/sub01/*.root Ds_${date}25_${Attempt}_ccbar/sub02/*.root Ds_${date}25_${Attempt}_ccbar/sub03/*.root
rm -r Ds_${date}25_${Attempt}_ccbar/ -f

echo -e "Merge Generic Charm Background Start"
hadd -f /group/belle2/users2022/amubarak/02-Grid/Sample_KEKCC/Ds2D0e-Generic_Ds_${date}25_${Attempt}_charged.root Ds_${date}25_${Attempt}_charged/sub00/*.root Ds_${date}25_${Attempt}_charged/sub01/*.root
rm -r Ds_${date}25_${Attempt}_charged/ -f

echo -e "Merge Generic Charm Background Start"
hadd -f /group/belle2/users2022/amubarak/02-Grid/Sample_KEKCC/Ds2D0e-Generic_Ds_${date}25_${Attempt}_ddbar.root Ds_${date}25_${Attempt}_ddbar/sub00/*.root Ds_${date}25_${Attempt}_ddbar/sub01/*.root
rm -r Ds_${date}25_${Attempt}_ddbar/ -f

echo -e "Merge Generic Charm Background Start"
hadd -f /group/belle2/users2022/amubarak/02-Grid/Sample_KEKCC/Ds2D0e-Generic_Ds_${date}25_${Attempt}_mixed.root Ds_${date}25_${Attempt}_mixed/sub00/*.root Ds_${date}25_${Attempt}_mixed/sub01/*.root
rm -r Ds_${date}25_${Attempt}_mixed/ -f

echo -e "Merge Generic Charm Background Start"
hadd -f /group/belle2/users2022/amubarak/02-Grid/Sample_KEKCC/Ds2D0e-Generic_Ds_${date}25_${Attempt}_ssbar.root Ds_${date}25_${Attempt}_ssbar/sub00/*.root Ds_${date}25_${Attempt}_ssbar/sub01/*.root
rm -r Ds_${date}25_${Attempt}_ssbar/ -f

echo -e "Merge Generic Charm Background Start"
hadd -f /group/belle2/users2022/amubarak/02-Grid/Sample_KEKCC/Ds2D0e-Generic_Ds_${date}25_${Attempt}_uubar.root Ds_${date}25_${Attempt}_uubar/sub00/*.root Ds_${date}25_${Attempt}_uubar/sub01/*.root Ds_${date}25_${Attempt}_uubar/sub02/*.root Ds_${date}25_${Attempt}_uubar/sub03/*.root
rm -r Ds_${date}25_${Attempt}_uubar/ -f

# hadd -f /group/belle2/users2022/amubarak/02-Grid/Sample_KEKCC/Ds2D0e-Generic_Ds_${date}25_${Attempt}_All.root /group/belle2/users2022/amubarak/02-Grid/Sample_KEKCC/Ds2D0e-Generic_Ds_${date}25_${Attempt}_ccbar.root /group/belle2/users2022/amubarak/02-Grid/Sample_KEKCC/Ds2D0e-Generic_Ds_${date}25_${Attempt}_charged.root /group/belle2/users2022/amubarak/02-Grid/Sample_KEKCC/Ds2D0e-Generic_Ds_${date}25_${Attempt}_ddbar.root /group/belle2/users2022/amubarak/02-Grid/Sample_KEKCC/Ds2D0e-Generic_Ds_${date}25_${Attempt}_mixed.root /group/belle2/users2022/amubarak/02-Grid/Sample_KEKCC/Ds2D0e-Generic_Ds_${date}25_${Attempt}_ssbar.root /group/belle2/users2022/amubarak/02-Grid/Sample_KEKCC/Ds2D0e-Generic_Ds_${date}25_${Attempt}_uubar.root
# hadd -f /group/belle2/users2022/amubarak/02-Grid/Sample_KEKCC/Ds2D0e-Generic_Ds_${date}25_${Attempt}_uds.root /group/belle2/users2022/amubarak/02-Grid/Sample_KEKCC/Ds2D0e-Generic_Ds_${date}25_${Attempt}_uubar.root /group/belle2/users2022/amubarak/02-Grid/Sample_KEKCC/Ds2D0e-Generic_Ds_${date}25_${Attempt}_ddbar.root /group/belle2/users2022/amubarak/02-Grid/Sample_KEKCC/Ds2D0e-Generic_Ds_${date}25_${Attempt}_ssbar.root

fi

echo -e "\nDone\n"