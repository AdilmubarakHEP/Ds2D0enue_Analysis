#!/usr/bin/env bash
# Prompt the user to enter a value for date
read -p "Submit or Merge (s or m) : " op

echo -e "Option Choosen: $op\n"

echo "Reconstruction Start"

if [ "$op" = "s" ]; then

rm C01-Simulated_Events/Ds2D0enu-Signal.root -f
mkdir C01-Simulated_Events/Ds

fnames_1=(/group/belle/users/amubarak/00-Generation/Signal/Dsp/*.root)
for ((i=1;i<=5;i++))
do
        bsub basf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -i ${fnames_1[i]} -o C01-Simulated_Events/Ds/Ds2D0enu-Signal_p_$i.root
done

fnames_2=(/group/belle/users/amubarak/00-Generation/Signal/Dsm/*.root)
for ((i=1;i<=5;i++))
do
        bsub basf2 /home/belle2/amubarak/Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -i ${fnames_2[i]} -o C01-Simulated_Events/Ds/Ds2D0enu-Signal_m_$i.root
done

# bsub -q h "basf2 Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -i C00-Generation/Signal/Dsp/mdst_000001_prod00046485_task10020000001.root -o C01-Simulated_Events/Ds/Ds2D0enu-Signal_p_1.root"
# bsub -q h "basf2 Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -i C00-Generation/Signal/Dsp/mdst_000002_prod00046485_task10020000002.root -o C01-Simulated_Events/Ds/Ds2D0enu-Signal_p_2.root"
# bsub -q h "basf2 Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -i C00-Generation/Signal/Dsp/mdst_000003_prod00046485_task10020000003.root -o C01-Simulated_Events/Ds/Ds2D0enu-Signal_p_3.root"
# bsub -q h "basf2 Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -i C00-Generation/Signal/Dsp/mdst_000004_prod00046485_task10020000004.root -o C01-Simulated_Events/Ds/Ds2D0enu-Signal_p_4.root"
# bsub -q h "basf2 Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -i C00-Generation/Signal/Dsp/mdst_000005_prod00046485_task10020000005.root -o C01-Simulated_Events/Ds/Ds2D0enu-Signal_p_5.root"
# # Ds2D0enue_Analysis/02-Reconstruction_Scripts/Signal_Grid/Dsp-Reconstruction.py

# bsub -q h "basf2 Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -i C00-Generation/Signal/Dsm/mdst_000001_prod00046486_task10020000001.root -o C01-Simulated_Events/Ds/Ds2D0enu-Signal_m_1.root"
# bsub -q h "basf2 Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -i C00-Generation/Signal/Dsm/mdst_000002_prod00046486_task10020000002.root -o C01-Simulated_Events/Ds/Ds2D0enu-Signal_m_2.root"
# bsub -q h "basf2 Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -i C00-Generation/Signal/Dsm/mdst_000003_prod00046486_task10020000003.root -o C01-Simulated_Events/Ds/Ds2D0enu-Signal_m_3.root"
# bsub -q h "basf2 Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -i C00-Generation/Signal/Dsm/mdst_000004_prod00046486_task10020000004.root -o C01-Simulated_Events/Ds/Ds2D0enu-Signal_m_4.root"
# bsub -q h "basf2 Ds2D0enue_Analysis/02-Reconstruction_Scripts/Ds2D0e-Reconstruction.py -i C00-Generation/Signal/Dsm/mdst_000005_prod00046486_task10020000005.root -o C01-Simulated_Events/Ds/Ds2D0enu-Signal_m_5.root"
# # Ds2D0enue_Analysis/02-Reconstruction_Scripts/Signal_Grid/Dsm-Reconstruction.py

echo -e "\nAll Signal Submitted \n"
elif [ "$op" = "m" ]; then

hadd -f C01-Simulated_Events/Ds2D0enu-Signal.root C01-Simulated_Events/Ds/*.root
rm -r C01-Simulated_Events/Ds -f

echo -e "\nAll Signal Merged \n"
fi

echo "Done"
# rm /home/belle2/amubarak/Ds2D0enue_Analysis/09-Sources/Signal/Reconstruction.sh -f
# cp Reconstruction.sh /home/belle2/amubarak/Ds2D0enue_Analysis/09-Sources/Signal
# echo "Updated Source File"