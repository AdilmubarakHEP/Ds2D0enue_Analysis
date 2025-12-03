echo "BDT Start"

# echo -n "Enter the Date: "
# read date
# echo -n "Attempt: "
# read Attempt

# echo -e "\nDate Choosen: $date"
# echo -e "Attempt: $Attempt"

# rm /home/belle2/amubarak/Ds2D0enue_Analysis/08-Sources/submit_BDT.sh -f
rm C02-MVA/Completed/*.root /home/belle2/amubarak/C02-MVA/Splits/*.root C02-MVA/BackgroundSuppressionEvaluation.pdf -f

basf2_mva_merge_mc -s /home/belle2/amubarak/C01-Simulated_Events/Ds2D0enu-Signal.root \
                   -b /group/belle2/users2022/amubarak/02-Grid/Completed_TopoAna/Ds2D0e-Generic_Ds_021325_2_All.root \
                   -t Dstree \
                   --columns Ds_extraInfo_FakeD0BDT Ds_Ds_starminusDs_M_Correction Ds_chiProb_noIP Ds_chiProb Ds_gammaveto_M_Correction e_pt K_kaonID K_abs_dr pi_pionID pi_abs_dr D0_cos_decayAngle_1 D0_chiProb D0_flightDistance D0_useCMSFrame_p D0_M Ds_isSignal \
                   --fsig 0.5 \
                   --ftest 0.7 \
                   -o /home/belle2/amubarak/C02-MVA/Splits/Ds2D0e.root
                #    --cut "-0.02 <= D0_dM <= 0.02"

# basf2 /home/belle2/amubarak/Ds2D0enue_Analysis/04-ML/basf2/Ds2D0e-XGBoost.py
basf2 /home/belle2/amubarak/Ds2D0enue_Analysis/04-ML/basf2/BackgroundSuppression/Ds2D0e-FastBDT.py
# basf2 /home/belle2/amubarak/Ds2D0enue_Analysis/04-MVA/basf2/Ds2D0e-Sideband_Subtraction.py
# basf2 /home/belle2/amubarak/Ds2D0enue_Analysis/04-MVA/basf2/Ds2D0e-sPlot.py

# General:
# basf2_mva_expert --identifiers C02-MVA/Completed/MVAFastBDT.xml \
#                 --datafiles /home/belle2/amubarak/C03-Grid/Ds2D0e-Generic_Ds_${date}25_${Attempt}_ccbar_C.root \
#                 --treename Dstree \
#                 --outputfile C02-MVA/Completed/Ds2D0enu_MVA.root

basf2_mva_expert --identifiers C02-MVA/Completed/MVAFastBDT.xml \
                --datafiles /home/belle2/amubarak/C02-MVA/Splits/Ds2D0e_test.root \
                --treename Dstree \
                --outputfile C02-MVA/Completed/Ds2D0enu_FoM.root

# # Signal:
# basf2_mva_expert --identifiers C02-MVA/Completed/MVAFastBDT.root \
#                 --datafiles /home/belle2/amubarak/C03-Grid/Completed/Ds2D0e-Generic_Ds_${date}24_${Attempt}_ccbar.root \
#                 --treename Dstree \
#                 --outputfile C02-MVA/Completed/Ds2D0enu-Signal_MVA.root

# # Background:
# basf2_mva_expert --identifiers C02-MVA/Completed/MVAFastBDT.root \
#                 --datafiles /home/belle2/amubarak/C03-Grid/Completed/Ds2D0e-Generic_Ds_${date}24_${Attempt}_ccbar.root \
#                 --treename Dstree \
#                 --outputfile C02-MVA/Completed/Ds2D0enu-Background_MVA.root

echo "BDT Completed"
rm /home/belle2/amubarak/Ds2D0enue_Analysis/11-Sources/submit_BDT.sh -f
cp submit_BDT.sh /home/belle2/amubarak/Ds2D0enue_Analysis/11-Sources
echo "Updated Source File"