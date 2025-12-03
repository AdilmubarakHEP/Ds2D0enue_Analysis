#!/usr/bin/env python3

##########################################################################
# basf2 (Belle II Analysis Software Framework)                           #
# Author: The Belle II Collaboration                                     #
#                                                                        #
# See git log for contributors and copyright holders.                    #
# This file is licensed under LGPL-3.0, see LICENSE.md.                  #
##########################################################################

import basf2_mva
import basf2_mva_util
import multiprocessing
import copy
import subprocess

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

train_data = '/home/belle2/amubarak/C02-MVA/Splits/Ds2D0e_train.root'
test_data = '/home/belle2/amubarak/C02-MVA/Splits/Ds2D0e_test.root'

training_data = basf2_mva.vector(train_data)
testing_data = basf2_mva.vector(test_data)

# Define the variables for training.
variables = ['Ds_extraInfo_FakeD0BDT',
             'Ds_Ds_starminusDs_M_Correction',
             'Ds_chiProb_noIP',
            #  'Ds_chiProb',
             'Ds_gammaveto_M_Correction',
            #  'e_pt'
             ]
# variables = ['D0_decayAngle_0','D0_M']

# variables.remove('Ds_massDifference_0')

general_options = basf2_mva.GeneralOptions()
general_options.m_datafiles = training_data
general_options.m_treename = "Dstree"
general_options.m_identifier = "/home/belle2/amubarak/C02-MVA/Completed/MVAFastBDT.xml"
general_options.m_variables = basf2_mva.vector(*variables)
general_options.m_target_variable = "Ds_isSignal"

xgboost_options = basf2_mva.PythonOptions()
xgboost_options.m_framework = "xgboost"
# param = ('{"max_depth": 3, "eta": 0.1, "silent": 1, "objective": "binary:logistic",'
#             '"subsample": 0.5, "nthread": 1, "nTrees": 400}')
# xgboost_options.m_config = param

# Train a MVA method and store the weightfile (MVAFastBDT.root) locally.
basf2_mva.teacher(general_options, xgboost_options)

# # Evaluate training:
# Check my results.
subprocess.call('basf2_mva_evaluate.py ' + 
                '-id /home/belle2/amubarak/C02-MVA/Completed/MVAFastBDT.xml ' + 
                '-train '+train_data+' ' +
                '-data '+test_data+' ' +
                '-tree Dstree ' +
                '-out /home/belle2/amubarak/C02-MVA/BackgroundSuppressionEvaluation.pdf ' +
                '-c ',
                # '-w MVA/Photos',
                shell=True
                )

# def roc_for_variable_set(variables):
#     method = basf2_mva_util.Method(general_options.m_identifier)
#     options = copy.copy(general_options)
#     options.m_variables = basf2_mva.vector(*variables)
#     m = method.train_teacher(training_data, general_options.m_treename, general_options=options)
#     p, t = m.apply_expert(testing_data, general_options.m_treename)
#     return basf2_mva_util.calculate_auc_efficiency_vs_background_retention(p, t)

method = basf2_mva_util.Method(general_options.m_identifier)
# p, t = method.apply_expert(testing_data, general_options.m_treename)
# global_auc = basf2_mva_util.calculate_auc_efficiency_vs_background_retention(p, t)

# Creating Empty DataFrame and Storing it in variable df
df = pd.DataFrame()

var = []
score = []

# Feature Importance:
#------------------------
# Approach 1: Read out the importance calculated by the method itself
print("Variable importances returned my method")
for variable in method.variables:
    print(variable, method.importances.get(variable, 0.0))
    var.append(variable)
    score.append(method.importances.get(variable, 0.0))

df["Variables"] = var
df["Feature Importance"] = score

# Printing Empty DataFrame
print(df)

df.to_pickle("/home/belle2/amubarak/C02-MVA/Features/Feature_Importance.pkl")

x = [df['Variables'][ind] for ind in df.index]
y = [df['Feature Importance'][ind] for ind in df.index]

# plt.figure(figsize=(10,8))
plt.ylabel(r'Feature Importance Variables')
plt.title(r'Feature Importance')
plt.tick_params(axis='both', which='major')
plt.barh(x,y)

plt.show()
plt.savefig("/home/belle2/amubarak/C02-MVA/Features/Feature_Importance.png")