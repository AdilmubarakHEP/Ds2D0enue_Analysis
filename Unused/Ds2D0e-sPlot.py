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
variables = ['K_kaonID',
             'pi_pionID',
             'D0_decayAngle_0','D0_chiProb','D0_flightDistance','D0_useCMSFrame_p','D0_M']

# variables.remove('Ds_massDifference_0')

# Perform an sPlot training
general_options = basf2_mva.GeneralOptions()
general_options.m_datafiles = training_data
general_options.m_identifier = "MVAFull"
general_options.m_treename = "tree"
general_options.m_variables = basf2_mva.vector(*variables)
general_options.m_target_variable = "isSignal"

fastbdt_options = basf2_mva.FastBDTOptions()
# SPlot is more stable if one doesn't use the randRatio
# FastBDT has a special sPlot mode, but which isn't implemented yet in the mva package
fastbdt_options.m_nTrees = 100
fastbdt_options.m_randRatio = 1.0
basf2_mva.teacher(general_options, fastbdt_options)

general_options.m_identifier = "MVAOrdinary"
general_options.m_variables = basf2_mva.vector(*variables[1:])
basf2_mva.teacher(general_options, fastbdt_options)

meta_options = basf2_mva.MetaOptions()
meta_options.m_use_splot = True
meta_options.m_splot_variable = "M"
# SPlot training assumes that the datafile given to the general options contains only data
# It requires an additional file with MC information from which it can extract the distribution
# of the discriminating variable (in this case M).
# Here we use the same file
general_options.m_datafiles = basf2_mva.vector("train_data.root")
meta_options.m_splot_mc_files = basf2_mva.vector("train_mc.root")

# First we do an ordinary sPlot training
general_options.m_identifier = "MVASPlot"
meta_options.m_splot_combined = False
meta_options.m_splot_boosted = False
basf2_mva.teacher(general_options, fastbdt_options, meta_options)

# Now we combine the sPlot training with a PDF classifier for M, in one step
general_options.m_identifier = "MVASPlotCombined"
meta_options.m_splot_combined = True
meta_options.m_splot_boosted = False
basf2_mva.teacher(general_options, fastbdt_options, meta_options)

# Now we use a boosted sPlot training
general_options.m_identifier = "MVASPlotBoosted"
meta_options.m_splot_combined = False
meta_options.m_splot_boosted = True
basf2_mva.teacher(general_options, fastbdt_options, meta_options)

# And finally a boosted and combined training
general_options.m_identifier = "MVASPlotCombinedBoosted"
meta_options.m_splot_combined = True
meta_options.m_splot_boosted = True
basf2_mva.teacher(general_options, fastbdt_options, meta_options)

# Also do a training of only the pdf classifier
pdf_options = basf2_mva.PDFOptions()
general_options.m_method = 'PDF'
general_options.m_identifier = "MVAPdf"
general_options.m_variables = basf2_mva.vector('M')
basf2_mva.teacher(general_options, pdf_options)

# Apply the trained methods on data
basf2_mva.expert(basf2_mva.vector('MVAPdf', 'MVAFull', 'MVAOrdinary', 'MVASPlot',
                                    'MVASPlotCombined', 'MVASPlotBoosted', 'MVASPlotCombinedBoosted'),
                 basf2_mva.vector('train.root'), 'tree', 'expert.root')
