#!/usr/bin/env python3

##########################################################################
# basf2 (Belle II Analysis Software Framework)                           #
# Author: The Belle II Collaboration                                     #
#                                                                        #
# See git log for contributors and copyright holders.                    #
# This file is licensed under LGPL-3.0, see LICENSE.md.                  #
##########################################################################

import sys
import basf2 as b2
import modularAnalysis as ma
import variables.collections as vc
import variables.utils as vu
import vertex as vx
# from variables.MCGenTopo import mc_gen_topo
from variables import variables as vm

# create path
my_path = b2.create_path()

# load input ROOT file
ma.inputMdst(environmentType='default',
             filename="C00-Generation/ccbarDs+EventGeneration.root",
             path=my_path)
#/==========================================================================================================================/
# Gamma
Gamma_Cut = 'daughter(0, p) <= 0.6 and daughter(1, p) <= 0.6 and daughter(0, abs(z0)) <= 3 and daughter(1, abs(z0)) <= 3'
Gamma_Cut += ' and daughter(0,electronID)>=0.5 and daughter(1,electronID)>=0.5'
ma.fillConvertedPhotonsList('gamma:converted -> e+ e-',"", path=my_path)
#
# Vertex Fitting D_s+
#----------------------
vx.treeFit('gamma:converted',conf_level=-1, path=my_path)
ma.applyCuts('gamma:converted', Gamma_Cut, path=my_path)
# #/==========================================================================================================================/
# perform MC matching (MC truth asociation)
ma.matchMCTruth(list_name='gamma:converted', path=my_path) 
#/==========================================================================================================================/
# Gamma Variables

vm.addAlias("e1_p", "daughter(0,p)")
vm.addAlias("e2_p", 'daughter(1,p)')
vm.addAlias("e1_mcPDG", "daughter(0,mcPDG)")
vm.addAlias("e2_mcPDG", 'daughter(1,mcPDG)')

var = ['p','M','isSignal','mcPDG','genMotherPDG','mcErrors','e1_p','e2_p','e1_mcPDG','e2_mcPDG','daughterAngle(0,1)',
       'nMCDaughters','nMCMatches','daughter(0,electronID)','daughter(1,electronID)','chiProb']
#/==========================================================================================================================/
# Saving variables to ntuple
output_file = 'C01-Simulated_Events/Conversion.root'

ma.variablesToNtuple('gamma:converted', var,
                     filename=output_file, treename='gtree', path=my_path)
#/==========================================================================================================================/

# Process the events
b2.process(my_path)

# print out the summary
print(b2.statistics)