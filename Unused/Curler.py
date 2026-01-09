#!/usr/bin/env python3

##########################################################################
# basf2 (Belle II Analysis Software Framework)                           #
# Author: The Belle II Collaboration                                     #
#                                                                        #
# See git log for contributors and copyright holders.                    #
# This file is licensed under LGPL-3.0, see LICENSE.md.                  #
##########################################################################

##########################################################################
#                                                                        #
# Stuck? Ask for help at questions.belle2.org                            #
#                                                                        #
# This tutorial demonstrates how to reconstruct the                      #
# following  decay chain:                                                #
#                                                                        #
# D_s+* -> D_s+ gamma                                                    #
#           |                                                            #
#           +->  D_s+ -> D0 e+ nu_e                                      #
#                        |                                               #
#                        +-> K- pi+                                      #
#                                                                        #
##########################################################################

import basf2 as b2
import modularAnalysis as ma
import variables.collections as vc
import variables.utils as vu
import vertex as vx
import variables as va
# from variables.MCGenTopo import mc_gen_topo
from variables import variables as vm
from stdPi0s import stdPi0s
from stdPhotons import stdPhotons

b2.conditions.prepend_globaltag(ma.getAnalysisGlobaltag()) # wait 180ms, doesnt work when confl down, 
b2.conditions.prepend_globaltag('user_adilmub_FakeD0Suppression')

#/================================================================================================================================/
# Curlers
#--------------
input_file_name = '/group/belle2/dataprod/MC/MC15ri/ccbar/sub00/*.root'

training_path = b2.Path()

ma.inputMdst(environmentType='default',
             filename=input_file_name,
             path=training_path)

ma.fillParticleList("e+:all", cut='', path=training_path)

ma.tagCurlTracks("e+:all", selectorType='mva', expert_train=True, expert_filename='test.root', ptCut=0.6, path=training_path)

b2.process(training_path, int(2e5))
#/================================================================================================================================/
