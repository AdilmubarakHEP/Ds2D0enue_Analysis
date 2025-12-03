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
# The number of Kaons is being counted                                   #
#                                                                        #
##########################################################################

import basf2 as b2
import modularAnalysis as ma
import variables.collections as vc
import variables.utils as vu

# create path
my_path = b2.create_path()
#/================================================================================================================================/
# --I/O----------------------------------------------------------------------------------------
# load input ROOT file
ma.inputMdst(environmentType='default',
             filename="",
             path=my_path)
#/================================================================================================================================/
# Create Particle Lists
#---------------------------------
# Kaon
# Creates Kaon Particle List (and c.c.)
ma.fillParticleListFromMC('K-:gen', 'abs(mcPDG) == 321', path=my_path)
#
# K Short
# Creates Kaon Particle List (and c.c.)
ma.fillParticleListFromMC('K_S0:gen', 'abs(mcPDG) == 310', path=my_path)
#/================================================================================================================================/
# Variables
#-------------------
# Select variables that we want to store to ntuple
vars = ['mcP','mcPhi','mcTheta','M','genMotherPDG',
        'nMCDaughters','mcPDG','genMotherPDG(0)','genMotherPDG(1)']
#/================================================================================================================================/
# Saving variables to ntuple
output_file = 'Truth_Info/Finished_Truth/test.root'
ma.variablesToNtuple('K-:gen', vars,
                     filename=output_file, treename='KMCtree', path=my_path)
ma.variablesToNtuple('K_S0:gen', vars,
                     filename=output_file, treename='KSMCtree', path=my_path)
#/================================================================================================================================/

# Process the events
b2.process(my_path)

# print out the summary
print(b2.statistics)