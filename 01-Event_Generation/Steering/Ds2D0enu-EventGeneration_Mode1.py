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
# e+e- -> ccbar -> D_s+ anything event generation                        #
#                                                                        #
# This tutorial demonstrates how to generate                             #
#                                                                        #
# e+e- -> ccbar -> D_s+ anything                                         #
#                                                                        #
# events with EvtGen in BASF2, where the decay of D*+                    #
# is specified by the given .dec file.                                   #
#                                                                        #
# The generated events are saved to the output ROOT file.                #
# In each event the generated particles (MCParticle objects)             #
# are stored in the StoreArray<MCParticle>.                              #
#                                                                        #
##########################################################################

"""
This script saves e+ or e- from photon conversions into a pair in MCParticles.

<header>
  <contact>dorisykim@ssu.ac.kr</contact>
  <description>
      Saves 100 generic BBbar events with EvtGen + secondary e+ or e- from pair conversions created by Geant4 in MCParticles.
      The corresponding secondaryPhysicsProcess ID is 14, which is defined as fGammaConversion in G4EmProcessSubType.h.
      The detector simulation mixed with background, trigger simulation, and standard reconstruction is done.
  </description>
</header>
"""

import basf2 as b2
import modularAnalysis as ma
import mdst as mdst
import simulation as si
import reconstruction as re
import generators as ge
import glob
import os

# Defining one path
my_path = b2.create_path()

ma.setupEventInfo(10000, path=my_path)

# # background files
# # location of the files is obtained from a shell variable - check first if it is set
# if 'BELLE2_BACKGROUND_DIR' not in os.environ:
#     b2.B2FATAL(
#         'BELLE2_BACKGROUND_DIR variable is not set. \n'
#         'Please export (setenv) the variable to the location of BG overlay sample. \n'
#         'Check https://xwiki.desy.de/xwiki/rest/p/90869 to find them')
# # get list of files and check the list length
# bg = glob.glob(os.environ['BELLE2_BACKGROUND_DIR'] + '/*.root')
# if len(bg) == 0:
#     b2.B2FATAL('No files found in ', os.environ['BELLE2_BACKGROUND_DIR'])

# generate ccbar events
ge.add_inclusive_continuum_generator(finalstate="ccbar",
                                     particles=["D_s+"],
                                     userdecfile=b2.find_file('/home/belle2/amubarak/Ds2D0enue_Analysis/01-Event_Generation/dec/Ds2D0enu-EventGeneration_Mode1.dec'),
                                     include_conjugates=True,
                                     path=my_path) 

# simulation
si.add_simulation(path=my_path)
# si.add_simulation(path=my_path, bkgfiles=bg)

# # saving e+ or e- from pair conversions with kinetic energy > 10.0 MeV.
# b2.set_module_parameters(my_path, "FullSim", StorePairConversions=True, PairConversionsEnergyCut=10.0)

# reconstruction
re.add_reconstruction(path=my_path)

# dump in MDST format
mdst.add_mdst_output(mc=True, filename='/home/belle2/amubarak/C00-Generation/ccbarDs+EventGeneration_Mode1.root', path=my_path)

# Show progress of processing
my_path.add_module('ProgressBar')

# Process the events
b2.process(my_path)

# print out the summary
print(b2.statistics)