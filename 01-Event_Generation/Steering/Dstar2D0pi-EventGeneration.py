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
# e+e- -> ccbar -> D*+ antyhing event generation                         #
#                                                                        #
# This tutorial demonstrates how to generate                             #
#                                                                        #
# e+e- -> ccbar -> D*+ anything                                          #
#                                                                        #
# events with EvtGen in BASF2, where the decay of D*+                    #
# is specified by the given .dec file.                                   #
#                                                                        #
# The generated events are saved to the output ROOT file.                #
# In each event the generated particles (MCParticle objects)             #
# are stored in the StoreArray<MCParticle>.                              #
#                                                                        #
##########################################################################

import basf2 as b2
import modularAnalysis as ma
import mdst as mdst
import simulation as si
import reconstruction as re
import generators as ge

# Defining one path
my_path = b2.create_path()

ma.setupEventInfo(10000, path=my_path)

# generate ccbar events
ge.add_inclusive_continuum_generator(finalstate="ccbar",
                                     particles=["D*+"],
                                     userdecfile=b2.find_file('/home/belle2/amubarak/Ds2D0enue_Analysis/01-Event_Generation/dec/Dstar2D0pi-EventGeneration.dec'),
                                     include_conjugates=True,
                                     path=my_path) 

# simulation
si.add_simulation(path=my_path)

# reconstruction
re.add_reconstruction(path=my_path)

# dump in MDST format
mdst.add_mdst_output(mc=True, filename='Event/Dstar+/ccbarDstarEventGeneration.root', path=my_path)

# Show progress of processing
my_path.add_module('ProgressBar')

# Process the events
b2.process(my_path)

# print out the summary
print(b2.statistics)