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
# D_s+ -> D0 e+ nu_e                                                     #
#         |                                                              #
#         +-> K- pi+                                                     #
#                                                                        #
#  This code includes both vetoes                                        #
#                                                                        #
##########################################################################

import basf2 as b2
import modularAnalysis as ma
import variables.utils as vu

# b2.conditions.prepend_globaltag(ma.getAnalysisGlobaltag()) # wait 180ms, doesnt work when confl down, 

#/================================================================================================================================/
# Create Path
#--------------
my_path = b2.Path()
#/================================================================================================================================/
# --I/O----------------------------------------------------------------------------------------
# load input ROOT file
ma.inputMdst(environmentType='default',
             filename="",
             path=my_path)
#/================================================================================================================================/
# Find MC Decay
#-----------------
ma.findMCDecay('D_s+:MC', 'D_s+ -> [D0:kpi -> K- pi+ ] e+:uncorrected ?nu',
               skipNonPrimaryDaughters=False, 
               path=my_path)

ma.printMCParticles(suppressPrint=True,
                    showMomenta=False, 
                    showVertices=True, 
                    showStatus=True, 
                    path=my_path)
#/================================================================================================================================/
# Save Variables
#-------------------
Ds_vars=[]
Ds_vars += vu.create_aliases_for_selected(
    list_of_variables=["seenInARICH","seenInCDC","seenInECL","seenInKLM","seenInPXD","seenInSVD","seenInTOP"],
    decay_string='D_s+ -> [D0:kpi -> K- pi+ ] ^e+:corrected ?nu',
    prefix=['e'])

Ds_vars += vu.create_aliases_for_selected(
    list_of_variables=["M","mcPDG",'nMCDaughters'],
    decay_string='^D_s+ -> [D0:kpi -> K- pi+ ] e+:corrected ?nu',
    prefix=['Ds'])
#/================================================================================================================================/
# Save
#--------
# Saving variables to ntuple
output_file = 'Reconstruction_D_s+.root'
ma.variablesToNtuple('D_s+:MC',Ds_vars,
                     filename=output_file, treename='Dstree', path=my_path)
#/================================================================================================================================/

b2.process(my_path)