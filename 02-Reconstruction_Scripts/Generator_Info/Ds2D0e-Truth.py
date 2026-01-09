#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
##########################################################################

import basf2 as b2
import modularAnalysis as ma
import variables.collections as vc
import variables.utils as vu

# create path
my_path = b2.create_path()

# load input ROOT file
ma.inputMdst(environmentType='default',
             filename=b2.find_file(''),
             path=my_path)

# Creating a MC list
# Daughters
ma.fillParticleListFromMC('K-:gen',  'abs(genMotherPDG) == 421 and abs(genMotherPDG(1)) == 431 and abs(mcPDG) == 321', path=my_path)
ma.fillParticleListFromMC('pi+:gen', 'abs(genMotherPDG) == 421 and abs(genMotherPDG(1)) == 431 and abs(mcPDG) == 211', path=my_path)
ma.fillParticleListFromMC('e+:gen',  'abs(genMotherPDG) == 431 and abs(mcPDG) == 11', path=my_path)
ma.fillParticleListFromMC('nu_e:gen','abs(genMotherPDG) == 431 and abs(mcPDG) == 12', path=my_path)
# Parents
ma.fillParticleListFromMC('D0:gen',  'abs(genMotherPDG) == 431 and abs(mcPDG) == 421', path=my_path)
ma.fillParticleListFromMC('D_s+:gen','abs(mcPDG) == 431', path=my_path)
#/==========================================================================================================================/
# Select variables that we want to store to ntuple
vars = ["mcE","mcP","mcPDG",'genMotherPDG','genMotherPDG(1)']
#/==========================================================================================================================/

# Saving variables to ntuple
output_file = '/home/belle2/amubarak/C04-Generator_Level/D_s+_Truth.root'
ma.variablesToNtuple('K-:gen', vars,
                     filename=output_file, treename='KMCtree', path=my_path)
ma.variablesToNtuple('pi+:gen', vars,
                     filename=output_file, treename='PiMCtree', path=my_path)
ma.variablesToNtuple('e+:gen', vars,
                     filename=output_file, treename='eMCtree', path=my_path)
ma.variablesToNtuple('nu_e:gen', vars,
                     filename=output_file, treename='nuMCtree', path=my_path)
# #
ma.variablesToNtuple('D_s+:gen', vars + ["nMCDaughters"],
                     filename=output_file, treename='D_s+MCtree', path=my_path)
ma.variablesToNtuple('D0:gen', vars + ["nMCDaughters","D0Mode"],
                     filename=output_file, treename='D0MCtree', path=my_path)
#/==========================================================================================================================/

# Process the events
b2.process(my_path)

# print out the summary
print(b2.statistics)

