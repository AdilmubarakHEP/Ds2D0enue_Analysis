#!/usr/bin/env python3

##########################################################################
# basf2 (Belle II Analysis Software Framework)                           #
# Author: The Belle II Collaboration                                     #
#                                                                        #
# See git log for contributors and copyright holders.                    #
# This file is licensed under LGPL-3.0, see LICENSE.md.                  #
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

b2.conditions.prepend_globaltag(ma.getAnalysisGlobaltag()) # wait 180ms, doesnt work when confl down, 

# if sys.argv[1]==0:
#     Charge = True
# else:
#     Charge = False
#/================================================================================================================================/
# Charge Conjugation:
#-----------------------
ChargeC = True
#/================================================================================================================================/
# Cuts
#-----------
Electron_Cut = 'electronID > 0.5 and abs(d0) < 1 and abs(z0) < 3'
#/================================================================================================================================/
# Create Path
#--------------
my_path = b2.create_path()
#/================================================================================================================================/
# --I/O----------------------------------------------------------------------------------------
# load input ROOT file
ma.inputMdst(environmentType='default',
             filename="/home/belle2/amubarak/C00-Generation/ccbarDs+EventGeneration.root",
             path=my_path)
#/================================================================================================================================/
# Create Particle Lists
#---------------------------------
# Electron
# Creates Electron Particle List (and c.c.)
# Cut on Variables: electronID, abs(d0), and abs(z0)
# ma.fillParticleList("e+:uncorrected",'', path=my_path)
ma.fillParticleList("e+:uncorrected", Electron_Cut, path=my_path)
#/================================================================================================================================/
# Bremsstrahlung Correction
#-----------------------------
# apply Bremsstrahlung correction to electrons [S10|S20]
vm.addAlias(
    "goodFWDGamma", "passesCut(clusterReg == 1 and clusterE > 0.075)"
)  # [E10]
vm.addAlias(
    "goodBRLGamma", "passesCut(clusterReg == 2 and clusterE > 0.05)"
)
vm.addAlias(
    "goodBWDGamma", "passesCut(clusterReg == 3 and clusterE > 0.1)"
)
vm.addAlias(
    "goodGamma", "passesCut(goodFWDGamma or goodBRLGamma or goodBWDGamma)"
)

ma.fillParticleList("gamma:brems", "goodGamma", path=my_path)

ma.correctBrems("e+:corrected", "e+:uncorrected", "gamma:brems", path=my_path)
vm.addAlias("isBremsCorrected", "extraInfo(bremsCorrected)")  # [E30]

ma.applyCuts('e+:corrected',Electron_Cut, path=my_path)
#/================================================================================================================================/
# Photon Conversion
#--------------------------

# Reconstruct gamma:V0_array -> e+ e-
#-----------------------------------------------
ma.fillConvertedPhotonsList(decayString='gamma:V0_array -> e+ e-',
                            cut='', 
                            path=my_path)
ma.copyParticles('gamma:V0_full', 'gamma:V0', writeOut=False, path=my_path)
#
# Vertex Fitting gamma:V0_array 
#----------------------------------
vx.treeFit('gamma:V0_full', 
           conf_level=-1,
           updateAllDaughters=True, 
           path=my_path)
ma.applyCuts('gamma:V0_full', 'M<=0.1', 
             path=my_path)

# Reconstruct gamma:v0_roe -> e+:sig e-:roe
#-----------------------------------------------
# make combinations of signal electron candidates with all electron from ROE
ma.reconstructDecay(decayString='gamma:RD -> e+:corrected e-:corrected', 
                    cut='',
                    chargeConjugation=ChargeC,
                    path=my_path)
vx.treeFit('gamma:RD', 
           conf_level=-1,
           updateAllDaughters=True, 
           path=my_path)
ma.applyCuts('gamma:RD', 'M<=0.1', path=my_path)

# Merge List
ma.mergeListsWithBestDuplicate('gamma:merged', ['gamma:V0_full', 'gamma:RD'],
                                variable='particleSource', preferLowest=True, path=my_path)
# #/==========================================================================================================================/
# perform MC matching (MC truth asociation)
ma.matchMCTruth(list_name='gamma:V0_full', path=my_path)
ma.matchMCTruth(list_name='gamma:RD', path=my_path)

# Perform MC Matching (MC truth asociation)
ma.matchMCTruth(list_name='gamma:merged', 
                path=my_path)
#/==========================================================================================================================/
# Photon Conversion Variables
#------------------------------
tracks = ['abs(dr)','abs(dz)','z0','d0','pValue',
          'firstCDCLayer','firstPXDLayer','firstSVDLayer',
          'nPXDHits','nVXDHits','nSVDHits','nCDCHits',
          'seenInCDC','seenInPXD','seenInSVD']

V0_full = []
V0_full += vu.create_aliases_for_selected(
    list_of_variables = ['InvM','M','convertedPhotonInvariantMass(0, 1)',
                         'convertedPhotonDelR(0,1)',
                         'isFromV0','mdstIndex','particleSource',
                         'distance','mcDecayVertexFromIPDistance','flightDistance',
                         'daughterAngle(0,1)',
                         'isSignal','mcPDG','genMotherPDG','mcErrors',
                         'genNMissingDaughter(11)','genNMissingDaughter(22)',
                         'chiProb',
                         'nMCDaughters','nMCMatches'
                         ],
    decay_string='^gamma:V0_full -> e+ e-',
    prefix=['gammaV0'])
V0_full += vu.create_aliases_for_selected(
    list_of_variables = ['electronID',
                         'isFromTrack','mdstIndex',
                         'p','pt',
                         'mcPDG','genMotherPDG',
                         'sourceObjectIsInList(e+:corrected)'
                         ] + tracks,
    decay_string='gamma:V0_full -> ^e+ e-',
    prefix=['eplusV0'])
V0_full += vu.create_aliases_for_selected(
    list_of_variables = ['electronID',
                         'isFromTrack','mdstIndex',
                         'p','pt',
                         'mcPDG','genMotherPDG',
                         'sourceObjectIsInList(e+:corrected)'
                         ] + tracks,
    decay_string='gamma:V0_full -> e+ ^e-',
    prefix=['eminusV0'])

RD = []
vm.addAlias("Lost", "passesCut(genNMissingDaughter(22)==1)")
RD += vu.create_aliases_for_selected(
    list_of_variables = ['InvM','M','convertedPhotonInvariantMass(0, 1)',
                         'convertedPhotonDelR(0,1)',
                         'isFromV0','mdstIndex',
                         'distance','mcDecayVertexFromIPDistance','flightDistance',
                         'daughterAngle(0,1)',
                         'isSignal','mcPDG','genMotherPDG','mcErrors',
                         'genNMissingDaughter(11)','genNMissingDaughter(22)','Lost',
                         'chiProb',
                         'nMCDaughters','nMCMatches'
                         ],
    decay_string='^gamma:RD -> e+:corrected e-:corrected',
    prefix=['gammaRD'])
RD += vu.create_aliases_for_selected(
    list_of_variables = ['electronID',
                         'isFromTrack','mdstIndex',
                         'p','pt',
                         'mcPDG','genMotherPDG',
                         ] + tracks,
    decay_string='gamma:RD -> ^e+:corrected e-:corrected',
    prefix=['eplusRD'])
RD += vu.create_aliases_for_selected(
    list_of_variables = ['electronID',
                         'isFromTrack','mdstIndex',
                         'p','pt',
                         'mcPDG','genMotherPDG',
                         ] + tracks,
    decay_string='gamma:RD -> e+:corrected ^e-:corrected',
    prefix=['eminusRD'])
#/==========================================================================================================================/
# Overall Cut:
#------------------
# ma.applyCuts('gamma:V0_full', "gammaV0_distance <= 20", path=my_path)
# ma.applyCuts('gamma:RD', "gammaRD_distance <= 10", path=my_path)
#/==========================================================================================================================/
# Saving variables to ntuple
output_file = 'C01-Simulated_Events/Conversion.root'

ma.variablesToNtuple('gamma:RD', RD,
                     filename=output_file, treename='CPtreeRD', path=my_path)
ma.variablesToNtuple('gamma:V0_full', V0_full,
                     filename=output_file, treename='CPtreeV0', path=my_path)

# ma.variablesToNtuple('e+:test', ['M','particleSource','mcMother(mcPDG)','isDaughterOfList(gamma:test)'],
#                      filename=output_file, treename='test2', path=my_path)
# ma.variablesToNtuple('gamma:test', ['M','particleSource','daughter(0,isDescendantOfList(gamma:RD,1))'],
#                      filename=output_file, treename='test1', path=my_path)
#/==========================================================================================================================/

# Process the events
b2.process(my_path)

# print out the summary
print(b2.statistics)