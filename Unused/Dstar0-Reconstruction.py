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
# D*0 -> D0 (gamma or pi0)                                               #
#         |                                                              #
#         +-> K- pi+                                                     #
#                                                                        #
#  This code includes both vetoes                                        #
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
Pion_Cut = 'pionID > 0.1 and abs(dr) < 1 and abs(dz) < 3'
Kaon_Cut = 'kaonID > 0.5 and abs(dr) < 1 and abs(dz) < 3'
Electron_Cut = 'electronID > 0.5 and abs(dr) < 1 and abs(dz) < 3'
Gamma_Cut = 'E >= 0.05'

D0_Cut = "-0.02 <= dM <= 0.02 and useCMSFrame(p) > 2.5"
Dstar0_Cut = "massDifference(0) <= 0.25 and -0.02 <= daughter(0,dM) <= 0.02"
#/================================================================================================================================/
# Create Path
#--------------
my_path = b2.create_path()
#/================================================================================================================================/
# --I/O----------------------------------------------------------------------------------------
# load input ROOT file
ma.inputMdst(environmentType='default',
             filename="/home/belle2/amubarak/C00-Generation/ccbarDstar0EventGeneration.root",
             path=my_path)
#/================================================================================================================================/
# Create Particle Lists
#---------------------------------
# Pion
# Creates Pion Particle List (and c.c.)
# Cut on Variables: pionID, abs(d0), and abs(z0)
# ma.fillParticleList('pi+','', path=my_path)
ma.fillParticleList('pi+:loose', cut=Pion_Cut, path=my_path)
# stdc.stdPi(listtype='loose', path=my_path)
#
# Kaon
# Creates Kaon Particle List (and c.c.)
# Cut on Variables: kaonID, abs(d0), and abs(z0)
# ma.fillParticleList('K+','', path=my_path)f
ma.fillParticleList('K-:loose', cut=Kaon_Cut, path=my_path)
# stdc.stdK(listtype='loose', path=my_path)
#
# Electron
# Creates Electron Particle List (and c.c.)
# Cut on Variables: electronID, abs(d0), and abs(z0)
# ma.fillParticleList("e+:uncorrected",'', path=my_path)
ma.fillParticleList("e+:uncorrected", Electron_Cut, path=my_path)
#
# Gamma
# Creates gamma Particle List (and c.c.)
ma.fillParticleList("gamma:reco", Gamma_Cut, path=my_path)
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
ma.fillParticleList("e+:V0_full",'isDescendantOfList(gamma:V0_full,1)==1', path=my_path)
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

# Partial Reconstruction
#-------------------------
# ma.fillParticleList("e+:Partial",'isDescendantOfList(gamma:merged,1)==1', path=my_path)
ma.cutAndCopyList('e+:Partial', 'e+:corrected', 'isDescendantOfList(gamma:RD,1)==1', path=my_path)
# ma.reconstructDecay(decayString='gamma:RD_partial -> e+:Partial ...',
#                     cut='',
#                     chargeConjugation=ChargeC,
#                     allowChargeViolation=True,
#                     path=my_path)
ma.reconstructDecay(decayString='gamma:RD_partial =direct=> e+:corrected ...',
                    cut='',
                    chargeConjugation=ChargeC,
                    allowChargeViolation=True,
                    path=my_path)
#/================================================================================================================================/
# Reconstruct D0 -> K- pi+ decay
#------------------------------------
ma.reconstructDecay(decayString='D0:kpi -> K-:loose pi+:loose', 
                    cut=D0_Cut, 
                    chargeConjugation=ChargeC, 
                    path=my_path)
#chargeConjugation=False: The option is true by default
#
# Vertex Fitting D0
#-------------------
vx.treeFit('D0:kpi',
           conf_level=0,
           updateAllDaughters=False,
           path=my_path)

# Reconstruct pi0-> gamma gamma
#--------------------------------------------------
ma.reconstructDecay(decayString='pi0:2photons =direct=> gamma:RD_partial ?gamma',
                    cut='',
                    chargeConjugation=ChargeC,
                    allowChargeViolation=True,
                    path=my_path)
# ma.applyCuts('pi0:ga', pi0_Cut, path=my_path)

# Reconstruct pi0:dalitz -> e+:corrected e-:corrected gamma:reco
#-----------------------------------------------------------------
ma.reconstructDecay(decayString='pi0:fulldalitz -> e+:corrected e-:corrected gamma:reco',
                    cut='0.080 < M < 0.200',
                    chargeConjugation=ChargeC,
                    path=my_path)
# ma.fillParticleList('e+:Dalitz','isDescendantOfList(pi0:fulldalitz,1)==1', path=my_path)
ma.cutAndCopyList('e+:Dalitz', 'e+:corrected', 'isDescendantOfList(pi0:fulldalitz,1)==1', path=my_path)
# ma.replaceMass("e+:Dalitz", particleLists="e+:Dalitz", pdgCode=211, path=my_path)
# ma.reconstructDecay(decayString='pi0:dalitz -> e+:Dalitz ... ?gamma',
#                     cut='',
#                     chargeConjugation=ChargeC,
#                     allowChargeViolation=True,
#                     path=my_path)
ma.reconstructDecay(decayString='pi0:dalitz =direct=> e+:corrected ... ?gamma',
                    cut='',
                    chargeConjugation=ChargeC,
                    allowChargeViolation=True,
                    path=my_path)
# ma.applyCuts('pi0:dalitz', pi0dalitz_Cut, path=my_path)

# Reconstruct D*0
#-----------------------------------------
ma.reconstructDecay(decayString='D*0:Ch1 =direct=> [D0:kpi -> K-:loose pi+:loose ] pi0:2photons', 
                    cut=Dstar0_Cut,
                    chargeConjugation=ChargeC,
                    allowChargeViolation=True,
                    path=my_path)
ma.reconstructDecay(decayString='D*0:Ch2 =direct=> [D0:kpi -> K-:loose pi+:loose ] pi0:dalitz', 
                    cut=Dstar0_Cut,
                    chargeConjugation=ChargeC,
                    allowChargeViolation=True, 
                    path=my_path) 
ma.reconstructDecay(decayString='D*0:Ch3 =direct=> [D0:kpi -> K-:loose pi+:loose ] gamma:RD_partial', 
                    cut=Dstar0_Cut,
                    chargeConjugation=ChargeC,
                    allowChargeViolation=True,
                    path=my_path)
# merge the D*0 lists together into one single list
ma.copyLists(outputListName='D*0:all',
             inputListNames=['D*0:Ch1', 'D*0:Ch2', 'D*0:Ch3'],
             path=my_path)
ma.applyCuts("D*0:all", 
             cut=Dstar0_Cut, 
             path=my_path)
#/================================================================================================================================/
# Perform MC Matching (MC truth asociation)
#----------------------------------------------
ma.matchMCTruth(list_name='D*0:Ch1', path=my_path)
ma.matchMCTruth(list_name='D*0:Ch2', path=my_path)
ma.matchMCTruth(list_name='D*0:Ch3', path=my_path)
ma.matchMCTruth(list_name='D*0:all', path=my_path)

# # Perform MC Matching (MC truth asociation)
# ma.matchMCTruth(list_name='gamma:merged', path=my_path)
#/================================================================================================================================/
# D_s+ Variables
#-------------------

# Channel 1
Dstar0CH1_vars=[]
Dstar0CH1_vars += vu.create_aliases_for_selected(
                list_of_variables=['M','InvM','massDifference(0)','mcPDG','isSignal','genNMissingDaughter(11)','genNMissingDaughter(22)'],
                decay_string='^D*0:Ch1 =direct=> [D0:kpi -> K-:loose pi+:loose ] pi0:2photons',
                prefix=['Dstar0_1'])
Dstar0CH1_vars += vu.create_aliases_for_selected(
                list_of_variables=['M','InvM','p','mcPDG'],
                decay_string='D*0:Ch1 =direct=> [^D0:kpi -> K-:loose pi+:loose ] pi0:2photons',
                prefix=['D0_1'])
Dstar0CH1_vars += vu.create_aliases_for_selected(
                list_of_variables=['M','InvM','p','mcPDG','daughter(0,genNMissingDaughter(11))','genNMissingDaughter(22)'],
                decay_string='D*0:Ch1 =direct=> [D0:kpi -> K-:loose pi+:loose ] ^pi0:2photons',
                prefix=['pi0_1'])

# Channel 2
Dstar0CH2_vars=[]
Dstar0CH2_vars += vu.create_aliases_for_selected(
                list_of_variables=['M','InvM','massDifference(0)','mcPDG','isSignal','genNMissingDaughter(11)','genNMissingDaughter(22)'],
                decay_string='^D*0:Ch2 =direct=> [D0:kpi -> K-:loose pi+:loose ] pi0:dalitz',
                prefix=['Dstar0_2'])
Dstar0CH2_vars += vu.create_aliases_for_selected(
                list_of_variables=['M','InvM','p','mcPDG'],
                decay_string='D*0:Ch2 =direct=> [^D0:kpi -> K-:loose pi+:loose ] pi0:dalitz',
                prefix=['D0_2'])
Dstar0CH2_vars += vu.create_aliases_for_selected(
                list_of_variables=['M','InvM','p','mcPDG','genNMissingDaughter(11)','genNMissingDaughter(22)'],
                decay_string='D*0:Ch2 =direct=> [D0:kpi -> K-:loose pi+:loose ] ^pi0:dalitz',
                prefix=['pi0_2'])

# Channel 3
Dstar0CH3_vars=[]
Dstar0CH3_vars += vu.create_aliases_for_selected(
                list_of_variables=['M','InvM','massDifference(0)','mcPDG','isSignal','nMCDaughters','genNMissingDaughter(11)','genNMissingDaughter(22)'],
                decay_string='^D*0:Ch3 =direct=> [D0:kpi -> K-:loose pi+:loose ] gamma:RD_partial',
                prefix=['Dstar0_3'])
Dstar0CH3_vars += vu.create_aliases_for_selected(
                list_of_variables=['M','InvM','p','mcPDG'],
                decay_string='D*0:Ch3 =direct=> [^D0:kpi -> K-:loose pi+:loose ] gamma:RD_partial',
                prefix=['D0_3'])
Dstar0CH3_vars += vu.create_aliases_for_selected(
                list_of_variables=['M','InvM','p','mcPDG','nMCDaughters','daughter(0,mcSecPhysProc)','genNMissingDaughter(11)','genNMissingDaughter(22)'],
                decay_string='D*0:Ch3 =direct=> [D0:kpi -> K-:loose pi+:loose ] ^gamma:RD_partial',
                prefix=['gamma_3'])
#/==========================================================================================================================/
# # Overall Cut:
# #------------------
# ma.applyCuts('D*0:Ch1', "pi0_1_daughter_0_genNMissingDaughter_11==1 and pi0_1_genNMissingDaughter_22==1", path=my_path)
# ma.applyCuts('D*0:Ch2', "pi0_2_genNMissingDaughter_11==1 and pi0_2_genNMissingDaughter_22==1", path=my_path)
# ma.applyCuts('D*0:Ch3', "gamma_3_genNMissingDaughter_11==1", path=my_path)
# ma.applyCuts('D*0:all', "Dstar0_isSignal==1", path=my_path)
#/==========================================================================================================================/
# Save
#--------
# Saving variables to ntuple
output_file = 'Reconstruction_Dstar0.root'

# D*0
ma.variablesToNtuple('D*0:all', ['M','isSignal','massDifference(0)','useAlternativeDaughterHypothesis(massDifference(0), 1:pi+)','mcPDG','dM'],
                     filename=output_file, treename='D*0tree', path=my_path)

# Channel 3
ma.variablesToNtuple('D*0:Ch3', Dstar0CH3_vars,
                     filename=output_file, treename='D*0treeCH3', path=my_path)
# Channel 2
ma.variablesToNtuple('D*0:Ch2', Dstar0CH2_vars,
                     filename=output_file, treename='D*0treeCH2', path=my_path)
# Channel 1
ma.variablesToNtuple('D*0:Ch1', Dstar0CH1_vars,
                     filename=output_file, treename='D*0treeCH1', path=my_path)
#/==========================================================================================================================/

# Process the events
b2.process(my_path)

# print out the summary
print(b2.statistics)