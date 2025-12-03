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
# D*0 -> D0 gamma                                                        #
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
b2.conditions.prepend_globaltag('user_adilmub_FakeD0Suppression') 

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

D0_Cut = "-0.02 <= dM <= 0.02 and useCMSFrame(p) > 2.5"
Dstar0_Cut = "0 <= massDifference(0) <= 0.25 and -0.02 <= daughter(0,dM) <= 0.02"

# Veto
#-------
ElectronROE_Cut = 'isInRestOfEvent == 1 and abs(dr) < 20.0 and abs(dz) < 20.0' # Cuts: abs(d0) < 20.0 and abs(z0) < 20.0
GammaROE_Cut = 'isInRestOfEvent == 1 and E >= 0.05' # No Energy cut (E >= 0.05)
pi0_Cut = ''
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
ma.applyCuts('D0:kpi', D0_Cut, path=my_path)

# Photon Conversion
#--------------------------
ma.reconstructDecay(decayString='gamma:RD_partial =direct=> e+:corrected ...',
                    cut='',
                    chargeConjugation=ChargeC,
                    allowChargeViolation=True,
                    path=my_path)
ma.applyCuts('gamma:RD_partial', '', path=my_path)

# Reconstruct D*0
#-----------------------------------------
ma.reconstructDecay(decayString='D*0 =direct=> [D0:kpi -> K-:loose pi+:loose ] [gamma:RD_partial =direct=> e+:corrected ...]', 
                    cut=Dstar0_Cut,
                    chargeConjugation=ChargeC,
                    allowChargeViolation=True,
                    path=my_path)
#chargeConjugation=False: The option is true by default
#
# Vertex Fitting D0
#-------------------
vx.treeFit('D*0',
           conf_level=0,
           ipConstraint=True,
           updateAllDaughters=False, 
           path=my_path)
ma.applyCuts('D*0', Dstar0_Cut, path=my_path)

# Best Candidate Selection
#---------------------------
# perform best candidate selection

# # D*0 Best Candidate
# vm.addAlias('electron_pt', 'daughter(1,pt)')
# ma.rankByHighest('D*0', 
#                 variable='electron_pt',
#                 # numBest=1,
#                 path=my_path)
# vm.addAlias('electron_pt_rank', 'extraInfo(electron_pt_rank)')

ma.rankByHighest('D*0', 
                variable='chiProb',
                outputVariable='chiProb_Ds_rank',
                numBest=1,
                path=my_path)
# vm.addAlias('chiProb_Ds_rank', 'extraInfo(chiProb_Ds_rank)')

# ma.rankByHighest('D*0', 
#                 variable='random',
#                 # numBest=1,
#                 path=my_path)
# vm.addAlias('random_rank', 'extraInfo(random_rank)')
#/================================================================================================================================/
# Perform MC Matching (MC truth asociation)
#----------------------------------------------
ma.matchMCTruth(list_name='D*0', path=my_path)

# # Perform MC Matching (MC truth asociation)
# ma.matchMCTruth(list_name='gamma:merged', path=my_path)
#/================================================================================================================================/
# D*0 VETO Starts Here
#------------------------

# ROE
#---------------------------------------------------------------------------------------------
# build RestOfEvent (ROE) object for each D*0 candidate
# ROE is required by the veto
ma.buildRestOfEvent(target_list_name='D*0',
                    fillWithMostLikely=False,
                    # chargedPIDPriors=[1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    path=my_path)
# Create a mask tuple:
track_based_cuts = "abs(d0) < 20.0 and abs(z0) < 20.0"
ecl_based_cuts = ""  # [E30]
roe_mask = ("my_mask", track_based_cuts, ecl_based_cuts)
ma.appendROEMasks("D*0", [roe_mask], path=my_path)  # [E40]

ma.updateROEUsingV0Lists(target_particle_list="D*0", 
                         mask_names="my_mask", 
                         default_cleanup=True, 
                         selection_cuts=None, 
                         apply_mass_fit=True, 
                         fitter='treefit', 
                         path=my_path)
#---------------------------------------------------------------------------------------------

# Create a new path (called ROE path) which will be executed for
# each ROE in an event.
# Note that ROE exists for each D*0 candidate, so when we loop
# over each ROE, we effectively loop over signal D*0 candidates

roe_path = b2.create_path()

# The ROE objects might in general be related to Particle from multiple
# particle lists therefore we need to check if the current ROE object
# is related to the Particle from our signal decay. If it is not
# the execution of roe_path will be finished (by starting empty,
# dead end path). Note that in this example this x-check is not
# necessary, but is anyway added for sake of completeness
deadEndPath = b2.create_path()

# Note again: all actions (modules) included in roe_path will be
# executed for each ROE in the event
# First we check that the current ROE is related to D*0 candidate
ma.signalSideParticleFilter(particleList='D*0',
                            selection='',
                            roe_path=roe_path,
                            deadEndPath=deadEndPath)

# Particle Lists
#---------------------------------------------------------------------------------------------
# Electron
# Creates Electron Particle List (and c.c.)
# all electrons found in ROE that are not used to reconstruct current D*0 candidate
# (one can add any cut)
ma.fillParticleList(decayString='e-:roe',
                    cut=ElectronROE_Cut,
                    path=roe_path)

# Gamma
# Creates gamma Particle List (and c.c.)
ma.fillParticleList("gamma:roe", 
                    cut=GammaROE_Cut,
                    path=roe_path)
ma.getFakePhotonProbability('gamma:roe',weight='MC15rd',path=roe_path) # needs prepend_globaltag for MC15ri/MC15rd
ma.getBeamBackgroundProbability('gamma:roe',weight='MC15rd',path=roe_path) # prepend_globaltag for MC15ri/MC15rd
#---------------------------------------------------------------------------------------------

# in order to be able to use modularAnalysis functions (reconstructDecay in particular)
# we need a ParticleList containing the photon candidate used to reconstruct the
# current D*0 as well
# The DecayString is used to specify the selected particle (^)
ma.fillSignalSideParticleList(outputListName='e+:sig',
                              decayString='D*0 =direct=> [D0:kpi -> K-:loose pi+:loose ] [gamma:RD_partial =direct=> ^e+:corrected ...]',
                              path=roe_path)

# Conversion Veto
#-------------------------------------------------------------------------------------------
# Reconstruct gamma:veto -> e+:sig e-:roe
#-----------------------------------------------
# make combinations of signal electron candidates with all electron from ROE
ma.reconstructDecay(decayString='gamma:veto -> e+:sig e-:roe', 
                    cut='',
                    chargeConjugation=ChargeC, 
                    path=roe_path)
#
# Vertex Fitting gamma:veto
#----------------------------------
vx.treeFit('gamma:veto', 
           conf_level=-1,
           updateAllDaughters=True,
           path=roe_path)
ma.applyCuts('gamma:veto', 'InvM <= 2 and M<=3', path=roe_path)

# Reconstruct gamma:singleT -> e+:sig
#-----------------------------------------------
# ma.reconstructDecay(decayString='gamma:reject -> e+:sig e-:roe', 
#                     cut='M>=0.1',
#                     chargeConjugation=ChargeC, 
#                     path=roe_path)
# ma.cutAndCopyList('e+:singleT', 'e+:sig', 'isDescendantOfList(gamma:reject,1)==1', path=roe_path)
ma.reconstructDecay(decayString='gamma:singleT =direct=> e+:sig', 
                    cut='',
                    chargeConjugation=ChargeC,
                    allowChargeViolation=True,
                    path=roe_path)
# #
# # Vertex Fitting gamma:singleT
# #----------------------------------
# vx.treeFit('gamma:singleT', 
#            conf_level=-1,
#         #    massConstraint=['gamma'],
#            updateAllDaughters=False, 
#            path=roe_path)
ma.applyCuts('gamma:singleT', '', path=roe_path)

# # Reconstruct Compton Scattering
# #-----------------------------------------------
# ma.reconstructDecay(decayString='gamma:CS -> e+:sig gamma:roe', 
#                     cut='',
#                     chargeConjugation=ChargeC,
#                     allowChargeViolation=True,
#                     path=roe_path)
# #
# # Vertex Fitting gamma:singleT
# #----------------------------------
# vx.treeFit('gamma:CS', 
#            conf_level=-1,
#         #    massConstraint=['gamma'],
#            updateAllDaughters=False, 
#            path=roe_path)
# ma.applyCuts('gamma:CS', '', path=roe_path)
#-------------------------------------------------------------------------------------------

# Dalitz Decay
#-------------------------------------------------------------------------------------------
# Reconstruct pi0:veto -> e+:sig e-:roe gamma:roe
#--------------------------------------------------
# make combinations of signal electron candidates with all electron from ROE and gamma
ma.reconstructDecay(decayString='pi0:veto -> e+:sig e-:roe gamma:roe',
                    cut='',
                    chargeConjugation=ChargeC,
                    path=roe_path)
vx.treeFit('pi0:veto', 
           conf_level=0,
           updateAllDaughters=False, 
           path=roe_path)
ma.applyCuts('pi0:veto', cut="0.080 < M < 0.200", path=roe_path)
#-------------------------------------------------------------------------------------------

# Perform MC Matching (MC truth asociation)
ma.matchMCTruth(list_name='gamma:veto', path=roe_path)
ma.matchMCTruth(list_name='gamma:singleT', path=roe_path)
# ma.matchMCTruth(list_name='gamma:CS', path=roe_path)
ma.matchMCTruth(list_name='pi0:veto', path=roe_path)

# Best Candidate Selection
#-------------------------------------------------------
# perform best candidate selection
vm.addAlias("gamma_abs_dM", "formula(abs(M - 0.0))")
ma.rankByLowest('gamma:veto', 
                variable='gamma_abs_dM',
                numBest=1, 
                path=roe_path)
ma.rankByLowest('gamma:singleT', 
                variable='gamma_abs_dM',
                numBest=1, 
                path=roe_path)
# ma.rankByHighest('gamma:CS', 
#                  variable='random',
#                  numBest=1, 
#                  path=roe_path)
ma.rankByLowest(particleList='pi0:veto',
                variable='abs(dM)',
                numBest=1,
                path=roe_path)
#-------------------------------------------------------

# Electron MVA Tool
#-----------------------------------------------------------------------------
# Electron
ma.applyChargedPidMVA(particleLists=["e-:roe"], 
                      path=roe_path, 
                      trainingMode=0, 
                      chargeIndependent=False, 
                      binaryHypoPDGCodes=(11, 211))
#-----------------------------------------------------------------------------

# Veto Variables
#-----------------------------------------------------------------------------------------------------------------------------------
# Extra Variables
#-----------------
vm.addAlias("pointangle", "formula(acos(cosAngleBetweenMomentumAndVertexVector))")

vm.addAlias("top", "formula((daughter(0,px)*daughter(1,px)) + (daughter(0,py)*daughter(1,py)) + (daughter(0,pz)*daughter(1,pz)))")
vm.addAlias("bottom", "formula(daughter(0,p)*daughter(1,p))")
vm.addAlias("openangle", "acos(formula(top/bottom))")

vm.addAlias("phi_diff", "daughterDiffOf(0, 1, phi)")
vm.addAlias("theta_diff", "daughterDiffOf(0, 1, theta)")
vm.addAlias("psi", "asin(formula(theta_diff/openangle))")

vm.addAlias("Goodgamma", "passesCut(mcPDG==22)")
vm.addAlias("GoodElectron", "passesCut(mcPDG==11)")

vm.addAlias("Truegamma", "passesCut(M<=0.1)")

# Saving Variables
#--------------------------------------------------------------------------------------------
# Gamma:veto Dictionary
gamma_dict = {'M': 'gammaveto_M',
              'InvM':'gammaveto_InvM',
              'E': 'gammaveto_E',
              'mcE': 'gammaveto_mcE',
              'p': 'gammaveto_p',
              'mcP': 'gammaveto_mcP',
              'pt': 'gammaveto_pt',
              'convertedPhotonInvariantMass(0, 1)':'gammaveto_ConvertM',
              'distance':'gammaveto_distance',
              'mcDecayVertexFromIPDistance':'gammaveto_mcDecayVertexFromIPDistance',
              'daughterAngle(0,1)':'gammaveto_daughterAngle',
              'decayAngle(0)':'gammaveto_decayAngle_0',
              'decayAngle(1)':'gammaveto_decayAngle_1',
              'cos(decayAngle(0))':'gammaveto_cos_decayAngle_0',
              'cos(decayAngle(1))':'gammaveto_cos_decayAngle_1',
              'isSignal':'gammaveto_isSignal',
              'isSignalAcceptWrongFSPs':'gammaveto_isSignalAcceptWrongFSPs',
              'mcPDG':'gammaveto_mcPDG',
              'Goodgamma':'gammaveto_Goodgamma',
              'Truegamma':'gammaveto_Truegamma',
              'genMotherPDG':'gammaveto_genMotherPDG',
              'genMotherPDG(1)':'gammaveto_genMotherPDG_1',
              'nMCDaughters':'gammaveto_nMCDaughters',
              'mcErrors':'gammaveto_mcErrors',
              'daughter(0,M)':'gammaveto_ep_M',
              'daughter(0,p)':'gammaveto_ep_p',
              'daughter(0,mcP)':'gammaveto_ep_mcP',
              'daughter(0,pt)':'gammaveto_ep_pt',
              'daughter(0,E)':'gammaveto_ep_E',
              'daughter(0,mcE)':'gammaveto_ep_mcE',
              'daughter(0,mcPDG)':'gammaveto_ep_mcPDG',
              'daughter(0,GoodElectron)':'gammaveto_ep_GoodElectron',
              'daughter(0,genMotherPDG)':'gammaveto_ep_genMotherPDG',
              'daughter(0,mcMother(E))':'gammaveto_ep_mcMother_E',
              'daughter(0,mcMother(p))':'gammaveto_ep_mcMother_p',
              'daughter(0,mcMother(nMCDaughters))':'gammaveto_ep_mcMother_nMCDaughters',
              'daughter(0,mcMother(mcDaughter(0, PDG)))':'gammaveto_ep_mcMother_daughter_0',
              'daughter(0,mcMother(mcDaughter(1, PDG)))':'gammaveto_ep_mcMother_daughter_1',
              'daughter(0,mcMother(mcDaughter(2, PDG)))':'gammaveto_ep_mcMother_daughter_2',
              'daughter(1,electronID)':'gammaveto_em_electronID',
              'daughter(1,d0)':'gammaveto_em_d0',
              'daughter(1,z0)':'gammaveto_em_z0',
              'daughter(1,pValue)':'gammaveto_em_pValue',
              'daughter(1,firstCDCLayer)':'gammaveto_em_firstCDCLayer',
              'daughter(1,firstPXDLayer)':'gammaveto_em_firstPXDLayer',
              'daughter(1,firstSVDLayer)':'gammaveto_em_firstSVDLayer',
              'daughter(1,M)':'gammaveto_em_M',
              'daughter(1,p)':'gammaveto_em_p',
              'daughter(1,mcP)':'gammaveto_em_mcP',
              'daughter(1,pt)':'gammaveto_em_pt',
              'daughter(1,E)':'gammaveto_em_E',
              'daughter(1,mcE)':'gammaveto_em_mcE',
              'daughter(1,mcPDG)':'gammaveto_em_mcPDG',
              'daughter(1,GoodElectron)':'gammaveto_em_GoodElectron',
              'daughter(1,genMotherPDG)':'gammaveto_em_genMotherPDG',
              'daughter(1,mcSecPhysProc)':'gammaveto_em_mcSecPhysProc',
              'daughter(1,seenInPXD)':'gammaveto_em_seenInPXD',
              'daughter(1,seenInSVD)':'gammaveto_em_seenInSVD',
              'daughter(1,seenInCDC)':'gammaveto_em_seenInCDC',
              'daughter(1,binaryPID(11,211))':'gammaveto_em_binaryPID',
              'daughter(1,pidPairChargedBDTScore(11, 211, All))':'gammaveto_em_pidPairChargedBDTScore',
              'pointangle':'gammaveto_pointangle',
              'chiProb':'gammaveto_chiProb',
              'openangle':'gammaveto_openangle',
              'psi':'gammaveto_psi',
              'genNMissingDaughter(11)':'gammaveto_genNMissingDaughter'
              }
ma.variableToSignalSideExtraInfo(particleList='gamma:veto', 
                                 varToExtraInfo=gamma_dict, 
                                 path=roe_path)
# Gamma:ROE Dictionary
gammaROE_dict = {'M': 'gammaROE_M',
                 'InvM':'gammaROE_InvM',
                 'E': 'gammaROE_E',
                 'decayAngle(0)':'gammaROE_decayAngle_0',
                 'cos(decayAngle(0))':'gammaROE_cos_decayAngle_0',
                 'isSignal':'gammaROE_isSignal',
                 'isSignalAcceptWrongFSPs':'gammaROE_isSignalAcceptWrongFSPs',
                 'mcPDG':'gammaROE_mcPDG',
                 'Goodgamma':'gammaROE_Goodgamma',
                 'genMotherPDG':'gammaROE_genMotherPDG',
                 'genMotherPDG(1)':'gammaROE_genMotherPDG_1',
                 'genMotherPDG(2)':'gammaROE_genMotherPDG_2',
                 'nMCDaughters':'gammaROE_nMCDaughters',
                 'mcErrors':'gammaROE_mcErrors',
                 'daughter(0,p)':'gammaROE_ep_p',
                 'daughter(0,pt)':'gammaROE_ep_pt',
                 'daughter(0,mcPDG)':'gammaROE_ep_mcPDG',
                 'daughter(0,genMotherPDG)':'gammaROE_ep_genMotherPDG',
                 'daughter(0,mcSecPhysProc)':'gammaROE_ep_mcSecPhysProc',
                 'pointangle':'gammaROE_pointangle',
                 'chiProb':'gammaROE_chiProb',
                 'genNMissingDaughter(11)':'gammaROE_genNMissingDaughter'
                }
ma.variableToSignalSideExtraInfo(particleList='gamma:singleT', 
                                 varToExtraInfo=gammaROE_dict, 
                                 path=roe_path)
# # Compton Scattering
# gammaCS_dict = {'M': 'gammaCS_M',
#                 'InvM':'gammaCS_InvM',
#                 'distance':'gammaCS_distance',
#                 'mcDecayVertexFromIPDistance':'gammaCS_mcDecayVertexFromIPDistance',
#                 'E': 'gammaCS_E',
#                 'daughterAngle(0,1)':'gammaCS_daughterAngle',
#                 'decayAngle(0)':'gammaCS_decayAngle_0',
#                 'isSignal':'gammaCS_isSignal',
#                 'isSignalAcceptWrongFSPs':'gammaCS_isSignalAcceptWrongFSPs',
#                 'mcPDG':'gammaCS_mcPDG',
#                 'Goodgamma':'gammaCS_Goodgamma',
#                 'genMotherPDG':'gammaCS_genMotherPDG',
#                 'genMotherPDG(1)':'gammaCS_genMotherPDG_1',
#                 'genMotherPDG(2)':'gammaCS_genMotherPDG_2',
#                 'nMCDaughters':'gammaCS_nMCDaughters',
#                 'mcErrors':'gammaCS_mcErrors',
#                 'daughter(0,p)':'gammaCS_ep_p',
#                 'daughter(0,pt)':'gammaCS_ep_pt',
#                 'daughter(0,mcPDG)':'gammaCS_ep_mcPDG',
#                 'daughter(0,genMotherPDG)':'gammaCS_ep_genMotherPDG',
#                 'daughter(0,mcSecPhysProc)':'gammaCS_ep_mcSecPhysProc',
#                 'daughter(1,p)':'gammaCS_g_p',
#                 'daughter(1,pt)':'gammaCS_g_pt',
#                 'daughter(1,mcPDG)':'gammaCS_g_mcPDG',
#                 'daughter(1,genMotherPDG)':'gammaCS_g_genMotherPDG',
#                 'daughter(1,mcSecPhysProc)':'gammaCS_g_mcSecPhysProc',
#                 'daughter(1,beamBackgroundSuppression)':'gammaCS_beamBackgroundSuppression',
#                 'daughter(1,fakePhotonSuppression)':'gammaCS_fakePhotonSuppression',
#                 'pointangle':'gammaCS_pointangle',
#                 'chiProb':'gammaCS_chiProb',
#                 }
# ma.variableToSignalSideExtraInfo(particleList='gamma:CS', 
#                                  varToExtraInfo=gammaCS_dict, 
#                                  path=roe_path)
# pi0 Dictionary
pi0_dict = {'M': 'pi0veto_M',
            'useAlternativeDaughterHypothesis(M, 1:pi+)':'pi0veto_M_epi',
            'useAlternativeDaughterHypothesis(M, 0:pi+, 1:pi+)':'pi0veto_M_pipi',
            'dM': 'pi0veto_dM',
            'E': 'pi0veto_E',
            'distance':'pi0veto_distance',
            'chiProb':'pi0veto_chiProb',
            'mcPDG':'pi0veto_mcPDG',
            'mcErrors':'pi0veto_mcErrors',
            'genMotherPDG':'pi0veto_genMotherPDG',
            'nMCDaughters':'pi0veto_nMCDaughters',
            'nDaughterPhotons':'pi0veto_nDaughterPhotons',
            'daughter(2,beamBackgroundSuppression)':'pi0veto_beamBackgroundSuppression',
            'daughter(2,fakePhotonSuppression)':'pi0veto_fakePhotonSuppression',
            }
ma.variableToSignalSideExtraInfo(particleList='pi0:veto', 
                                 varToExtraInfo=pi0_dict, 
                                 path=roe_path)
#--------------------------------------------------------------------------------------------

# execute roe_path for each RestOfEvent in the event
my_path.for_each('RestOfEvent', 'RestOfEvents', roe_path)

# Variable List
#--------------------------------------------------------------------------------------------
gamma_ROE = []
var_V0 = []
for key in gamma_dict:
    vm.addAlias(gamma_dict[key], "extraInfo("+ gamma_dict[key] +")")
    gamma_ROE.append(gamma_dict[key])
    var_V0.append(key)

gammaROE_ROE = []
for key in gammaROE_dict:
    vm.addAlias(gammaROE_dict[key], "extraInfo("+ gammaROE_dict[key] +")")
    gammaROE_ROE.append(gammaROE_dict[key])

# gammaCS_ROE = []
# for key in gammaCS_dict:
#     vm.addAlias(gammaCS_dict[key], "extraInfo("+ gammaCS_dict[key] +")")
#     gammaCS_ROE.append(gammaCS_dict[key])

pi0_ROE = []
for key in pi0_dict:
    vm.addAlias(pi0_dict[key], "extraInfo("+ pi0_dict[key] +")")
    pi0_ROE.append(pi0_dict[key])

vm.addAlias('gammaveto_M_Correction', 'ifNANgiveX(gammaveto_M,10)')
gamma_ROE.append('gammaveto_M_Correction')
# vm.addAlias('M_Correction', 'ifNANgiveX(M,10)')
# vm.addAlias('pidPairChargedBDTScore_Correction', 'ifNANgiveX(daughter(1,pidPairChargedBDTScore(11, 211, All)),10)')
# var_V0.append('M_Correction')
# var_V0.append('pidPairChargedBDTScore_Correction')

vm.addAlias('gammaROE_M_Correction', 'ifNANgiveX(gammaROE_M,10)')
gammaROE_ROE.append('gammaROE_M_Correction')

# vm.addAlias('gammaCS_M_Correction', 'ifNANgiveX(gammaCS_M,10)')
# gammaCS_ROE.append('gammaCS_M_Correction')

vm.addAlias('pi0veto_M_Correction', 'ifNANgiveX(pi0veto_M,10)')
pi0_ROE.append('pi0veto_M_Correction')
#--------------------------------------------------------------------------------------------
# D*0 Veto Ends Here
#/================================================================================================================================/
# Other Tools:
#-------------------

# # MVAExpert
# my_path.add_module('MVAExpert', 
#                    listNames=['D*0'], 
#                    extraInfoName='FastBDT', 
#                    identifier='user_adilmub_FakeD0Suppression')
# # Variables from MVAExpert.
# expertVars = ['extraInfo(FastBDT)'] #'transformedNetworkOutput(FastBDT,0.1,1.0)'
#/================================================================================================================================/
# D*0 Variables
#-------------------

tracks = ['abs(dr)','abs(dz)','z0','d0','pValue',
          'firstCDCLayer','firstPXDLayer','firstSVDLayer',
          'nPXDHits','nVXDHits','nSVDHits','nCDCHits',
          'seenInCDC','seenInPXD','seenInSVD']

Dstar0_vars=[]
Dstar0_vars += vu.create_aliases_for_selected(
                list_of_variables=['electronID','binaryPID(11,211)',
                                   'cos(theta)','phi','mcP','M','pt','p']
                                   + tracks,
                decay_string='D*0 =direct=> [D0:kpi -> K-:loose pi+:loose ] [gamma:RD_partial =direct=> ^e+:corrected ...]',
                prefix=['e'])

Dstar0_vars += vu.create_aliases_for_selected(
                list_of_variables=['M','InvM','p','E','mcPDG','genMotherPDG','nMCDaughters',
                                   'daughterMotherDiffOf(0,p)','daughterMotherDiffOf(0,E)',
                                   'daughter(0,mcSecPhysProc)',
                                   'genNMissingDaughter(11)','genNMissingDaughter(22)'],
                decay_string='D*0 =direct=> [D0:kpi -> K-:loose pi+:loose ] [^gamma:RD_partial =direct=> e+:corrected ...]',
                prefix=['gamma'])
Dstar0_vars += vu.create_aliases_for_selected(
                list_of_variables=['M','InvM','p','mcPDG','genMotherPDG'],
                decay_string='D*0 =direct=> [^D0:kpi -> K-:loose pi+:loose ] [gamma:RD_partial =direct=> e+:corrected ...]',
                prefix=['D0'])
Dstar0_vars += vu.create_aliases_for_selected(
                list_of_variables=['M','InvM','massDifference(0)',
                                   "chiProb",
                                   "cos(theta)",
                                   "decayAngle(0)","cos(decayAngle(0))",
                                   "decayAngle(1)","cos(decayAngle(1))",
                                   'useAlternativeDaughterHypothesis(massDifference(0), 1:pi+)',
                                   'mcPDG','genMotherPDG','isSignal','nMCDaughters',
                                   'genNMissingDaughter(11)','genNMissingDaughter(22)']
                                   + gamma_ROE + gammaROE_ROE + pi0_ROE,
                decay_string='^D*0 =direct=> [D0:kpi -> K-:loose pi+:loose ] [gamma:RD_partial =direct=> e+:corrected ...]',
                prefix=['Dstar0'])
#/==========================================================================================================================/
# # Overall Cut:
# #------------------
# ma.applyCuts('D*0', "Dstar0_isSignal==1", path=my_path)
#/==========================================================================================================================/
# Save
#--------
# Saving variables to ntuple
output_file = 'Reconstruction_Dstar0.root'

# Mode 3
ma.variablesToNtuple('D*0', Dstar0_vars,
                     filename=output_file, treename='D*0tree', path=my_path)
#/==========================================================================================================================/

# Process the events
b2.process(my_path)

# print out the summary
print(b2.statistics)