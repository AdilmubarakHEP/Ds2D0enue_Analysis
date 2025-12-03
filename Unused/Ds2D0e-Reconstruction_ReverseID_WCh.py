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
# D*0 -> D0 gamma/pi0                                                    #
#        |                                                               #
#        +-> K- pi+                                                      #
#                                                                        #
# The main purpose of this code is to obtain the D^*0 peaking background #
##########################################################################

import argparse
import basf2 as b2

import modularAnalysis as ma
import variables.collections as vc
import variables.utils as vu
import vertex as vx
import variables as va
from variables import variables as vm
from stdPi0s import stdPi0s
from stdPhotons import stdPhotons

from variables.MCGenTopo import mc_gen_topo

b2.conditions.prepend_globaltag(ma.getAnalysisGlobaltag()) # wait 180ms, doesnt work when confl down, 
# b2.conditions.prepend_globaltag('user_adilmub_BkgSuppressionBDT_v1_03082025')
# b2.conditions.prepend_globaltag('user_adilmub_BkgSuppression')
# b2.conditions.prepend_globaltag('user_adilmub_FakeD0Suppression')

#/================================================================================================================================/
# Charge Conjugation:
#-----------------------
ChargeC = True
#/================================================================================================================================/
# Cuts
#-----------
Pion_Cut = 'pionID > 0.2 and abs(dr) < 1 and abs(dz) < 3'
Kaon_Cut = 'kaonID > 0.5 and abs(dr) < 1 and abs(dz) < 3'
# Electron_Cut = 'abs(dr) < 1 and abs(dz) < 3'
Electron_Cut = 'electronID <= 0.5 and abs(dr) < 1 and abs(dz) < 3'
Gamma_Cut = 'E >= 0.1'

D0_Cut = "-0.04 <= dM <= 0.04 and useCMSFrame(p) > 2.5"
# D0_Cut = "[-0.23 <= dM <= -0.03 or 0.03 <= dM <= 0.23] and useCMSFrame(p) > 2.5"
Ds_Cut = "massDifference(0) <= 0.35 and -0.02 <= daughter(0,dM) <= 0.02"
# Ds_Cut = "useAlternativeDaughterHypothesis(massDifference(0), 1:pi+) <= 0.6 and [-0.2 <= daughter(0,dM) <= -0.04 or 0.04 <= daughter(0,dM) <= 0.2]"
# Ds_Cut = "massDifference(0) <= 0.25"
Ds_star_Cut = "massDifference(0) <= 0.4"

# Veto
#-------
ElectronROE_Cut = 'isInRestOfEvent == 1 and abs(dr) < 20.0 and abs(dz) < 20.0' # Cuts: abs(d0) < 20.0 and abs(z0) < 20.0
PionROE_Cut = 'isInRestOfEvent == 1 and abs(dr) < 20.0 and abs(dz) < 20.0' # Cuts: abs(d0) < 20.0 and abs(z0) < 20.0
GammaROE_Cut = 'isInRestOfEvent == 1 and E > 0.1'
#/================================================================================================================================/

def Reconstruction(path, ChargeC, Pion_Cut, Kaon_Cut, Electron_Cut, Gamma_Cut, D0_Cut, Ds_Cut, Ds_star_Cut):
    #/================================================================================================================================/
    # # Cuts
    # #-----------
    # Pion_Cut = 'pionID > 0.2 and abs(dr) < 1 and abs(dz) < 3'
    # Kaon_Cut = 'kaonID > 0.5 and abs(dr) < 1 and abs(dz) < 3'
    # # Electron_Cut = 'abs(dr) < 1 and abs(dz) < 3'
    # Electron_Cut = 'electronID >= 0.5 and abs(dr) < 1 and abs(dz) < 3'
    # Gamma_Cut = 'E >= 0.1'

    # D0_Cut = "-0.04 <= dM <= 0.04 and useCMSFrame(p) > 2.5"
    # # D0_Cut = "[-0.23 <= dM <= -0.03 or 0.03 <= dM <= 0.23] and useCMSFrame(p) > 2.5"
    # Ds_Cut = "useAlternativeDaughterHypothesis(massDifference(0), 1:pi+) <= 0.6 and -0.02 <= daughter(0,dM) <= 0.02"
    # # Ds_Cut = "useAlternativeDaughterHypothesis(massDifference(0), 1:pi+) <= 0.6 and [-0.2 <= daughter(0,dM) <= -0.04 or 0.04 <= daughter(0,dM) <= 0.2]"
    # # Ds_Cut = "massDifference(0) <= 0.25"
    # Ds_star_Cut = "massDifference(0) <= 0.4"

    # # Veto
    # #-------
    # ElectronROE_Cut = 'isInRestOfEvent == 1 and abs(dr) < 20.0 and abs(dz) < 20.0' # Cuts: abs(d0) < 20.0 and abs(z0) < 20.0
    # PionROE_Cut = 'isInRestOfEvent == 1 and abs(dr) < 20.0 and abs(dz) < 20.0' # Cuts: abs(d0) < 20.0 and abs(z0) < 20.0
    # GammaROE_Cut = 'isInRestOfEvent == 1 and E > 0.050' # No Energy cut (E >= 0.05)
    # pi0_Cut = ''
    # #/================================================================================================================================/
    # # Photon Conversion
    # #----------------------

    # # Reconstruct gamma:V0_array -> e+ e-
    # #-----------------------------------------------
    # ma.fillConvertedPhotonsList(decayString='gamma:V0_array -> e+ e-',
    #                             cut='', 
    #                             path=path)
    # #
    # # Vertex Fitting gamma:V0_array 
    # #----------------------------------
    # vx.treeFit('gamma:V0_array', 
    #            conf_level=-1,
    #            updateAllDaughters=True, 
    #            path=path)
    # ma.applyCuts('gamma:V0_array', '-0.5<=M<=0.2 and -0.5<=InvM<=0.2', 
    #              path=path)

    # ma.matchMCTruth(list_name='gamma:V0_array', path=path) 

    # ma.fillParticleList("e+:V0_array",'isDescendantOfList(gamma:V0_array,1)==1', path=path)
    #/================================================================================================================================/
    # Create Particle Lists
    #---------------------------------
    # Pion
    # Creates Pion Particle List (and c.c.)
    # Cut on Variables: pionID, abs(d0), and abs(z0)
    # ma.fillParticleList('pi+','', path=path)
    ma.fillParticleList('pi+:loose', cut=Pion_Cut, path=path)
    # ma.fillParticleListFromMC('pi+:gen','',path=path)
    # stdPi0s(listtype='eff50_May2020', path=path)
    #
    # Kaon
    # Creates Kaon Particle List (and c.c.)
    # Cut on Variables: kaonID, abs(d0), and abs(z0)
    # ma.fillParticleList('K+','', path=path)f
    ma.fillParticleList('K-:loose', cut=Kaon_Cut, path=path)
    # stdc.stdK(listtype='loose', path=path)
    #
    # Electron
    # Creates Electron Particle List (and c.c.)
    # Cut on Variables: electronID, abs(d0), and abs(z0)
    ma.fillParticleList("e+:uncorrected", cut='', path=path)
    # ma.fillParticleListFromMC('e+:gen','',path=path)
    # ma.tagCurlTracks("e+:uncorrected", selectorType='mva', ptCut=0.6, mcTruth=True, expert_filename='test.root', path=path)
    # vm.addAlias("isCurl", "extraInfo(isCurl)") 
    # ma.replaceMass("e+:uncorrected", particleLists="e+:uncorrected", pdgCode=211, path=path) # Does not recalculate anything, still e mass
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

    ma.fillParticleList("gamma:brems", "goodGamma", path=path)

    ma.correctBrems(outputList="e+:corrected", 
                    inputList="e+:uncorrected",
                    gammaList="gamma:brems",
                    path=path)
    vm.addAlias("isBremsCorrected", "extraInfo(bremsCorrected)")  # [E30]

    ma.applyCuts('e+:corrected',Electron_Cut, path=path)

    # ma.tagCurlTracks("e+:corrected", selectorType='mva', ptCut=0.6, mcTruth=True, path=path)
    # vm.addAlias('isCurl', 'extraInfo(isCurl)')
    # vm.addAlias('isTruthCurl', 'extraInfo(isTruthCurl)')
    Curler = []
    # Curler = ['isCurl','isTruthCurl']
    #/================================================================================================================================/
    # MVA
    #------------------------
    # Photon
    stdPhotons('cdc', beamBackgroundMVAWeight="MC15rd", fakePhotonMVAWeight="MC15rd", path=path)

    # apply additional cuts on photon list if needed
    ma.cutAndCopyList("gamma:recon", "gamma:cdc", cut=Gamma_Cut, path=path)
    
    # return the extra info - note this step is not mandatory, but B2WARNING messages will appear if the extra info is not explicitly returned
    vm.addAlias("beamBackgroundSuppressionScore", "extraInfo(beamBackgroundSuppression)")
    vm.addAlias("fakePhotonSuppressionScore", "extraInfo(fakePhotonSuppression)")

    # ma.getFakePhotonProbability('gamma:recon',weight='MC15rd',path=path) # needs prepend_globaltag for MC15ri/MC15rd
    # ma.getBeamBackgroundProbability('gamma:recon',weight='MC15rd',path=path) # prepend_globaltag for MC15ri/MC15rd

    # # Electron
    # ma.applyChargedPidMVA(particleLists=["e+:corrected"], 
    #                       path=path, 
    #                       trainingMode=1, 
    #                       chargeIndependent=False, 
    #                       binaryHypoPDGCodes=(0, 0))
    #/================================================================================================================================/
    # Reconstruct D0 decay
    #------------------------------------
    ma.reconstructDecay(decayString='D0:kmpip -> K-:loose pi+:loose', 
                        cut=D0_Cut, 
                        chargeConjugation=ChargeC,
                        path=path)
    # ma.reconstructDecay(decayString='D0:kmpippi0 -> K-:loose pi+:loose pi0:eff50_May2020', 
    #                     cut="-0.06 <= dM <= 0.06 and useCMSFrame(p) > 2", 
    #                     chargeConjugation=ChargeC,
    #                     path=path)
    # Vertex Fitting D0
    #-------------------
    vx.treeFit('D0:kmpip',
            conf_level=0,
            updateAllDaughters=True,
            path=path)
    ma.applyCuts("D0:kmpip", D0_Cut, path=path)

    # vx.treeFit('D0:kmpippi0',
    #         conf_level=0,
    #         updateAllDaughters=False,
    #         path=path)
    # ma.applyCuts("D0:kmpippi0", cut="-0.05 <= dM <= 0.05 and useCMSFrame(p) > 2.5", path=path)

    # # merge the D0 lists together into one single list
    # ma.copyLists(outputListName='D0:merged',
    #              inputListNames=['D0:kmpip','D0:kmpippi0'],
    #              path=path)

    # Reconstruct D_s+ -> D0 e+ nu_e decay
    #-----------------------------------------
    ma.reconstructDecay(decayString='D_s+ -> [D0:kmpip -> K-:loose pi+:loose] e-:corrected ?nu ?addbrems', 
                        cut=Ds_Cut, 
                        chargeConjugation=ChargeC, 
                        allowChargeViolation=True,
                        path=path) 
    #chargeConjugation=False: The option is true by default
    #
    # Vertex Fitting D_s+
    #-------------------
    # Before Vertex Fitting
    ma.variablesToExtraInfo('D_s+', variables={"diff_D0e": "diff_D0e_noVF", 
                                               "diff_D0pi": "diff_D0pi_noVF",
                                               "diff_D0K": "diff_D0K_noVF",
                                               "cos(daughterAngle(0, 1))":"cos_theta_D0e_noVF"}, 
                                               path=path)
    vm.addAlias('diff_D0e_noVF', 'extraInfo(diff_D0e_noVF)')
    vm.addAlias('diff_D0pi_noVF', 'extraInfo(diff_D0pi_noVF)')
    vm.addAlias('diff_D0K_noVF', 'extraInfo(diff_D0K_noVF)')
    vm.addAlias('cos_theta_D0e_noVF', 'extraInfo(cos_theta_D0e_noVF)')
    Before_VF = ["diff_D0e_noVF","diff_D0pi_noVF","diff_D0K_noVF","cos_theta_D0e_noVF"]

    # After Vertex Fit
    vx.treeFit('D_s+', 
            conf_level=0, 
            ipConstraint=True,
            updateAllDaughters=True, 
            path=path)
    ma.applyCuts("D_s+", Ds_Cut, path=path)

    # D_s+ Best Candidate
    ma.rankByHighest('D_s+', 
                    variable='chiProb',
                    outputVariable='chiProb_Ds_rank',
                    # numBest=1,
                    path=path)
    vm.addAlias('chiProb_Ds_rank', 'extraInfo(chiProb_Ds_rank)')
    ma.rankByHighest('D_s+', 
                    variable='daughter(1,p)',
                    outputVariable='Electron_p_Ds_rank',
                    # numBest=1,
                    path=path)
    vm.addAlias('Electron_p_Ds_rank', 'extraInfo(Electron_p_Ds_rank)')
    ma.rankByHighest('D_s+', 
                    variable='random',
                    outputVariable='random_Ds_rank',
                    # numBest=1,
                    path=path)
    vm.addAlias('random_Ds_rank', 'extraInfo(random_Ds_rank)')
    BCS = ['chiProb_Ds_rank','Electron_p_Ds_rank','random_Ds_rank']
    #/================================================================================================================================/
    # Reconstruct D_s+* -> D_s+ gamma
    #-----------------------------------------
    ma.reconstructDecay(decayString='D_s*+ -> [D_s+ -> [D0:kmpip -> K-:loose pi+:loose] e-:corrected ?nu ?addbrems] gamma:recon', 
                        cut=Ds_star_Cut, 
                        chargeConjugation=ChargeC, 
                        allowChargeViolation=True,
                        path=path)
    #chargeConjugation=False: The option is true by default
    #
    # Vertex Fitting D_s*+
    # ----------------------
    # vx.treeFit('D_s*+', 
    #            conf_level=-1,
    #            path=path)
    ma.applyCuts("D_s*+", Ds_star_Cut, path=path)

    # # D_s*+ Best Candidate
    # vm.addAlias('DeltaM_diff', 'formula(abs(massDifference(0) - 0.1438))')
    # ma.rankByLowest('D_s*+', 
    #                 variable='DeltaM_diff',
    #                 # numBest=1,
    #                 path=path)
    # vm.addAlias('DeltaM_diff_rank', 'extraInfo(DeltaM_diff_rank)')

    # ma.rankByLowest('D_s*+', 
    #                 variable='dM',
    #                 # numBest=1,
    #                 path=path)
    # vm.addAlias('dM_rank', 'extraInfo(dM_rank)')

    # ma.rankByHighest('D_s*+', 
    #                 variable='random',
    #                 # numBest=1,
    #                 path=path)
    # vm.addAlias('random_rank', 'extraInfo(random_rank)')
    #/================================================================================================================================/
    # Perform MC Matching (MC truth asociation)
    #----------------------------------------------
    # ma.matchMCTruth(list_name='D0:kmpip', path=path) 
    # ma.matchMCTruth(list_name='D0:kmpippi0', path=path) 

    # ma.matchMCTruth(list_name='D_s+', path=path) 
    # ma.matchMCTruth(list_name='D_s*+', path=path) 
    #/================================================================================================================================/
    return Curler, Before_VF, BCS

def Rest_of_Event(roe_path, ChargeC, ElectronROE_Cut, PionROE_Cut, GammaROE_Cut):
    #/================================================================================================================================/
    # D_s+ VETO Starts Here
    #------------------------

    # Particle Lists
    #---------------------------------------------------------------------------------------------
    # Electron
    # Creates Electron Particle List (and c.c.)
    # all electrons found in ROE that are not used to reconstruct current D_s+ candidate
    # (one can add any cut)
    ma.fillParticleList(decayString='e+:roe',
                        cut=ElectronROE_Cut,
                        path=roe_path)

    # Gamma
    # Creates gamma Particle List (and c.c.)
    ma.fillParticleList("gamma:roe", 
                        cut=GammaROE_Cut,
                        path=roe_path)
    ma.getFakePhotonProbability('gamma:roe',weight='MC15rd',path=roe_path) # needs prepend_globaltag for MC15ri/MC15rd
    ma.getBeamBackgroundProbability('gamma:roe',weight='MC15rd',path=roe_path) # prepend_globaltag for MC15ri/MC15rd
    ma.applyCuts('gamma:roe', cut="beamBackgroundSuppression >= 0.5 and fakePhotonSuppression >= 0.5", path=roe_path)

    # Pion
    ma.fillParticleList('pi+:ROE', cut=PionROE_Cut, path=roe_path)
    #---------------------------------------------------------------------------------------------

    # in order to be able to use modularAnalysis functions (reconstructDecay in particular)
    # we need a ParticleList containing the photon candidate used to reconstruct the
    # current D_s+ as well
    # The DecayString is used to specify the selected particle (^)
    ma.fillSignalSideParticleList(outputListName='e-:sig',
                                decayString='D_s+ -> [D0:kmpip -> K-:loose pi+:loose] ^e-:corrected ?nu',
                                path=roe_path)
    ma.fillSignalSideParticleList(outputListName='D0:sig',
                                decayString='D_s+ -> [^D0:kmpip -> K-:loose pi+:loose] e-:corrected ?nu',
                                path=roe_path)
    # ma.fillSignalSideParticleList(outputListName='D_s+:sig',
    #                             decayString='^D_s+ -> [D0:kmpip -> K-:loose pi+:loose] e+:corrected ?nu',
    #                             path=roe_path)

    # Conversion Veto
    #-------------------------------------------------------------------------------------------
    # Reconstruct gamma:veto -> e+:sig e-:roe
    #-----------------------------------------------
    # make combinations of signal electron candidates with all electron from ROE
    ma.reconstructDecay(decayString='gamma:veto -> e-:sig e+:roe', 
                        cut='',
                        chargeConjugation=ChargeC, 
                        path=roe_path)
    #
    # Vertex Fitting gamma:veto
    #----------------------------------
    vx.treeFit('gamma:veto', 
            conf_level=-1,
            updateAllDaughters=False,
            path=roe_path)
    ma.applyCuts('gamma:veto', 'InvM <= 2 and M<=3', path=roe_path)

    # # Reconstruct Compton Scattering
    # #-----------------------------------------------
    # ma.reconstructDecay(decayString='gamma:CS -> e+:sig pi-:ROE', 
    #                     cut='',
    #                     chargeConjugation=ChargeC,
    #                     allowChargeViolation=True,
    #                     path=roe_path)
    # #
    # # Vertex Fitting gamma:singleT
    # #----------------------------------
    # vx.treeFit('gamma:CS', 
    #         conf_level=-1,
    #         #    massConstraint=['gamma'],
    #         updateAllDaughters=False, 
    #         path=roe_path)
    # ma.applyCuts('gamma:CS', '', path=roe_path)

    # # Reconstruct gamma:singleT -> e+:sig
    # #-----------------------------------------------
    # ma.reconstructDecay(decayString='gamma:singleT -> e+:sig', 
    #                     cut='',
    #                     chargeConjugation=ChargeC,
    #                     allowChargeViolation=True,
    #                     path=roe_path)
    # ma.applyCuts('gamma:singleT', '', path=roe_path)
    #-------------------------------------------------------------------------------------------

    # pi0 Decay
    #-------------------------------------------------------------------------------------------
    # Reconstruct pi0:ROE -> gamma:roe gamma:roe
    #--------------------------------------------------
    ma.reconstructDecay(decayString='pi0:ROE -> gamma:roe gamma:roe',
                        cut='0.080 < M < 0.200',
                        path=roe_path)
    vx.treeFit('pi0:ROE', 
               conf_level=-1,
               updateAllDaughters=False, 
               path=roe_path)
    ma.applyCuts('pi0:ROE', cut="0.080 < M < 0.200", path=roe_path)
    #-------------------------------------------------------------------------------------------

    # D0 Veto
    #-------------------------------------------------------------------------------------------
    # Reconstruct D*0:veto -> D0:sig X
    #--------------------------------------------------
    # make combinations of signal electron candidates with all electron from ROE and gamma
    ma.reconstructDecay(decayString='D*0:Mode1 -> D0:sig gamma:roe',
                        cut='',
                        chargeConjugation=ChargeC,
                        path=roe_path)
    vx.treeFit('D*0:Mode1', 
            conf_level=-1,
            updateAllDaughters=False, 
            path=roe_path)
    ma.applyCuts('D*0:Mode1', cut="0 <= massDifference(0) <= 0.25", path=roe_path)

    ma.reconstructDecay(decayString='D*0:Mode2 -> D0:sig pi0:ROE',
                        cut='',
                        chargeConjugation=ChargeC,
                        path=roe_path)
    vx.treeFit('D*0:Mode2', 
            conf_level=-1,
            updateAllDaughters=False, 
            path=roe_path)
    ma.applyCuts('D*0:Mode2', cut="0 <= massDifference(0) <= 0.25", path=roe_path)

    ma.reconstructDecay(decayString='D*+:ROE -> D0:sig pi+:ROE',
                        cut='',
                        chargeConjugation=ChargeC,
                        path=roe_path)
    vx.treeFit('D*+:ROE', 
            conf_level=-1,
            updateAllDaughters=False, 
            path=roe_path)
    ma.applyCuts('D*+:ROE', cut="0 <= massDifference(0) <= 0.5", path=roe_path)
    #-------------------------------------------------------------------------------------------

    # # Angle Between D_s+ and Kaons
    # #-------------------------------------------------------------------------------------------
    # # make combinations of signal electron candidates with all electron from ROE and gamma
    # ma.reconstructDecay(decayString='vpho:1 -> D_s+:sig K+:ROE',
    #                     cut='',
    #                     chargeConjugation=ChargeC,
    #                     allowChargeViolation=True,
    #                     path=roe_path)
    # ma.applyCuts('vpho:1', cut="", path=roe_path)
    # ma.variablesToDaughterExtraInfo(particleList='vpho:1', 
    #                                 decayString='vpho:1 -> ^D_s+:sig K+:ROE', 
    #                                 variables={'M': 'vpho_M',
    #                                         'daughterAngle(0,1)':'vpho_daughterAngle',
    #                                         'useCMSFrame(daughterAngle(0,1))':'vpho_CMS_daughterAngle',
    #                                         'cos(daughterAngle(0,1))':'vpho_cos_daughterAngle',
    #                                         'useCMSFrame(cos(daughterAngle(0,1)))':'vpho_cos_CMS_daughterAngle',
    #                                         },
    #                                 option=0, 
    #                                 path=roe_path)

    # ma.reconstructDecay(decayString='vpho:2 -> D_s+:sig K_S0:legacyGoodKS',
    #                     cut='',
    #                     chargeConjugation=ChargeC,
    #                     allowChargeViolation=True,
    #                     path=roe_path)
    # ma.applyCuts('vpho:2', cut="", path=roe_path)
    # ma.variablesToDaughterExtraInfo(particleList='vpho:2', 
    #                                 decayString='vpho:2 -> ^D_s+:sig K+:ROE', 
    #                                 variables={'M': 'vpho_M',
    #                                         'daughterAngle(0,1)':'vpho_daughterAngle',
    #                                         'useCMSFrame(daughterAngle(0,1))':'vpho_CMS_daughterAngle',
    #                                         'cos(daughterAngle(0,1))':'vpho_cos_daughterAngle',
    #                                         'useCMSFrame(cos(daughterAngle(0,1)))':'vpho_cos_CMS_daughterAngle',
    #                                         }, 
    #                                 option=0, 
    #                                 path=roe_path)

    # ma.copyLists('vpho', ['vpho:1', 'vpho:2'], path=roe_path)                               
    # #-------------------------------------------------------------------------------------------

    # Perform MC Matching (MC truth asociation)
    # ma.matchMCTruth(list_name='K+:ROE', path=roe_path)

    # ma.matchMCTruth(list_name='gamma:veto', path=roe_path)
    # ma.matchMCTruth(list_name='gamma:CS', path=roe_path)
    # ma.matchMCTruth(list_name='gamma:singleT', path=roe_path)
    # ma.matchMCTruth(list_name='pi0:veto', path=roe_path)
    # ma.matchMCTruth(list_name='D*0:Mode1', path=roe_path)
    # ma.matchMCTruth(list_name='D*0:Mode2', path=roe_path)
    # ma.matchMCTruth(list_name='D*+:ROE',   path=roe_path)

    # Best Candidate Selection
    #-------------------------------------------------------
    # perform best candidate selection
    vm.addAlias("gamma_abs_dM", "formula(abs(M - 0.0))")
    ma.rankByLowest('gamma:veto', 
                    variable='gamma_abs_dM',
                    numBest=1, 
                    path=roe_path)
    # ma.rankByLowest('gamma:singleT', 
    #                 variable='gamma_abs_dM',
    #                 numBest=1, 
    #                 path=roe_path)
    # ma.rankByHighest('gamma:CS', 
    #                 variable='random',
    #                 numBest=1, 
    #                 path=roe_path)
    # ma.rankByLowest(particleList='pi0:veto',
    #                 variable='abs(dM)',
    #                 numBest=1,
    #                 path=roe_path)
    ma.rankByLowest(particleList='D*0:Mode1',
                    variable='abs(dM)',
                    numBest=1,
                    path=roe_path)
    ma.rankByLowest(particleList='D*0:Mode2',
                    variable='abs(dM)',
                    numBest=1,
                    path=roe_path)
    ma.rankByLowest(particleList='D*+:ROE',
                    variable='abs(dM)',
                    numBest=1,
                    path=roe_path)
    #-------------------------------------------------------

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
                'MomentumAsymmetry':'gammaveto_MomentumAsymmetry',
                'convertedPhotonInvariantMass(0, 1)':'gammaveto_ConvertM',
                'daughterAngle(0,1)':'gammaveto_daughterAngle',
                'decayAngle(0)':'gammaveto_decayAngle_0',
                'decayAngle(1)':'gammaveto_decayAngle_1',
                'cos(decayAngle(0))':'gammaveto_cos_decayAngle_0',
                'cos(decayAngle(1))':'gammaveto_cos_decayAngle_1',
                'daughter(0,M)':'gammaveto_ep_M',
                'daughter(0,p)':'gammaveto_ep_p',
                'daughter(0,mcP)':'gammaveto_ep_mcP',
                'daughter(0,pt)':'gammaveto_ep_pt',
                'daughter(0,E)':'gammaveto_ep_E',
                'daughter(1,electronID)':'gammaveto_em_electronID',
                'daughter(1,pionID)':'gammaveto_em_pionID',
                'daughter(1,d0)':'gammaveto_em_d0',
                'daughter(1,z0)':'gammaveto_em_z0',
                'daughter(1,M)':'gammaveto_em_M',
                'daughter(1,p)':'gammaveto_em_p',
                'daughter(1,pt)':'gammaveto_em_pt',
                'daughter(1,E)':'gammaveto_em_E',
                'daughter(1,binaryPID(11,211))':'gammaveto_em_binaryPID',
                'pointangle':'gammaveto_pointangle',
                'chiProb':'gammaveto_chiProb',
                'openangle':'gammaveto_openangle',
                'psi':'gammaveto_psi',
                'daughterMotherDiffOf(0,theta)':'gammaveto_daughterMotherDiffOf_theta',
                'daughterMotherDiffOf(0,phi)':'gammaveto_daughterMotherDiffOf_phi'
                }
    ma.variableToSignalSideExtraInfo(particleList='gamma:veto', 
                                    varToExtraInfo=gamma_dict, 
                                    path=roe_path)
    # # Gamma:ROE Dictionary
    # gammaROE_dict = {'M': 'gammaROE_M',
    #                  'InvM':'gammaROE_InvM',
    #                  'E': 'gammaROE_E',
    #                  'decayAngle(0)':'gammaROE_decayAngle_0',
    #                  'cos(decayAngle(0))':'gammaROE_cos_decayAngle_0',
    #                  'dr':'gammaROE_dr',
    #                  'distance':'gammaROE_distance',
    #                  'mcDecayVertexFromIPDistance':'gammaROE_mcDecayVertexFromIPDistance',
    #                  'isSignal':'gammaROE_isSignal',
    #                  'isSignalAcceptWrongFSPs':'gammaROE_isSignalAcceptWrongFSPs',
    #                  'mcPDG':'gammaROE_mcPDG',
    #                  'Goodgamma':'gammaROE_Goodgamma',
    #                  'genMotherPDG':'gammaROE_genMotherPDG',
    #                  'genMotherPDG(1)':'gammaROE_genMotherPDG_1',
    #                  'genMotherPDG(2)':'gammaROE_genMotherPDG_2',
    #                  'nMCDaughters':'gammaROE_nMCDaughters',
    #                  'mcErrors':'gammaROE_mcErrors',
    #                  'daughter(0,p)':'gammaROE_ep_p',
    #                  'daughter(0,pt)':'gammaROE_ep_pt',
    #                  'daughter(0,mcPDG)':'gammaROE_ep_mcPDG',
    #                  'daughter(0,genMotherPDG)':'gammaROE_ep_genMotherPDG',
    #                  'daughter(0,mcSecPhysProc)':'gammaROE_ep_mcSecPhysProc',
    #                  'pointangle':'gammaROE_pointangle',
    #                  'chiProb':'gammaROE_chiProb',
    #                  'genNMissingDaughter(11)':'gammaROE_genNMissingDaughter',
    #                  'daughterMotherDiffOf(0,theta)':'gammaROE_daughterMotherDiffOf_theta',
    #                  'daughterMotherDiffOf(0,phi)':'gammaROE_daughterMotherDiffOf_phi',
    #                 }
    # ma.variableToSignalSideExtraInfo(particleList='gamma:singleT', 
    #                                  varToExtraInfo=gammaROE_dict, 
    #                                  path=roe_path)
    # Compton Scattering
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
    #                 # 'Goodgamma':'gammaCS_Goodgamma',
    #                 'genMotherPDG':'gammaCS_genMotherPDG',
    #                 'genMotherPDG(1)':'gammaCS_genMotherPDG_1',
    #                 'genMotherPDG(2)':'gammaCS_genMotherPDG_2',
    #                 'nMCDaughters':'gammaCS_nMCDaughters',
    #                 'mcErrors':'gammaCS_mcErrors',
    #                 'daughter(0,electronID)':'gammaCS_ep_electronID',
    #                 'daughter(0,p)':'gammaCS_ep_p',
    #                 'daughter(0,pt)':'gammaCS_ep_pt',
    #                 'daughter(0,mcPDG)':'gammaCS_ep_mcPDG',
    #                 'daughter(0,genMotherPDG)':'gammaCS_ep_genMotherPDG',
    #                 'daughter(0,mcSecPhysProc)':'gammaCS_ep_mcSecPhysProc',
    #                 'daughter(1,electronID)':'gammaCS_pi_electronID',
    #                 'daughter(1,p)':'gammaCS_pi_p',
    #                 'daughter(1,pt)':'gammaCS_pi_pt',
    #                 'daughter(1,mcPDG)':'gammaCS_pi_mcPDG',
    #                 'daughter(1,genMotherPDG)':'gammaCS_pi_genMotherPDG',
    #                 'daughter(1,mcSecPhysProc)':'gammaCS_pi_mcSecPhysProc',
    #                 'daughter(1,GoodElectron)':'gammaCS_pi_GoodElectron',
    #                 'pointangle':'gammaCS_pointangle',
    #                 'chiProb':'gammaCS_chiProb',
    #                 }
    # ma.variableToSignalSideExtraInfo(particleList='gamma:CS', 
    #                                 varToExtraInfo=gammaCS_dict, 
    #                                 path=roe_path)
    # # pi0 Dictionary
    # pi0_dict = {'M': 'pi0veto_M',
    #             'dM': 'pi0veto_dM',
    #             'E': 'pi0veto_E',
    #             'distance':'pi0veto_distance',
    #             'chiProb':'pi0veto_chiProb',
    #             'mcPDG':'pi0veto_mcPDG',
    #             'mcErrors':'pi0veto_mcErrors',
    #             'genMotherPDG':'pi0veto_genMotherPDG',
    #             'nMCDaughters':'pi0veto_nMCDaughters',
    #             'nDaughterPhotons':'pi0veto_nDaughterPhotons',
    #             'daughter(2,beamBackgroundSuppression)':'pi0veto_beamBackgroundSuppression',
    #             'daughter(2,fakePhotonSuppression)':'pi0veto_fakePhotonSuppression',
    #             }
    # ma.variableToSignalSideExtraInfo(particleList='pi0:veto', 
    #                                  varToExtraInfo=pi0_dict, 
    #                                  path=roe_path)
    # Dstar0 Dictionary
    Dstar0Mode1_dict = {'M': 'Dstar0Mode1_M',
                'dM': 'Dstar0Mode1_dM',
                'massDifference(0)':'Dstar0Mode1_massDifference',
                'E': 'Dstar0Mode1_E',
                'chiProb':'Dstar0Mode1_chiProb',
                'decayAngle(0)':'Dstar0Mode1_decayAngle_0',
                'decayAngle(1)':'Dstar0Mode1_decayAngle_1',
                'cos(decayAngle(0))':'Dstar0Mode1_cos_decayAngle_0',
                'cos(decayAngle(1))':'Dstar0Mode1_cos_decayAngle_1',
                'daughter(1,M)':'Dstar0Mode1_g_M',
                'daughter(1,pt)':'Dstar0Mode1_g_pt',
                'daughter(1,p)':'Dstar0Mode1_g_p',
                'daughter(1,E)':'Dstar0Mode1_g_E',
                'daughter(1,beamBackgroundSuppression)':'Dstar0Mode1_beamBackgroundSuppression',
                'daughter(1,fakePhotonSuppression)':'Dstar0Mode1_fakePhotonSuppression',
                }
    ma.variableToSignalSideExtraInfo(particleList='D*0:Mode1', 
                                    varToExtraInfo=Dstar0Mode1_dict, 
                                    path=roe_path)
    Dstar0Mode2_dict = {'M': 'Dstar0Mode2_M',
                        'dM': 'Dstar0Mode2_dM',
                        'massDifference(0)':'Dstar0Mode2_massDifference',
                        'E': 'Dstar0Mode2_E',
                        'chiProb':'Dstar0Mode2_chiProb',
                        'decayAngle(0)':'Dstar0Mode2_decayAngle_0',
                        'decayAngle(1)':'Dstar0Mode2_decayAngle_1',
                        'cos(decayAngle(0))':'Dstar0Mode2_cos_decayAngle_0',
                        'cos(decayAngle(1))':'Dstar0Mode2_cos_decayAngle_1',
                        'daughter(1,M)':'Dstar0Mode2_pi0_M',
                        'daughter(1,pt)':'Dstar0Mode2_pi0_pt',
                        'daughter(1,p)':'Dstar0Mode2_pi0_p',
                        'daughter(1,E)':'Dstar0Mode2_pi0_E',
                        }
    ma.variableToSignalSideExtraInfo(particleList='D*0:Mode2', 
                                    varToExtraInfo=Dstar0Mode2_dict, 
                                    path=roe_path)
    # Dstar+
    Dstarplus_dict = {'M': 'Dstarplus_M',
                    'dM': 'Dstarplus_dM',
                    'massDifference(0)':'Dstarplus_massDifference',
                    'E': 'Dstarplus_E',
                    'chiProb':'Dstarplus_chiProb',
                    'decayAngle(0)':'Dstarplus_decayAngle_0',
                    'decayAngle(1)':'Dstarplus_decayAngle_1',
                    'cos(decayAngle(0))':'Dstarplus_cos_decayAngle_0',
                    'cos(decayAngle(1))':'Dstarplus_cos_decayAngle_1',
                    'daughter(1,M)':'Dstarplus_pi_M',
                    'daughter(1,pt)':'Dstarplus_pi_pt',
                    'daughter(1,p)':'Dstarplus_pi_p',
                    'daughter(1,E)':'Dstarplus_pi_E',
                    'daughter(1,electronID)':'Dstarplus_pi_electronID',
                    'daughter(1,pionID)':'Dstarplus_pi_pionID',
                    }
    ma.variableToSignalSideExtraInfo(particleList='D*+:ROE', 
                                    varToExtraInfo=Dstarplus_dict, 
                                    path=roe_path)
    # # vpho
    # vpho_dict = {'vpho_M': 'vpho_M',
    #             'vpho_daughterAngle':'vpho_daughterAngle',
    #             'vpho_CMS_daughterAngle':'vpho_CMS_daughterAngle',
    #             'vpho_cos_daughterAngle':'vpho_cos_daughterAngle',
    #             'vpho_cos_CMS_daughterAngle':'vpho_cos_CMS_daughterAngle',
    #             }
    # ma.variableToSignalSideExtraInfo(particleList='vpho', 
    #                                 varToExtraInfo=vpho_dict, 
    #                                 path=roe_path)
    #--------------------------------------------------------------------------------------------

    # execute roe_path for each RestOfEvent in the event
    path.for_each('RestOfEvent', 'RestOfEvents', roe_path)

    # Variable List
    #--------------------------------------------------------------------------------------------
    gamma_ROE = []
    var_V0 = []
    for key in gamma_dict:
        vm.addAlias(gamma_dict[key], "extraInfo("+ gamma_dict[key] +")")
        gamma_ROE.append(gamma_dict[key])
        var_V0.append(key)

    # gammaROE_ROE = []
    # for key in gammaROE_dict:
    #     vm.addAlias(gammaROE_dict[key], "extraInfo("+ gammaROE_dict[key] +")")
    #     gammaROE_ROE.append(gammaROE_dict[key])

    # gammaCS_ROE = []
    # for key in gammaCS_dict:
    #     vm.addAlias(gammaCS_dict[key], "extraInfo("+ gammaCS_dict[key] +")")
    #     gammaCS_ROE.append(gammaCS_dict[key])

    # pi0_ROE = []
    # for key in pi0_dict:
    #     vm.addAlias(pi0_dict[key], "extraInfo("+ pi0_dict[key] +")")
    #     pi0_ROE.append(pi0_dict[key])

    Dstar0Mode1_ROE = []
    Dstar0Mode2_ROE = []
    for key in Dstar0Mode1_dict:
        vm.addAlias(Dstar0Mode1_dict[key], "extraInfo("+ Dstar0Mode1_dict[key] +")")
        Dstar0Mode1_ROE.append(Dstar0Mode1_dict[key])
    for key in Dstar0Mode2_dict:
        vm.addAlias(Dstar0Mode2_dict[key], "extraInfo("+ Dstar0Mode2_dict[key] +")")
        Dstar0Mode2_ROE.append(Dstar0Mode2_dict[key])

    Dstarplus_ROE = []
    for key in Dstarplus_dict:
        vm.addAlias(Dstarplus_dict[key], "extraInfo("+ Dstarplus_dict[key] +")")
        Dstarplus_ROE.append(Dstarplus_dict[key])

    # vpho_ROE = []
    # vpho_Base_ROE = []
    # for key in vpho_dict:
    #     vm.addAlias(vpho_dict[key], "extraInfo("+ vpho_dict[key] +")")
    #     vpho_ROE.append(vpho_dict[key])
    #     vpho_Base_ROE.append(vpho_dict[key])

    vm.addAlias('gammaveto_M_Correction', 'ifNANgiveX(gammaveto_M,10)')
    gamma_ROE.append('gammaveto_M_Correction')
    # vm.addAlias('M_Correction', 'ifNANgiveX(M,10)')
    # vm.addAlias('pidPairChargedBDTScore_Correction', 'ifNANgiveX(daughter(1,pidPairChargedBDTScore(11, 211, All)),10)')
    # var_V0.append('M_Correction')
    # var_V0.append('pidPairChargedBDTScore_Correction')

    # vm.addAlias('gammaROE_M_Correction', 'ifNANgiveX(gammaROE_M,10)')
    # gammaROE_ROE.append('gammaROE_M_Correction')

    # vm.addAlias('gammaCS_M_Correction', 'ifNANgiveX(gammaCS_M,10)')
    # gammaCS_ROE.append('gammaCS_M_Correction')

    # vm.addAlias('pi0veto_M_Correction', 'ifNANgiveX(pi0veto_M,10)')
    # pi0_ROE.append('pi0veto_M_Correction')

    vm.addAlias('Dstar0Mode1_M_Correction', 'ifNANgiveX(Dstar0Mode1_M,10)')
    Dstar0Mode1_ROE.append('Dstar0Mode1_M_Correction')

    vm.addAlias('Dstar0Mode2_M_Correction', 'ifNANgiveX(Dstar0Mode2_M,10)')
    Dstar0Mode2_ROE.append('Dstar0Mode2_M_Correction')

    vm.addAlias('Dstarplus_massDifference_Correction', 'ifNANgiveX(Dstarplus_massDifference,10)')
    Dstarplus_ROE.append('Dstarplus_massDifference_Correction')
    vm.addAlias('Dstarplus_M_Correction', 'ifNANgiveX(Dstarplus_M,10)')
    Dstarplus_ROE.append('Dstarplus_M_Correction')
    #--------------------------------------------------------------------------------------------
    # D_s+ Veto Ends Here
    #/================================================================================================================================/
    # # D_s*+ VETO Starts Here
    # #------------------------
    # # build RestOfEvent (ROE) object for each D_s*+ candidate
    # # ROE is required by the veto
    # ma.buildRestOfEvent(target_list_name='D_s*+',
    #                     path=path)

    # # Create a new path (called ROE path) which will be executed for
    # # each ROE in an event.
    # # Note that ROE exists for each B0 candidate, so when we loop
    # # over each ROE, we effectively loop over signal D_s*+ candidates

    # roe_path_Ds_star = b2.create_path()

    # # The ROE objects might in general be related to Particle from multiple
    # # particle lists therefore we need to check if the current ROE object
    # # is related to the Particle from our signal decay. If it is not
    # # the execution of roe_path will be finished (by starting empty,
    # # dead end path). Note that in this example this x-check is not
    # # necessary, but is anyway added for sake of completeness
    # deadEndPath_Ds_star = b2.create_path()

    # # Note again: all actions (modules) included in roe_path will be
    # # executed for each ROE in the event
    # # First we check that the current ROE is related to D_s*+ candidate
    # ma.signalSideParticleFilter(particleList='D_s*+',
    #                             selection='',
    #                             roe_path=roe_path_Ds_star,
    #                             deadEndPath=deadEndPath_Ds_star)

    # # create and fill gamma ParticleList that will contain
    # # all photons found in ROE (not used to reconstruct current B0 candidate)
    # # The photons need to have energy above 50 MeV to be considered
    # # (one can add any cut)
    # ma.fillParticleList(decayString='gamma:ROE',
    #                     cut='isInRestOfEvent == 1 and E > 0.050',
    #                     path=roe_path_Ds_star)

    # # in order to be able to use modularAnalysis functions (reconstructDecay in particular)
    # # we need a ParticleList containing the photon candidate used to reconstruct the
    # # current D_s*+ meson as well
    # # The DecayString is used to specify the selected particle (^)
    # ma.fillSignalSideParticleList(outputListName='gamma:sig',
    #                               decayString='D_s*+ -> D_s+ ^gamma:recon',
    #                               path=roe_path_Ds_star)

    # # make combinations of signal photon candidates with all photons from ROE
    # # keep only combinations in given invariant mass range
    # ma.reconstructDecay(decayString='pi0:vetoDsstar -> gamma:sig gamma:ROE',
    #                     cut='0.080 < M < 0.200',
    #                     path=roe_path_Ds_star)

    # # at this point one could use all features provided by the analysis software
    # # to make the veto as effective as possible. For example, one can perform truth
    # # matching, training/applying TMVA classifier, save pi0 candidates with ntuple
    # # maker for offline analysis/study.
    # ma.matchMCTruth(list_name='pi0:vetoDsstar', path=roe_path_Ds_star)

    # # in this example the variable, which is used to veto pi0 is very simple:
    # # invariant mass of pi0 that is closest to the pi0's nominal mass
    # # Therefore, we just simply rank pi0 candidates according to their distance
    # # from nominal mass (dM variable) and keep only the best candidate
    # ma.rankByLowest(particleList='pi0:vetoDsstar',
    #                 variable='abs(dM)',
    #                 numBest=1,
    #                 path=roe_path_Ds_star)

    # # write the invariant mass of the best pi0 candidate to the D_s*+ current
    # # candidate as the 'pi0veto' extraInfo
    # Dsstar_dict = {'M': 'pi0vetoDsstar_M',
    #                'decayAngle(0)':'pi0vetoDsstar_decayAngle_0',
    #                'decayAngle(1)':'pi0vetoDsstar_decayAngle_1',
    #                'cos(decayAngle(0))':'pi0vetoDsstar_cos_decayAngle_0',
    #                'cos(decayAngle(1))':'pi0vetoDsstar_cos_decayAngle_1',
    #                'isSignal':'pi0vetoDsstar_isSignal',
    #                'mcErrors':'pi0vetoDsstar_mcErrors',
    #                'mcPDG':'pi0vetoDsstar_mcPDG',
    #                'genMotherPDG':'pi0vetoDsstar_genMotherPDG',
    #                'daughter(1,p)':'pi0vetoDsstar_g_p',
    #                'daughter(1,mcP)':'pi0vetoDsstar_g_mcP',
    #                'daughter(1,pt)':'pi0vetoDsstar_g_pt',
    #                'daughter(1,E)':'pi0vetoDsstar_g_E',
    #                }
    # ma.variableToSignalSideExtraInfo(particleList='pi0:vetoDsstar', 
    #                                  varToExtraInfo=Dsstar_dict, 
    #                                  path=roe_path_Ds_star)

    # # execute roe_path for each RestOfEvent in the event
    # path.for_each('RestOfEvent', 'RestOfEvents', roe_path_Ds_star)

    # Dsstar_ROE = []
    # for key in Dsstar_dict:
    #     vm.addAlias(Dsstar_dict[key], "extraInfo("+ Dsstar_dict[key] +")")
    #     Dsstar_ROE.append(Dsstar_dict[key])

    # vm.addAlias('pi0vetoDsstar_M_Correction', 'ifNANgiveX(pi0vetoDsstar_M,10)')
    # Dsstar_ROE.append('pi0vetoDsstar_M_Correction')
    # #--------------------------------------------------------------------------------------------
    # # D_s*+ Veto Ends Here
    #/================================================================================================================================/
    return gamma_ROE, Dstar0Mode1_ROE, Dstar0Mode2_ROE, Dstarplus_ROE

def Suggestion(path):
    #/================================================================================================================================/
    # D_s*+ Info Saved to D_s+
    ma.variablesToDaughterExtraInfo(particleList='D_s*+', 
                                    decayString='D_s*+ -> [^D_s+ -> [D0:kmpip -> K-:loose pi+:loose] e-:corrected ?nu ?addbrems] gamma:recon', 
                                    variables={'massDifference(0)':'Ds_starminusDs',
                                            # 'DeltaM_diff_rank':'Ds_star_rank',
                                            #    'pi0vetoDsstar_M':'pi0veto_Ds_star',
                                            #    'pi0vetoDsstar_M_Correction':'pi0veto_Correction_Ds_star',
                                            }, 
                                    option=0, 
                                    path=path)
    vm.addAlias("Ds_starminusDs", "extraInfo(Ds_starminusDs)")
    vm.addAlias('Ds_starminusDs_M_Correction', 'ifNANgiveX(Ds_starminusDs,10)')

    # vm.addAlias("Ds_star_rank", "extraInfo(Ds_star_rank)")
    # vm.addAlias('Ds_star_rank_Correction', 'ifNANgiveX(Ds_star_rank,10)')

    # vm.addAlias("pi0veto_Ds_star", "extraInfo(pi0veto_Ds_star)")
    # vm.addAlias('pi0veto_Ds_star_Correction', 'ifNANgiveX(pi0veto_Ds_star,10)')

    # vm.addAlias("pi0veto_Correction_Ds_star", "extraInfo(pi0veto_Correction_Ds_star)")
    # vm.addAlias('pi0veto_Correction_Ds_star_Correction', 'ifNANgiveX(pi0veto_Correction_Ds_star,10)')

    vm.addAlias("goodDsplus", "passesCut(Ds_starminusDs_M_Correction >= 0.12 and Ds_starminusDs_M_Correction <= 0.165)")
    #/================================================================================================================================/
    # Kaons:
    #-------------
    # ma.fillParticleListFromROE('K+:ROE', '', maskName='cleanMask',
    #   sourceParticleListName='D_s+', useMissing = True, path=path)
    #/================================================================================================================================/
    # Other Tools:
    #-------------------
    # ma.calculateDistance('e+:corrected', 'D_s+ -> [^D0:kmpip -> K- pi+ ] e+:corrected ?nu', "trackvertex", path=path)
    # vm.addAlias("CalculatedDistance", "extraInfo(CalculatedDistance)")  # [E30]

    # ma.estimateAndAttachTrackFitResult('e+:corrected', path=path)

    # # BDT
    # #------------
    # # Fake D^0
    # path.add_module('MVAExpert', 
    #                 listNames=['D_s+'], 
    #                 extraInfoName='FakeD0BDT', 
    #                 identifier='user_adilmub_FakeD0Suppression')
    # # Bkg Suppression
    # path.add_module('MVAExpert', 
    #                 listNames=['D_s+'], 
    #                 extraInfoName='BkgBDT', 
    #                 identifier='user_adilmub_BkgSuppression'
    #                 # identifier='/home/belle2/amubarak/C02-MVA/Completed/MVAFastBDT.xml'
    #                 )
    # Variables from MVAExpert.
    expertVars = []
    # expertVars += ['extraInfo(FakeD0BDT_wdM)','extraInfo(FakeD0BDT_wodM)','extraInfo(FakeD0BDT_sigside)']
    # expertVars += ['extraInfo(FakeD0BDT)','extraInfo(BkgBDT)'] #'transformedNetworkOutput(FastBDT,0.1,1.0)'
    #/================================================================================================================================/
    # Overall Cut:
    #------------------
    # ma.applyCuts('D_s+', "e_genMotherPDG==22", path=path)
    # ma.applyCuts('D_s+', "Ds_chiProb_rank==1 and Ds_D0_sideband==1", path=path)
    # ma.applyCuts('D_s+', "e_mcSecPhysProc>=10 and e_mcSecPhysProc<=20 and gammaveto_M>=0.1", path=path)
    # ma.applyCuts('D_s*+', "gamma_beamBackgroundSuppressionScore>=0.5 and gamma_fakePhotonSuppressionScore>=0.5", path=path)
    #/================================================================================================================================/
    return expertVars

def Reconstruction_Variable(Curler, Before_VF, BCS, gamma_ROE, Dstar0Mode1_ROE, Dstar0Mode2_ROE, Dstarplus_ROE, expertVars):
    #/================================================================================================================================/
    # Save Variables
    #-------------------

    # D0 Variables
    #-------------------
    truth_D0 = ["mcSecPhysProc",'mcPDG','genMotherPDG','genMotherID','mcErrors']

    D0_kmpip_vars = []
    D0_kmpip_vars += vu.create_aliases_for_selected(
        list_of_variables=["M","dM","useCMSFrame(E)","useCMSFrame(p)",
                           "chiProb","cos(theta)","daughterAngle(0, 1)",
                           'daughterDiffOf(0,1,theta)','daughterDiffOf(0,1,cos(theta))','daughterDiffOf(0,1,phi)',
                           'daughterMotherDiffOf(0,x)','daughterMotherDiffOf(0,y)','daughterMotherDiffOf(0,z)',
                           "decayAngle(0)","cos(decayAngle(0))",
                           "decayAngle(1)","cos(decayAngle(1))",
                           "significanceOfDistance","flightDistance",
                           "isSignal",'ifNANgiveX(isSignal,0)','D0Mode','Dbar0Mode'] + truth_D0,
        decay_string='^D0:kmpip -> K-:loose pi+:loose',
        prefix=['D0_kmpip'])
    
    D0_kmpip_vars += vu.create_aliases_for_selected(
        list_of_variables=["M",'kaonID','nSVDHits','nCDCHits','abs(dr)','abs(dz)','z0','d0','pValue','chi2','ndf',
                           'theta','phi','omega'] + truth_D0,
        decay_string='D0:kmpip -> ^K-:loose pi+:loose',
        prefix=['K_kmpip'])
    D0_kmpip_vars += vu.create_aliases_for_selected(
        list_of_variables=["M",'pionID','nSVDHits','nCDCHits','abs(dr)','abs(dz)','z0','d0','pValue','chi2','ndf',
                           'theta','phi','omega'] + truth_D0,
        decay_string='D0:kmpip -> K-:loose ^pi+:loose',
        prefix=['pi_kmpip'])

    # D0_kmpippi0_vars = []
    # D0_kmpippi0_vars += vc.dalitz_3body
    # D0_kmpippi0_vars += vu.create_aliases_for_selected(
    #     list_of_variables=["M","dM","useCMSFrame(E)","useCMSFrame(p)",
    #                        "cosAngleBetweenMomentumAndVertexVector",
    #                        "chiProb",
    #                        "cos(theta)",
    #                        "decayAngle(0)","cos(decayAngle(0))",
    #                        "decayAngle(1)","cos(decayAngle(1))",
    #                        "flightDistance",
    #                        "isSignal",'ifNANgiveX(isSignal,0)','D0Mode','Dbar0Mode'] + truth_D0,
    #     decay_string='^D0:kmpippi0 -> K-:loose pi+:loose pi0:eff50_May2020',
    #     prefix=['D0_kmpippi0'])
    
    # D0_kmpippi0_vars += vu.create_aliases_for_selected(
    #     list_of_variables=["M",'kaonID','nSVDHits','nCDCHits','abs(dr)','abs(dz)','z0','d0','pValue','chi2','ndf',
    #                        'theta','phi','omega'] + truth_D0,
    #     decay_string='D0:kmpippi0 -> ^K-:loose pi+ pi0:eff50_May2020',
    #     prefix=['K_kmpippi0'])
    # D0_kmpippi0_vars += vu.create_aliases_for_selected(
    #     list_of_variables=["M",'pionID','nSVDHits','nCDCHits','abs(dr)','abs(dz)','z0','d0','pValue','chi2','ndf',
    #                        'theta','phi','omega'] + truth_D0,
    #     decay_string='D0:kmpippi0 -> K- ^pi+ pi0:eff50_May2020',
    #     prefix=['pi_kmpippi0'])
    
    # vm.addAlias('E1', 'daughter(0, E)')
    # vm.addAlias('E2', 'daughter(1, E)')
    # vm.addAlias('energyAsymmetry', 'formula(abs(E1 - E2) / (E1 + E2))')
    # D0_kmpippi0_vars += vu.create_aliases_for_selected(
    #     list_of_variables=["M","InvM","E",
    #                        "E1","E2","max(E1, E2)","energyAsymmetry"] + truth_D0,
    #     decay_string='D0:kmpippi0 -> K- pi+ ^pi0:eff50_May2020',
    #     prefix=['pi0_kmpippi0'])

    # D_s+ Variables
    #-------------------
    tracks = ['isCloneTrack',
            'dr','dz','abs(dr)','abs(dz)','z0','d0','pValue',
            'firstCDCLayer','firstPXDLayer','firstSVDLayer',
            'nPXDHits','nVXDHits','nSVDHits','nCDCHits',
            'seenInCDC','seenInPXD','seenInSVD','seenInTOP',
            'inARICHAcceptance','inCDCAcceptance','inTOPAcceptance']

    Ds_vars=[]
    var_1 = ['ImpactXY','cos(theta)','phi','mcP','M','pt','E','p','px','py','pz','abs(pz)','isOrHasCloneTrack',"charge"]

    Ds_vars += vu.create_aliases_for_selected(
        list_of_variables=Curler + [
                           'clusterE1E9','clusterE9E21',
                           "trackFitCovariance(0, 2)","trackFitCovariance(2, 3)","trackFitCovariance(0, 3)",
                           "trackFitHypothesisPDG",
                           'chi2','ndf','trackTime',
                           'pionID','electronID','binaryPID(11,211)',
                           'omega',
                           'daughter(0,isCloneTrack)',
                           "flightTime","formula(E/p)"]
                           +  var_1 + tracks,
        decay_string='D_s+ -> [D0:kmpip -> K-:loose pi+:loose] ^e-:corrected ?nu',
        prefix=['e'])

    Ds_vars += vu.create_aliases_for_selected(
        list_of_variables=['kaonID','nSVDHits','nCDCHits','abs(dr)','abs(dz)','z0','d0','pValue','chi2','ndf',
                           'theta','phi','omega'] + var_1 + tracks,
        decay_string='D_s+ -> [D0:kmpip -> ^K-:loose pi+:loose] e-:corrected ?nu',
        prefix=['K'])

    Ds_vars += vu.create_aliases_for_selected(
        list_of_variables=['pionID','nSVDHits','nCDCHits','abs(dr)','abs(dz)','z0','d0','pValue','chi2','ndf',
                           'theta','phi','omega'] + var_1 + tracks,
        decay_string='D_s+ -> [D0:kmpip -> K-:loose ^pi+:loose] e-:corrected ?nu',
        prefix=['pi'])   

    vm.addAlias('p1', 'daughter(0, p)')
    vm.addAlias('p2', 'daughter(1, p)')
    vm.addAlias('MomentumAsymmetry', 'formula((p1 - p2) / (p1 + p2))')
    vm.addAlias("charged_product", "formula(daughter(0,charge)*daughter(1,charge))")
    vm.addAlias("D0orD0bar", "passesCut(daughter(0,charge)==-1 and daughter(1,charge)==1)")
    Ds_vars += vu.create_aliases_for_selected(
        list_of_variables=["charged_product","MomentumAsymmetry",
                           "useCMSFrame(E)","useCMSFrame(p)",
                           "dM","useAlternativeDaughterHypothesis(M, 1:K+)",
                           "chiProb","cos(theta)","daughterAngle(0, 1)",
                           'daughterDiffOf(0,1,theta)','daughterDiffOf(0,1,cos(theta))','daughterDiffOf(0,1,phi)',
                           'daughterMotherDiffOf(0,theta)','daughterMotherDiffOf(0,cos(theta))','daughterMotherDiffOf(0,phi)',
                           "decayAngle(0)","cos(decayAngle(0))",
                           "decayAngle(1)","cos(decayAngle(1))",
                           "significanceOfDistance","flightDistance",
                           "useRestFrame(daughterAngle(0, 1))",
                           "formula(daughter(0,dr) - daughter(1,dr))","formula(daughter(0,dz) - daughter(1,dz))"]
                           + var_1,
        decay_string='D_s+ -> [^D0:kmpip -> K-:loose pi+:loose] e-:corrected ?nu',
        prefix=['D0'])

    vm.addAlias("Mode1Veto", "passesCut(1.98<=Ds_Dstar0Mode1_M_Correction<=2.04)")
    vm.addAlias("Mode2Veto", "passesCut(1.98<=Ds_Dstar0Mode2_M_Correction<=2.03)")
    vm.addAlias("DstarplusVeto", "passesCut(1.99<=Ds_Dstarplus_M_Correction<=2.025)")

    vm.addAlias("L_diff", "formula(((x - daughter(0,x))**2 + (y - daughter(0,y))**2 + (z - daughter(0,z))**2)**(1/2))")

    # vm.addAlias("nROE_Kplus", "nROE_ParticlesInList(K+:plus)")
    # vm.addAlias("nROE_Kminus", "nROE_ParticlesInList(K-:minus)")
    # vm.addAlias("nROE_KS0_1", "nROE_ParticlesInList(K_S0:merged)")
    # vm.addAlias("nROE_KS0_2", "nROE_ParticlesInList(K_S0:legacyGoodKS)")

    # vm.addAlias("KaonCount_1", "formula(nROE_Kminus - nROE_Kplus + nROE_KS0_1)")
    # vm.addAlias("KaonCount_2", "formula(nROE_Kminus + nROE_Kplus + nROE_KS0_1)")
    # vm.addAlias("KaonCount_3", "formula(nROE_Kminus - nROE_Kplus + nROE_KS0_2)")
    # vm.addAlias("KaonCount_4", "formula(nROE_Kminus + nROE_Kplus + nROE_KS0_2)")
    Ds_vars += vu.create_aliases_for_selected(
        list_of_variables=[
                        "M_uncorrected",'M_pi','massDifference(0)',"diff_D0pi","diff_D0K"]
                        + Before_VF + [
                        "formula(diff_D0pi - massDifference(0))",
                        'InvMLambda',
                        "Mode1Veto","Mode2Veto","DstarplusVeto",
                        # "chiProb_noIP","flightDistance_noIP",
                        'chiProb',
                        'useCMSFrame(p)','useCMSFrame(E)',
                        "phi_diff",
                        "mcM","mcM_D0e_emass",'mcM_D0e_pimass',"MminusMtrue_D0e_emass","MminusMtrue_D0e_pimass",
                        "decayAngle(0)","cos(decayAngle(0))",
                        "decayAngle(1)","cos(decayAngle(1))",
                        'pointingAngle(0)',
                        'daughterDiffOfPhi(0, 1)',
                        'pointangle',
                        'daughterAngle(0, 1)',
                        'cos(daughterAngle(0, 1))',
                        "Angle_D0e","Angle_Ke","Angle_pie",
                        "useCMSFrame(daughterAngle(0, 1))",
                        'useRestFrame(daughterAngle(0, 1))',
                        'useDaughterRestFrame(daughterAngle(0, 1),0:1)',
                        'useDaughterRestFrame(daughterAngle(0, 1),0:1,1)',"psi",
                        'azimuthalAngleInDecayPlane(0,1)',
                        'daughterDiffOf(0, 1, cos(theta))',
                        'cosAngleBetweenMomentumAndVertexVector',
                        'cosAngleBetweenMomentumAndVertexVectorInXYPlane',
                        'useRestFrame(daughterDiffOf(0, 1, p))',
                        'useRestFrame(daughterMotherDiffOf(0, p))',
                        'flightDistance','distance',
                        "flightDistanceOfDaughter(0)","flightDistanceOfDaughterErr(0)",
                        "flightTimeOfDaughter(0)","flightTimeOfDaughterErr(0)",
                        "vertexDistanceOfDaughter(0)","vertexDistanceOfDaughterErr(0)",
                        'daughterDiffOf(0,1,p)','daughterDiffOf(0,1,E)',
                        'daughterDiffOf(0,1,px)','daughterDiffOf(0,1,py)','daughterDiffOf(0,1,pz)',
                        'daughterDiffOf(0,1,x)','daughterDiffOf(0,1,y)','daughterDiffOf(0,1,z)',
                        'daughterMotherDiffOf(0,p)','daughterMotherDiffOf(0,E)',
                        'daughterMotherDiffOf(0,px)','daughterMotherDiffOf(0,py)','daughterMotherDiffOf(0,pz)',
                        'daughterMotherDiffOf(0,x)','daughterMotherDiffOf(0,y)','daughterMotherDiffOf(0,z)',
                        'abs(daughterMotherDiffOf(0,distance))',"L_diff",
                        'daughterMotherDiffOf(0,flightDistance)','daughterMotherDiffOf(0,vertexDistance)',
                        'flightDistanceOfDaughter(0)']
                        #  + ["nROE_Kplus","nROE_Kminus","nROE_KS0_1","nROE_KS0_2",
                            #   "nROE_Charged(roe_mask, 211)","nROE_ParticlesInList(pi+:loose, roe_mask)",
                            #   "nROE_Charged(roe_mask, 321)",
                            #   "nROE_ParticlesInList(K+:loose, roe_mask)",
                            #   "nROE_ParticlesInList(K+:plus, roe_mask)","nROE_ParticlesInList(K+:minus, roe_mask)",
                            # "KaonCount_1","KaonCount_2","KaonCount_3","KaonCount_4"]
                        + gamma_ROE + Dstar0Mode1_ROE + Dstar0Mode2_ROE + Dstarplus_ROE
                        + ["Ds_starminusDs","Ds_starminusDs_M_Correction",
                            # "Ds_star_rank","Ds_star_rank_Correction",
                            #   "pi0veto_Ds_star",'pi0veto_Ds_star_Correction',
                            #   "pi0veto_Correction_Ds_star",'pi0veto_Correction_Ds_star_Correction',
                            "goodDsplus"] 
                        + var_1 + truth 
                        + ['genNStepsToDaughter(0)','genNStepsToDaughter(1)',
                        'genNMissingDaughter(11)','genNMissingDaughter(22)',
                        'ifNANgiveX(mcPDG,9999)',
                        'nMCDaughters','genParticleID'] 
                        + ["Comb","Failed","Dstar0","Dstarplus","NoDstarplus","Signal"] 
                        + ["D0_Dstarplus","D0_Dstar0","D0_NoDstarplusDstar0","D0_Other"]
                        + ["mcDaughter(0, PDG)"]
                        + expertVars + BCS,
        decay_string='^D_s+ -> [D0:kmpip -> K-:loose pi+:loose] e-:corrected ?nu',
        prefix=['Ds'])

    # D_s*+ Variables
    #-------------------
    var_2 = ['cos(theta)','phi','mcP','M','InvM','pt','p','isOrHasCloneTrack',"charge"]

    Ds_star_vars=[]
    Ds_star_vars += vu.create_aliases_for_selected(
        list_of_variables= ['formula(clusterTotalMCMatchWeight/clusterE)',
                            'phi','genMotherID',"E",
                            "inCDCAcceptance",
                            "clusterSecondMoment",
                            "clusterReg","clusterE",
                            "clusterTiming","clusterPulseShapeDiscriminationMVA","clusterTheta","clusterZernikeMVA",
                            "minC2TDist","clusterZernikeMVA","clusterTheta",
                            "clusterErrorTiming","clusterE1E9","clusterE9E21",
                            'mcSecPhysProc','isMisidentified',
                            "beamBackgroundSuppressionScore","fakePhotonSuppressionScore"] + var_2,
        decay_string='D_s*+ -> [D_s+ -> [D0:kmpip -> K-:loose pi+:loose] e-:corrected ?nu ?addbrems] ^gamma:recon',
        prefix=['gamma'])

    Ds_star_vars += vu.create_aliases_for_selected(
        list_of_variables= expertVars + ['gammaveto_M','gammaveto_M_Correction'] + var_2,
        decay_string='D_s*+ -> [^D_s+ -> [D0:kmpip -> K-:loose pi+:loose] e-:corrected ?nu ?addbrems] gamma:recon',
        prefix=['Ds_c'])

    Ds_star_vars += vu.create_aliases_for_selected(
        list_of_variables= ['massDifference(0)',
                            'isPrimarySignal','mcSecPhysProc','ifNANgiveX(isSignal,5)',
                            'genNStepsToDaughter(0)','genNStepsToDaughter(1)',
                            'daughterAngle(0, 1)',
                            'M_pi',"mass_diff_pi",
                            # 'DeltaM_diff','dM','random',
                            ] + var_2,
        decay_string='^D_s*+ -> [D_s+ -> [D0:kmpip -> K-:loose pi+:loose] e-:corrected ?nu ?addbrems] gamma:recon',
        prefix=['Ds_star'])

    # # Photon Conversion Variables
    # #------------------------------
    # gammaV0 = []
    # gammaV0 += vu.create_aliases_for_selected(
    #     list_of_variables = ['ArmenterosLongitudinalMomentumAsymmetry',
    #                          'M','InvM',
    #                         #  'convertedPhotonInvariantMass(0, 1)','convertedPhotonDelR(0, 1)',
    #                          'convertedPhotonX(0, 1)','convertedPhotonY(0, 1)','convertedPhotonZ(0, 1)',
    #                         #  'v0DaughterD0PullWithOriginAsPivot(0)','v0DaughterD0PullWithOriginAsPivot(1)',
    #                         #  'v0DaughterFirstPXDLayer(0)','v0DaughterFirstPXDLayer(1)',
    #                         #  'v0DaughtersShare1stHit','v0DaughtersShare1stUHit','v0DaughtersShare1stVHit',
    #                          'isFromV0',
    #                          'distance','mcDecayVertexFromIPDistance','flightDistance',
    #                          'daughterAngle(0,1)',
    #                          'isSignal','mcPDG','genMotherPDG','mcErrors',
    #                          "chiProb"
    #                          ],
    #     decay_string='^gamma:V0_array -> e+ e-',
    #     prefix=['gammaV0'])
    # gammaV0 += vu.create_aliases_for_selected(
    #     list_of_variables = ['electronID',
    #                          'InvM','p','pt',
    #                          'd0','z0',
    #                          'mcPDG','genMotherPDG',"mcSecPhysProc"
    #                          ],
    #     decay_string='gamma:V0_array -> ^e+ e-',
    #     prefix=['eplusV0'])
    # gammaV0 += vu.create_aliases_for_selected(
    #     list_of_variables = ['electronID',
    #                          'InvM','p','pt',
    #                          'd0','z0',
    #                          'mcPDG','genMotherPDG',"mcSecPhysProc"
    #                          ],
    #     decay_string='gamma:V0_array -> e+ ^e-',
    #     prefix=['eminusV0'])

    # Extra Variables
    #--------------------------------------------------------------------------------------------
    # Some extra D_s+ Variables
    vm.addAlias("diff_D0e", "massDifference(0)")
    vm.addAlias("diff_D0pi", "useAlternativeDaughterHypothesis(massDifference(0), 1:pi+)")
    vm.addAlias("diff_D0K", "useAlternativeDaughterHypothesis(massDifference(0), 1:K+)")
    vm.addAlias("Angle_D0e", "daughterAngle(0, 1)")

    vm.addAlias("Angle_Ke", "daughterAngle(0:0, 1)")
    vm.addAlias("Angle_pie", "daughterAngle(0:1, 1)")

    vm.addAlias("M_uncorrected", "daughterCombination(M,0,1:0)")
    #--------------------------------------------------------------------------------------------
    # Curlers (Won't be needed too long)
    vm.addAlias("M_pi", "useAlternativeDaughterHypothesis(M, 1:pi+)")
    vm.addAlias("mcM", "formula(((mcE)**2 - (mcP)**2)**(1/2))")

    vm.addAlias("mcP_D0e_mag", "formula((daughter(0,mcP))**2 + (daughter(1,mcP))**2)")
    vm.addAlias("mcP_D0e_comp", "formula((daughter(0,mcPX)*daughter(1,mcPX)) + (daughter(0,mcPY)*daughter(1,mcPY)) + (daughter(0,mcPZ)*daughter(1,mcPZ)))")
    # Positron has electron mass
    vm.addAlias("mcE_D0e_emass", "formula(daughter(0,mcE) + daughter(1,mcE))")

    vm.addAlias("mcM_D0e_emass", "formula(((mcE_D0e_emass)**2 - mcP_D0e_mag - (2*mcP_D0e_comp))**(1/2))")

    vm.addAlias("MminusMtrue_D0e_emass", "formula(M - mcM_D0e_emass)")
    # Positron has pion mass
    vm.addAlias("mcE_e", "formula(((0.13957039)**2 + (daughter(1,mcP))**2)**(1/2))")
    vm.addAlias("mcE_D0e_pimass", "formula(daughter(0,mcE) + mcE_e)")

    vm.addAlias("mcM_D0e_pimass", "formula(((mcE_D0e_pimass)**2 - mcP_D0e_mag - (2*mcP_D0e_comp))**(1/2))")

    vm.addAlias("MminusMtrue_D0e_pimass", "formula(M_pi - mcM_D0e_pimass)")
    #--------------------------------------------------------------------------------------------
    # Mass Difference with positron given a pion mass
    vm.addAlias("Ds_E_D0pi", 'daughter(0,useAlternativeDaughterHypothesis(E, 1:pi+))')
    vm.addAlias("Ds_star_E", "formula(Ds_E_D0pi + daughter(1,E))")

    vm.addAlias("Ds_star_P_mag", "formula((daughter(0,p))**2 + (daughter(1,p))**2)")
    vm.addAlias("Ds_star_P_comp", "formula((daughter(0,px)*daughter(1,px)) + (daughter(0,py)*daughter(1,py)) + (daughter(0,pz)*daughter(1,pz)))")

    vm.addAlias("M_Dsph_pi", "formula(((Ds_star_E)**2 - Ds_star_P_mag - (2*Ds_star_P_comp))**(1/2))")
    vm.addAlias("Ds_M_pi_star", "formula(daughter(0,useAlternativeDaughterHypothesis(M, 1:pi+)))")

    vm.addAlias("mass_diff_pi", "formula(M_Dsph_pi - Ds_M_pi_star)")
    #--------------------------------------------------------------------------------------------
    Event = ['IPX','IPY','IPZ']
    #/==========================================================================================================================/
    return D0_kmpip_vars, Ds_vars, Ds_star_vars, Event

if __name__ == '__main__':

    #/==========================================================================================================================/
    # Create  Main Path
    #--------------------
    path = b2.create_path()
    #/================================================================================================================================/
    # --I/O----------------------------------------------------------------------------------------
    # load input ROOT file
    ma.inputMdst(environmentType='default',
                filename="/home/belle2/amubarak/C00-Generation/ccbarDs+EventGeneration.root",
                path=path)
    #/================================================================================================================================/
    # Reconstruction
    #-----------------------------
    Curler, Before_VF, BCS = Reconstruction(path, ChargeC, Pion_Cut, Kaon_Cut, Electron_Cut, Gamma_Cut, D0_Cut, Ds_Cut, Ds_star_Cut)

    # ROE
    #---------------------------------------------------------------------------------------------
    # build RestOfEvent (ROE) object for each D_s+ candidate
    # ROE is required by the veto
    ma.buildRestOfEvent(target_list_name='D_s+',
                        fillWithMostLikely=False,
                        # chargedPIDPriors=[1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        path=path)
    # Create a mask tuple:
    track_based_cuts = "abs(d0) < 20.0 and abs(z0) < 20.0"
    # track_based_cuts = ""
    ecl_based_cuts = ""  # [E30]
    roe_mask = ("roe_mask", track_based_cuts, ecl_based_cuts)
    ma.appendROEMasks("D_s+", [roe_mask], path=path)  # [E40]

    ma.updateROEUsingV0Lists(target_particle_list="D_s+", 
                            mask_names="roe_mask", 
                            default_cleanup=True, 
                            selection_cuts=None, 
                            apply_mass_fit=True, 
                            fitter='treefit', 
                            path=path)

    # Create a new path (called ROE path) which will be executed for
    # each ROE in an event.
    # Note that ROE exists for each D_s+ candidate, so when we loop
    # over each ROE, we effectively loop over signal D_s+ candidates

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
    # First we check that the current ROE is related to D_s+ candidate
    ma.signalSideParticleFilter(particleList='D_s+',
                                selection='',
                                roe_path=roe_path,
                                deadEndPath=deadEndPath)
    #---------------------------------------------------------------------------------------------
    gamma_ROE, Dstar0Mode1_ROE, Dstar0Mode2_ROE, Dstarplus_ROE = Rest_of_Event(roe_path, ChargeC, ElectronROE_Cut, PionROE_Cut, GammaROE_Cut)
    expertVars = Suggestion(path)

    D0_kmpip_vars, Ds_vars, Ds_star_vars, Event = Reconstruction_Variable(Curler, Before_VF, BCS, gamma_ROE, Dstar0Mode1_ROE, Dstar0Mode2_ROE, Dstarplus_ROE, expertVars)
    #/================================================================================================================================/
    # Save
    #--------
    # Saving variables to ntuple

    output_file = 'output_test.root'

    # D_s*+ -> [D_s+ -> [D0:kmpip -> K- pi+ ] e+:corrected ?nu] gamma:recon
    ma.variablesToNtuple('D_s*+', Ds_star_vars,
                filename=output_file, treename='Ds*tree', path=path)
    # D_s+ -> D0 e+ nu_e decay
    ma.variablesToNtuple('D_s+', Ds_vars + Event,
                        filename=output_file, treename='Dstree', path=path)
    # D0
    ma.variablesToNtuple('D0:kmpip', D0_kmpip_vars,
                        filename=output_file, treename='D02kmpiptree', path=path)

    # # Converted Photons
    # ma.variablesToNtuple('gamma:veto', ["M"],# + mc_gen_topo(n=200),
    #                     filename=output_file, treename='gammaV0tree', path=roe_path)

    # ma.variablesToNtuple('D*0:Mode1', ['M'],
    #                     filename=output_file, treename='Dstar01tree', path=roe_path)
    # ma.variablesToNtuple('D*0:Mode2', ['M'],
    #                     filename=output_file, treename='Dstar02tree', path=roe_path)
    # ma.variablesToNtuple('D*+:ROE', ['M'],
    #                     filename=output_file, treename='Dstarplustree', path=roe_path)

    # # Generator Level Particle
    # ma.variablesToNtuple('e+:gen', ["mcSecPhysProc","mcTheta","mcPhi","mcPX","mcPY","mcPZ","mcPDG","genMotherPDG"],
    #             filename=output_file, treename='eMCtree', path=path)
    #/==========================================================================================================================/
    # -----------------
    # Monitor progress.
    # -----------------

    progress = b2.register_module("Progress")
    path.add_module(progress)

    # ---------------
    # Process events.
    # ---------------

    # Start processing of modules.
    b2.process(path)

    # Print basf2 call statistics.
    print(b2.statistics)