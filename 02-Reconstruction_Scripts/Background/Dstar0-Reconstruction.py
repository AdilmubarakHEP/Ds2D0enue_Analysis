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
# D*0 -> D0 pi0/gamma                                                    #
#         |                                                              #
#         +-> K- pi+                                                     #
#                                                                        #
##########################################################################

import argparse

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
b2.conditions.prepend_globaltag('user_adilmub_BkgSuppression')
b2.conditions.prepend_globaltag('user_adilmub_FakeD0Suppression')

def argparser():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("--Mode",
                        type=int,
                        # nargs="+",
                        default=1,
                        choices=list(range(1, 4)),
                        help="Different D*0 Modes:\n"
                        "First Mode is where pi0 decays to two photons.\n"
                        "Second Mode is the dalitz decay.\n"
                        "Third Mode is where there is a single photon.")


    return parser

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
Electron_Cut = 'electronID >= 0.5 and abs(dr) < 1 and abs(dz) < 3'
Gamma_Cut = 'E >= 0.1'

D0_Cut = "-0.04 <= dM <= 0.04 and useCMSFrame(p) > 2.5"
# D0_Cut = "[-0.23 <= dM <= -0.03 or 0.03 <= dM <= 0.23] and useCMSFrame(p) > 2.5"
Ds_Cut = "useAlternativeDaughterHypothesis(massDifference(0), 1:pi+) <= 0.6 and -0.02 <= daughter(0,dM) <= 0.02"
# Ds_Cut = "useAlternativeDaughterHypothesis(massDifference(0), 1:pi+) <= 0.6 and [-0.2 <= daughter(0,dM) <= -0.04 or 0.04 <= daughter(0,dM) <= 0.2]"
# Ds_Cut = "massDifference(0) <= 0.25"
Ds_star_Cut = "massDifference(0) <= 0.4"

# Veto
#-------
ElectronROE_Cut = 'isInRestOfEvent == 1 and abs(dr) < 20.0 and abs(dz) < 20.0' # Cuts: abs(d0) < 20.0 and abs(z0) < 20.0
PionROE_Cut = 'isInRestOfEvent == 1 and abs(dr) < 20.0 and abs(dz) < 20.0' # Cuts: abs(d0) < 20.0 and abs(z0) < 20.0
GammaROE_Cut = 'isInRestOfEvent == 1 and E > 0.1'
#/================================================================================================================================/

def Mode_1(path):
    #/================================================================================================================================/
    # Reconstruct D0 -> K- pi+ decay
    #------------------------------------
    ma.reconstructDecay(decayString='D0:kmpip -> K-:loose pi+:loose', 
                        cut=D0_Cut, 
                        chargeConjugation=ChargeC,
                        path=path)
    #chargeConjugation=False: The option is true by default
    #
    # Vertex Fitting D0
    #-------------------
    vx.treeFit('D0:kmpip',
            conf_level=0,
            updateAllDaughters=False,
            path=path)
    ma.applyCuts("D0:kmpip", D0_Cut, path=path)

    # Photon Conversion
    #--------------------------
    ma.reconstructDecay(decayString='gamma:RD_partial -> e+:corrected ...',
                        cut='',
                        chargeConjugation=ChargeC,
                        allowChargeViolation=True,
                        path=path)
    ma.applyCuts('gamma:RD_partial', '', path=path)

    # Reconstruct pi0-> gamma gamma
    #--------------------------------------------------
    ma.reconstructDecay(decayString='pi0:2photons -> [gamma:RD_partial -> e+:corrected ...] ?gamma',
                        cut='',
                        chargeConjugation=ChargeC,
                        allowChargeViolation=True,
                        path=path)
    ma.applyCuts('pi0:2photons', '', path=path)

    # Reconstruct D*0
    #-----------------------------------------
    ma.reconstructDecay(decayString='D*0:Mode1 -> [D0:kmpip -> K-:loose pi+:loose ] [pi0:2photons -> [gamma:RD_partial -> e+:corrected ...] ?gamma]', 
                        cut=Ds_Cut, 
                        chargeConjugation=ChargeC, 
                        path=path) 
    #chargeConjugation=False: The option is true by default
    #
    # Vertex Fitting D*0
    #-------------------
    vx.treeFit('D*0:Mode1', 
            conf_level=0, 
            ipConstraint=True,
            updateAllDaughters=False, 
            path=path)
    ma.applyCuts("D*0:Mode1", Ds_Cut, path=path)
    # Save some Variables Before IP Constraint
    ma.variablesToExtraInfo('D*0:Mode1', variables={"chiProb": "chiProb_noIP", "flightDistance": "flightDistance_noIP"}, path=path)
    vm.addAlias('chiProb_noIP', 'extraInfo(chiProb_noIP)')
    vm.addAlias('flightDistance_noIP', 'extraInfo(flightDistance_noIP)')
    # Save Variable After IP Constraint
    vx.treeFit('D*0:Mode1', 
                conf_level=0, 
                ipConstraint=True,
                updateAllDaughters=False, 
                path=path)

    # D*0 Best Candidate
    ma.rankByHighest('D*0:Mode1', 
                    variable='chiProb',
                    outputVariable='chiProb_Ds_rank',
                    # numBest=1,
                    path=path)
    vm.addAlias('chiProb_Ds_rank', 'extraInfo(chiProb_Ds_rank)')
    ma.rankByHighest('D*0:Mode1', 
                    variable='daughter(1,p)',
                    outputVariable='Electron_p_Ds_rank',
                    # numBest=1,
                    path=path)
    vm.addAlias('Electron_p_Ds_rank', 'extraInfo(Electron_p_Ds_rank)')
    ma.rankByHighest('D*0:Mode1', 
                    variable='random',
                    outputVariable='random_Ds_rank',
                    # numBest=1,
                    path=path)
    vm.addAlias('random_Ds_rank', 'extraInfo(random_Ds_rank)')
    BCS = ['chiProb_Ds_rank','Electron_p_Ds_rank','random_Ds_rank']
    #/================================================================================================================================/
    # Reconstruct vpho -> D*0 gamma
    #-----------------------------------------
    ma.reconstructDecay(decayString='vpho -> [D*0:Mode1 -> [D0:kmpip -> K-:loose pi+:loose ] [pi0:2photons -> [gamma:RD_partial -> e+:corrected ...] ?gamma]] gamma:recon', 
                        cut=Ds_star_Cut, 
                        chargeConjugation=ChargeC, 
                        allowChargeViolation=True,
                        path=path)
    #chargeConjugation=False: The option is true by default
    #
    # Vertex Fitting vpho
    # ----------------------
    # vx.treeFit('vpho', 
    #            conf_level=-1,
    #            path=my_path)
    ma.applyCuts("vpho", Ds_star_Cut, path=path)

    # vpho Best Candidate
    vm.addAlias('DeltaM_diff', 'formula(abs(massDifference(0) - 0.1438))')
    ma.rankByLowest('vpho', 
                    variable='DeltaM_diff',
                    # numBest=1,
                    path=path)
    vm.addAlias('DeltaM_diff_rank', 'extraInfo(DeltaM_diff_rank)')

    ma.rankByLowest('vpho', 
                    variable='dM',
                    # numBest=1,
                    path=path)
    vm.addAlias('dM_rank', 'extraInfo(dM_rank)')

    ma.rankByHighest('vpho', 
                    variable='random',
                    # numBest=1,
                    path=path)
    vm.addAlias('random_rank', 'extraInfo(random_rank)')
    #/================================================================================================================================/
    # Perform MC Matching (MC truth asociation)
    #----------------------------------------------
    ma.matchMCTruth(list_name='D*0:Mode1', path=path) 
    ma.matchMCTruth(list_name='vpho', path=path) 
    #/================================================================================================================================/
    # D*0 VETO Starts Here
    #------------------------

    # ROE
    #---------------------------------------------------------------------------------------------
    # build RestOfEvent (ROE) object for each D*0 candidate
    # ROE is required by the veto
    ma.buildRestOfEvent(target_list_name='D*0:Mode1',
                        fillWithMostLikely=False,
                        # chargedPIDPriors=[1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        path=path)
    # Create a mask tuple:
    track_based_cuts = "abs(d0) < 20.0 and abs(z0) < 20.0"
    # track_based_cuts = ""
    ecl_based_cuts = ""  # [E30]
    roe_mask = ("roe_mask", track_based_cuts, ecl_based_cuts)
    ma.appendROEMasks("D*0:Mode1", [roe_mask], path=path)  # [E40]

    ma.updateROEUsingV0Lists(target_particle_list="D*0:Mode1", 
                            mask_names="roe_mask", 
                            default_cleanup=True, 
                            selection_cuts=None, 
                            apply_mass_fit=True, 
                            fitter='treefit', 
                            path=path)
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
    ma.signalSideParticleFilter(particleList='D*0:Mode1',
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
    #---------------------------------------------------------------------------------------------

    # in order to be able to use modularAnalysis functions (reconstructDecay in particular)
    # we need a ParticleList containing the photon candidate used to reconstruct the
    # current D*0 as well
    # The DecayString is used to specify the selected particle (^)
    ma.fillSignalSideParticleList(outputListName='e+:sig',
                                decayString='D*0:Mode1 -> [D0:kmpip -> K-:loose pi+:loose ] [pi0:2photons -> [gamma:RD_partial -> ^e+:corrected ...] ?gamma]',
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
            updateAllDaughters=False,
            path=roe_path)
    ma.applyCuts('gamma:veto', 'InvM <= 2 and M<=3', path=roe_path)
    #-------------------------------------------------------------------------------------------

    # Perform MC Matching (MC truth asociation)
    ma.matchMCTruth(list_name='gamma:veto', path=roe_path)

    # Best Candidate Selection
    #-------------------------------------------------------
    # perform best candidate selection
    vm.addAlias("gamma_abs_dM", "formula(abs(M - 0.0))")
    ma.rankByLowest('gamma:veto', 
                    variable='gamma_abs_dM',
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

    # vm.addAlias("Goodgamma", "passesCut(mcPDG==22)")
    vm.addAlias("GoodElectron", "passesCut(abs(mcPDG)==11)")

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
                'daughterAngle(0,1)':'gammaveto_daughterAngle',
                'decayAngle(0)':'gammaveto_decayAngle_0',
                'decayAngle(1)':'gammaveto_decayAngle_1',
                'cos(decayAngle(0))':'gammaveto_cos_decayAngle_0',
                'cos(decayAngle(1))':'gammaveto_cos_decayAngle_1',
                'isSignal':'gammaveto_isSignal',
                'isSignalAcceptWrongFSPs':'gammaveto_isSignalAcceptWrongFSPs',
                'mcPDG':'gammaveto_mcPDG',
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
                'daughter(0,genMotherPDG)':'gammaveto_ep_genMotherPDG',
                'daughter(1,electronID)':'gammaveto_em_electronID',
                'daughter(1,d0)':'gammaveto_em_d0',
                'daughter(1,z0)':'gammaveto_em_z0',
                'daughter(1,M)':'gammaveto_em_M',
                'daughter(1,p)':'gammaveto_em_p',
                'daughter(1,mcP)':'gammaveto_em_mcP',
                'daughter(1,pt)':'gammaveto_em_pt',
                'daughter(1,E)':'gammaveto_em_E',
                'daughter(1,mcE)':'gammaveto_em_mcE',
                'daughter(1,mcPDG)':'gammaveto_em_mcPDG',
                'daughter(1,genMotherPDG)':'gammaveto_em_genMotherPDG',
                'daughter(1,mcSecPhysProc)':'gammaveto_em_mcSecPhysProc',
                'daughter(1,binaryPID(11,211))':'gammaveto_em_binaryPID',
                'daughter(1,GoodElectron)':'gammaveto_em_GoodElectron',
                'pointangle':'gammaveto_pointangle',
                'chiProb':'gammaveto_chiProb',
                'openangle':'gammaveto_openangle',
                'psi':'gammaveto_psi',
                'genNMissingDaughter(11)':'gammaveto_genNMissingDaughter',
                'daughterMotherDiffOf(0,theta)':'gammaveto_daughterMotherDiffOf_theta',
                'daughterMotherDiffOf(0,phi)':'gammaveto_daughterMotherDiffOf_phi'
                }
    ma.variableToSignalSideExtraInfo(particleList='gamma:veto', 
                                    varToExtraInfo=gamma_dict, 
                                    path=roe_path)
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

    vm.addAlias('gammaveto_M_Correction', 'ifNANgiveX(gammaveto_M,10)')
    gamma_ROE.append('gammaveto_M_Correction')
    #--------------------------------------------------------------------------------------------
    # D*0 Veto Ends Here
    #/================================================================================================================================/
    # vpho Info Saved to D*0
    ma.variablesToDaughterExtraInfo(particleList='vpho', 
                                    decayString='vpho -> [^D*0:Mode1 -> [D0:kmpip -> K-:loose pi+:loose ] [pi0:2photons -> [gamma:RD_partial -> e+:corrected ...] ?gamma]] gamma:recon', 
                                    variables={'massDifference(0)':'Ds_starminusDs',
                                            'DeltaM_diff_rank':'Ds_star_rank',
                                            }, 
                                    option=0, 
                                    path=path)
    vm.addAlias("Ds_starminusDs", "extraInfo(Ds_starminusDs)")
    vm.addAlias('Ds_starminusDs_M_Correction', 'ifNANgiveX(Ds_starminusDs,10)')

    vm.addAlias("Ds_star_rank", "extraInfo(Ds_star_rank)")
    vm.addAlias('Ds_star_rank_Correction', 'ifNANgiveX(Ds_star_rank,10)')

    vm.addAlias("goodDsplus", "passesCut(Ds_starminusDs_M_Correction >= 0.12 and Ds_starminusDs_M_Correction <= 0.165)")
    #/================================================================================================================================/
    # Other Tools:
    #-------------------
    # BDT
    #------------
    # MVAExpert
    path.add_module('MVAExpert', 
                    listNames=['D*0:Mode1'], 
                    extraInfoName='FakeD0BDT', 
                    identifier='user_adilmub_FakeD0Suppression')
    path.add_module('MVAExpert', 
                    listNames=['D*0:Mode1'], 
                    extraInfoName='BkgBDT', 
                    identifier='BkgMVA')
    # Variables from MVAExpert.
    expertVars = ['extraInfo(FakeD0BDT)','extraInfo(BkgBDT)'] #'transformedNetworkOutput(FastBDT,0.1,1.0)'
    #/================================================================================================================================/
    return BCS, gamma_ROE, expertVars

def Mode_2(path):
    #/================================================================================================================================/
    # Reconstruct D0 -> K- pi+ decay
    #------------------------------------
    ma.reconstructDecay(decayString='D0:kmpip -> K-:loose pi+:loose', 
                        cut=D0_Cut, 
                        chargeConjugation=ChargeC,
                        path=path)
    #chargeConjugation=False: The option is true by default
    #
    # Vertex Fitting D0
    #-------------------
    vx.treeFit('D0:kmpip',
            conf_level=0,
            updateAllDaughters=False,
            path=path)
    ma.applyCuts("D0:kmpip", D0_Cut, path=path)

    # Dalitz Decay
    #--------------------------
    ma.reconstructDecay(decayString='pi0:dalitz -> e+:corrected ... ?gamma',
                        cut='',
                        chargeConjugation=ChargeC,
                        allowChargeViolation=True,
                        path=path)
    ma.applyCuts('pi0:dalitz', '', path=path)

    # Reconstruct D*0
    #-----------------------------------------
    ma.reconstructDecay(decayString='D*0:Mode2 -> [D0:kmpip -> K-:loose pi+:loose ] [pi0:dalitz -> e+:corrected ... ?gamma]', 
                        cut=Ds_Cut, 
                        chargeConjugation=ChargeC, 
                        path=path) 
    #chargeConjugation=False: The option is true by default
    #
    # Vertex Fitting D*0
    #-------------------
    vx.treeFit('D*0:Mode2', 
            conf_level=0, 
            ipConstraint=True,
            updateAllDaughters=False, 
            path=path)
    ma.applyCuts("D*0:Mode2", Ds_Cut, path=path)
    # Save some Variables Before IP Constraint
    ma.variablesToExtraInfo('D*0:Mode2', variables={"chiProb": "chiProb_noIP", "flightDistance": "flightDistance_noIP"}, path=path)
    vm.addAlias('chiProb_noIP', 'extraInfo(chiProb_noIP)')
    vm.addAlias('flightDistance_noIP', 'extraInfo(flightDistance_noIP)')
    # Save Variable After IP Constraint
    vx.treeFit('D*0:Mode2', 
                conf_level=0, 
                ipConstraint=True,
                updateAllDaughters=False, 
                path=path)

    # D*0 Best Candidate
    ma.rankByHighest('D*0:Mode2', 
                    variable='chiProb',
                    outputVariable='chiProb_Ds_rank',
                    # numBest=1,
                    path=path)
    vm.addAlias('chiProb_Ds_rank', 'extraInfo(chiProb_Ds_rank)')
    ma.rankByHighest('D*0:Mode2', 
                    variable='daughter(1,p)',
                    outputVariable='Electron_p_Ds_rank',
                    # numBest=1,
                    path=path)
    vm.addAlias('Electron_p_Ds_rank', 'extraInfo(Electron_p_Ds_rank)')
    ma.rankByHighest('D*0:Mode2', 
                    variable='random',
                    outputVariable='random_Ds_rank',
                    # numBest=1,
                    path=path)
    vm.addAlias('random_Ds_rank', 'extraInfo(random_Ds_rank)')
    BCS = ['chiProb_Ds_rank','Electron_p_Ds_rank','random_Ds_rank']
    #/================================================================================================================================/
    # Reconstruct vpho -> D*0 gamma
    #-----------------------------------------
    ma.reconstructDecay(decayString='vpho -> [D*0:Mode2 -> [D0:kmpip -> K-:loose pi+:loose ] [pi0:dalitz -> e+:corrected ... ?gamma]] gamma:recon', 
                        cut=Ds_star_Cut, 
                        chargeConjugation=ChargeC, 
                        allowChargeViolation=True,
                        path=path)
    #chargeConjugation=False: The option is true by default
    #
    # Vertex Fitting vpho
    # ----------------------
    # vx.treeFit('vpho', 
    #            conf_level=-1,
    #            path=my_path)
    ma.applyCuts("vpho", Ds_star_Cut, path=path)

    # vpho Best Candidate
    vm.addAlias('DeltaM_diff', 'formula(abs(massDifference(0) - 0.1438))')
    ma.rankByLowest('vpho', 
                    variable='DeltaM_diff',
                    # numBest=1,
                    path=path)
    vm.addAlias('DeltaM_diff_rank', 'extraInfo(DeltaM_diff_rank)')

    ma.rankByLowest('vpho', 
                    variable='dM',
                    # numBest=1,
                    path=path)
    vm.addAlias('dM_rank', 'extraInfo(dM_rank)')

    ma.rankByHighest('vpho', 
                    variable='random',
                    # numBest=1,
                    path=path)
    vm.addAlias('random_rank', 'extraInfo(random_rank)')
    #/================================================================================================================================/
    # Perform MC Matching (MC truth asociation)
    #----------------------------------------------
    ma.matchMCTruth(list_name='D*0:Mode2', path=path) 
    ma.matchMCTruth(list_name='vpho', path=path) 
    #/================================================================================================================================/
    # D*0 VETO Starts Here
    #------------------------

    # ROE
    #---------------------------------------------------------------------------------------------
    # build RestOfEvent (ROE) object for each D*0 candidate
    # ROE is required by the veto
    ma.buildRestOfEvent(target_list_name='D*0:Mode2',
                        fillWithMostLikely=False,
                        # chargedPIDPriors=[1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        path=path)
    # Create a mask tuple:
    track_based_cuts = "abs(d0) < 20.0 and abs(z0) < 20.0"
    # track_based_cuts = ""
    ecl_based_cuts = ""  # [E30]
    roe_mask = ("roe_mask", track_based_cuts, ecl_based_cuts)
    ma.appendROEMasks("D*0:Mode2", [roe_mask], path=path)  # [E40]

    ma.updateROEUsingV0Lists(target_particle_list="D*0:Mode2", 
                            mask_names="roe_mask", 
                            default_cleanup=True, 
                            selection_cuts=None, 
                            apply_mass_fit=True, 
                            fitter='treefit', 
                            path=path)
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
    ma.signalSideParticleFilter(particleList='D*0:Mode2',
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
    #---------------------------------------------------------------------------------------------

    # in order to be able to use modularAnalysis functions (reconstructDecay in particular)
    # we need a ParticleList containing the photon candidate used to reconstruct the
    # current D*0 as well
    # The DecayString is used to specify the selected particle (^)
    ma.fillSignalSideParticleList(outputListName='e+:sig',
                                decayString='D*0:Mode2 -> [D0:kmpip -> K-:loose pi+:loose ] [pi0:dalitz -> ^e+:corrected ... ?gamma]',
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
            updateAllDaughters=False,
            path=roe_path)
    ma.applyCuts('gamma:veto', 'InvM <= 2 and M<=3', path=roe_path)
    #-------------------------------------------------------------------------------------------

    # Perform MC Matching (MC truth asociation)
    ma.matchMCTruth(list_name='gamma:veto', path=roe_path)

    # Best Candidate Selection
    #-------------------------------------------------------
    # perform best candidate selection
    vm.addAlias("gamma_abs_dM", "formula(abs(M - 0.0))")
    ma.rankByLowest('gamma:veto', 
                    variable='gamma_abs_dM',
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

    # vm.addAlias("Goodgamma", "passesCut(mcPDG==22)")
    vm.addAlias("GoodElectron", "passesCut(abs(mcPDG)==11)")

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
                'daughterAngle(0,1)':'gammaveto_daughterAngle',
                'decayAngle(0)':'gammaveto_decayAngle_0',
                'decayAngle(1)':'gammaveto_decayAngle_1',
                'cos(decayAngle(0))':'gammaveto_cos_decayAngle_0',
                'cos(decayAngle(1))':'gammaveto_cos_decayAngle_1',
                'isSignal':'gammaveto_isSignal',
                'isSignalAcceptWrongFSPs':'gammaveto_isSignalAcceptWrongFSPs',
                'mcPDG':'gammaveto_mcPDG',
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
                'daughter(0,genMotherPDG)':'gammaveto_ep_genMotherPDG',
                'daughter(1,electronID)':'gammaveto_em_electronID',
                'daughter(1,d0)':'gammaveto_em_d0',
                'daughter(1,z0)':'gammaveto_em_z0',
                'daughter(1,M)':'gammaveto_em_M',
                'daughter(1,p)':'gammaveto_em_p',
                'daughter(1,mcP)':'gammaveto_em_mcP',
                'daughter(1,pt)':'gammaveto_em_pt',
                'daughter(1,E)':'gammaveto_em_E',
                'daughter(1,mcE)':'gammaveto_em_mcE',
                'daughter(1,mcPDG)':'gammaveto_em_mcPDG',
                'daughter(1,genMotherPDG)':'gammaveto_em_genMotherPDG',
                'daughter(1,mcSecPhysProc)':'gammaveto_em_mcSecPhysProc',
                'daughter(1,binaryPID(11,211))':'gammaveto_em_binaryPID',
                'daughter(1,GoodElectron)':'gammaveto_em_GoodElectron',
                'pointangle':'gammaveto_pointangle',
                'chiProb':'gammaveto_chiProb',
                'openangle':'gammaveto_openangle',
                'psi':'gammaveto_psi',
                'genNMissingDaughter(11)':'gammaveto_genNMissingDaughter',
                'daughterMotherDiffOf(0,theta)':'gammaveto_daughterMotherDiffOf_theta',
                'daughterMotherDiffOf(0,phi)':'gammaveto_daughterMotherDiffOf_phi'
                }
    ma.variableToSignalSideExtraInfo(particleList='gamma:veto', 
                                    varToExtraInfo=gamma_dict, 
                                    path=roe_path)
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

    vm.addAlias('gammaveto_M_Correction', 'ifNANgiveX(gammaveto_M,10)')
    gamma_ROE.append('gammaveto_M_Correction')
    #--------------------------------------------------------------------------------------------
    # D*0 Veto Ends Here
    #/================================================================================================================================/
    # vpho Info Saved to D*0
    ma.variablesToDaughterExtraInfo(particleList='vpho', 
                                    decayString='vpho -> [^D*0:Mode2 -> [D0:kmpip -> K-:loose pi+:loose ] [pi0:dalitz -> e+:corrected ... ?gamma]] gamma:recon', 
                                    variables={'massDifference(0)':'Ds_starminusDs',
                                            'DeltaM_diff_rank':'Ds_star_rank',
                                            }, 
                                    option=0, 
                                    path=path)
    vm.addAlias("Ds_starminusDs", "extraInfo(Ds_starminusDs)")
    vm.addAlias('Ds_starminusDs_M_Correction', 'ifNANgiveX(Ds_starminusDs,10)')

    vm.addAlias("Ds_star_rank", "extraInfo(Ds_star_rank)")
    vm.addAlias('Ds_star_rank_Correction', 'ifNANgiveX(Ds_star_rank,10)')

    vm.addAlias("goodDsplus", "passesCut(Ds_starminusDs_M_Correction >= 0.12 and Ds_starminusDs_M_Correction <= 0.165)")
    #/================================================================================================================================/
    # Other Tools:
    #-------------------
    # BDT
    #------------
    # MVAExpert
    path.add_module('MVAExpert', 
                    listNames=['D*0:Mode2'], 
                    extraInfoName='FakeD0BDT', 
                    identifier='user_adilmub_FakeD0Suppression')
    path.add_module('MVAExpert', 
                    listNames=['D*0:Mode2'], 
                    extraInfoName='BkgBDT', 
                    identifier='BkgMVA')
    # Variables from MVAExpert.
    expertVars = ['extraInfo(FakeD0BDT)','extraInfo(BkgBDT)'] #'transformedNetworkOutput(FastBDT,0.1,1.0)'
    #/================================================================================================================================/
    return BCS, gamma_ROE, expertVars

def Mode_3(path):
    #/================================================================================================================================/
    # Reconstruct D0 -> K- pi+ decay
    #------------------------------------
    ma.reconstructDecay(decayString='D0:kmpip -> K-:loose pi+:loose', 
                        cut=D0_Cut, 
                        chargeConjugation=ChargeC,
                        path=path)
    #chargeConjugation=False: The option is true by default
    #
    # Vertex Fitting D0
    #-------------------
    vx.treeFit('D0:kmpip',
            conf_level=0,
            updateAllDaughters=False,
            path=path)
    ma.applyCuts("D0:kmpip", D0_Cut, path=path)

    # Photon Conversion
    #--------------------------
    ma.reconstructDecay(decayString='gamma:RD_partial -> e+:corrected ...',
                        cut='',
                        chargeConjugation=ChargeC,
                        allowChargeViolation=True,
                        path=path)
    ma.applyCuts('gamma:RD_partial', '', path=path)

    # Reconstruct D*0
    #-----------------------------------------
    ma.reconstructDecay(decayString='D*0:Mode3 -> [D0:kmpip -> K-:loose pi+:loose ] [gamma:RD_partial -> e+:corrected ...]', 
                        cut=Ds_Cut, 
                        chargeConjugation=ChargeC, 
                        path=path) 
    #chargeConjugation=False: The option is true by default
    #
    # Vertex Fitting D*0
    #-------------------
    vx.treeFit('D*0:Mode3', 
            conf_level=0, 
            ipConstraint=True,
            updateAllDaughters=False, 
            path=path)
    ma.applyCuts("D*0:Mode3", Ds_Cut, path=path)
    # Save some Variables Before IP Constraint
    ma.variablesToExtraInfo('D*0:Mode3', variables={"chiProb": "chiProb_noIP", "flightDistance": "flightDistance_noIP"}, path=path)
    vm.addAlias('chiProb_noIP', 'extraInfo(chiProb_noIP)')
    vm.addAlias('flightDistance_noIP', 'extraInfo(flightDistance_noIP)')
    # Save Variable After IP Constraint
    vx.treeFit('D*0:Mode3', 
                conf_level=0, 
                ipConstraint=True,
                updateAllDaughters=False, 
                path=path)

    # D*0 Best Candidate
    ma.rankByHighest('D*0:Mode3', 
                    variable='chiProb',
                    outputVariable='chiProb_Ds_rank',
                    # numBest=1,
                    path=path)
    vm.addAlias('chiProb_Ds_rank', 'extraInfo(chiProb_Ds_rank)')
    ma.rankByHighest('D*0:Mode3', 
                    variable='daughter(1,p)',
                    outputVariable='Electron_p_Ds_rank',
                    # numBest=1,
                    path=path)
    vm.addAlias('Electron_p_Ds_rank', 'extraInfo(Electron_p_Ds_rank)')
    ma.rankByHighest('D*0:Mode3', 
                    variable='random',
                    outputVariable='random_Ds_rank',
                    # numBest=1,
                    path=path)
    vm.addAlias('random_Ds_rank', 'extraInfo(random_Ds_rank)')
    BCS = ['chiProb_Ds_rank','Electron_p_Ds_rank','random_Ds_rank']
    #/================================================================================================================================/
    # Reconstruct vpho -> D*0 gamma
    #-----------------------------------------
    ma.reconstructDecay(decayString='vpho -> [D*0:Mode3 -> [D0:kmpip -> K-:loose pi+:loose ] [gamma:RD_partial -> e+:corrected ...]] gamma:recon', 
                        cut=Ds_star_Cut, 
                        chargeConjugation=ChargeC, 
                        allowChargeViolation=True,
                        path=path)
    #chargeConjugation=False: The option is true by default
    #
    # Vertex Fitting vpho
    # ----------------------
    # vx.treeFit('vpho', 
    #            conf_level=-1,
    #            path=my_path)
    ma.applyCuts("vpho", Ds_star_Cut, path=path)

    # vpho Best Candidate
    vm.addAlias('DeltaM_diff', 'formula(abs(massDifference(0) - 0.1438))')
    ma.rankByLowest('vpho', 
                    variable='DeltaM_diff',
                    # numBest=1,
                    path=path)
    vm.addAlias('DeltaM_diff_rank', 'extraInfo(DeltaM_diff_rank)')

    ma.rankByLowest('vpho', 
                    variable='dM',
                    # numBest=1,
                    path=path)
    vm.addAlias('dM_rank', 'extraInfo(dM_rank)')

    ma.rankByHighest('vpho', 
                    variable='random',
                    # numBest=1,
                    path=path)
    vm.addAlias('random_rank', 'extraInfo(random_rank)')
    #/================================================================================================================================/
    # Perform MC Matching (MC truth asociation)
    #----------------------------------------------
    ma.matchMCTruth(list_name='D*0:Mode3', path=path) 
    ma.matchMCTruth(list_name='vpho', path=path) 
    #/================================================================================================================================/
    # D*0 VETO Starts Here
    #------------------------

    # ROE
    #---------------------------------------------------------------------------------------------
    # build RestOfEvent (ROE) object for each D*0 candidate
    # ROE is required by the veto
    ma.buildRestOfEvent(target_list_name='D*0:Mode3',
                        fillWithMostLikely=False,
                        # chargedPIDPriors=[1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        path=path)
    # Create a mask tuple:
    track_based_cuts = "abs(d0) < 20.0 and abs(z0) < 20.0"
    # track_based_cuts = ""
    ecl_based_cuts = ""  # [E30]
    roe_mask = ("roe_mask", track_based_cuts, ecl_based_cuts)
    ma.appendROEMasks("D*0:Mode3", [roe_mask], path=path)  # [E40]

    ma.updateROEUsingV0Lists(target_particle_list="D*0:Mode3", 
                            mask_names="roe_mask", 
                            default_cleanup=True, 
                            selection_cuts=None, 
                            apply_mass_fit=True, 
                            fitter='treefit', 
                            path=path)
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
    ma.signalSideParticleFilter(particleList='D*0:Mode3',
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
    #---------------------------------------------------------------------------------------------

    # in order to be able to use modularAnalysis functions (reconstructDecay in particular)
    # we need a ParticleList containing the photon candidate used to reconstruct the
    # current D*0 as well
    # The DecayString is used to specify the selected particle (^)
    ma.fillSignalSideParticleList(outputListName='e+:sig',
                                decayString='D*0:Mode3 -> [D0:kmpip -> K-:loose pi+:loose ] [gamma:RD_partial -> ^e+:corrected ...]',
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
            updateAllDaughters=False,
            path=roe_path)
    ma.applyCuts('gamma:veto', 'InvM <= 2 and M<=3', path=roe_path)
    #-------------------------------------------------------------------------------------------

    # Perform MC Matching (MC truth asociation)
    ma.matchMCTruth(list_name='gamma:veto', path=roe_path)

    # Best Candidate Selection
    #-------------------------------------------------------
    # perform best candidate selection
    vm.addAlias("gamma_abs_dM", "formula(abs(M - 0.0))")
    ma.rankByLowest('gamma:veto', 
                    variable='gamma_abs_dM',
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

    # vm.addAlias("Goodgamma", "passesCut(mcPDG==22)")
    vm.addAlias("GoodElectron", "passesCut(abs(mcPDG)==11)")

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
                'daughterAngle(0,1)':'gammaveto_daughterAngle',
                'decayAngle(0)':'gammaveto_decayAngle_0',
                'decayAngle(1)':'gammaveto_decayAngle_1',
                'cos(decayAngle(0))':'gammaveto_cos_decayAngle_0',
                'cos(decayAngle(1))':'gammaveto_cos_decayAngle_1',
                'isSignal':'gammaveto_isSignal',
                'isSignalAcceptWrongFSPs':'gammaveto_isSignalAcceptWrongFSPs',
                'mcPDG':'gammaveto_mcPDG',
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
                'daughter(0,genMotherPDG)':'gammaveto_ep_genMotherPDG',
                'daughter(1,electronID)':'gammaveto_em_electronID',
                'daughter(1,d0)':'gammaveto_em_d0',
                'daughter(1,z0)':'gammaveto_em_z0',
                'daughter(1,M)':'gammaveto_em_M',
                'daughter(1,p)':'gammaveto_em_p',
                'daughter(1,mcP)':'gammaveto_em_mcP',
                'daughter(1,pt)':'gammaveto_em_pt',
                'daughter(1,E)':'gammaveto_em_E',
                'daughter(1,mcE)':'gammaveto_em_mcE',
                'daughter(1,mcPDG)':'gammaveto_em_mcPDG',
                'daughter(1,genMotherPDG)':'gammaveto_em_genMotherPDG',
                'daughter(1,mcSecPhysProc)':'gammaveto_em_mcSecPhysProc',
                'daughter(1,binaryPID(11,211))':'gammaveto_em_binaryPID',
                'daughter(1,GoodElectron)':'gammaveto_em_GoodElectron',
                'pointangle':'gammaveto_pointangle',
                'chiProb':'gammaveto_chiProb',
                'openangle':'gammaveto_openangle',
                'psi':'gammaveto_psi',
                'genNMissingDaughter(11)':'gammaveto_genNMissingDaughter',
                'daughterMotherDiffOf(0,theta)':'gammaveto_daughterMotherDiffOf_theta',
                'daughterMotherDiffOf(0,phi)':'gammaveto_daughterMotherDiffOf_phi'
                }
    ma.variableToSignalSideExtraInfo(particleList='gamma:veto', 
                                    varToExtraInfo=gamma_dict, 
                                    path=roe_path)
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

    vm.addAlias('gammaveto_M_Correction', 'ifNANgiveX(gammaveto_M,10)')
    gamma_ROE.append('gammaveto_M_Correction')
    #--------------------------------------------------------------------------------------------
    # D*0 Veto Ends Here
    #/================================================================================================================================/
    # vpho Info Saved to D*0
    ma.variablesToDaughterExtraInfo(particleList='vpho', 
                                    decayString='vpho -> [^D*0:Mode3 -> [D0:kmpip -> K-:loose pi+:loose ] [gamma:RD_partial -> e+:corrected ...]] gamma:recon', 
                                    variables={'massDifference(0)':'Ds_starminusDs',
                                            'DeltaM_diff_rank':'Ds_star_rank',
                                            }, 
                                    option=0, 
                                    path=path)
    vm.addAlias("Ds_starminusDs", "extraInfo(Ds_starminusDs)")
    vm.addAlias('Ds_starminusDs_M_Correction', 'ifNANgiveX(Ds_starminusDs,10)')

    vm.addAlias("Ds_star_rank", "extraInfo(Ds_star_rank)")
    vm.addAlias('Ds_star_rank_Correction', 'ifNANgiveX(Ds_star_rank,10)')

    vm.addAlias("goodDsplus", "passesCut(Ds_starminusDs_M_Correction >= 0.12 and Ds_starminusDs_M_Correction <= 0.165)")
    #/================================================================================================================================/
    # Other Tools:
    #-------------------
    # BDT
    #------------
    # MVAExpert
    path.add_module('MVAExpert', 
                    listNames=['D*0:Mode3'], 
                    extraInfoName='FakeD0BDT', 
                    identifier='user_adilmub_FakeD0Suppression')
    path.add_module('MVAExpert', 
                    listNames=['D*0:Mode3'], 
                    extraInfoName='BkgBDT', 
                    identifier='BkgMVA')
    # Variables from MVAExpert.
    expertVars = ['extraInfo(FakeD0BDT)','extraInfo(BkgBDT)'] #'transformedNetworkOutput(FastBDT,0.1,1.0)'
    #/================================================================================================================================/
    return BCS, gamma_ROE, expertVars

def Mode_1_Variable(BCS, gamma_ROE, expertVars):
    #/================================================================================================================================/
    # Save Variables
    #-------------------

    # D_s+ Variables
    #-------------------
    tracks = ['isCloneTrack',
            'dr','dz','abs(dr)','abs(dz)','z0','d0','pValue',
            'firstCDCLayer','firstPXDLayer','firstSVDLayer',
            'nPXDHits','nVXDHits','nSVDHits','nCDCHits',
            'seenInCDC','seenInPXD','seenInSVD',
            'inARICHAcceptance','inCDCAcceptance','inTOPAcceptance']

    Ds_vars=[]
    var_1 = ['ImpactXY','cos(theta)','phi','mcP','M','pt','p','px','py','pz','abs(pz)','isOrHasCloneTrack',"charge"]
    truth = ["mcSecPhysProc",'nMCMatches','mcPrimary','isSignal','mcPDG','genMotherPDG','genMotherID','mcErrors',
            'mcMatchWeight']

    vm.addAlias("pminusptrue_p", "formula(p - mcP)")
    vm.addAlias("pminusptrue_old", "formula(daughter(0,p) - mcP)")

    vm.addAlias("Motherptrueminusp_p", "formula(mcMother(p) - p)")
    vm.addAlias("Motherptrueminusp_E", "formula(mcMother(E) - E)")
    Ds_vars += vu.create_aliases_for_selected(
        list_of_variables=['chi2','ndf',
                        'electronID','binaryPID(11,211)',
                        'omega',
                        'daughter(0,isCloneTrack)',
                        #    'sourceObjectIsInList(e+:V0_array)',
                        'mcVirtual',
                        #    'pidMostLikelyPDG(ePrior=1/6, muPrior=1/6, piPrior=1/6, KPrior=1/6, pPrior=1/6, dPrior=1/6)',"isBremsCorrected",
                        'isMisidentified',
                        "pminusptrue_p","pminusptrue_old",
                        "Motherptrueminusp_p","Motherptrueminusp_E",
                        'genMotherPDG(1)','genMotherPDG(2)',
                        'genMotherID(1)','genMotherID(2)']
                        +  var_1 + tracks + truth ,
        decay_string='D*0:Mode1 -> [D0:kmpip -> K-:loose pi+:loose ] [pi0:2photons -> [gamma:RD_partial -> ^e+:corrected ...] ?gamma]',
        prefix=['e'])

    Ds_vars += vu.create_aliases_for_selected(
        list_of_variables=['theta','phi','omega','kaonID'] + var_1 + tracks + truth,
        decay_string='D*0:Mode1 -> [D0:kmpip -> ^K-:loose pi+:loose ] [pi0:2photons -> [gamma:RD_partial -> e+:corrected ...] ?gamma]',
        prefix=['K'])

    Ds_vars += vu.create_aliases_for_selected(
        list_of_variables=['theta','phi','omega','pionID'] + var_1 + tracks + truth,
        decay_string='D*0:Mode1 -> [D0:kmpip -> K-:loose ^pi+:loose ] [pi0:2photons -> [gamma:RD_partial -> e+:corrected ...] ?gamma]',
        prefix=['pi'])   

    Ds_vars += vu.create_aliases_for_selected(
        list_of_variables=['useCMSFrame(p)',
                        "dM","useAlternativeDaughterHypothesis(M, 1:K+)",
                        "chiProb",
                        "decayAngle(0)","cos(decayAngle(0))",
                        "decayAngle(1)","cos(decayAngle(1))",
                        "flightDistance","distance",
                        "useRestFrame(daughterAngle(0, 1))"]
                        + var_1 
                        + ['genMotherPDG(1)','genMotherPDG(2)','genMotherID(1)','genMotherID(2)'] 
                        + truth + ['ifNANgiveX(isSignal,5)','D0Mode','Dbar0Mode','nMCDaughters'],
        decay_string='D*0:Mode1 -> [^D0:kmpip -> K-:loose pi+:loose ] [pi0:2photons -> [gamma:RD_partial -> e+:corrected ...] ?gamma]',
        prefix=['D0'])

    vm.addAlias("D0_sideband", "passesCut(-0.02 <= D0_dM <= 0.02)")

    vm.addAlias("D0_Dstarplus", "passesCut(abs(D0_mcPDG)==421 and abs(D0_genMotherPDG)==413)")
    vm.addAlias("D0_Dstar0", "passesCut(abs(D0_mcPDG)==421 and abs(D0_genMotherPDG)==423)")
    vm.addAlias("D0_NoDstarplusDstar0", "passesCut(abs(D0_mcPDG)==421 and abs(D0_genMotherPDG)!=413 and abs(D0_genMotherPDG)!=423)")
    vm.addAlias("D0_Other", "passesCut(abs(D0_mcPDG)!=421)")

    vm.addAlias("Comb", "passesCut(abs(mcPDG)==23 or abs(mcPDG)==300553)")
    vm.addAlias("Failed", "passesCut(ifNANgiveX(isSignal,5)==5)")
    vm.addAlias("Dstar0", "passesCut(abs(mcPDG)==423)")
    vm.addAlias("Dstarplus", "passesCut(abs(mcPDG)==413)")
    vm.addAlias("Other", "passesCut(Comb!=1 and Failed!=1 and Dstar0!=1 and Dstarplus!=1)")
    vm.addAlias("Signal", "passesCut(abs(mcPDG)==431)")

    vm.addAlias("L_diff", "formula(((x - daughter(0,x))**2 + (y - daughter(0,y))**2 + (z - daughter(0,z))**2)**(1/2))")

    vm.addAlias("TrueDstarplus", "passesCut(DstplusMode==1001 or DstminusMode==-1001)")
    Ds_vars += vu.create_aliases_for_selected(
        list_of_variables=["chiProb_noIP","flightDistance_noIP",
                        'chiProb',
                        #    'originalDaughter(1, cos(theta))',
                        'InvMLambda',
                        'massDifference(0)',"phi_diff",
                        "M_uncorrected",
                        "diff_D0pi",
                        'azimuthalAngleInDecayPlane(0,1)',
                        'daughterDiffOf(0, 1, cos(theta))',
                        'useRestFrame(daughterDiffOf(0, 1, p))',
                        'useRestFrame(daughterMotherDiffOf(0, p))',
                        'flightDistance','distance',
                        'daughterDiffOf(0,1,x)','daughterDiffOf(0,1,y)','daughterDiffOf(0,1,z)',
                        'daughterMotherDiffOf(0,x)','daughterMotherDiffOf(0,y)','daughterMotherDiffOf(0,z)',
                        'abs(daughterMotherDiffOf(0,distance))',"L_diff",
                        'daughterMotherDiffOf(0,flightDistance)','daughterMotherDiffOf(0,vertexDistance)',
                        'flightDistanceOfDaughter(0)']
                        + BCS
                        + gamma_ROE
                        + ["Ds_starminusDs","Ds_starminusDs_M_Correction",
                            "Ds_star_rank","Ds_star_rank_Correction",
                            "goodDsplus"] 
                        + var_1 + truth 
                        + ['genNStepsToDaughter(0)','genNStepsToDaughter(1)',
                        'genNMissingDaughter(11)','genNMissingDaughter(22)',
                        'ifNANgiveX(isSignal,5)',
                        'nMCDaughters','genParticleID'] 
                        + ["Comb","Failed","Dstar0","Dstarplus","Other","Signal","D0_sideband"] 
                        + ["D0_Dstarplus","D0_Dstar0","D0_NoDstarplusDstar0","D0_Other"] 
                        + ["DstplusMode","DstminusMode","TrueDstarplus"]
                        + expertVars,
        decay_string='^D*0:Mode1 -> [D0:kmpip -> K-:loose pi+:loose ] [pi0:2photons -> [gamma:RD_partial -> e+:corrected ...] ?gamma]',
        prefix=['Ds'])

    # Extra Variables
    #--------------------------------------------------------------------------------------------
    # Some extra D_s+ Variables
    vm.addAlias("diff_D0pi", "useAlternativeDaughterHypothesis(massDifference(0), 1:pi+)")
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
    return Ds_vars, Event

def Mode_2_Variable(BCS, gamma_ROE, expertVars):
    #/================================================================================================================================/
    # Save Variables
    #-------------------

    # D_s+ Variables
    #-------------------
    tracks = ['isCloneTrack',
            'dr','dz','abs(dr)','abs(dz)','z0','d0','pValue',
            'firstCDCLayer','firstPXDLayer','firstSVDLayer',
            'nPXDHits','nVXDHits','nSVDHits','nCDCHits',
            'seenInCDC','seenInPXD','seenInSVD',
            'inARICHAcceptance','inCDCAcceptance','inTOPAcceptance']

    Ds_vars=[]
    var_1 = ['ImpactXY','cos(theta)','phi','mcP','M','pt','p','px','py','pz','abs(pz)','isOrHasCloneTrack',"charge"]
    truth = ["mcSecPhysProc",'nMCMatches','mcPrimary','isSignal','mcPDG','genMotherPDG','genMotherID','mcErrors',
            'mcMatchWeight']

    vm.addAlias("pminusptrue_p", "formula(p - mcP)")
    vm.addAlias("pminusptrue_old", "formula(daughter(0,p) - mcP)")

    vm.addAlias("Motherptrueminusp_p", "formula(mcMother(p) - p)")
    vm.addAlias("Motherptrueminusp_E", "formula(mcMother(E) - E)")
    Ds_vars += vu.create_aliases_for_selected(
        list_of_variables=['chi2','ndf',
                        'electronID','binaryPID(11,211)',
                        'omega',
                        'daughter(0,isCloneTrack)',
                        #    'sourceObjectIsInList(e+:V0_array)',
                        'mcVirtual',
                        #    'pidMostLikelyPDG(ePrior=1/6, muPrior=1/6, piPrior=1/6, KPrior=1/6, pPrior=1/6, dPrior=1/6)',"isBremsCorrected",
                        'isMisidentified',
                        "pminusptrue_p","pminusptrue_old",
                        "Motherptrueminusp_p","Motherptrueminusp_E",
                        'genMotherPDG(1)','genMotherPDG(2)',
                        'genMotherID(1)','genMotherID(2)']
                        +  var_1 + tracks + truth ,
        decay_string='D*0:Mode2 -> [D0:kmpip -> K-:loose pi+:loose ] [pi0:dalitz -> ^e+:corrected ... ?gamma]',
        prefix=['e'])

    Ds_vars += vu.create_aliases_for_selected(
        list_of_variables=['theta','phi','omega','kaonID'] + var_1 + tracks + truth,
        decay_string='D*0:Mode2 -> [D0:kmpip -> ^K-:loose pi+:loose ] [pi0:dalitz -> e+:corrected ... ?gamma]',
        prefix=['K'])

    Ds_vars += vu.create_aliases_for_selected(
        list_of_variables=['theta','phi','omega','pionID'] + var_1 + tracks + truth,
        decay_string='D*0:Mode2 -> [D0:kmpip -> K-:loose ^pi+:loose ] [pi0:dalitz -> e+:corrected ... ?gamma]',
        prefix=['pi'])   

    Ds_vars += vu.create_aliases_for_selected(
        list_of_variables=['useCMSFrame(p)',
                        "dM","useAlternativeDaughterHypothesis(M, 1:K+)",
                        "chiProb",
                        "decayAngle(0)","cos(decayAngle(0))",
                        "decayAngle(1)","cos(decayAngle(1))",
                        "flightDistance","distance",
                        "useRestFrame(daughterAngle(0, 1))"]
                        + var_1 
                        + ['genMotherPDG(1)','genMotherPDG(2)','genMotherID(1)','genMotherID(2)'] 
                        + truth + ['ifNANgiveX(isSignal,5)','D0Mode','Dbar0Mode','nMCDaughters'],
        decay_string='D*0:Mode2 -> [^D0:kmpip -> K-:loose pi+:loose ] [pi0:dalitz -> e+:corrected ... ?gamma]',
        prefix=['D0'])

    vm.addAlias("D0_sideband", "passesCut(-0.02 <= D0_dM <= 0.02)")

    vm.addAlias("D0_Dstarplus", "passesCut(abs(D0_mcPDG)==421 and abs(D0_genMotherPDG)==413)")
    vm.addAlias("D0_Dstar0", "passesCut(abs(D0_mcPDG)==421 and abs(D0_genMotherPDG)==423)")
    vm.addAlias("D0_NoDstarplusDstar0", "passesCut(abs(D0_mcPDG)==421 and abs(D0_genMotherPDG)!=413 and abs(D0_genMotherPDG)!=423)")
    vm.addAlias("D0_Other", "passesCut(abs(D0_mcPDG)!=421)")

    vm.addAlias("Comb", "passesCut(abs(mcPDG)==23 or abs(mcPDG)==300553)")
    vm.addAlias("Failed", "passesCut(ifNANgiveX(isSignal,5)==5)")
    vm.addAlias("Dstar0", "passesCut(abs(mcPDG)==423)")
    vm.addAlias("Dstarplus", "passesCut(abs(mcPDG)==413)")
    vm.addAlias("Other", "passesCut(Comb!=1 and Failed!=1 and Dstar0!=1 and Dstarplus!=1)")
    vm.addAlias("Signal", "passesCut(abs(mcPDG)==431)")

    vm.addAlias("L_diff", "formula(((x - daughter(0,x))**2 + (y - daughter(0,y))**2 + (z - daughter(0,z))**2)**(1/2))")

    vm.addAlias("TrueDstarplus", "passesCut(DstplusMode==1001 or DstminusMode==-1001)")
    Ds_vars += vu.create_aliases_for_selected(
        list_of_variables=["chiProb_noIP","flightDistance_noIP",
                        'chiProb',
                        #    'originalDaughter(1, cos(theta))',
                        'InvMLambda',
                        'massDifference(0)',"phi_diff",
                        "M_uncorrected",
                        "diff_D0pi",
                        'azimuthalAngleInDecayPlane(0,1)',
                        'daughterDiffOf(0, 1, cos(theta))',
                        'useRestFrame(daughterDiffOf(0, 1, p))',
                        'useRestFrame(daughterMotherDiffOf(0, p))',
                        'flightDistance','distance',
                        'daughterDiffOf(0,1,x)','daughterDiffOf(0,1,y)','daughterDiffOf(0,1,z)',
                        'daughterMotherDiffOf(0,x)','daughterMotherDiffOf(0,y)','daughterMotherDiffOf(0,z)',
                        'abs(daughterMotherDiffOf(0,distance))',"L_diff",
                        'daughterMotherDiffOf(0,flightDistance)','daughterMotherDiffOf(0,vertexDistance)',
                        'flightDistanceOfDaughter(0)']
                        + BCS
                        + gamma_ROE
                        + ["Ds_starminusDs","Ds_starminusDs_M_Correction",
                            "Ds_star_rank","Ds_star_rank_Correction",
                            "goodDsplus"] 
                        + var_1 + truth 
                        + ['genNStepsToDaughter(0)','genNStepsToDaughter(1)',
                        'genNMissingDaughter(11)','genNMissingDaughter(22)',
                        'ifNANgiveX(isSignal,5)',
                        'nMCDaughters','genParticleID'] 
                        + ["Comb","Failed","Dstar0","Dstarplus","Other","Signal","D0_sideband"] 
                        + ["D0_Dstarplus","D0_Dstar0","D0_NoDstarplusDstar0","D0_Other"] 
                        + ["DstplusMode","DstminusMode","TrueDstarplus"]
                        + expertVars,
        decay_string='^D*0:Mode2 -> [D0:kmpip -> K-:loose pi+:loose ] [pi0:dalitz -> e+:corrected ... ?gamma]',
        prefix=['Ds'])

    # Extra Variables
    #--------------------------------------------------------------------------------------------
    # Some extra D_s+ Variables
    vm.addAlias("diff_D0pi", "useAlternativeDaughterHypothesis(massDifference(0), 1:pi+)")
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
    return Ds_vars, Event

def Mode_3_Variable(BCS, gamma_ROE, expertVars):
    #/================================================================================================================================/
    # Save Variables
    #-------------------

    # D_s+ Variables
    #-------------------
    tracks = ['isCloneTrack',
            'dr','dz','abs(dr)','abs(dz)','z0','d0','pValue',
            'firstCDCLayer','firstPXDLayer','firstSVDLayer',
            'nPXDHits','nVXDHits','nSVDHits','nCDCHits',
            'seenInCDC','seenInPXD','seenInSVD',
            'inARICHAcceptance','inCDCAcceptance','inTOPAcceptance']

    Ds_vars=[]
    var_1 = ['ImpactXY','cos(theta)','phi','mcP','M','pt','p','px','py','pz','abs(pz)','isOrHasCloneTrack',"charge"]
    truth = ["mcSecPhysProc",'nMCMatches','mcPrimary','isSignal','mcPDG','genMotherPDG','genMotherID','mcErrors',
            'mcMatchWeight']

    vm.addAlias("pminusptrue_p", "formula(p - mcP)")
    vm.addAlias("pminusptrue_old", "formula(daughter(0,p) - mcP)")

    vm.addAlias("Motherptrueminusp_p", "formula(mcMother(p) - p)")
    vm.addAlias("Motherptrueminusp_E", "formula(mcMother(E) - E)")
    Ds_vars += vu.create_aliases_for_selected(
        list_of_variables=['chi2','ndf',
                        'electronID','binaryPID(11,211)',
                        'omega',
                        'daughter(0,isCloneTrack)',
                        #    'sourceObjectIsInList(e+:V0_array)',
                        'mcVirtual',
                        #    'pidMostLikelyPDG(ePrior=1/6, muPrior=1/6, piPrior=1/6, KPrior=1/6, pPrior=1/6, dPrior=1/6)',"isBremsCorrected",
                        'isMisidentified',
                        "pminusptrue_p","pminusptrue_old",
                        "Motherptrueminusp_p","Motherptrueminusp_E",
                        'genMotherPDG(1)','genMotherPDG(2)',
                        'genMotherID(1)','genMotherID(2)']
                        +  var_1 + tracks + truth ,
        decay_string='D*0 -> [D0:kmpip -> K-:loose pi+:loose ] [gamma:RD_partial -> ^e+:corrected ...]',
        prefix=['e'])

    Ds_vars += vu.create_aliases_for_selected(
        list_of_variables=['theta','phi','omega','kaonID'] + var_1 + tracks + truth,
        decay_string='D*0 -> [D0:kmpip -> ^K-:loose pi+:loose ] [gamma:RD_partial -> e+:corrected ...]',
        prefix=['K'])

    Ds_vars += vu.create_aliases_for_selected(
        list_of_variables=['theta','phi','omega','pionID'] + var_1 + tracks + truth,
        decay_string='D*0 -> [D0:kmpip -> K-:loose ^pi+:loose ] [gamma:RD_partial -> e+:corrected ...]',
        prefix=['pi'])   

    Ds_vars += vu.create_aliases_for_selected(
        list_of_variables=['useCMSFrame(p)',
                        "dM","useAlternativeDaughterHypothesis(M, 1:K+)",
                        "chiProb",
                        "decayAngle(0)","cos(decayAngle(0))",
                        "decayAngle(1)","cos(decayAngle(1))",
                        "flightDistance","distance",
                        "useRestFrame(daughterAngle(0, 1))"]
                        + var_1 
                        + ['genMotherPDG(1)','genMotherPDG(2)','genMotherID(1)','genMotherID(2)'] 
                        + truth + ['ifNANgiveX(isSignal,5)','D0Mode','Dbar0Mode','nMCDaughters'],
        decay_string='D*0 -> [^D0:kmpip -> K-:loose pi+:loose ] [gamma:RD_partial -> e+:corrected ...]',
        prefix=['D0'])

    vm.addAlias("D0_sideband", "passesCut(-0.02 <= D0_dM <= 0.02)")

    vm.addAlias("D0_Dstarplus", "passesCut(abs(D0_mcPDG)==421 and abs(D0_genMotherPDG)==413)")
    vm.addAlias("D0_Dstar0", "passesCut(abs(D0_mcPDG)==421 and abs(D0_genMotherPDG)==423)")
    vm.addAlias("D0_NoDstarplusDstar0", "passesCut(abs(D0_mcPDG)==421 and abs(D0_genMotherPDG)!=413 and abs(D0_genMotherPDG)!=423)")
    vm.addAlias("D0_Other", "passesCut(abs(D0_mcPDG)!=421)")

    vm.addAlias("Comb", "passesCut(abs(mcPDG)==23 or abs(mcPDG)==300553)")
    vm.addAlias("Failed", "passesCut(ifNANgiveX(isSignal,5)==5)")
    vm.addAlias("Dstar0", "passesCut(abs(mcPDG)==423)")
    vm.addAlias("Dstarplus", "passesCut(abs(mcPDG)==413)")
    vm.addAlias("Other", "passesCut(Comb!=1 and Failed!=1 and Dstar0!=1 and Dstarplus!=1)")
    vm.addAlias("Signal", "passesCut(abs(mcPDG)==431)")

    vm.addAlias("L_diff", "formula(((x - daughter(0,x))**2 + (y - daughter(0,y))**2 + (z - daughter(0,z))**2)**(1/2))")

    vm.addAlias("TrueDstarplus", "passesCut(DstplusMode==1001 or DstminusMode==-1001)")
    Ds_vars += vu.create_aliases_for_selected(
        list_of_variables=["chiProb_noIP","flightDistance_noIP",
                        'chiProb',
                        #    'originalDaughter(1, cos(theta))',
                        'InvMLambda',
                        'massDifference(0)',"phi_diff",
                        "M_uncorrected",
                        "diff_D0pi",
                        'azimuthalAngleInDecayPlane(0,1)',
                        'daughterDiffOf(0, 1, cos(theta))',
                        'useRestFrame(daughterDiffOf(0, 1, p))',
                        'useRestFrame(daughterMotherDiffOf(0, p))',
                        'flightDistance','distance',
                        'daughterDiffOf(0,1,x)','daughterDiffOf(0,1,y)','daughterDiffOf(0,1,z)',
                        'daughterMotherDiffOf(0,x)','daughterMotherDiffOf(0,y)','daughterMotherDiffOf(0,z)',
                        'abs(daughterMotherDiffOf(0,distance))',"L_diff",
                        'daughterMotherDiffOf(0,flightDistance)','daughterMotherDiffOf(0,vertexDistance)',
                        'flightDistanceOfDaughter(0)']
                        + BCS
                        + gamma_ROE
                        + ["Ds_starminusDs","Ds_starminusDs_M_Correction",
                            "Ds_star_rank","Ds_star_rank_Correction",
                            "goodDsplus"] 
                        + var_1 + truth 
                        + ['genNStepsToDaughter(0)','genNStepsToDaughter(1)',
                        'genNMissingDaughter(11)','genNMissingDaughter(22)',
                        'ifNANgiveX(isSignal,5)',
                        'nMCDaughters','genParticleID'] 
                        + ["Comb","Failed","Dstar0","Dstarplus","Other","Signal","D0_sideband"] 
                        + ["D0_Dstarplus","D0_Dstar0","D0_NoDstarplusDstar0","D0_Other"] 
                        + ["DstplusMode","DstminusMode","TrueDstarplus"]
                        + expertVars,
        decay_string='^D*0 -> [D0:kmpip -> K-:loose pi+:loose ] [gamma:RD_partial -> e+:corrected ...]',
        prefix=['Ds'])

    # Extra Variables
    #--------------------------------------------------------------------------------------------
    # Some extra D_s+ Variables
    vm.addAlias("diff_D0pi", "useAlternativeDaughterHypothesis(massDifference(0), 1:pi+)")
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
    return Ds_vars, Event

if __name__ == '__main__':

    args = argparser().parse_args()

    import basf2
    import modularAnalysis as ma
    from ROOT import Belle2
    import pdg

    #/==========================================================================================================================/
    # Create  Main Path
    #--------------------
    path = basf2.create_path()
    #/================================================================================================================================/
    # --I/O----------------------------------------------------------------------------------------
    # load input ROOT file
    ma.inputMdst(environmentType='default',
                filename="/home/belle2/amubarak/C00-Generation/ccbarDs+EventGeneration.root",
                path=path)
    #/================================================================================================================================/
    # Create Particle Lists
    #---------------------------------
    # Pion
    # Creates Pion Particle List (and c.c.)
    # Cut on Variables: pionID, abs(d0), and abs(z0)
    # ma.fillParticleList('pi+','', path=path)
    ma.fillParticleList('pi+:loose', cut=Pion_Cut, path=path)
    # stdc.stdPi(listtype='loose', path=path)
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

    # ma.tagCurlTracks("e+:corrected", selectorType='mva', ptCut=0.6, mcTruth=True, expert_filename='test.root', path=path)
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
    # Reconstruction
    #-----------------------------
    # if args.Mode==0:
    #     Mode_All(path)
    if args.Mode==1:
        gamma_ROE, expertVars = Mode_1(path)
        Ds_vars, Event = Mode_1_Variable(gamma_ROE, expertVars)
    elif args.Mode==2:
        gamma_ROE, expertVars = Mode_2(path)
        Ds_vars, Event = Mode_2_Variable(gamma_ROE, expertVars)
    elif args.Mode==3:
        gamma_ROE, expertVars = Mode_3(path)
        Ds_vars, Event = Mode_3_Variable(gamma_ROE, expertVars)
    #/================================================================================================================================/
    # Save
    #--------
    # Saving variables to ntuple

    # if args.Mode==0:
    #     output_file = 'C01-Simulated_Events/Test/output_ModeAll.root'
    #     ma.variablesToNtuple('D*0:All', ['M','isSignal','massDifference(0)',
    #                                      'useAlternativeDaughterHypothesis(massDifference(0), 1:pi+)','mcPDG','dM'],
    #                          filename=output_file, treename='Dstree', path=path)
    #     ma.variablesToNtuple('D*0:All_Peak', ['M','isSignal','massDifference(0)',
    #                                           'useAlternativeDaughterHypothesis(massDifference(0), 1:pi+)','mcPDG','dM'],
    #                          filename=output_file, treename='Dstree_Peak', path=path)
    if args.Mode==1:
        output_file = 'C01-Simulated_Events/Test/output_Mode1.root'
        ma.variablesToNtuple('D*0:Mode1', Ds_vars + Event,
                             filename=output_file, treename='Dstree', path=path)
    elif args.Mode==2:
        output_file = 'C01-Simulated_Events/Test/output_Mode2.root'
        ma.variablesToNtuple('D*0:Mode2', Ds_vars + Event,
                             filename=output_file, treename='Dstree', path=path)
    elif args.Mode==3:
        output_file = 'C01-Simulated_Events/Test/output_Mode3.root'
        ma.variablesToNtuple('D*0:Mode3', Ds_vars + Event,
                             filename=output_file, treename='Dstree', path=path)
    #/==========================================================================================================================/
    # -----------------
    # Monitor progress.
    # -----------------

    progress = basf2.register_module("Progress")
    path.add_module(progress)

    # ---------------
    # Process events.
    # ---------------

    # Start processing of modules.
    basf2.process(path)

    # Print basf2 call statistics.
    print(basf2.statistics)