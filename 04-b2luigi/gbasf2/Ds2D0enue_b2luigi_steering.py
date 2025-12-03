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
# D_s+* -> D_s+ gamma                                                    #
#           |                                                            #
#           +->  D_s+ -> D0 e+ nu_e                                      #
#                        |                                               #
#                        +-> K- pi+                                      #
#                                                                        #
# This script is compatible with both:                                   #
#   1. Standalone testing: basf2 Ds2D0enue_b2luigi_steering.py --        #
#   2. b2luigi integration: via create_analysis_path()                   #
#                                                                        #
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

b2.conditions.prepend_globaltag(ma.getAnalysisGlobaltag())

#/================================================================================================================================/
# D0 Modes:
#-----------------------
MODE_CFG = {
    "kmpip":    {"d0list": "D0:kmpip",    "dslist": "D_s+:Ch1",
                 "infile": "/home/belle2/amubarak/C00-Generation/ccbarDs+EventGeneration_Mode1.root"},
    "kmpippi0": {"d0list": "D0:kmpippi0", "dslist": "D_s+:Ch2",
                 "infile": "/home/belle2/amubarak/C00-Generation/ccbarDs+EventGeneration_Mode2.root"},
    "km3pi":    {"d0list": "D0:km3pi",    "dslist": "D_s+:Ch3",
                 "infile": "/home/belle2/amubarak/C00-Generation/ccbarDs+EventGeneration_Mode3.root"},
}

MODE_ALIASES = {
    "1": "kmpip",      "kmpip": "kmpip",
    "2": "kmpippi0",   "kmpippi0": "kmpippi0",
    "3": "km3pi",      "km3pi": "km3pi",
}

def parse_mode(value: str) -> str:
    key = value.strip().lower()
    if key in MODE_ALIASES:
        return MODE_ALIASES[key]
    raise argparse.ArgumentTypeError(
        "Mode must be one of 1|2|3 or kmpip|kmpippi0|km3pi"
    )
#/================================================================================================================================/
# Charge Conjugation:
#-----------------------
ChargeC = True
#/================================================================================================================================/
# Cuts
#-----------
Pion_Cut = 'pionID > 0.2 and abs(dr) < 1 and abs(dz) < 3'
Kaon_Cut = 'kaonID > 0.5 and abs(dr) < 1 and abs(dz) < 3'
Electron_Cut = 'electronID >= 0.5 and abs(dr) < 1 and abs(dz) < 3'
Gamma_Cut = 'E >= 0.1'

D0kmpip_Cut = "-0.04 <= dM <= 0.04 and useCMSFrame(p) > 2.5"
D0kmpippi0_Cut = "-0.15 <= dM <= 0.10 and useCMSFrame(p) > 2.5"
D0km3pi_Cut = "-0.04 <= dM <= 0.04 and useCMSFrame(p) > 2.5"

Ds_D0kmpip_Cut = "massDifference(0) <= 0.5 and -0.03 <= daughter(0,dM) <= 0.03"
Ds_D0kmpippi0_Cut = "massDifference(0) <= 0.5 and -0.15 <= daughter(0,dM) <= 0.10"
Ds_D0km3pi_Cut = "massDifference(0) <= 0.5 and -0.03 <= daughter(0,dM) <= 0.03"

Ds_star_Cut = "massDifference(0) <= 0.4"

# Veto
#-------
ElectronROE_Cut = 'isInRestOfEvent == 1 and abs(dr) < 20.0 and abs(dz) < 20.0'
PionROE_Cut = 'isInRestOfEvent == 1 and abs(dr) < 20.0 and abs(dz) < 20.0'
GammaROE_Cut = 'isInRestOfEvent == 1 and E > 0.1'
#/================================================================================================================================/

def argparser():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "--truth", "-t",
        action="store_true", 
        default=False,
        help="Apply truth-matching on particles.")
    parser.add_argument(
        "--mode", "-m",
        default="kmpip",
        type=parse_mode,
        metavar="{1|2|3|kmpip|kmpippi0|km3pi}",
        help="D0 decay mode: 1=kmpip, 2=kmpippi0, 3=km3pi. Default: 1 (kmpip).",
    )
    parser.add_argument("--data", action="store_true", default=False,
                        help="Running on data (disable MC-only lists/vars).")
    parser.add_argument("--infile", default="", help="Override input file path.")
    parser.add_argument(
        "--pi0-list",
        default="eff50_May2020",
        choices=[
            "eff10_May2020", "eff20_May2020", "eff30_May2020", "eff40_May2020",
            "eff50_May2020", "eff60_May2020",
            "eff10_May2020Fit", "eff20_May2020Fit", "eff30_May2020Fit",
            "eff40_May2020Fit", "eff50_May2020Fit", "eff60_May2020Fit",
        ],
        help="Standard π0 list to use when mode=kmpippi0. Default: eff50_May2020.",
    )
    parser.add_argument(
        "--outfile",
        default="",
        help="Output ROOT filename for variablesToNtuple."
    )
    return parser

def build_truth_mode(mode, path):
    #/================================================================================================================================/
    """Create mode-dependent MC lists for Ds and D0 using findMCDecay."""
    ch_tag = {"kmpip": "Ch1", "kmpippi0": "Ch2", "km3pi": "Ch3"}[mode]

    d0_decay = {
        "kmpip":    "D0:kmpip -> K- pi+",
        "kmpippi0": "D0:kmpippi0 -> K- pi+ pi0",
        "km3pi":    "D0:km3pi -> K- pi+ pi- pi+",
    }[mode]

    ds_decay = f"D_s+ -> [{d0_decay}] e+ ?nu"

    ds_mc_list = f"D_s+:MC_{ch_tag}"
    d0_mc_list = f"D0:MC_{mode}"

    # Build the MC lists directly from truth
    ma.findMCDecay(ds_mc_list, ds_decay,
                   skipNonPrimaryDaughters=False, path=path)
    ma.findMCDecay(d0_mc_list, d0_decay,
                   skipNonPrimaryDaughters=False, path=path)
    #/================================================================================================================================/
    return ds_mc_list, d0_mc_list, ch_tag

def Truth_Info_Variable(mode):
    #/================================================================================================================================/
    """Variables to save for the MC truth lists, per mode."""
    base_truth = ["isSignal","mcPDG", "genMotherPDG", "nMCDaughters", "M", "mcP", "mcE"]

    # D0 block per mode (candidate + daughters)
    if mode == "kmpip":
        D0_vars = []
        D0_vars += vu.create_aliases_for_selected(
            base_truth + ["D0Mode","Dbar0Mode"],
            decay_string="^D0:kmpip -> K- pi+",
            prefix=["D0"])
        D0_vars += vu.create_aliases_for_selected(
            ["mcPDG","M","mcP","mcE"],
            decay_string="D0:kmpip -> ^K- pi+",
            prefix=["K"])
        D0_vars += vu.create_aliases_for_selected(
            ["mcPDG","M","mcP","mcE"],
            decay_string="D0:kmpip -> K- ^pi+",
            prefix=["pi"])

    elif mode == "kmpippi0":
        D0_vars = []
        D0_vars += vu.create_aliases_for_selected(
            base_truth + ["D0Mode","Dbar0Mode",
                          "passesCut(abs(D0_mcPDG)==421 and abs(D0_D0Mode)==1036)"],
            decay_string="^D0:kmpippi0 -> K- pi+ pi0",
            prefix=["D0"])
        D0_vars += vu.create_aliases_for_selected(
            ["mcPDG","M","mcP","mcE"],
            decay_string="D0:kmpippi0 -> ^K- pi+ pi0",
            prefix=["K"])
        D0_vars += vu.create_aliases_for_selected(
            ["mcPDG","M","mcP","mcE"],
            decay_string="D0:kmpippi0 -> K- ^pi+ pi0",
            prefix=["pi"])
        D0_vars += vu.create_aliases_for_selected(
            ["mcPDG","genMotherPDG","M","mcP","mcE","daughter(0, E)","daughter(1, E)",
             "passesCut(abs(pi0_daughter_0_E)==421)"],
            decay_string="D0:kmpippi0 -> K- pi+ ^pi0",
            prefix=["pi0"])

    else:  # km3pi
        D0_vars = []
        D0_vars += vu.create_aliases_for_selected(
            base_truth + ["D0Mode","Dbar0Mode"],
            decay_string="^D0:km3pi -> K- pi+ pi- pi+",
            prefix=["D0"])
        D0_vars += vu.create_aliases_for_selected(
            ["mcPDG","M","mcP","mcE"],
            decay_string="D0:km3pi -> ^K- pi+ pi- pi+",
            prefix=["K"])
        D0_vars += vu.create_aliases_for_selected(
            ["mcPDG","M","mcP","mcE"],
            decay_string="D0:km3pi -> K- ^pi+ ^pi- ^pi+",
            prefix=["pi1","pi2","pi3"])

    # Ds block (candidate + e+)
    Ds_vars = []
    Ds_vars += vu.create_aliases_for_selected(
        base_truth + ["DsplusMode","DsminusMode"],
        decay_string="^D_s+ -> D0 e+ ?nu",
        prefix=["Ds"])
    Ds_vars += vu.create_aliases_for_selected(
        base_truth + ["D0Mode","Dbar0Mode","daughter(2, genMotherPDG)","daughter(2, passesCut(abs(pi0_genMotherPDG)==421))"],
        decay_string="D_s+ -> ^D0 e+ ?nu",
        prefix=["D0"])
    Ds_vars += vu.create_aliases_for_selected(
        ["mcPDG","M","mcP","mcE"],
        decay_string="D_s+ -> D0 ^e+ ?nu",
        prefix=["e"])
    #/================================================================================================================================/
    return Ds_vars, D0_vars

def build_selected_mode(
    mode, pi0_list, path, *,
    ChargeC,
    Pion_Cut, Kaon_Cut, Electron_Cut, Gamma_Cut,
    D0kmpip_Cut, D0kmpippi0_Cut, D0km3pi_Cut,
    Ds_D0kmpip_Cut, Ds_D0kmpippi0_Cut, Ds_D0km3pi_Cut, 
    Ds_star_Cut,
    is_data=False
):
    cfg = MODE_CFG[mode]
    d0list = cfg["d0list"]
    dslist = cfg["dslist"]
    ch_tag = dslist.split(":")[1]
    #/================================================================================================================================/
    # Create Particle Lists
    #---------------------------------

    # Gamma
    stdPhotons('cdc', beamBackgroundMVAWeight="MC15rd", fakePhotonMVAWeight="MC15rd", path=path)
    ma.cutAndCopyList("gamma:recon", "gamma:cdc", cut=Gamma_Cut, path=path)
    vm.addAlias("beamBackgroundSuppressionScore", "extraInfo(beamBackgroundSuppression)")
    vm.addAlias("fakePhotonSuppressionScore", "extraInfo(fakePhotonSuppression)")

    # Pion
    ma.fillParticleList('pi+:loose', cut=Pion_Cut, path=path)
    if not is_data:
        ma.fillParticleListFromMC('pi+:gen', '', path=path)
    if mode == "kmpippi0":
        stdPi0s(listtype=pi0_list, path=path, beamBackgroundMVAWeight="MC15rd", fakePhotonMVAWeight="MC15rd")

    # Kaon
    ma.fillParticleList('K-:loose', cut=Kaon_Cut, path=path)

    # Electron
    ma.fillParticleList("e+:uncorrected", cut='', path=path)
    if not is_data:
        ma.fillParticleListFromMC('e+:gen', '', path=path)
    
    # Bremsstrahlung Correction
    vm.addAlias("goodFWDGamma", "passesCut(clusterReg == 1 and clusterE > 0.01)")
    vm.addAlias("goodBRLGamma", "passesCut(clusterReg == 2 and clusterE > 0.01)")
    vm.addAlias("goodBWDGamma", "passesCut(clusterReg == 3 and clusterE > 0.01)")
    vm.addAlias("goodGamma", "passesCut(goodFWDGamma or goodBRLGamma or goodBWDGamma)")

    ma.fillParticleList("gamma:brems", "goodGamma", path=path)

    ma.correctBrems(outputList="e+:corrected", 
                    inputList="e+:uncorrected",
                    multiplePhotons=True,
                    gammaList="gamma:brems",
                    path=path)
    vm.addAlias("isBremsCorrected", "extraInfo(bremsCorrected)")

    ma.applyCuts('e+:corrected',Electron_Cut, path=path)
    
    # Curler
    ma.tagCurlTracks("e+:corrected", selectorType='mva', ptCut=0.6, mcTruth=True, path=path)
    vm.addAlias('isCurl', 'extraInfo(isCurl)')
    vm.addAlias('isTruthCurl', 'extraInfo(isTruthCurl)')

    Electron = ["isBremsCorrected",'isCurl','isTruthCurl']
    #/================================================================================================================================/
    # Reconstruct D0 decay
    #------------------------------------
    if mode == "kmpip":
        ma.reconstructDecay('D0:kmpip -> K-:loose pi+:loose',
                            cut=D0kmpip_Cut, chargeConjugation=ChargeC, path=path)
        vx.treeFit('D0:kmpip', conf_level=0, updateAllDaughters=True, path=path)
        ma.applyCuts("D0:kmpip", D0kmpip_Cut, path=path)

    elif mode == "kmpippi0":
        ma.reconstructDecay(f'D0:kmpippi0 -> K-:loose pi+:loose pi0:{pi0_list}',
                            cut=D0kmpippi0_Cut, chargeConjugation=ChargeC, path=path)
        vx.treeFit('D0:kmpippi0', conf_level=0, updateAllDaughters=True, path=path)
        ma.applyCuts("D0:kmpippi0", D0kmpippi0_Cut, path=path)

    else:  # km3pi
        ma.reconstructDecay('D0:km3pi -> K-:loose pi+:loose pi-:loose pi+:loose',
                            cut=D0km3pi_Cut, chargeConjugation=ChargeC, path=path)
        vx.treeFit('D0:km3pi', conf_level=0, updateAllDaughters=True, path=path)
        ma.applyCuts("D0:km3pi", D0km3pi_Cut, path=path)

    # Ds -> D0 e nu per mode
    if mode == "kmpip":
        ma.reconstructDecay('D_s+:Ch1 -> [D0:kmpip -> K-:loose pi+:loose] e+:corrected ?nu ?addbrems',
                            cut=Ds_D0kmpip_Cut, chargeConjugation=ChargeC, path=path)
        vx.treeFit("D_s+:Ch1", conf_level=0, ipConstraint=True, updateAllDaughters=True, path=path)
        ma.applyCuts("D_s+:Ch1", Ds_D0kmpip_Cut, path=path)
    elif mode == "kmpippi0":
        ma.reconstructDecay(f'D_s+:Ch2 -> [D0:kmpippi0 -> K-:loose pi+:loose pi0:{pi0_list}] e+:corrected ?nu ?addbrems',
                            cut=Ds_D0kmpippi0_Cut, chargeConjugation=ChargeC, path=path)
        vx.treeFit("D_s+:Ch2", conf_level=0, ipConstraint=True, updateAllDaughters=True, path=path)
        ma.applyCuts("D_s+:Ch2", Ds_D0kmpippi0_Cut, path=path)
    else:
        ma.reconstructDecay('D_s+:Ch3 -> [D0:km3pi -> K-:loose pi+:loose pi-:loose pi+:loose] e+:corrected ?nu ?addbrems',
                            cut=Ds_D0km3pi_Cut, chargeConjugation=ChargeC, path=path)
        vx.treeFit("D_s+:Ch3", conf_level=0, ipConstraint=True, updateAllDaughters=True, path=path)
        ma.applyCuts("D_s+:Ch3", Ds_D0km3pi_Cut, path=path)

    # Ds* -> Ds gamma per mode
    if mode == "kmpip":
        ma.reconstructDecay('D_s*+:Ch1 -> [D_s+:Ch1 -> [D0:kmpip -> K-:loose pi+:loose] e+:corrected ?nu ?addbrems] gamma:recon',
                            cut=Ds_star_Cut, chargeConjugation=ChargeC, path=path)
        ma.applyCuts("D_s*+:Ch1", Ds_star_Cut, path=path)
    elif mode == "kmpippi0":
        ma.reconstructDecay(f'D_s*+:Ch2 -> [D_s+:Ch2 -> [D0:kmpippi0 -> K-:loose pi+:loose pi0:{pi0_list}] e+:corrected ?nu ?addbrems] gamma:recon',
                            cut=Ds_star_Cut, chargeConjugation=ChargeC, path=path)
        ma.applyCuts("D_s*+:Ch2", Ds_star_Cut, path=path)
    else:
        ma.reconstructDecay('D_s*+:Ch3 -> [D_s+:Ch3 -> [D0:km3pi -> K-:loose pi+:loose pi-:loose pi+:loose] e+:corrected ?nu ?addbrems] gamma:recon',
                            cut=Ds_star_Cut, chargeConjugation=ChargeC, path=path)
        ma.applyCuts("D_s*+:Ch3", Ds_star_Cut, path=path)

    # Best-candidate ranks
    ma.rankByHighest(dslist, variable='chiProb',
                     outputVariable=f'chiProb_Ds_rank_{ch_tag}', path=path)
    vm.addAlias(f'chiProb_Ds_rank_{ch_tag}', f'extraInfo(chiProb_Ds_rank_{ch_tag})')

    ma.rankByHighest(dslist, variable='daughter(1,p)',
                     outputVariable=f'Electron_p_Ds_rank_{ch_tag}', path=path)
    vm.addAlias(f'Electron_p_Ds_rank_{ch_tag}', f'extraInfo(Electron_p_Ds_rank_{ch_tag})')

    ma.rankByHighest(dslist, variable='random',
                     outputVariable=f'random_Ds_rank_{ch_tag}', path=path)
    vm.addAlias(f'random_Ds_rank_{ch_tag}', f'extraInfo(random_Ds_rank_{ch_tag})')

    BCS = [f'chiProb_Ds_rank_{ch_tag}',
           f'Electron_p_Ds_rank_{ch_tag}',
           f'random_Ds_rank_{ch_tag}']
    #/================================================================================================================================/
    # Perform MC Matching
    if not is_data:
        ma.matchMCTruth(list_name=d0list, path=path)
        ma.matchMCTruth(list_name=dslist, path=path)
        dsstar = {"kmpip": "D_s*+:Ch1", "kmpippi0": "D_s*+:Ch2", "km3pi": "D_s*+:Ch3"}[mode]
        ma.matchMCTruth(list_name=dsstar, path=path)
    #/================================================================================================================================/
    return d0list, dslist, Electron, BCS

def Rest_of_Event_Ch(main_path, roe_path_Ch, ChargeC, ElectronROE_Cut, PionROE_Cut, GammaROE_Cut, *, is_mc):
    #/================================================================================================================================/
    # D_s+ VETO Starts Here

    # Particle Lists
    ma.fillParticleList(decayString='e-:roe', cut=ElectronROE_Cut, path=roe_path_Ch)

    ma.fillParticleList("gamma:roe", cut=GammaROE_Cut, path=roe_path_Ch)
    ma.getFakePhotonProbability('gamma:roe', weight='MC15rd', path=roe_path_Ch)
    ma.getBeamBackgroundProbability('gamma:roe', weight='MC15rd', path=roe_path_Ch)
    ma.applyCuts('gamma:roe', "beamBackgroundSuppression >= 0.5 and fakePhotonSuppression >= 0.5", path=roe_path_Ch)

    ma.fillParticleList('pi+:ROE', cut=PionROE_Cut, path=roe_path_Ch)

    ma.fillSignalSideParticleList(outputListName='e+:sig', decayString='D_s+ -> D0 ^e+:corrected ?nu', path=roe_path_Ch)
    ma.fillSignalSideParticleList(outputListName='D0:sig', decayString='D_s+ -> ^D0 e+:corrected ?nu', path=roe_path_Ch)

    # Conversion Veto
    ma.reconstructDecay(decayString='gamma:veto -> e+:sig e-:roe', 
                        cut='',
                        chargeConjugation=ChargeC, 
                        path=roe_path_Ch)
    vx.treeFit('gamma:veto', 
            conf_level=-1,
            updateAllDaughters=False,
            path=roe_path_Ch)
    ma.applyCuts('gamma:veto', 'InvM <= 3 and M<=3', path=roe_path_Ch)

    # pi0 Decay
    ma.reconstructDecay(decayString='pi0:ROE -> gamma:roe gamma:roe',
                        cut='0.080 < M < 0.200',
                        path=roe_path_Ch)
    vx.treeFit('pi0:ROE', 
               conf_level=-1,
               updateAllDaughters=False, 
               path=roe_path_Ch)
    ma.applyCuts('pi0:ROE', cut="0.080 < M < 0.200", path=roe_path_Ch)

    # D0 Veto
    ma.reconstructDecay(decayString='D*0:Mode1 -> D0:sig gamma:roe',
                        cut='',
                        chargeConjugation=ChargeC,
                        path=roe_path_Ch)
    vx.treeFit('D*0:Mode1', 
            conf_level=-1,
            updateAllDaughters=False, 
            path=roe_path_Ch)
    ma.applyCuts('D*0:Mode1', cut="0 <= massDifference(0) <= 0.25", path=roe_path_Ch)

    ma.reconstructDecay(decayString='D*0:Mode2 -> D0:sig pi0:ROE',
                        cut='',
                        chargeConjugation=ChargeC,
                        path=roe_path_Ch)
    vx.treeFit('D*0:Mode2', 
            conf_level=-1,
            updateAllDaughters=False, 
            path=roe_path_Ch)
    ma.applyCuts('D*0:Mode2', cut="0 <= massDifference(0) <= 0.25", path=roe_path_Ch)

    ma.reconstructDecay(decayString='D*+:ROE -> D0:sig pi+:ROE',
                        cut='',
                        chargeConjugation=ChargeC,
                        path=roe_path_Ch)
    vx.treeFit('D*+:ROE', 
            conf_level=-1,
            updateAllDaughters=False, 
            path=roe_path_Ch)
    ma.applyCuts('D*+:ROE', cut="0 <= massDifference(0) <= 0.5", path=roe_path_Ch)

    # MC truth only on MC
    if is_mc:
        ma.matchMCTruth('gamma:veto', path=roe_path_Ch)
        ma.matchMCTruth('D*0:Mode1', path=roe_path_Ch)
        ma.matchMCTruth('D*0:Mode2', path=roe_path_Ch)
        ma.matchMCTruth('D*+:ROE',   path=roe_path_Ch)

    # Best Candidate Selection
    vm.addAlias("gamma_abs_dM", "formula(abs(M - 0.0))")
    ma.rankByLowest(particleList='gamma:veto', 
                    variable='gamma_abs_dM',
                    numBest=1, 
                    path=roe_path_Ch)
    ma.rankByLowest(particleList='D*0:Mode1',
                    variable='abs(dM)',
                    numBest=1,
                    path=roe_path_Ch)
    ma.rankByLowest(particleList='D*0:Mode2',
                    variable='abs(dM)',
                    numBest=1,
                    path=roe_path_Ch)
    ma.rankByLowest(particleList='D*+:ROE',
                    variable='abs(dM)',
                    numBest=1,
                    path=roe_path_Ch)

    # Veto Variables
    vm.addAlias("pointangle", "formula(acos(cosAngleBetweenMomentumAndVertexVector))")
    vm.addAlias("top", "formula((daughter(0,px)*daughter(1,px)) + (daughter(0,py)*daughter(1,py)) + (daughter(0,pz)*daughter(1,pz)))")
    vm.addAlias("bottom", "formula(daughter(0,p)*daughter(1,p))")
    vm.addAlias("openangle", "acos(formula(top/bottom))")
    vm.addAlias("phi_diff", "daughterDiffOf(0, 1, phi)")
    vm.addAlias("theta_diff", "daughterDiffOf(0, 1, theta)")
    vm.addAlias("psi", "asin(formula(theta_diff/openangle))")

    # Saving Variables
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
                'daughter(1,pionID)':'gammaveto_em_pionID',
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
                                    path=roe_path_Ch)
    
    Dstar0Mode1_dict = {'M': 'Dstar0Mode1_M',
                'dM': 'Dstar0Mode1_dM',
                'massDifference(0)':'Dstar0Mode1_massDifference',
                'E': 'Dstar0Mode1_E',
                'chiProb':'Dstar0Mode1_chiProb',
                'decayAngle(0)':'Dstar0Mode1_decayAngle_0',
                'decayAngle(1)':'Dstar0Mode1_decayAngle_1',
                'cos(decayAngle(0))':'Dstar0Mode1_cos_decayAngle_0',
                'cos(decayAngle(1))':'Dstar0Mode1_cos_decayAngle_1',
                'isSignal':'Dstar0Mode1_isSignal',
                'mcPDG':'Dstar0Mode1_mcPDG',
                'mcErrors':'Dstar0Mode1_mcErrors',
                'genMotherPDG':'Dstar0Mode1_genMotherPDG',
                'nMCDaughters':'Dstar0Mode1_nMCDaughters',
                'daughter(1,M)':'Dstar0Mode1_g_M',
                'daughter(1,pt)':'Dstar0Mode1_g_pt',
                'daughter(1,p)':'Dstar0Mode1_g_p',
                'daughter(1,E)':'Dstar0Mode1_g_E',
                'daughter(1,isSignal)':'Dstar0Mode1_g_isSignal',
                'daughter(1,mcPDG)':'Dstar0Mode1_g_mcPDG',
                'daughter(1,mcErrors)':'Dstar0Mode1_g_mcErrors',
                'daughter(1,genMotherPDG)':'Dstar0Mode1_g_genMotherPDG',
                'daughter(1,genMotherPDG(1))':'Dstar0Mode1_g_genMotherPDG_1',
                'daughter(1,beamBackgroundSuppression)':'Dstar0Mode1_beamBackgroundSuppression',
                'daughter(1,fakePhotonSuppression)':'Dstar0Mode1_fakePhotonSuppression',
                }
    ma.variableToSignalSideExtraInfo(particleList='D*0:Mode1', 
                                    varToExtraInfo=Dstar0Mode1_dict, 
                                    path=roe_path_Ch)
    
    Dstar0Mode2_dict = {'M': 'Dstar0Mode2_M',
                        'dM': 'Dstar0Mode2_dM',
                        'massDifference(0)':'Dstar0Mode2_massDifference',
                        'E': 'Dstar0Mode2_E',
                        'chiProb':'Dstar0Mode2_chiProb',
                        'decayAngle(0)':'Dstar0Mode2_decayAngle_0',
                        'decayAngle(1)':'Dstar0Mode2_decayAngle_1',
                        'cos(decayAngle(0))':'Dstar0Mode2_cos_decayAngle_0',
                        'cos(decayAngle(1))':'Dstar0Mode2_cos_decayAngle_1',
                        'mcPDG':'Dstar0Mode2_mcPDG',
                        'mcErrors':'Dstar0Mode2_mcErrors',
                        'genMotherPDG':'Dstar0Mode2_genMotherPDG',
                        'nMCDaughters':'Dstar0Mode2_nMCDaughters',
                        'daughter(1,M)':'Dstar0Mode2_pi0_M',
                        'daughter(1,pt)':'Dstar0Mode2_pi0_pt',
                        'daughter(1,p)':'Dstar0Mode2_pi0_p',
                        'daughter(1,E)':'Dstar0Mode2_pi0_E',
                        'daughter(1,isSignal)':'Dstar0Mode2_pi0_isSignal',
                        'daughter(1,mcPDG)':'Dstar0Mode2_pi0_mcPDG',
                        'daughter(1,mcErrors)':'Dstar0Mode2_pi0_mcErrors',
                        'daughter(1,genMotherPDG)':'Dstar0Mode2_pi0_genMotherPDG',
                        'daughter(1,genMotherPDG(1))':'Dstar0Mode2_pi0_genMotherPDG_1',
                        }
    ma.variableToSignalSideExtraInfo(particleList='D*0:Mode2', 
                                    varToExtraInfo=Dstar0Mode2_dict, 
                                    path=roe_path_Ch)
    
    Dstarplus_dict = {'M': 'Dstarplus_M',
                    'dM': 'Dstarplus_dM',
                    'massDifference(0)':'Dstarplus_massDifference',
                    'E': 'Dstarplus_E',
                    'chiProb':'Dstarplus_chiProb',
                    'decayAngle(0)':'Dstarplus_decayAngle_0',
                    'decayAngle(1)':'Dstarplus_decayAngle_1',
                    'cos(decayAngle(0))':'Dstarplus_cos_decayAngle_0',
                    'cos(decayAngle(1))':'Dstarplus_cos_decayAngle_1',
                    'mcPDG':'Dstarplus_mcPDG',
                    'mcErrors':'Dstarplus_mcErrors',
                    'genMotherPDG':'Dstarplus_genMotherPDG',
                    'nMCDaughters':'Dstarplus_nMCDaughters',
                    'daughter(1,M)':'Dstarplus_pi_M',
                    'daughter(1,pt)':'Dstarplus_pi_pt',
                    'daughter(1,p)':'Dstarplus_pi_p',
                    'daughter(1,E)':'Dstarplus_pi_E',
                    'daughter(1,electronID)':'Dstarplus_pi_electronID',
                    'daughter(1,pionID)':'Dstarplus_pi_pionID',
                    'daughter(1,isSignal)':'Dstarplus_pi_isSignal',
                    'daughter(1,mcPDG)':'Dstarplus_pi_mcPDG',
                    'daughter(1,mcErrors)':'Dstarplus_pi_mcErrors',
                    'daughter(1,genMotherPDG)':'Dstarplus_pi_genMotherPDG',
                    'daughter(1,genMotherPDG(1))':'Dstarplus_pi_genMotherPDG_1',
                    }
    ma.variableToSignalSideExtraInfo(particleList='D*+:ROE', 
                                    varToExtraInfo=Dstarplus_dict, 
                                    path=roe_path_Ch)

    # execute roe_path for each RestOfEvent in the event
    main_path.for_each('RestOfEvent', 'RestOfEvents', roe_path_Ch)

    # Variable List
    gamma_ROE_Ch, Dstar0Mode1_ROE_Ch, Dstar0Mode2_ROE_Ch, Dstarplus_ROE_Ch = [], [], [], []
    for key in gamma_dict:
        vm.addAlias(gamma_dict[key], f"extraInfo({gamma_dict[key]})")
        gamma_ROE_Ch.append(gamma_dict[key])
    vm.addAlias('gammaveto_M_Correction', 'ifNANgiveX(gammaveto_M,10)')
    gamma_ROE_Ch.append('gammaveto_M_Correction')

    for key in Dstar0Mode1_dict:
        vm.addAlias(Dstar0Mode1_dict[key], f"extraInfo({Dstar0Mode1_dict[key]})")
        Dstar0Mode1_ROE_Ch.append(Dstar0Mode1_dict[key])
    vm.addAlias('Dstar0Mode1_M_Correction', 'ifNANgiveX(Dstar0Mode1_M,10)')
    Dstar0Mode1_ROE_Ch.append('Dstar0Mode1_M_Correction')

    for key in Dstar0Mode2_dict:
        vm.addAlias(Dstar0Mode2_dict[key], f"extraInfo({Dstar0Mode2_dict[key]})")
        Dstar0Mode2_ROE_Ch.append(Dstar0Mode2_dict[key])
    vm.addAlias('Dstar0Mode2_M_Correction', 'ifNANgiveX(Dstar0Mode2_M,10)')
    Dstar0Mode2_ROE_Ch.append('Dstar0Mode2_M_Correction')

    for key in Dstarplus_dict:
        vm.addAlias(Dstarplus_dict[key], f"extraInfo({Dstarplus_dict[key]})")
        Dstarplus_ROE_Ch.append(Dstarplus_dict[key])
    vm.addAlias('Dstarplus_massDifference_Correction', 'ifNANgiveX(Dstarplus_massDifference,10)')
    Dstarplus_ROE_Ch.append('Dstarplus_massDifference_Correction')
    vm.addAlias('Dstarplus_M_Correction', 'ifNANgiveX(Dstarplus_M,10)')
    Dstarplus_ROE_Ch.append('Dstarplus_M_Correction')
    #/================================================================================================================================/
    return gamma_ROE_Ch, Dstar0Mode1_ROE_Ch, Dstar0Mode2_ROE_Ch, Dstarplus_ROE_Ch

def build_roe_for_selected(main_path, dslist, ChargeC, ElectronROE_Cut, PionROE_Cut, GammaROE_Cut, *, is_mc):
    #/================================================================================================================================/
    # ROE
    ma.buildRestOfEvent(target_list_name=dslist, fillWithMostLikely=False, path=main_path)
    
    track_based_cuts = "abs(d0) < 20.0 and abs(z0) < 20.0"
    ecl_based_cuts = ""
    roe_mask = ("roe_mask", track_based_cuts, ecl_based_cuts)
    ma.appendROEMasks(dslist, [roe_mask], path=main_path)

    ma.updateROEUsingV0Lists(target_particle_list=dslist,
                             mask_names="roe_mask",
                             default_cleanup=True,
                             selection_cuts=None,
                             apply_mass_fit=True,
                             fitter='treefit',
                             path=main_path)

    roe_path = b2.create_path()
    dead_end = b2.create_path()
    ma.signalSideParticleFilter(particleList=dslist, selection='', roe_path=roe_path, deadEndPath=dead_end)
    #/================================================================================================================================/
    return Rest_of_Event_Ch(main_path, roe_path, ChargeC, ElectronROE_Cut, PionROE_Cut, GammaROE_Cut, is_mc=is_mc)

def Suggestion(path, mode, pi0_list):
    #/================================================================================================================================/
    """Attach Ds* -> Ds gamma ΔM"""
    ds_by_mode = {"kmpip": "D_s+:Ch1", "kmpippi0": "D_s+:Ch2", "km3pi": "D_s+:Ch3"}
    dsstar_by_mode = {"kmpip": "D_s*+:Ch1", "kmpippi0": "D_s*+:Ch2", "km3pi": "D_s*+:Ch3"}
    d0_decay_by_mode = {
        "kmpip":    "[D0:kmpip -> K-:loose pi+:loose]",
        "kmpippi0": f"[D0:kmpippi0 -> K-:loose pi+:loose pi0:{pi0_list}]",
        "km3pi":    "[D0:km3pi -> K-:loose pi+:loose pi-:loose pi+:loose]",
    }

    ds     = ds_by_mode[mode]
    dsstar = dsstar_by_mode[mode]
    ch_tag = ds.split(":")[1]
    d0_blk = d0_decay_by_mode[mode]

    decay_string = f"{dsstar} -> [^{ds} -> {d0_blk} e+:corrected ?nu ?addbrems] gamma:recon"
    ma.variablesToDaughterExtraInfo(
        particleList=dsstar,
        decayString=decay_string,
        variables={"massDifference(0)": f"Ds_starminusDs_{ch_tag}"},
        option=0,
        path=path,
    )

    vm.addAlias(f"Ds_starminusDs_{ch_tag}", f"extraInfo(Ds_starminusDs_{ch_tag})")
    vm.addAlias(f"Ds_starminusDs_M_Correction_{ch_tag}", f"ifNANgiveX(Ds_starminusDs_{ch_tag},10)")
    vm.addAlias(
        f"goodDsplus_{ch_tag}",
        f"passesCut(Ds_starminusDs_M_Correction_{ch_tag} >= 0.12 and "
        f"Ds_starminusDs_M_Correction_{ch_tag} <= 0.165)"
    )

    vm.addAlias("Ds_starminusDs", f"extraInfo(Ds_starminusDs_{ch_tag})")
    vm.addAlias("Ds_starminusDs_M_Correction", "ifNANgiveX(Ds_starminusDs,10)")
    vm.addAlias("goodDsplus", "passesCut(Ds_starminusDs_M_Correction >= 0.12 and Ds_starminusDs_M_Correction <= 0.165)")
    
    expertVars = []
    #/================================================================================================================================/
    return ch_tag, expertVars

def Reconstruction_Variable(
    ch_tag,
    mode, pi0_list,
    Electron, Before_VF, BCS, expertVars,
    roe_vars,
    *, is_data
):
    gamma_ROE, Dstar0Mode1_ROE, Dstar0Mode2_ROE, Dstarplus_ROE = roe_vars
    ROE_flat = list(gamma_ROE) + list(Dstar0Mode1_ROE) + list(Dstar0Mode2_ROE) + list(Dstarplus_ROE)
    #/================================================================================================================================/
    # Save Variables

    truth_D0 = ["mcSecPhysProc","mcPDG","genMotherPDG","genMotherID","mcErrors"]
    if is_data:
        truth_D0 = []

    D0_kmpip_vars, D0_kmpippi0_vars, D0_km3pi_vars = [], [], []

    if mode == "kmpip":
        D0_kmpip_vars += vu.create_aliases_for_selected(
            ["M","dM","useCMSFrame(E)","useCMSFrame(p)","chiProb","cos(theta)","daughterAngle(0, 1)",
             'daughterDiffOf(0,1,theta)','daughterDiffOf(0,1,cos(theta))','daughterDiffOf(0,1,phi)',
             'daughterMotherDiffOf(0,x)','daughterMotherDiffOf(0,y)','daughterMotherDiffOf(0,z)',
             "decayAngle(0)","cos(decayAngle(0))","decayAngle(1)","cos(decayAngle(1))",
             "significanceOfDistance","flightDistance","isSignal",'ifNANgiveX(isSignal,0)',
             'D0Mode','Dbar0Mode'] + truth_D0,
            decay_string='^D0:kmpip -> K-:loose pi+:loose', 
            prefix=['D0_kmpip'])
        D0_kmpip_vars += vu.create_aliases_for_selected(
            ["M",'kaonID','nSVDHits','nCDCHits','abs(dr)','abs(dz)','z0','d0','pValue','chi2','ndf',
             'theta','phi','omega'] + truth_D0,
            decay_string='D0:kmpip -> ^K-:loose pi+:loose', 
            prefix=['K_kmpip'])
        D0_kmpip_vars += vu.create_aliases_for_selected(
            ["M",'pionID','nSVDHits','nCDCHits','abs(dr)','abs(dz)','z0','d0','pValue','chi2','ndf',
             'theta','phi','omega'] + truth_D0,
            decay_string='D0:kmpip -> K-:loose ^pi+:loose', 
            prefix=['pi_kmpip'])

    if mode == "kmpippi0":
        D0_kmpippi0_vars += vc.dalitz_3body
        D0_kmpippi0_vars += vu.create_aliases_for_selected(
            ["M","dM","useCMSFrame(E)","useCMSFrame(p)","cosAngleBetweenMomentumAndVertexVector",
             "chiProb","cos(theta)","decayAngle(0)","cos(decayAngle(0))","decayAngle(1)","cos(decayAngle(1))",
             "flightDistance","isSignal",'ifNANgiveX(isSignal,0)',"nMCDaughters",
             "passesCut(abs(D0_kmpippi0_nMCDaughters)==2)","passesCut(abs(D0_kmpippi0_nMCDaughters)==3)",
             "passesCut(abs(D0_kmpippi0_mcPDG)==421)",
             'D0Mode','Dbar0Mode'] + truth_D0,
            decay_string=f'^D0:kmpippi0 -> K-:loose pi+:loose pi0:{pi0_list}', 
            prefix=['D0_kmpippi0'])
        D0_kmpippi0_vars += vu.create_aliases_for_selected(
            ["M",'kaonID','nSVDHits','nCDCHits','abs(dr)','abs(dz)','z0','d0','pValue','chi2','ndf',
             'theta','phi','omega'] + truth_D0,
            decay_string=f'D0:kmpippi0 -> ^K-:loose pi+ pi0:{pi0_list}', 
            prefix=['K_kmpippi0'])
        D0_kmpippi0_vars += vu.create_aliases_for_selected(
            ["M",'pionID','nSVDHits','nCDCHits','abs(dr)','abs(dz)','z0','d0','pValue','chi2','ndf',
             'theta','phi','omega'] + truth_D0,
            decay_string=f'D0:kmpippi0 -> K- ^pi+ pi0:{pi0_list}', 
            prefix=['pi_kmpippi0'])
        vm.addAlias('E1', 'daughter(0, E)'); vm.addAlias('E2', 'daughter(1, E)')
        vm.addAlias('energyAsymmetry', 'formula(abs(E1 - E2) / (E1 + E2))')
        D0_kmpippi0_vars += vu.create_aliases_for_selected(
            ["M","InvM","E","E1","E2","max(E1, E2)","energyAsymmetry",
             "daughter(0, isSignal)","daughter(0, mcPDG)","daughter(0, genMotherPDG)",
             "daughter(1, isSignal)","daughter(1, mcPDG)","daughter(1, genMotherPDG)",
             "passesCut(pi0_kmpippi0_mcPDG==111)","isSignal"] + truth_D0,
            decay_string=f'D0:kmpippi0 -> K- pi+ ^pi0:{pi0_list}', 
            prefix=['pi0_kmpippi0'])

    if mode == "km3pi":
        D0_km3pi_vars += vu.create_aliases_for_selected(
            ["M","dM","useCMSFrame(E)","useCMSFrame(p)","chiProb","cos(theta)","daughterAngle(0, 1)",
             'daughterDiffOf(0,1,theta)','daughterDiffOf(0,1,cos(theta))','daughterDiffOf(0,1,phi)',
             'daughterMotherDiffOf(0,x)','daughterMotherDiffOf(0,y)','daughterMotherDiffOf(0,z)',
             "decayAngle(0)","cos(decayAngle(0))","decayAngle(1)","cos(decayAngle(1))",
             "significanceOfDistance","flightDistance","isSignal",'ifNANgiveX(isSignal,0)',
             'D0Mode','Dbar0Mode'] + truth_D0,
            decay_string='^D0:km3pi -> K-:loose pi+:loose pi-:loose pi+:loose', 
            prefix=['D0_km3pi'])
        D0_km3pi_vars += vu.create_aliases_for_selected(
            ["M",'kaonID','nSVDHits','nCDCHits','abs(dr)','abs(dz)','z0','d0','pValue','chi2','ndf',
             'theta','phi','omega'] + truth_D0,
            decay_string='D0:km3pi -> ^K-:loose pi+:loose pi-:loose pi+:loose', 
            prefix=['K_km3pi'])
        D0_km3pi_vars += vu.create_aliases_for_selected(
            ["M",'pionID','nSVDHits','nCDCHits','abs(dr)','abs(dz)','z0','d0','pValue','chi2','ndf',
             'theta','phi','omega'] + truth_D0,
            decay_string='D0:km3pi -> K-:loose ^pi+:loose ^pi-:loose ^pi+:loose', 
            prefix=['pi1_km3pi','pi2_km3pi','pi3_km3pi'])

    # D_s+ Variables
    tracks = ['isCloneTrack',
            'dr','dz','abs(dr)','abs(dz)','z0','d0','pValue',
            'firstCDCLayer','firstPXDLayer','firstSVDLayer',
            'nPXDHits','nVXDHits','nSVDHits','nCDCHits',
            'seenInCDC','seenInPXD','seenInSVD','seenInTOP',
            'inARICHAcceptance','inCDCAcceptance','inTOPAcceptance']

    Ds_vars=[]
    var_1 = ['ImpactXY','cos(theta)','phi','mcP','M','pt','E','p','px','py','pz','abs(pz)','isOrHasCloneTrack',"charge"]
    truth = ["mcP","mcE","mcSecPhysProc",'nMCMatches','mcPrimary','isSignal','ifNANgiveX(isSignal,5)','mcPDG','genMotherPDG',
             'genMotherID','mcErrors','mcMatchWeight']

    vm.addAlias("pminusptrue_p", "formula(p - mcP)")
    vm.addAlias("pminusptrue_old", "formula(daughter(0,p) - mcP)")
    vm.addAlias("Motherptrueminusp_p", "formula(mcMother(p) - p)")
    vm.addAlias("Motherptrueminusp_E", "formula(mcMother(E) - E)")
    
    mc_assoc_e = [
        'angleToClosestInList(e+:gen)',
        'closestInList(e+:gen, mcPDG)','closestInList(e+:gen, genMotherPDG)',
        'closestInList(e+:gen, pt)','closestInList(e+:gen, px)','closestInList(e+:gen, py)','closestInList(e+:gen, pz)',
        'closestInList(pi+:gen, mcPDG)','closestInList(pi+:gen, genMotherPDG)',
        'closestInList(pi+:gen, pt)','closestInList(pi+:gen, px)','closestInList(pi+:gen, py)','closestInList(pi+:gen, pz)',
    ]
    e_base = Electron + [
                        'clusterE1E9','clusterE9E21',
                        "trackFitCovariance(0, 2)","trackFitCovariance(2, 3)","trackFitCovariance(0, 3)",
                        "trackFitHypothesisPDG",
                        'chi2','ndf','trackTime',
                        'pionID','electronID','binaryPID(11,211)',
                        'omega',
                        'daughter(0,p)',
                        'daughter(0,E)',
                        'daughter(0,electronID)',
                        'daughter(0,isCloneTrack)',
                        "flightTime","formula(E/p)",
                        'mcVirtual',
                        'isMisidentified',
                        "pminusptrue_p","pminusptrue_old",
                        'mcMother(E)','mcMother(p)','mcMother(mcDecayVertexFromIPDistance)','mcMother(nMCDaughters)',
                        'mcMother(mcDaughter(0, p))','mcMother(mcDaughter(1, p))',
                        'mcMother(mcDaughter(0, E))','mcMother(mcDaughter(1, E))',
                        'mcMother(mcDaughter(0, PDG))','mcMother(mcDaughter(1, PDG))','mcMother(mcDaughter(2, PDG))',
                        "Motherptrueminusp_p","Motherptrueminusp_E",
                        'genMotherPDG(1)','genMotherPDG(2)','genMotherID(1)','genMotherID(2)'
                        ] +  var_1 + tracks + truth
    if is_data:
        e_use = e_base
    else:
        e_use = e_base + mc_assoc_e
    e_vars = vu.create_aliases_for_selected(
        list_of_variables= e_use,
        decay_string='D_s+ -> D0 ^e+:corrected ?nu',
        prefix=['e'])

    km_vars = []
    pip_vars = []
    pi0_vars = []

    if mode == "kmpip":
        km_vars  = vu.create_aliases_for_selected(
            ['kaonID','nSVDHits','nCDCHits','abs(dr)','abs(dz)','z0','d0','pValue','chi2','ndf','theta','phi','omega']
            + var_1 + tracks + truth,
            decay_string='D_s+ -> [D0:kmpip -> ^K-:loose pi+:loose] e+:corrected ?nu', 
            prefix=['K_Ch1'])
        pip_vars = vu.create_aliases_for_selected(
            ['pionID','nSVDHits','nCDCHits','abs(dr)','abs(dz)','z0','d0','pValue','chi2','ndf','theta','phi','omega']
            + var_1 + tracks + truth,
            'D_s+ -> [D0:kmpip -> K-:loose ^pi+:loose] e+:corrected ?nu', 
            prefix=['pi_Ch1'])

    if mode == "kmpippi0":
        km_vars  = vu.create_aliases_for_selected(
            ['kaonID','nSVDHits','nCDCHits','abs(dr)','abs(dz)','z0','d0','pValue','chi2','ndf','theta','phi','omega']
            + var_1 + tracks + truth,
            decay_string=f'D_s+ -> [D0:kmpippi0 -> ^K-:loose pi+:loose pi0:{pi0_list}] e+:corrected ?nu', 
            prefix=['K_Ch2'])
        pip_vars = vu.create_aliases_for_selected(
            ['pionID','nSVDHits','nCDCHits','abs(dr)','abs(dz)','z0','d0','pValue','chi2','ndf','theta','phi','omega']
            + var_1 + tracks + truth,
            decay_string=f'D_s+ -> [D0:kmpippi0 -> K-:loose ^pi+:loose pi0:{pi0_list}] e+:corrected ?nu', 
            prefix=['pi_Ch2'])
        pi0_vars = vu.create_aliases_for_selected(
            ["M","InvM","E","E1","E2","max(E1, E2)","energyAsymmetry"] + var_1 + truth,
            decay_string=f'D_s+ -> [D0:kmpippi0 -> K-:loose pi+:loose ^pi0:{pi0_list}] e+:corrected ?nu', 
            prefix=['pi0_Ch2'])

    if mode == "km3pi":
        km_vars  = vu.create_aliases_for_selected(
            ['kaonID','nSVDHits','nCDCHits','abs(dr)','abs(dz)','z0','d0','pValue','chi2','ndf','theta','phi','omega']
            + var_1 + tracks + truth,
            decay_string='D_s+ -> [D0:km3pi -> ^K-:loose pi+:loose pi-:loose pi+:loose] e+:corrected ?nu', 
            prefix=['K_Ch3'])
        pip_vars = vu.create_aliases_for_selected(
            ['pionID','nSVDHits','nCDCHits','abs(dr)','abs(dz)','z0','d0','pValue','chi2','ndf','theta','phi','omega']
            + var_1 + tracks + truth,
            decay_string='D_s+ -> [D0:km3pi -> K-:loose ^pi+:loose ^pi-:loose ^pi+:loose] e+:corrected ?nu',
            prefix=['pi1_Ch3','pi2_Ch3','pi3_Ch3'])

    vm.addAlias('p1', 'daughter(0, p)')
    vm.addAlias('p2', 'daughter(1, p)')
    vm.addAlias('MomentumAsymmetry', 'formula((p1 - p2) / (p1 + p2))')
    vm.addAlias("charged_product", "formula(daughter(0,charge)*daughter(1,charge))")
    vm.addAlias("D0orD0bar", "passesCut(daughter(0,charge)==-1 and daughter(1,charge)==1)")
    D0_vars = vu.create_aliases_for_selected(
        list_of_variables=vc.dalitz_3body +
                          ["charged_product","MomentumAsymmetry",
                           "useCMSFrame(E)","useCMSFrame(p)",
                           "dM","useAlternativeDaughterHypothesis(M, 1:K+)",
                           "chiProb","cos(theta)","daughterAngle(0, 1)",
                           'daughterDiffOf(0,1,theta)','daughterDiffOf(0,1,cos(theta))','daughterDiffOf(0,1,phi)',
                           'daughterMotherDiffOf(0,theta)','daughterMotherDiffOf(0,cos(theta))','daughterMotherDiffOf(0,phi)',
                           "decayAngle(0)","cos(decayAngle(0))",
                           "decayAngle(1)","cos(decayAngle(1))",
                           "mcDecayTime","mcLifeTime","mcFlightTime",
                           "significanceOfDistance",
                           "flightDistance","flightDistanceErr","flightTime","flightTimeErr",
                           "vertexDistance","vertexDistanceErr",
                           "useRestFrame(daughterAngle(0, 1))",
                           "formula(daughter(0,dr) - daughter(1,dr))","formula(daughter(0,dz) - daughter(1,dz))"]
                           + var_1 
                           + [
                              'mcMother(nMCDaughters)',
                              'mcMother(mcDaughter(1, PDG))','mcMother(mcDaughter(1, nMCDaughters))',
                              'mcMother(mcDaughter(1, pt))','mcMother(mcDaughter(1, pz))','mcMother(mcDaughter(1, cos(theta)))',
                              'mcMother(mcDaughter(1, mcDaughter(0,nMCDaughters)))','mcMother(mcDaughter(1, mcDaughter(1,nMCDaughters)))'
                             ]
                           + [
                            'isSignalAcceptBremsPhotons','isSignalAcceptMissing','isSignalAcceptMissingGamma',
                            'isSignalAcceptMissingMassive','isSignalAcceptMissingNeutrino','isSignalAcceptWrongFSPs',
                            "passesCut(abs(D0_mcPDG)==421)"
                           ]
                           + ['genMotherPDG(1)','genMotherPDG(2)','genMotherID(1)','genMotherID(2)'] 
                           + truth + ['D0Mode','Dbar0Mode','nMCDaughters','D0orD0bar'],
        decay_string='D_s+ -> ^D0 e+:corrected ?nu',
        prefix=['D0'])

    vm.addAlias("D0_Dstarplus", "passesCut(abs(D0_mcPDG)==421 and abs(D0_genMotherPDG)==413)")
    vm.addAlias("D0_Dstar0", "passesCut(abs(D0_mcPDG)==421 and abs(D0_genMotherPDG)==423)")
    vm.addAlias("D0_NoDstarplusDstar0", "passesCut(abs(D0_mcPDG)==421 and abs(D0_genMotherPDG)!=413 and abs(D0_genMotherPDG)!=423)")
    vm.addAlias("D0_Other", "passesCut(abs(D0_mcPDG)!=421)")

    vm.addAlias("Dstarplus", "passesCut(abs(Ds_ifNANgiveX_mcPDG_9999)==413)")
    vm.addAlias("NoDstarplus", "passesCut(abs(Ds_ifNANgiveX_mcPDG_9999)!=413)")
    vm.addAlias("Dstar0", "passesCut(abs(Ds_ifNANgiveX_mcPDG_9999)==423)")
    vm.addAlias("Failed", "passesCut(abs(Ds_ifNANgiveX_mcPDG_9999)==9999)")
    vm.addAlias("Comb", "passesCut(abs(Ds_ifNANgiveX_mcPDG_9999)!=413 and abs(Ds_ifNANgiveX_mcPDG_9999)!=423 and abs(Ds_ifNANgiveX_mcPDG_9999)!=9999)")
    vm.addAlias("Signal", "passesCut(abs(Ds_ifNANgiveX_mcPDG_9999)==431)")

    vm.addAlias("L_diff", "formula(((x - daughter(0,x))**2 + (y - daughter(0,y))**2 + (z - daughter(0,z))**2)**(1/2))")

    vm.addAlias("diff_D0e_noBrems", "formula(daughterCombination(M,0,1:0) - daughter(0, M))")
    vm.addAlias("diff_D0pi_noBrems", "formula(M_D0pi_noBrems - daughter(0, M))")
    Ds_vars += vu.create_aliases_for_selected(
        list_of_variables=[
                        "M_uncorrected",'M_pi','massDifference(0)',"diff_D0pi","diff_D0K",
                        "diff_D0e_noBrems","diff_D0pi_noBrems"]
                        + Before_VF + [
                        "formula(diff_D0pi - massDifference(0))",
                        'InvMLambda',
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
                        + [f"Ds_starminusDs_{ch_tag}",
                           f"Ds_starminusDs_M_Correction_{ch_tag}",
                           f"goodDsplus_{ch_tag}"] 
                        + var_1 + truth 
                        + ['genNStepsToDaughter(0)','genNStepsToDaughter(1)',
                        'genNMissingDaughter(11)','genNMissingDaughter(22)',
                        'isSignalAcceptBremsPhotons','isSignalAcceptMissing','isSignalAcceptMissingGamma',
                        'isSignalAcceptMissingMassive','isSignalAcceptMissingNeutrino','isSignalAcceptWrongFSPs',
                        'ifNANgiveX(mcPDG,9999)',
                        'nMCDaughters','genParticleID'] 
                        + ["Comb","Failed","Dstar0","Dstarplus","NoDstarplus","Signal"] 
                        + ["D0_Dstarplus","D0_Dstar0","D0_NoDstarplusDstar0","D0_Other"]
                        + ["mcDaughter(0, PDG)"]
                        + expertVars + BCS,
        decay_string='^D_s+ -> D0 e+:corrected ?nu',
        prefix=['Ds'])

    Ds_star_vars=[]

    # Extra Variables
    vm.addAlias("diff_D0e", "massDifference(0)")
    vm.addAlias("diff_D0pi", "useAlternativeDaughterHypothesis(massDifference(0), 1:pi+)")
    vm.addAlias("diff_D0K", "useAlternativeDaughterHypothesis(massDifference(0), 1:K+)")
    vm.addAlias("Angle_D0e", "daughterAngle(0, 1)")
    vm.addAlias("Angle_Ke", "daughterAngle(0:0, 1)")
    vm.addAlias("Angle_pie", "daughterAngle(0:1, 1)")
    vm.addAlias("M_uncorrected", "daughterCombination(M,0,1:0)")
    vm.addAlias("M_pi", "useAlternativeDaughterHypothesis(M, 1:pi+)")
    vm.addAlias("mcM", "formula(((mcE)**2 - (mcP)**2)**(1/2))")

    vm.addAlias("mcP_D0e_mag", "formula((daughter(0,mcP))**2 + (daughter(1,mcP))**2)")
    vm.addAlias("mcP_D0e_comp", "formula((daughter(0,mcPX)*daughter(1,mcPX)) + (daughter(0,mcPY)*daughter(1,mcPY)) + (daughter(0,mcPZ)*daughter(1,mcPZ)))")
    vm.addAlias("mcE_D0e_emass", "formula(daughter(0,mcE) + daughter(1,mcE))")
    vm.addAlias("mcM_D0e_emass", "formula(((mcE_D0e_emass)**2 - mcP_D0e_mag - (2*mcP_D0e_comp))**(1/2))")
    vm.addAlias("MminusMtrue_D0e_emass", "formula(M - mcM_D0e_emass)")
    vm.addAlias("mcE_e", "formula(((0.13957039)**2 + (daughter(1,mcP))**2)**(1/2))")
    vm.addAlias("mcE_D0e_pimass", "formula(daughter(0,mcE) + mcE_e)")
    vm.addAlias("mcM_D0e_pimass", "formula(((mcE_D0e_pimass)**2 - mcP_D0e_mag - (2*mcP_D0e_comp))**(1/2))")
    vm.addAlias("MminusMtrue_D0e_pimass", "formula(M_pi - mcM_D0e_pimass)")
    vm.addAlias("e_E_noBrems", 'formula(((0.13957039)**2 + (daughter(1,daughter(0,p)))**2)**(1/2))')
    vm.addAlias("Ds_E_noBrems", "formula(daughter(0,E) + e_E_noBrems)")
    vm.addAlias("Ds_P_mag_noBrems", "formula((daughter(0,p))**2 + (daughter(1,daughter(0,p)))**2)")
    vm.addAlias("Ds_P_comp_noBrems", "formula((daughter(0,px)*daughter(1,daughter(0,px))) + (daughter(0,py)*daughter(1,daughter(0,py))) + (daughter(0,pz)*daughter(1,daughter(0,pz))))")
    vm.addAlias("M_D0pi_noBrems", "formula(((Ds_E_noBrems)**2 - Ds_P_mag_noBrems - (2*Ds_P_comp_noBrems))**(1/2))")

    if mode == "kmpip":
        D0_mode_vars = D0_kmpip_vars
    elif mode == "kmpippi0":
        D0_mode_vars = D0_kmpippi0_vars
    else:
        D0_mode_vars = D0_km3pi_vars

    Ds_vars += ROE_flat

    Event = ['IPX','IPY','IPZ']
    #/==========================================================================================================================/
    return D0_mode_vars, e_vars, km_vars, pip_vars, pi0_vars, D0_vars, Ds_vars, Event

#/================================================================================================================================/
# b2luigi ENTRY POINT
#/================================================================================================================================/

def create_analysis_path(mode="kmpip", pi0_list="eff50_May2020", output_filename="output.root"):
    """
    b2luigi entry point for gbasf2 integration.
    
    Args:
        mode: D0 decay mode ("kmpip", "kmpippi0", or "km3pi")
        pi0_list: Standard pi0 list for kmpippi0 mode
        output_filename: Output ROOT file name (relative path for gbasf2)
    
    Returns:
        basf2.Path: Complete analysis path
    """
    # Create path
    path = b2.Path()  # ← FIX THIS LINE (was b2.create_path())
    
    # Input - gbasf2 overrides this on the grid
    cfg = MODE_CFG[mode]
    ma.inputMdst(environmentType="default", filename=cfg["infile"], path=path)
    
    # Run reconstruction (same as standalone)
    d0list, dslist, Electron, BCS = build_selected_mode(
        mode, pi0_list, path,
        ChargeC=ChargeC,
        Pion_Cut=Pion_Cut, Kaon_Cut=Kaon_Cut, Electron_Cut=Electron_Cut, Gamma_Cut=Gamma_Cut,
        D0kmpip_Cut=D0kmpip_Cut, D0kmpippi0_Cut=D0kmpippi0_Cut, D0km3pi_Cut=D0km3pi_Cut,
        Ds_D0kmpip_Cut=Ds_D0kmpip_Cut, Ds_D0kmpippi0_Cut=Ds_D0kmpippi0_Cut, Ds_D0km3pi_Cut=Ds_D0km3pi_Cut,
        Ds_star_Cut=Ds_star_Cut,
        is_data=False
    )
    
    gamma_ROE, Dst0M1_ROE, Dst0M2_ROE, Dstp_ROE = build_roe_for_selected(
        main_path=path,
        dslist=dslist,
        ChargeC=ChargeC,
        ElectronROE_Cut=ElectronROE_Cut,
        PionROE_Cut=PionROE_Cut,
        GammaROE_Cut=GammaROE_Cut,
        is_mc=True
    )
    
    ch_tag, expertVars = Suggestion(path, mode, pi0_list)
    
    (D0_mode_vars, e_vars, km_vars, pip_vars, pi0_vars,
     D0_vars, Ds_vars, Event) = Reconstruction_Variable(
        ch_tag=ch_tag,
        mode=mode,
        pi0_list=pi0_list,
        Electron=Electron, Before_VF=[], BCS=BCS,
        expertVars=expertVars,
        roe_vars=(gamma_ROE, Dst0M1_ROE, Dst0M2_ROE, Dstp_ROE),
        is_data=False
    )
    
    # Save
    ds_list, ds_tree = {
        "kmpip":    ("D_s+:Ch1", "DstreeCh1"),
        "kmpippi0": ("D_s+:Ch2", "DstreeCh2"),
        "km3pi":    ("D_s+:Ch3", "DstreeCh3"),
    }[mode]
    
    ds_vars_for_save = e_vars + km_vars + pip_vars + D0_vars + Ds_vars + Event
    if mode == "kmpippi0":
        ds_vars_for_save = e_vars + km_vars + pip_vars + pi0_vars + D0_vars + Ds_vars + Event
    
    ma.variablesToNtuple(ds_list, ds_vars_for_save,
                        filename=output_filename, treename=ds_tree, path=path)
    
    if mode == "kmpip":
        ma.variablesToNtuple('D0:kmpip', D0_mode_vars,
                            filename=output_filename, treename='D02kmpiptree', path=path)
    elif mode == "kmpippi0":
        ma.variablesToNtuple('D0:kmpippi0', D0_mode_vars,
                            filename=output_filename, treename='D02kmpippi0tree', path=path)
    else:
        ma.variablesToNtuple('D0:km3pi', D0_mode_vars,
                            filename=output_filename, treename='D02km3pitree', path=path)
    
    return path

#/================================================================================================================================/
# STANDALONE TESTING
#/================================================================================================================================/

if __name__ == '__main__':

    args = argparser().parse_args()

    path = b2.create_path()
    
    # Load input
    cfg = MODE_CFG[args.mode]
    infile = args.infile or cfg["infile"]
    
    if args.data and not args.infile:
        raise ValueError("Data mode requires --infile. No default data path is set.")
    
    import os
    if not os.path.exists(infile):
        raise FileNotFoundError(f"Input not found: {infile}")
    
    ma.inputMdst(environmentType="default", filename=infile, path=path)
    
    # Reconstruction
    if args.truth:
        ds_mc_list, d0_mc_list, ch_tag = build_truth_mode(args.mode, path)
        Ds_mc_vars, D0_mc_vars = Truth_Info_Variable(args.mode)
    else:
        d0list, dslist, Electron, BCS = build_selected_mode(
            args.mode, args.pi0_list, path,
            ChargeC=ChargeC,
            Pion_Cut=Pion_Cut, Kaon_Cut=Kaon_Cut, Electron_Cut=Electron_Cut, Gamma_Cut=Gamma_Cut,
            D0kmpip_Cut=D0kmpip_Cut, D0kmpippi0_Cut=D0kmpippi0_Cut, D0km3pi_Cut=D0km3pi_Cut,
            Ds_D0kmpip_Cut=Ds_D0kmpip_Cut, Ds_D0kmpippi0_Cut=Ds_D0kmpippi0_Cut, Ds_D0km3pi_Cut=Ds_D0km3pi_Cut,
            Ds_star_Cut=Ds_star_Cut,
            is_data=args.data
        )
        
        gamma_ROE, Dst0M1_ROE, Dst0M2_ROE, Dstp_ROE = build_roe_for_selected(
            main_path=path,
            dslist=dslist,
            ChargeC=ChargeC,
            ElectronROE_Cut=ElectronROE_Cut,
            PionROE_Cut=PionROE_Cut,
            GammaROE_Cut=GammaROE_Cut,
            is_mc=not args.data
        )
    
        ch_tag, expertVars = Suggestion(path, args.mode, args.pi0_list)
    
        (D0_mode_vars, e_vars, km_vars, pip_vars, pi0_vars,
        D0_vars, Ds_vars, Event) = Reconstruction_Variable(
            ch_tag=ch_tag,
            mode=args.mode,
            pi0_list=args.pi0_list,
            Electron=Electron, Before_VF=[], BCS=BCS,
            expertVars=expertVars,
            roe_vars=(gamma_ROE, Dst0M1_ROE, Dst0M2_ROE, Dstp_ROE),
            is_data=args.data
        )
    
    # Save
    outfile = f"/home/belle2/amubarak/output_test_{args.mode}.root"
    if args.mode == "kmpippi0":
        outfile = f"/home/belle2/amubarak/output_test_{args.mode}_{args.pi0_list}.root"
    outfileMC = f"/home/belle2/amubarak/output_MC_test_{args.mode}.root"
    
    if args.outfile:
        outfile = args.outfile
    
    if args.truth:
        ma.variablesToNtuple(ds_mc_list, Ds_mc_vars,
                            filename=outfileMC, treename=f"DsMCtree_{ch_tag}", path=path)
        ma.variablesToNtuple(d0_mc_list, D0_mc_vars,
                            filename=outfileMC, treename=f"D0MCtree_{ch_tag}", path=path)
    else:
        ds_list, ds_tree = {
            "kmpip":    ("D_s+:Ch1", "DstreeCh1"),
            "kmpippi0": ("D_s+:Ch2", "DstreeCh2"),
            "km3pi":    ("D_s+:Ch3", "DstreeCh3"),
        }[args.mode]
        
        ds_vars_for_save = e_vars + km_vars + pip_vars + D0_vars + Ds_vars + Event
        if args.mode == "kmpippi0":
            ds_vars_for_save = e_vars + km_vars + pip_vars + pi0_vars + D0_vars + Ds_vars + Event
        
        ma.variablesToNtuple(ds_list, ds_vars_for_save,
                            filename=outfile, treename=ds_tree, path=path)
        
        if args.mode == "kmpip":
            ma.variablesToNtuple('D0:kmpip', D0_mode_vars,
                                filename=outfile, treename='D02kmpiptree', path=path)
        elif args.mode == "kmpippi0":
            ma.variablesToNtuple('D0:kmpippi0', D0_mode_vars,
                                filename=outfile, treename='D02kmpippi0tree', path=path)
        else:
            ma.variablesToNtuple('D0:km3pi', D0_mode_vars,
                                filename=outfile, treename='D02km3pitree', path=path)
    
    # Monitor progress
    progress = b2.register_module("Progress")
    path.add_module(progress)
    
    # Process events
    b2.process(path)
    print(b2.statistics)