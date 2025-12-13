"""
Final variables for BDT training after applying all drop rules
Generated from ML00-Fake_D0.ipynb

Usage:
    from final_variables import VARIABLES
    
    # Get all variables for kmpip mode
    vars_kmpip = VARIABLES["kmpip"]["all_vars"]
    
    # Get only D0 variables for km3pi mode
    d0_vars_km3pi = VARIABLES["km3pi"]["D0_vars"]
"""

VARIABLES = {
    "kmpip": {
        "D0_vars": [
            "D0_useCMSFrame_p",
            "D0_chiProb",
            "D0_daughterAngle_0_1",
            "D0_cos_decayAngle_0",
            "D0_flightDistance",
            "D0_flightDistanceErr",
        ],
        "K_vars": [
            "K_Ch1_kaonID",
            "K_Ch1_pt",
            "K_Ch1_dr",
        ],
        "pi_vars": [
            "pi_Ch1_pionID",
            "pi_Ch1_pt",
            "pi_Ch1_dr",
        ],
        "pi0_vars": [
        ],
        "all_vars": [
            "D0_useCMSFrame_p",
            "D0_chiProb",
            "D0_daughterAngle_0_1",
            "D0_cos_decayAngle_0",
            "D0_flightDistance",
            "D0_flightDistanceErr",
            "K_Ch1_kaonID",
            "K_Ch1_pt",
            "K_Ch1_dr",
            "pi_Ch1_pionID",
            "pi_Ch1_pt",
            "pi_Ch1_dr",
        ],
    },
    "km3pi": {
        "D0_vars": [
            "D0_cos_decayAngle_2",
            "D0_cos_decayAngle_3",
            "D0_daughterInvM_0_1",
            "D0_daughterInvM_0_3",
            "D0_daughterInvM_1_2",
            "D0_daughterInvM_2_3",
            "D0_useCMSFrame_p",
            "D0_chiProb",
            "D0_daughterAngle_0_1",
            "D0_cos_decayAngle_0",
            "D0_cos_decayAngle_1",
            "D0_flightDistance",
            "D0_flightDistanceErr",
        ],
        "K_vars": [
            "K_Ch3_kaonID",
            "K_Ch3_pt",
            "K_Ch3_dr",
        ],
        "pi_vars": [
            "pi1_Ch3_pionID",
            "pi1_Ch3_pt",
            "pi1_Ch3_dr",
            "pi2_Ch3_pionID",
            "pi2_Ch3_pt",
            "pi2_Ch3_dr",
            "pi3_Ch3_pionID",
            "pi3_Ch3_pt",
            "pi3_Ch3_dr",
        ],
        "pi0_vars": [
        ],
        "all_vars": [
            "D0_cos_decayAngle_2",
            "D0_cos_decayAngle_3",
            "D0_daughterInvM_0_1",
            "D0_daughterInvM_0_3",
            "D0_daughterInvM_1_2",
            "D0_daughterInvM_2_3",
            "D0_useCMSFrame_p",
            "D0_chiProb",
            "D0_daughterAngle_0_1",
            "D0_cos_decayAngle_0",
            "D0_cos_decayAngle_1",
            "D0_flightDistance",
            "D0_flightDistanceErr",
            "K_Ch3_kaonID",
            "K_Ch3_pt",
            "K_Ch3_dr",
            "pi1_Ch3_pionID",
            "pi1_Ch3_pt",
            "pi1_Ch3_dr",
            "pi2_Ch3_pionID",
            "pi2_Ch3_pt",
            "pi2_Ch3_dr",
            "pi3_Ch3_pionID",
            "pi3_Ch3_pt",
            "pi3_Ch3_dr",
        ],
    },
    "kmpippi0_eff20_May2020": {
        "D0_vars": [
            "D0_momentaTripleProduct_0_1_2",
            "D0_daughterInvM_0_1",
            "D0_daughterInvM_0_2",
            "D0_daughterInvM_1_2",
            "D0_useCMSFrame_p",
            "D0_chiProb",
            "D0_daughterAngle_0_1",
            "D0_cos_decayAngle_0",
            "D0_cos_decayAngle_1",
            "D0_flightDistance",
            "D0_flightDistanceErr",
        ],
        "K_vars": [
            "K_Ch2_kaonID",
            "K_Ch2_pt",
            "K_Ch2_dr",
        ],
        "pi_vars": [
            "pi_Ch2_pionID",
            "pi_Ch2_pt",
            "pi_Ch2_dr",
        ],
        "pi0_vars": [
            "pi0_Ch2_pt",
        ],
        "all_vars": [
            "D0_momentaTripleProduct_0_1_2",
            "D0_daughterInvM_0_1",
            "D0_daughterInvM_0_2",
            "D0_daughterInvM_1_2",
            "D0_useCMSFrame_p",
            "D0_chiProb",
            "D0_daughterAngle_0_1",
            "D0_cos_decayAngle_0",
            "D0_cos_decayAngle_1",
            "D0_flightDistance",
            "D0_flightDistanceErr",
            "K_Ch2_kaonID",
            "K_Ch2_pt",
            "K_Ch2_dr",
            "pi_Ch2_pionID",
            "pi_Ch2_pt",
            "pi_Ch2_dr",
            "pi0_Ch2_pt",
        ],
    },
}
