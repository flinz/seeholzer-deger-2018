from numpy import pi


class MetaParameters(object):
    pass


# static definitions
class MeanfieldParameters(object):
    """Parameters used for meanfield"""
    G_0 = "g_0"
    G_1 = "g_1"
    G_SIGMA = "g_sigma"
    G_R = "g_r"

    T_EXT = "T_EXT"
    T_I = "T_I"

    E = "e"
    I = "i"

    V_AVG_MF = "v_avg"
    PHI_FIRING = "phi_firing"
    CV = "cv"

    TAU_MEMB_EFF = "tau_memb_eff"
    MU_MF = "mu_mf"
    SIGMA_SQUARE = "sigma_square"
    W_J = "w_j"
    S_MF = "S_mf"

    G_EE = "g_ee"
    G_IE = "g_ie"
    G_EI = "g_ei"
    G_II = "g_ii"

    NU_I = "nu_i"
    JACOBIAN = "jacobian"
    JACOBIAN_VARS = [G_0, G_1, G_SIGMA, G_R, NU_I]
    JACOBIAN_VARS_P = {G_0: 0, G_1: 1, G_SIGMA: 2, G_R: 3, NU_I: 4}

    ALPHA = "alpha"
    BETA = "beta"

    P_EE = "p_ee"
    W_NOISE = "w_noise"
    EL_NOISE = "el_noise"
    EL_NOISE_ARRAY = "el_noise_array"
    EL_SEED = "el_seed"

    CORES = "cores"
    BASE_SEED = "base_seed"


class NeuronParameters(object):
    """Parameters used in network models and neurons"""

    GAMMA = "gamma"
    BETA = "beta"
    VE = "VE"
    VL = "VL"
    VI = "VI"
    VTHR = "Vthr"
    VRES = "Vres"
    TAU_NMDA = "tau_NMDA"
    TAU_NMDA_RISE = "tau_NMDA_RISE"
    TAU_GABA = "tau_GABA"
    TAU_AMPA = "tau_AMPA"
    TAU_RP = "tau_rp"
    TAU_M = "tau_m"
    CM = "Cm"
    GM = "gm"
    C_E = "CE"
    C_I = "CI"
    C_EXT = "Cext"
    NU_EXT = "nu_ext"
    G_AMPA = "gAMPA"
    G_GABA = "gGABA"
    G_NMDA = "gNMDA"

    W_SIGMA = "w_sigma"
    W_R = "w_r"
    W_1 = "w_1"
    W_0 = "w_0"

    STP_U = "stp_u"
    STP_TAU_R = "stp_tau_r"
    STP_TAU_F = "stp_tau_f"

    PHI_F_RATESCALE = "phi_f_ratescale"
    PHI_F_SCALE = "phi_f_scale"
    PHI_F_OFFSET = "phi_f_offset"


# These are deprecated
STANDARD_PARAMS = {
    "E": {
        NeuronParameters.GAMMA: 0.2801120448179272,
        NeuronParameters.BETA: 0.062,
        NeuronParameters.VE: 0.,
        NeuronParameters.VL: -70.,
        NeuronParameters.VI: -70.,
        NeuronParameters.VTHR: -50.,
        NeuronParameters.VRES: -60.,

        NeuronParameters.TAU_NMDA: 100.,
        NeuronParameters.TAU_NMDA_RISE: 2.,

        NeuronParameters.TAU_GABA: 10.,
        NeuronParameters.TAU_AMPA: 2.,
        NeuronParameters.TAU_RP: 2.,

        NeuronParameters.CM: 500.,
        NeuronParameters.GM: 25.,
        NeuronParameters.TAU_M: 20.,
        NeuronParameters.C_E: 800.,
        NeuronParameters.C_I: 200.,
        NeuronParameters.C_EXT: 1000,
        NeuronParameters.NU_EXT: 0.0024,

        NeuronParameters.G_AMPA: 2.08,
        NeuronParameters.G_NMDA: 0.327,
        NeuronParameters.G_GABA: 1.25,

        NeuronParameters.W_SIGMA: 18. / 360. * 2. * pi,
        NeuronParameters.W_1: 0.,
        NeuronParameters.W_0: 1.,
    },

    "I": {
        NeuronParameters.GAMMA: 0.2801120448179272,
        NeuronParameters.BETA: 0.062,
        NeuronParameters.VE: 0.,
        NeuronParameters.VL: -70.,
        NeuronParameters.VI: -70.,
        NeuronParameters.VTHR: -50.,
        NeuronParameters.VRES: -60.,

        NeuronParameters.TAU_NMDA: 100.,
        NeuronParameters.TAU_NMDA_RISE: 2.,

        NeuronParameters.TAU_GABA: 10.,
        NeuronParameters.TAU_AMPA: 2.,
        NeuronParameters.TAU_RP: 1.,

        NeuronParameters.CM: 200.,
        NeuronParameters.GM: 20.,
        NeuronParameters.TAU_M: 10.,
        NeuronParameters.C_E: 800.,
        NeuronParameters.C_I: 200.,
        NeuronParameters.C_EXT: 1000,
        NeuronParameters.NU_EXT: 0.0024,

        NeuronParameters.G_AMPA: 1.62,
        NeuronParameters.G_NMDA: 0.258,
        NeuronParameters.G_GABA: 0.973
    }
}


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
