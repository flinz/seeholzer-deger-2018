from classes.static import NeuronParameters as npr
from classes.static import MeanfieldParameters as mpr
from scipy.special import erf
from scipy.integrate import quad
from meanfield.shared_functions import w_function, g_function
import pylab as pl
import numpy as np
import warnings

DEBUG = False

# constants
PI = np.pi
SQRT_PI = np.sqrt(PI)
PI_2 = 2. * PI
SQRT_PI_2 = np.sqrt(PI_2)
E_1 = np.exp(-1.)

# parameters that, if changed, cause updates of dependent variables
update_params_list = [
    npr.W_SIGMA,
    npr.W_0,
    npr.W_1,
    npr.G_AMPA,
    npr.G_GABA,
    npr.C_I,
    npr.C_EXT,
    npr.TAU_GABA,
    npr.TAU_AMPA,
    npr.GM,
]

update_states_list = [
    mpr.W_J,
]

meanfield_params = [
    mpr.G_1,
    mpr.G_0,
    mpr.G_SIGMA,
    mpr.G_R,
    mpr.W_J,
    mpr.PHI_FIRING,
    mpr.V_AVG_MF,
    mpr.CV,
    mpr.TAU_MEMB_EFF,
    mpr.MU_MF,
    mpr.SIGMA_SQUARE,
]


# population keys
pop_keys = [mpr.E, mpr.I]


class MeanfieldHandler(object):

    def __init__(self, network):
        # do a deep copy - to be sure we don't mess up the db.
        print "[MeanfieldHandler] initializing for: %s" % network
        self.network = network

        self.p_mf = {
            # specific parameters
            mpr.E: {},
            mpr.I: {},

            # general parameters
            mpr.G_1: 0.,
            mpr.G_0: 1.,
            mpr.G_SIGMA: 1.,
            mpr.G_R: 2.,

            mpr.W_J: 1.,
        }
        self.update_params()

    def p_(self, pop_key, key):
        return self.network.p_(pop_key, key)

    @property
    def p_e(self):
        return self.network.p_e

    @property
    def p_i(self):
        return self.network.p_i

    # calculate connectivity as in Compte (2000), gaussian normalized to 1.
    def get_gauss_params(self):

        aa = -2. * np.sqrt(PI)
        bb = self.p_e[npr.W_SIGMA] * \
            np.sqrt(2) * erf(PI / np.sqrt(2) / self.p_e[npr.W_SIGMA])
        j0 = (aa + self.p_mf[mpr.W_J] * bb) / (aa + bb)
        assert j0 >= 0, "E-E connectivity could not be normalized to 1 (sigma too broad?)"

        out = {
            npr.W_0: j0,
            npr.W_1: (self.p_mf[mpr.W_J] - j0),
            npr.W_SIGMA: self.p_e[npr.W_SIGMA]
        }
        if DEBUG:
            print "[get_gauss_params] setting params %s (W_J was %.2f)" % (out, self.p_mf[mpr.W_J])

        return out

    # get derived meanfield constants
    def get_dependent_params(self, pop_key):
        return {
            mpr.T_I: self.p_(pop_key, npr.G_GABA) * self.p_(pop_key, npr.C_I) * self.p_(pop_key, npr.TAU_GABA) / self.p_(pop_key, npr.GM),
            mpr.T_EXT: self.p_(pop_key, npr.G_AMPA) * self.p_(pop_key, npr.C_EXT) * self.p_(pop_key, npr.TAU_AMPA) / self.p_(pop_key, npr.GM)
        }

    def update_params(self):

        # Recalculating derived state W_J
        self.p_mf[mpr.W_J] = self.p_(mpr.E, npr.W_0) + self.p_(mpr.E, npr.W_1)
        if DEBUG:
            print "[update_params] W_J changed to", self.p_mf[mpr.W_J]

        # set the gauss parameters to the E population
        self.network.__getattribute__(
            mpr.E).paramset.parameters.update(self.get_gauss_params())

        # set the population parameters
        for key in pop_keys:
            self.p_mf[key].update(self.get_dependent_params(key))

    def set_params(self, pop_key, dic):

        if DEBUG:
            print "[set_params: %s] %s" % (pop_key, dic)

        if dic.has_key(npr.W_0) and dic.has_key(npr.W_1):
            warnings.warn(
                "Setting w_0 or w_1 does only influence the value of W_J, which then rescales both again.", RuntimeWarning)

        update = False
        for key, value in dic.iteritems():
            assert self.network.__getattribute__(pop_key).paramset.parameters.has_key(
                key), "Parameter %s in %s does not exist" % (key, pop_key)
            self.network.__getattribute__(pop_key).paramset[key] = value
            if DEBUG:
                print "[set_params] %s changed to %s" % (key, self.network.__getattribute__(pop_key).paramset[key])
            if key in update_params_list:
                update = True

        if update:
            if DEBUG:
                print "[set_params] updating derived states"
            self.update_params()

    def update_states(self):
        if DEBUG:
            print "[update_states]"

        # set the gauss parameters to the E population
        self.network.__getattribute__(
            mpr.E).paramset.parameters.update(self.get_gauss_params())

    def set_states(self, dic):
        if DEBUG:
            print "[set_states] %s" % (dic)
        update = False
        for key, value in dic.iteritems():
            assert self.p_mf.has_key(
                key), "parameter %s in p_mf does not exist" % (key)
            self.p_mf[key] = value
            if DEBUG:
                print "[set_states] %s changed to %s" % (key, self.p_mf[key])
            if key in update_states_list:
                update = True
        if update:
            self.update_states()

    def set_conductances(self, dict):
        for key, value in dict.iteritems():
            if key == "g_ee":
                self.set_params(mpr.E, {npr.G_NMDA: value})
            elif key == "g_ie":
                self.set_params(mpr.I, {npr.G_NMDA: value})
            elif key == "g_ei":
                self.set_params(mpr.E, {npr.G_GABA: value})
            elif key == "g_ii":
                self.set_params(mpr.I, {npr.G_GABA: value})
            else:
                print "[set_conductances] discarding key %s" % key

    def get_conductances(self):
        return {
            "g_ee": self.p_(mpr.E, npr.G_NMDA),
            "g_ie": self.p_(mpr.I, npr.G_NMDA),
            "g_ei": self.p_(mpr.E, npr.G_GABA),
            "g_ii": self.p_(mpr.I, npr.G_GABA)
        }

    def w_func(self, x):
        return w_function(x, self.p_(mpr.E, npr.W_0), self.p_(mpr.E, npr.W_1), self.p_(mpr.E, npr.W_SIGMA))

    def g_func(self, x):
        if self.p_mf[mpr.G_1] == 0:
            return self.p_mf[mpr.G_0]
        return g_function(x, self.p_mf[mpr.G_0], self.p_mf[mpr.G_1], self.p_mf[mpr.G_SIGMA], self.p_mf[mpr.G_R])

    def plot_rad_func(self, func, **kwargs):
        xv = np.arange(0., 2. * np.pi, .01) - np.pi
        if func == self.g_func:
            yv = [func(x) for x in xv]
        else:
            yv = [func(x) for x in xv]
        pl.plot(xv, yv, **kwargs)

    # calculates the position dependent input from recurrent excitation by convolution of the
    # firing rate distribution with the EE nonlinearity

    def n_pos_func(self, phi):
        def integrand(x):
            # circular boundary!
            tmp = x - phi
            tmp = tmp - round(tmp / PI_2) * PI_2
            return self.w_func(tmp) * self.network.e.nmda_synapse(self.g_func(x))

        return 1. / PI_2 * quad(integrand, -PI, PI)[0]

        # full integral from -pi to pi neccessary.
        #
        # Q: but DO WE NORMALIZE BY 1/2PI OR NOT? Input should not be normalized, but the weights should.
        # A: By reproducing the brunel 2000 network as a limiting case, we get
        # the input normalization by 1/2pi

    def n_i_func(self):
        """
        Calculates from firing rate distribution the global input to inhibitory population.
        """
        def integrand(x):
            # here the inhibitory NMDA nonlinearity is used
            return self.network.i.nmda_synapse(self.g_func(x))
        return 1. / PI * (quad(integrand, 0, PI)[0])

    def beta(self, mu, sigma, tau_memb_eff, pop_key):
        tmp = (self.p_(pop_key, npr.VRES) -
               self.p_(pop_key, npr.VL) - mu) / sigma
        return tmp

    def alpha(self, mu, sigma, tau_memb_eff, pop_key):
        tmp =   -0.5 * self.p_(pop_key, npr.TAU_AMPA) / tau_memb_eff \
            + 1.03 * np.sqrt(self.p_(pop_key, npr.TAU_AMPA) / tau_memb_eff) \
            + (-mu - self.p_(pop_key, npr.VL) + self.p_(pop_key, npr.VTHR)) * \
            (1. + (0.5 * self.p_(pop_key, npr.TAU_AMPA) / tau_memb_eff)) / sigma
        return tmp

    def phi_firing_func_limits(self, alpha, beta, tau_memb_eff, pop_key):
        def integrand(x):
            if x < -7.:
                return np.exp(7.**2) * (1. + erf(7.))
            if x > 7.:
                return 0.
            return np.exp(x**2) * (1. + erf(x))
        return 1. / (self.p_(pop_key, npr.TAU_RP) + tau_memb_eff * SQRT_PI * quad(integrand, beta, alpha)[0])

    def phi_firing_func(self, mu, sigma, tau_memb_eff, pop_key, return_all=0):
        alpha = self.alpha(mu, sigma, tau_memb_eff, pop_key)
        beta = self.beta(mu, sigma, tau_memb_eff, pop_key)
        out = self.phi_firing_func_limits(alpha, beta, tau_memb_eff, pop_key)
        if return_all:
            return {mpr.PHI_FIRING: out,
                    mpr.ALPHA: alpha,
                    mpr.BETA: beta
                    }
        return out

    def cv_func(self, phi_firing, mu, sigma, tau_memb_eff, pop_key):

        def integrand1(x):
            return np.exp(x**2) * (1. + erf(x))**2

        def integrand2(x):
            return np.exp(x**2) * quad(integrand1, -10., x)[0]
        return (phi_firing**2 * 2. * PI * tau_memb_eff**2 * quad(integrand2, self.beta(mu, sigma, tau_memb_eff, pop_key), self.alpha(mu, sigma, tau_memb_eff, pop_key))[0])**(.5)


    def v_avg_mf_func(self, mu, nu, tau_memb_eff, pop_key):
        return self.p_(pop_key, npr.VL) + mu - (self.p_(pop_key, npr.VTHR) - self.p_(pop_key, npr.VRES)) * nu * tau_memb_eff


    def sigma_square_func(self, tau_memb_eff, v_avg, pop_key):

        return (
            self.p_(pop_key, npr.C_EXT) * (

                self.p_(pop_key, npr.G_AMPA) * self.p_(pop_key, npr.TAU_AMPA) * (

                    v_avg - self.p_(pop_key, npr.VE)

                ) / (self.p_(pop_key, npr.GM) * self.p_(pop_key, npr.TAU_M))

            )**2 * self.p_(pop_key, npr.NU_EXT) * tau_memb_eff
        )


    def S_mf_func(self, n_subpop, nu_i, v_avg, pop_key):
        # this is a mathematica simplification, the + sign in the exp is
        # correct
        tmp = np.exp(v_avg * self.p_(pop_key, npr.BETA))
        return 1. + (
            (
                    (
                        self.p_(pop_key, npr.C_E) * tmp * self.p_(pop_key, npr.G_NMDA) * n_subpop * (
                            tmp + self.p_(pop_key, npr.GAMMA) + v_avg * self.p_(pop_key, npr.BETA) * self.p_(
                                pop_key, npr.GAMMA) - self.p_(pop_key, npr.VE) * self.p_(pop_key, npr.BETA) * self.p_(pop_key, npr.GAMMA)
                        )
                    )
                / (tmp + self.p_(pop_key, npr.GAMMA))**2 / self.p_(pop_key, npr.GM)
            ) + self.p_mf[pop_key][mpr.T_EXT] * self.p_(pop_key, npr.NU_EXT) + self.p_mf[pop_key][mpr.T_I] * nu_i
        )

    def mu_mf_func(self, S_mf, n_subpop, nu_i, v_avg, pop_key):
        # this is a mathematica simplification, the + sign in the exp is
        # correct
        tmp = np.exp(v_avg * self.p_(pop_key, npr.BETA))
        tmp2 = (
            (self.p_(pop_key, npr.C_E) * self.p_(pop_key, npr.G_NMDA) * tmp * n_subpop * v_avg * (v_avg - self.p_(pop_key, npr.VE)) *
             self.p_(pop_key, npr.BETA) * self.p_(pop_key, npr.GAMMA)) / (self.p_(pop_key, npr.GM) * (tmp + self.p_(pop_key, npr.GAMMA))**2)
            + self.p_(pop_key, npr.VE) * ((self.p_(pop_key, npr.C_E) * self.p_(pop_key, npr.G_NMDA) * n_subpop) / (self.p_(pop_key,
                                                                                                                           npr.GM) * (1 + self.p_(pop_key, npr.GAMMA) / tmp)) + self.p_mf[pop_key][mpr.T_EXT] * self.p_(pop_key, npr.NU_EXT))
            + self.p_mf[pop_key][mpr.T_I] * self.p_(pop_key, npr.VI) * nu_i
        )
        return (tmp2 - self.p_(pop_key, npr.VL) * (S_mf - 1)) / S_mf

    def tau_memb_eff_func(self, S_mf, pop_key):
        return self.p_(pop_key, npr.TAU_M) / S_mf

    # JACOBIAN ------------------------------------------------------
    # derivative of the siegert formula, here to the distribution parameters in first iteration,
    # although the earlier derivatives can be reused later for the full model
    # probably.

    def d_g_dg0(self, x):
        return 1.

    def d_g_dg1(self, x):
        return np.exp(-(abs(x) / self.p_mf[mpr.G_SIGMA])**self.p_mf[mpr.G_R])

    def d_g_dgsig(self, x):
        return (self.g_func(x) - self.p_mf[mpr.G_0]) * (np.abs(x) / self.p_mf[mpr.G_SIGMA])**(self.p_mf[mpr.G_R]) * self.p_mf[mpr.G_R] / self.p_mf[mpr.G_SIGMA]

    def d_g_dgr(self, x):
        if x == 0:
            return 0.
        arg_x = (np.abs(x) / self.p_mf[mpr.G_SIGMA])
        return -(self.g_func(x) - self.p_mf[mpr.G_0]) * arg_x**(self.p_mf[mpr.G_R]) * np.log(arg_x)

    def d_simple_nEE_dg0(self, phi):
        assert self.network.e.nmda_synapse.is_nmda == False, "only implemented for non-nmda synapses"

        def integrand(x):
            # circular boundary!
            tmp = x - phi
            tmp = tmp - round(tmp / PI_2) * PI_2
            g_x = self.g_func(x)
            return (
                self.w_func(tmp) *
                (
                    self.network.e.nmda_synapse.stp_ur(g_x)
                    + g_x * self.network.e.nmda_synapse.stp_dur_dx(g_x)
                )
            )
        return self.network.e.nmda_synapse.tau_post / PI_2 * quad(integrand, -PI, PI)[0]

    def d_simple_nEE_dg1(self, phi):
        assert self.network.e.nmda_synapse.is_nmda == False, "only implemented for non-nmda synapses"

        def integrand(x):
            # circular boundary!
            tmp = x - phi
            tmp = tmp - round(tmp / PI_2) * PI_2
            g_x = self.g_func(x)
            exp_x = np.exp(-(abs(x) /
                             self.p_mf[mpr.G_SIGMA])**self.p_mf[mpr.G_R])
            return (
                self.w_func(tmp) *
                exp_x *
                (
                    self.network.e.nmda_synapse.stp_dur_dx(g_x) * g_x
                    + self.network.e.nmda_synapse.stp_ur(g_x)
                )
            )
        return self.network.e.nmda_synapse.tau_post / PI_2 * quad(integrand, -PI, PI)[0]

    def d_simple_nEE_dgsig(self, phi):
        assert self.network.e.nmda_synapse.is_nmda == False, "only implemented for non-nmda synapses"

        def integrand(x):

            # circular boundary!
            tmp = x - phi
            tmp = tmp - round(tmp / PI_2) * PI_2
            g_x = self.g_func(x)
            exp_x = (g_x - self.p_mf[mpr.G_0]) * (np.abs(x) / self.p_mf[mpr.G_SIGMA])**(
                self.p_mf[mpr.G_R]) * self.p_mf[mpr.G_R] / self.p_mf[mpr.G_SIGMA]
            return (
                self.w_func(tmp) *
                exp_x *
                (
                    self.network.e.nmda_synapse.stp_dur_dx(g_x) * g_x
                    + self.network.e.nmda_synapse.stp_ur(g_x)
                )
            )
        return self.network.e.nmda_synapse.tau_post / PI_2 * quad(integrand, -PI, PI)[0]

    def d_simple_nEE_dgr(self, phi):
        assert self.network.e.nmda_synapse.is_nmda == False, "only implemented for non-nmda synapses"

        def integrand(x):

            if x == 0:
                return 0.
            # circular boundary!
            tmp = x - phi
            tmp = tmp - round(tmp / PI_2) * PI_2
            g_x = self.g_func(x)
            arg_x = (np.abs(x) / self.p_mf[mpr.G_SIGMA])
            exp_x = - (g_x - self.p_mf[mpr.G_0]) * \
                arg_x**(self.p_mf[mpr.G_R]) * np.log(arg_x)
            return (
                self.w_func(tmp) *
                exp_x *
                (
                    self.network.e.nmda_synapse.stp_dur_dx(g_x) * g_x
                    + self.network.e.nmda_synapse.stp_ur(g_x)
                )
            )
        return self.network.e.nmda_synapse.tau_post / PI_2 * quad(integrand, -PI, PI)[0]

    def d_simple_nIE_dg0(self):
        assert self.network.i.nmda_synapse.is_nmda == False, "only implemented for non-nmda synapses"

        def integrand(x):
            # circular boundary!
            g_x = self.g_func(x)
            return (
                self.network.i.nmda_synapse.stp_ur(g_x)
                + g_x * self.network.i.nmda_synapse.stp_dur_dx(g_x)
            )
        return self.network.i.nmda_synapse.tau_post / PI * quad(integrand, 0, PI)[0]

    def d_simple_nIE_dg1(self):
        assert self.network.i.nmda_synapse.is_nmda == False, "only implemented for non-nmda synapses"

        def integrand(x):
            # circular boundary!
            g_x = self.g_func(x)
            exp_x = np.exp(-(abs(x) /
                             self.p_mf[mpr.G_SIGMA])**self.p_mf[mpr.G_R])
            return (
                exp_x *
                (
                    self.network.i.nmda_synapse.stp_dur_dx(g_x) * g_x
                    + self.network.i.nmda_synapse.stp_ur(g_x)
                )
            )
        return self.network.i.nmda_synapse.tau_post / PI * quad(integrand, 0, PI)[0]

    def d_simple_nIE_dgsig(self):
        assert self.network.i.nmda_synapse.is_nmda == False, "only implemented for non-nmda synapses"

        def integrand(x):
            # circular boundary!
            g_x = self.g_func(x)
            exp_x = (g_x - self.p_mf[mpr.G_0]) * (np.abs(x) / self.p_mf[mpr.G_SIGMA])**(
                self.p_mf[mpr.G_R]) * self.p_mf[mpr.G_R] / self.p_mf[mpr.G_SIGMA]
            return (
                exp_x *
                (
                    self.network.i.nmda_synapse.stp_dur_dx(g_x) * g_x
                    + self.network.i.nmda_synapse.stp_ur(g_x)
                )
            )
        return self.network.i.nmda_synapse.tau_post / PI * quad(integrand, 0, PI)[0]

    def d_simple_nIE_dgr(self):
        assert self.network.i.nmda_synapse.is_nmda == False, "only implemented for non-nmda synapses"

        def integrand(x):
            if x == 0:
                return 0.
            # circular boundary!
            g_x = self.g_func(x)
            arg_x = (np.abs(x) / self.p_mf[mpr.G_SIGMA])
            exp_x = - (g_x - self.p_mf[mpr.G_0]) * \
                arg_x**(self.p_mf[mpr.G_R]) * np.log(arg_x)
            return (
                exp_x *
                (
                    self.network.i.nmda_synapse.stp_dur_dx(g_x) * g_x
                    + self.network.i.nmda_synapse.stp_ur(g_x)
                )
            )
        return self.network.i.nmda_synapse.tau_post / PI * quad(integrand, 0, PI)[0]

    def d_S_mf_dnui(self, d_n_ddist, pop_key):
        return self.p_mf[pop_key][mpr.T_I]

    def d_S_mf_ddist(self, d_n_ddist, pop_key):
        return (self.p_(pop_key, npr.C_E) * self.p_(pop_key, npr.G_NMDA) * d_n_ddist) / self.p_(pop_key, npr.GM)

    def d_mu_mf_dnui(self, mu_mf, S_mf, d_S_mf_dnui, d_n_ddist, pop_key):
        return (
            - mu_mf * d_S_mf_dnui +
            self.p_mf[pop_key][mpr.T_I] * self.p_(pop_key, npr.VI)
            - self.p_(pop_key, npr.VL) * d_S_mf_dnui
        ) / S_mf

    def d_mu_mf_ddist(self, mu_mf, S_mf, d_S_mf_ddist, d_n_ddist, pop_key):
        return (
            - mu_mf * d_S_mf_ddist
            + self.p_(pop_key, npr.C_E) * self.p_(pop_key, npr.G_NMDA) /
            self.p_(pop_key, npr.GM) * d_n_ddist * self.p_(pop_key, npr.VE)
            - self.p_(pop_key, npr.VL) * d_S_mf_ddist
        ) / S_mf

    def d_alpha_deps(self, tau_memb_eff, d_tau_memb_eff_deps, mu_mf, d_mu_mf_deps, sigma_square, d_sigma_square_deps, pop_key):

        return (
            -1. * (
                ((1 + (0.5 * self.p_(pop_key, npr.TAU_AMPA)) / tau_memb_eff)
                 * d_mu_mf_deps) / np.sqrt(sigma_square)
            )
            - (
                (-self.p_(pop_key, npr.VL) + self.p_(pop_key, npr.VTHR) - mu_mf) * (1 +
                                                                                    (0.5 * self.p_(pop_key, npr.TAU_AMPA)) / tau_memb_eff) * d_sigma_square_deps
            ) / (2. * sigma_square**(1.5))

            + 0.5 * self.p_(pop_key, npr.TAU_AMPA) *
            d_tau_memb_eff_deps / tau_memb_eff**2
            - (0.5 * self.p_(pop_key, npr.TAU_AMPA) * (-self.p_(pop_key, npr.VL) + self.p_(pop_key, npr.VTHR) - mu_mf) * d_tau_memb_eff_deps) /
            (np.sqrt(sigma_square) * tau_memb_eff**2)
            - (0.515 * self.p_(pop_key, npr.TAU_AMPA) * d_tau_memb_eff_deps) /
            (np.sqrt(self.p_(pop_key, npr.TAU_AMPA) / tau_memb_eff) * tau_memb_eff**2)
        )

    def d_beta_deps(self, mu_mf, d_mu_mf_deps, sigma_square, d_sigma_square_deps, pop_key):
        return (
            - (d_mu_mf_deps / np.sqrt(sigma_square))
            - (
                (-self.p_(pop_key, npr.VL) + self.p_(pop_key,
                                                     npr.VRES) - mu_mf) * d_sigma_square_deps
                / (2. * sigma_square**1.5)
            )
        )

    def d_phi_firing_deps(self, phi_firing, tau_memb_eff, d_tau_memb_eff_deps, alpha, d_alpha_deps, beta, d_beta_deps, pop_key):
        return (
            phi_firing * (
                SQRT_PI * phi_firing * tau_memb_eff**2 * (
                    np.exp(beta**2) * (1 + erf(beta)) * d_beta_deps -
                    np.exp(alpha**2) * (1 + erf(alpha)) * d_alpha_deps
                )
                + (-1 + phi_firing * self.p_(pop_key, npr.TAU_RP)) *
                d_tau_memb_eff_deps
            ) / tau_memb_eff
        )

    def d_tau_memb_eff_deps(self, S_mf, d_S_mf_deps, pop_key):
        return -self.p_(pop_key, npr.CM) * d_S_mf_deps / (self.p_(pop_key, npr.GM) * S_mf**2)

    def d_v_avg_mf_deps(self, d_mu_mf_deps, nu, d_nu_deps, tau_memb_eff, d_tau_memb_eff_deps, pop_key):
        tmpdiff = (self.p_(pop_key, npr.VTHR) - self.p_(pop_key, npr.VRES))
        return (
            d_mu_mf_deps - tmpdiff *
            (nu * d_tau_memb_eff_deps + tau_memb_eff * d_nu_deps)
        )

    def d_sigma_square_deps(self, sigma_square, v_avg_mf, d_v_avg_mf_deps, tau_memb_eff, d_tau_memb_eff_deps, pop_key):
        return sigma_square * (
            2. * d_v_avg_mf_deps / (v_avg_mf - self.p_(pop_key, npr.VE))
            + d_tau_memb_eff_deps / tau_memb_eff
        )

    # MEANFIELD CALCULATION --------------------------------------------------

    def ret_zero(*args, **kwargs):
        return 0.

    def ret_one(*args, **kwargs):
        return 1.

    d_n_funcs_i = [d_simple_nIE_dg0, d_simple_nIE_dg1,
                   d_simple_nIE_dgsig, d_simple_nIE_dgr, ret_zero]
    d_n_funcs_e = [d_simple_nEE_dg0, d_simple_nEE_dg1,
                   d_simple_nEE_dgsig, d_simple_nEE_dgr, ret_zero]
    d_S_funcs = 4 * [d_S_mf_ddist] + [d_S_mf_dnui]
    d_mu_funcs = 4 * [d_mu_mf_ddist] + [d_mu_mf_dnui]

    d_tau_funcs = 5 * [d_tau_memb_eff_deps]

    d_nu_funcs = [d_g_dg0, d_g_dg1, d_g_dgsig, d_g_dgr, ret_zero]
    d_nu_in_funcs = 4 * [ret_zero] + [ret_one]

    d_v_avg_funcs = 5 * [d_v_avg_mf_deps]
    d_sigma_funcs = 5 * [d_sigma_square_deps]

    d_alpha_funcs = 5 * [d_alpha_deps]
    d_beta_funcs = 5 * [d_beta_deps]
    d_phi_funcs = 5 * [d_phi_firing_deps]

    def calc_simple_jacobian(self, nu, n, nu_i, pos, out, pop_key):
        is_inhib = pop_key == mpr.I

        if is_inhib:
            d_n = [f(self) for f in self.d_n_funcs_i]
        else:
            d_n = [f(self, pos) for f in self.d_n_funcs_e]

        d_S = [f(self, d_n[i], pop_key) for i, f in enumerate(self.d_S_funcs)]
        d_mu = [f(self, out[mpr.MU_MF], out[mpr.S_MF], d_S[i], d_n[i], pop_key)
                for i, f in enumerate(self.d_mu_funcs)]

        d_tau = [f(self, out[mpr.S_MF], d_S[i], pop_key)
                 for i, f in enumerate(self.d_tau_funcs)]

        if is_inhib:
            d_nu = [f(self, pos) for i, f in enumerate(self.d_nu_in_funcs)]
        else:
            d_nu = [f(self, pos) for i, f in enumerate(self.d_nu_funcs)]

        d_v_avg = [f(self, d_mu[i], nu, d_nu[i], out[mpr.TAU_MEMB_EFF], d_tau[
                     i], pop_key) for i, f in enumerate(self.d_v_avg_funcs)]

        d_sigma = [f(self, out[mpr.SIGMA_SQUARE], out[mpr.V_AVG_MF], d_v_avg[i], out[
                     mpr.TAU_MEMB_EFF], d_tau[i], pop_key) for i, f in enumerate(self.d_sigma_funcs)]

        d_alpha = [f(self, out[mpr.TAU_MEMB_EFF], d_tau[i], out[mpr.MU_MF], d_mu[i], out[
                     mpr.SIGMA_SQUARE], d_sigma[i], pop_key) for i, f in enumerate(self.d_alpha_funcs)]
        d_beta = [f(self, out[mpr.MU_MF], d_mu[i], out[mpr.SIGMA_SQUARE], d_sigma[
                    i], pop_key) for i, f in enumerate(self.d_beta_funcs)]
        d_phi = [f(self, out[mpr.PHI_FIRING], out[mpr.TAU_MEMB_EFF], d_tau[i], out[mpr.ALPHA], d_alpha[
                   i], out[mpr.BETA], d_beta[i], pop_key) for i, f in enumerate(self.d_phi_funcs)]

        return d_phi

    def calc_mf_simple_for_pos(self, pos, nu_i):

        nu = self.g_func(pos)
        n = self.n_pos_func(pos)

        # this is independent of Vavg for the simple case!
        S_mf = self.S_mf_func(n, nu_i, 0., mpr.E)
        tau_memb_eff = self.tau_memb_eff_func(S_mf, mpr.E)
        # this is independent of Vavg for the simple case!
        mu_mf = self.mu_mf_func(S_mf, n, nu_i, 0., mpr.E)

        # we can calculate the V_acg simply from the meanfield eqs.
        v_avg = self.v_avg_mf_func(mu_mf, nu, tau_memb_eff, mpr.E)
        sigma_square = self.sigma_square_func(tau_memb_eff, v_avg, mpr.E)

        phi_firing = self.phi_firing_func(mu_mf, np.sqrt(
            sigma_square), tau_memb_eff, mpr.E, return_all=1)
        out = {
            mpr.PHI_FIRING: phi_firing[mpr.PHI_FIRING],
            mpr.V_AVG_MF: v_avg
        }
        out[mpr.TAU_MEMB_EFF] = tau_memb_eff
        out[mpr.MU_MF] = mu_mf
        out[mpr.SIGMA_SQUARE] = sigma_square
        out[mpr.S_MF] = S_mf
        out[mpr.ALPHA] = phi_firing[mpr.ALPHA]
        out[mpr.BETA] = phi_firing[mpr.BETA]
        out["nu"] = nu
        out["n"] = n
        return out

    def calc_dphi_dDL(self, x, nu_i=None, out_mf=None):

        if not out_mf:
            out_mf = self.calc_mf_simple_for_pos(x, nu_i)

        nu = out_mf["nu"]
        n = out_mf["n"]

        d_S = 0.

        S = out_mf[mpr.S_MF]
        d_mu = 1. / S

        d_tau = 0.
        d_nu = 0.

        d_v_avg = self.d_v_avg_mf_deps(
            d_mu, nu, d_nu, out_mf[mpr.TAU_MEMB_EFF], d_tau, mpr.E)
        d_sigma = self.d_sigma_square_deps(out_mf[mpr.SIGMA_SQUARE], out_mf[
                                           mpr.V_AVG_MF], d_v_avg, out_mf[mpr.TAU_MEMB_EFF], d_tau, mpr.E)

        d_alpha = self.d_alpha_deps(out_mf[mpr.TAU_MEMB_EFF], d_tau, out_mf[
                                    mpr.MU_MF], d_mu, out_mf[mpr.SIGMA_SQUARE], d_sigma, mpr.E)
        d_beta = self.d_beta_deps(out_mf[mpr.MU_MF], d_mu, out_mf[
                                  mpr.SIGMA_SQUARE], d_sigma, mpr.E)
        d_phi = self.d_phi_firing_deps(out_mf[mpr.PHI_FIRING], out_mf[mpr.TAU_MEMB_EFF], d_tau, out_mf[
                                       mpr.ALPHA], d_alpha, out_mf[mpr.BETA], d_beta, mpr.E)

        return {
            "d_phi": d_phi,
            "d_alpha": d_alpha,
            "d_beta": d_beta,
            "d_sigma": d_sigma,

            "d_v_avg": d_v_avg,
            "d_mu": d_mu,
        }

    def calc_dphi_dJstruct(self, x, nu_i=None, out_mf=None, pop_key=mpr.E):

        if not out_mf:
            out_mf = self.calc_mf_simple_for_pos(x, nu_i)

        nu = out_mf["nu"]
        n = out_mf["n"]

        # from theory, J struct directly acts on n
        d_n = 1.

        d_S = self.d_S_mf_ddist(d_n, pop_key)
        d_mu = self.d_mu_mf_ddist(out_mf[mpr.MU_MF], out_mf[
                                  mpr.S_MF], d_S, d_n, pop_key)

        d_tau = self.d_tau_memb_eff_deps(out_mf[mpr.S_MF], d_S, pop_key)

        d_nu = 0.

        d_v_avg = self.d_v_avg_mf_deps(
            d_mu, nu, d_nu, out_mf[mpr.TAU_MEMB_EFF], d_tau, pop_key)
        d_sigma = self.d_sigma_square_deps(out_mf[mpr.SIGMA_SQUARE], out_mf[
                                           mpr.V_AVG_MF], d_v_avg, out_mf[mpr.TAU_MEMB_EFF], d_tau, pop_key)

        d_alpha = self.d_alpha_deps(out_mf[mpr.TAU_MEMB_EFF], d_tau, out_mf[
                                    mpr.MU_MF], d_mu, out_mf[mpr.SIGMA_SQUARE], d_sigma, pop_key)
        d_beta = self.d_beta_deps(out_mf[mpr.MU_MF], d_mu, out_mf[
                                  mpr.SIGMA_SQUARE], d_sigma, pop_key)
        d_phi = self.d_phi_firing_deps(out_mf[mpr.PHI_FIRING], out_mf[mpr.TAU_MEMB_EFF], d_tau, out_mf[
                                       mpr.ALPHA], d_alpha, out_mf[mpr.BETA], d_beta, pop_key)

        return {
            "d_phi": d_phi,
            "d_alpha": d_alpha,
            "d_beta": d_beta,
            "d_sigma": d_sigma,

            "d_v_avg": d_v_avg,
            "d_mu": d_mu,
            "d_S": d_S
        }

    def calc_dphi_dn_i(self, x, nu_i=None, out_mf=None, pop_key=mpr.E):

        if not out_mf:
            out_mf = self.calc_mf_simple_for_pos(x, nu_i)

        nu = out_mf["nu"]

        # from theory, J struct directly acts on n
        d_n = 1.

        d_S = self.d_S_mf_dnui(d_n, pop_key)
        d_mu = self.d_mu_mf_dnui(out_mf[mpr.MU_MF], out_mf[
                                 mpr.S_MF], d_S, d_n, pop_key)

        d_tau = self.d_tau_memb_eff_deps(out_mf[mpr.S_MF], d_S, pop_key)

        # from theory? shouldnt change...
        d_nu = 0.

        d_v_avg = self.d_v_avg_mf_deps(
            d_mu, nu, d_nu, out_mf[mpr.TAU_MEMB_EFF], d_tau, pop_key)
        d_sigma = self.d_sigma_square_deps(out_mf[mpr.SIGMA_SQUARE], out_mf[
                                           mpr.V_AVG_MF], d_v_avg, out_mf[mpr.TAU_MEMB_EFF], d_tau, pop_key)

        d_alpha = self.d_alpha_deps(out_mf[mpr.TAU_MEMB_EFF], d_tau, out_mf[
                                    mpr.MU_MF], d_mu, out_mf[mpr.SIGMA_SQUARE], d_sigma, pop_key)
        d_beta = self.d_beta_deps(out_mf[mpr.MU_MF], d_mu, out_mf[
                                  mpr.SIGMA_SQUARE], d_sigma, pop_key)
        d_phi = self.d_phi_firing_deps(out_mf[mpr.PHI_FIRING], out_mf[mpr.TAU_MEMB_EFF], d_tau, out_mf[
                                       mpr.ALPHA], d_alpha, out_mf[mpr.BETA], d_beta, pop_key)

        return {
            "d_phi": d_phi,
            "d_alpha": d_alpha,
            "d_beta": d_beta,
            "d_sigma": d_sigma,

            "d_v_avg": d_v_avg,
            "d_mu": d_mu,
            "d_S": d_S,
        }

    def calc_mf_simple(self, nu, n, nu_i, pop_key, pos=0., cv=0, return_all=0, jacobian=0):

        # this is independent of Vavg for the simple case!
        S_mf = self.S_mf_func(n, nu_i, 0., pop_key)
        tau_memb_eff = self.tau_memb_eff_func(S_mf, pop_key)
        # this is independent of Vavg for the simple case!
        mu_mf = self.mu_mf_func(S_mf, n, nu_i, 0., pop_key)

        # we can calculate the V_acg simply from the meanfield eqs.
        v_avg = self.v_avg_mf_func(mu_mf, nu, tau_memb_eff, pop_key)
        sigma_square = self.sigma_square_func(tau_memb_eff, v_avg, pop_key)

        phi_firing = self.phi_firing_func_spline(mu_mf, np.sqrt(
            sigma_square), tau_memb_eff, pop_key, return_all=1)
        out = {
            mpr.PHI_FIRING: phi_firing[mpr.PHI_FIRING],
            "nu": phi_firing[mpr.PHI_FIRING],
            "n": n,
            mpr.V_AVG_MF: v_avg
        }

        if return_all == 1 or jacobian:
            out[mpr.TAU_MEMB_EFF] = tau_memb_eff
            out[mpr.MU_MF] = mu_mf
            out[mpr.SIGMA_SQUARE] = sigma_square
            out[mpr.S_MF] = S_mf
            out[mpr.ALPHA] = phi_firing[mpr.ALPHA]
            out[mpr.BETA] = phi_firing[mpr.BETA]

        if jacobian:
            out[mpr.JACOBIAN] = self.calc_simple_jacobian(
                nu, n, nu_i, pos, out, pop_key)

        if cv == 1:
            cv = self.cv_func(phi_firing, mu_mf, np.sqrt(
                sigma_square), tau_memb_eff, pop_key)
            out[mpr.CV] = cv

        return out
