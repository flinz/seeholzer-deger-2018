import warnings

import numpy as np
import pylab as pl
from scipy.special import erf

from classes.static import MeanfieldParameters as mpr
from classes.static import NeuronParameters as npr
from meanfield.shared_functions import w_function, g_function

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
    """Handles parameters of the networks, and is used for theoretical predictions"""

    def __init__(self, network):
        """
        :param network: Network instance to attach to handler
        """
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
        """Return parameter for a given population key and parameter key

        See classes.static for keys

        :param pop_key: population key, commonly "e" or "i", but use classes.static.MeanfieldParameters.E/I
        :param key: classes.static.NeuronParameters keys
        :return: network parameter value
        """
        return self.network.p_(pop_key, key)

    @property
    def p_e(self):
        """Return excitatory parameters"""
        return self.network.p_e

    @property
    def p_i(self):
        """Return inhibitory parameters"""
        return self.network.p_i

    def get_gauss_params(self):
        """calculate connectivity as in Compte (2000), gaussian normalized to 1.

        We only need to calculate W_0 from W_1 and W_SIGMA and return all three."""

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


    def get_dependent_params(self, pop_key):
        # get derived meanfield constants
        return {
            mpr.T_I: self.p_(pop_key, npr.G_GABA) * self.p_(pop_key, npr.C_I) * self.p_(pop_key,
                                                                                        npr.TAU_GABA) / self.p_(pop_key,
                                                                                                                npr.GM),
            mpr.T_EXT: self.p_(pop_key, npr.G_AMPA) * self.p_(pop_key, npr.C_EXT) * self.p_(pop_key,
                                                                                            npr.TAU_AMPA) / self.p_(
                pop_key, npr.GM)
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
        """Set parameters to population from a dictionary

        See classes.static.NeuronParameters for usable dictionary keys.

        :param pop_key: population key, commonly "e" or "i", but use classes.static.MeanfieldParameters.E/I
        :param dic: dictionary with classes.static.NeuronParameters keys
        :return:
        """
        if DEBUG:
            print "[set_params: %s] %s" % (pop_key, dic)

        if dic.has_key(npr.W_0) and dic.has_key(npr.W_1):
            warnings.warn(
                "Setting w_0 or w_1 does only influence the value of W_J, which then rescales both again.",
                RuntimeWarning)

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
        """Sets the network connectivity parameters from the gauss parameters"""
        if DEBUG:
            print "[update_states]"

        # set the gauss parameters to the E population
        self.network.__getattribute__(
            mpr.E).paramset.parameters.update(self.get_gauss_params())

    def set_states(self, dic):
        """sets meanfield states to dictionary values"""
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
        """Sets network recurrent confuctances (NMDA and GABA)"""
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
        """Returns network recurrent conductances (NMDA, GABA)"""
        return {
            "g_ee": self.p_(mpr.E, npr.G_NMDA),
            "g_ie": self.p_(mpr.I, npr.G_NMDA),
            "g_ei": self.p_(mpr.E, npr.G_GABA),
            "g_ii": self.p_(mpr.I, npr.G_GABA)
        }

    def w_func(self, x):
        """Calls the W function (gauss) with network parameters"""
        return w_function(x, self.p_(mpr.E, npr.W_0), self.p_(mpr.E, npr.W_1), self.p_(mpr.E, npr.W_SIGMA))

    def g_func(self, x):
        """Calls the g function (generalized gaussian) with meanfield states"""
        if self.p_mf[mpr.G_1] == 0:
            return self.p_mf[mpr.G_0]
        return g_function(x, self.p_mf[mpr.G_0], self.p_mf[mpr.G_1], self.p_mf[mpr.G_SIGMA], self.p_mf[mpr.G_R])

    def plot_rad_func(self, func, **kwargs):
        """Plots a function on [-pi,pi)"""
        xv = np.arange(0., 2. * np.pi, .01) - np.pi
        if func == self.g_func:
            yv = [func(x) for x in xv]
        else:
            yv = [func(x) for x in xv]
        pl.plot(xv, yv, **kwargs)
