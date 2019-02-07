import time

import nest
import nest.raster_plot
import nest.topology as tp
import numpy as np
import pylab as pl
from NeuroTools import analysis, signals

from classes.static import MeanfieldParameters as mpr
from classes.static import NeuronParameters as npr
from tools.analysis import process_spikes

# GENERAL PREPARATIONS ----------------------- #
nest.ResetKernel()
nest.SetKernelStatus({'print_time': True, 'local_num_threads': 1})
nest.set_verbosity('M_ERROR')
try:
    nest.Install('mymodule')
except Exception, e:
    pass

# switch the tsodyks synapse, used for debugging
tsodyks_model = "tsodyks_mongillo_synapse"

gen_params_baserate = {
    mpr.W_J: 1.,  # bifurcation parameter for weights; 1. is a flat weight profile
    "sig": 18. / 360. * np.pi * 2.,  # standard deviation of weight

    # ANALYSIS / RUNNING
    "tmax": 3000.,

    "show_weights": False,
    "show_results": False,
    "show_spikes": False,

    "bump_shape_enddist": 100,  # end of bump averaging (reverse from tmax)
    "bump_i_rates_enddist": 300,  # end of bump averaging (reverse from tmax)
    "bump_i_rates_time": 750,  # end of bump averaging (reverse from tmax)

    "tstart_spikes_base": 1500.,  # when to start averaging base rate (from there, until signal_start)

    # GENERAL MODEL PARAMETERS
    mpr.W_NOISE: 0.,
    mpr.EL_NOISE: 0.,
    mpr.P_EE: 1.,
    mpr.EL_NOISE_ARRAY: None,
    mpr.EL_SEED: None,

    # start signal
    "sig_rate": 3000.,  # rate of external input Poisson spikes
    "sig_start": 500.,  # start time of external input
    "sig_len": 1000.,  # duration of external input
    "sig_fade": True,  # halves the external input rate during the second half of sig_len
    "sig_width": .2,  # width of external input (proportion of excitatory neurons)
    "sig_center": .5,  # center of external input (proportion of excitatory neurons)
    "sig_weight": .5,  # weight of connections onto AMPA synapses

    "return_w_mat": False,  # returns the weights without running simulation
    "do_bistab_pert": False,  # bistability perturbation to allow easier settling into base-state
    "base_seed_run": None,  # base seed for run
    "base_seed": 142152521,  # base seed for connectivity
}

# params for simbump
gen_params_bump = gen_params_baserate.copy()
gen_params_bump.update({
    "tmax": 15000.,
    # bump average
    "tstart_spikes_base": 300.  # when to start averaging base rate (until signal),
})

# params for simdrift
gen_params_drift = gen_params_baserate.copy()
gen_params_drift.update({

    "tmax": 8000.,
    # bump average
    "tstart_spikes_base": 300.,  # when to start averaging base rate (until signal),
})


def get_cue_neurons(width, center, n_neurons):
    """
    gets stimulus neuron numbers given width and center of stimulus

    :param width: width of stimulus (proportion of n_neurons)
    :param center: center of stimulus (proportion of n_neurons)
    :param n_neurons: number of neurons
    :return: numbers of neurons to receive stimulus
    """
    # get left and right boundaries modulo the number of neurons
    left = int((center - width / 2) * n_neurons % n_neurons)
    right = int((center + width / 2) * n_neurons % n_neurons)
    out = []

    if left > right:
        out = out + range(0, right)
        out = out + range(left, n_neurons)
    else:
        out = range(left, right)
    return out


class NestSimulator:
    """ Handles all nest simulations """

    def __init__(self, handler, cores=2):
        """
        :param handler: meanfield.error_functions.MeanfieldHandler instance
        :param cores: cores to use for simulation
        """
        self.params = {}
        self.params["gen"] = None
        self.handler = handler

        self.spikes = {}
        self.ex_cells = None
        self.in_cells = None
        self.cores = cores

    def reset(self):
        self.spikes = {}

    def run(self, do_only_setup=False):
        """Run a simulation in nest

        :param do_only_setup: only do network setup, do not simulation.
            This can be useul for debugging or to run theory predictions
        :return:
        """
        if self.params["gen"] == None:
            self.params["gen"] = gen_params_baserate

        print "Using network: %s (E: %s, I: %s)" % (
        self.handler.network, self.handler.network.e.neuron_model.name, self.handler.network.i.neuron_model.name)

        self.reset()
        nest.ResetKernel()

        n_e = int(self.handler.p_e[npr.C_E])
        n_i = int(self.handler.p_i[npr.C_I])

        n_threads = self.cores
        nest.SetKernelStatus({'print_time': True, 'local_num_threads': n_threads})
        N_vp = nest.GetKernelStatus(['total_num_virtual_procs'])[0]

        # create models     

        # excitatory neurons

        neuron_params = {
            'V_reset': self.handler.p_e[npr.VRES],
            # hardcoded in the model 'V_th': -50mV,
            't_ref': self.handler.p_e[npr.TAU_RP],
            'E_L': self.handler.p_e[npr.VL],
            'g_L': self.handler.p_e[npr.GM],
            'C_m': self.handler.p_e[npr.CM],

            'AMPA_g_peak': self.handler.p_e[npr.G_AMPA],
            'AMPA_Tau_1': self.handler.p_e[npr.TAU_AMPA],

            'GABA_g_peak': self.handler.p_e[npr.G_GABA],
            'GABA_Tau_1': self.handler.p_e[npr.TAU_GABA],

            # WE RESCALE THE NMDA CONDUCTANCE BY THE CONNECTIVITY, multiplies everything and leaves weights as they are
            'NMDA_g_peak': self.handler.p_e[npr.G_NMDA] / self.params["gen"][mpr.P_EE],
            'NMDA_Tau_1': self.handler.p_e[npr.TAU_NMDA]
        }

        nest.CopyModel(self.handler.network.e.neuron_model.name, 'ht_neuron_ex', params=neuron_params)

        # inhibitory neurons

        neuron_params = {
            'V_reset': self.handler.p_i[npr.VRES],
            # hardcoded in the model 'V_th': -50mV,
            't_ref': self.handler.p_i[npr.TAU_RP],
            'E_L': self.handler.p_i[npr.VL],
            'g_L': self.handler.p_i[npr.GM],
            'C_m': self.handler.p_i[npr.CM],

            'AMPA_g_peak': self.handler.p_i[npr.G_AMPA],
            'AMPA_Tau_1': self.handler.p_i[npr.TAU_AMPA],

            'GABA_Tau_1': self.handler.p_i[npr.TAU_GABA],
            'GABA_g_peak': self.handler.p_i[npr.G_GABA],

            'NMDA_g_peak': self.handler.p_i[npr.G_NMDA],
            'NMDA_Tau_1': self.handler.p_i[npr.TAU_NMDA],

        }

        nest.CopyModel(self.handler.network.i.neuron_model.name, 'ht_neuron_in', params=neuron_params)

        # get receptor ID information, so we can connect to the
        # different synapses
        ex_receptors = nest.GetDefaults('ht_neuron_ex')['receptor_types']
        in_receptors = nest.GetDefaults('ht_neuron_in')['receptor_types']

        # CREATE NEURON POPULATIONS ----------------------- #
        ex_cells = tp.CreateLayer({
            'rows': 1,
            'columns': n_e,
            'elements': 'ht_neuron_ex',
            'extent': [2. * np.pi, 1.],
            'edge_wrap': True
        })

        in_cells = tp.CreateLayer({
            'rows': 1,
            'columns': n_i,
            'elements': 'ht_neuron_in',
            'extent': [2. * np.pi, 1.],
            'edge_wrap': True
        })

        nu_ext_e = self.handler.p_e[npr.NU_EXT] * 1e3
        nu_ext_i = self.handler.p_i[npr.NU_EXT] * 1e3

        # create Poisson generators
        ex_noise = nest.Create('poisson_generator', 1, params={'rate': nu_ext_e * self.handler.p_e[npr.C_EXT]})
        in_noise = nest.Create('poisson_generator', 1, params={'rate': nu_ext_i * self.handler.p_i[npr.C_EXT]})

        nest.OneToOneConnect(ex_noise * n_e, nest.GetNodes(ex_cells)[0],
                             params={'weight': 1., 'receptor_type': ex_receptors['AMPA']}, model='static_synapse')
        nest.OneToOneConnect(in_noise * n_i, nest.GetNodes(in_cells)[0],
                             params={'weight': 1., 'receptor_type': in_receptors['AMPA']}, model='static_synapse')

        # signal initiators
        if self.params["gen"]["sig_len"] > 0:

            start = self.params["gen"]["sig_start"]
            stop = self.params["gen"]["sig_start"] + self.params["gen"]["sig_len"]

            print "Creating poisson signal from t=%i to t=%i" % (
            self.params["gen"]["sig_start"], self.params["gen"]["sig_start"] + self.params["gen"]["sig_len"])

            sig_len = int(self.params["gen"]["sig_width"] * n_e)
            sig_range = get_cue_neurons(self.params["gen"]["sig_width"], self.params["gen"]["sig_center"], n_e)
            print "Signal (width %g) will be sent to %i neurons centered at %g" % (
            self.params["gen"]["sig_width"], len(sig_range), self.params["gen"]["sig_center"])
            ex_neurons = nest.GetNodes(ex_cells)[0]

            # WATCH OUT, setting range to length. sometimes due to rounding this
            # was not equal
            sig_range = sig_range[:sig_len]

            if self.params["gen"]["sig_fade"]:
                print "Using fading poisson signal."
                ex_signal = nest.Create('poisson_generator', 1,
                                        params={'rate': self.params["gen"]["sig_rate"], 'start': start,
                                                'stop': start + self.params["gen"]["sig_len"] / 2.})
                nest.OneToOneConnect(ex_signal * sig_len, [ex_neurons[i] for i in sig_range],
                                     params={'weight': self.params["gen"]["sig_weight"],
                                             'receptor_type': ex_receptors['AMPA']}, model='static_synapse')
                ex_signal = nest.Create('poisson_generator', 1, params={'rate': 0.5 * self.params["gen"]["sig_rate"],
                                                                        'start': start + self.params["gen"][
                                                                            "sig_len"] / 2., 'stop': stop})
                nest.OneToOneConnect(ex_signal * sig_len, [ex_neurons[i] for i in sig_range],
                                     params={'weight': self.params["gen"]["sig_weight"],
                                             'receptor_type': ex_receptors['AMPA']}, model='static_synapse')
            else:
                print "Using non-fading poisson signal."
                ex_signal = nest.Create('poisson_generator', 1,
                                        params={'rate': self.params["gen"]["sig_rate"], 'start': start,
                                                'stop': start + self.params["gen"]["sig_len"]})
                nest.OneToOneConnect(ex_signal * sig_len, [ex_neurons[i] for i in sig_range],
                                     params={'weight': self.params["gen"]["sig_weight"],
                                             'receptor_type': ex_receptors['AMPA']}, model='static_synapse')

        # STP CONNECTIONS
        if self.handler.network.i.nmda_synapse.is_stp:
            print "Using E->I stp (%s) with params: U=%f tau_r=%f rau_f=%f" % (
            tsodyks_model, self.handler.network.i.nmda_synapse.U, self.handler.network.i.nmda_synapse.tau_r,
            self.handler.network.i.nmda_synapse.tau_f)
            nest.CopyModel(tsodyks_model, 'NMDA_EI', {
                'U': self.handler.network.i.nmda_synapse.U,
                'tau_fac': self.handler.network.i.nmda_synapse.tau_f,
                'tau_rec': self.handler.network.i.nmda_synapse.tau_r,
                'receptor_type': in_receptors['NMDA']
            })
        # REGULAR CONNECTIONS
        else:
            print "Using E->I static synapses"
            nest.CopyModel('static_synapse', 'NMDA_EI', {'receptor_type': in_receptors['NMDA']})

        # E->I all to all connected
        conndict_EI = {"connection_type": "divergent",
                       "synapse_model": "NMDA_EI",
                       "mask": {"grid": {"rows": 1, "columns": n_i}},
                       "weights": 1.,
                       "kernel": 1.,
                       "allow_autapses": False}
        tp.ConnectLayers(ex_cells, in_cells, conndict_EI)

        # I->I all to all connected
        nest.CopyModel('static_synapse', 'GABA_II', {'receptor_type': in_receptors['GABA']})
        conndict_II = {"connection_type": "divergent",
                       "synapse_model": "GABA_II",
                       "mask": {"grid": {"rows": 1, "columns": n_i}},
                       "weights": 1.,
                       "kernel": 1.,
                       "allow_autapses": False}
        tp.ConnectLayers(in_cells, in_cells, conndict_II)

        # I->E all to all connected
        nest.CopyModel('static_synapse', 'GABA_IE', {'receptor_type': ex_receptors['GABA']})
        conndict_IE = {"connection_type": "divergent",
                       "synapse_model": "GABA_IE",
                       "mask": {"grid": {"rows": 1, "columns": n_e}},
                       "weights": 1.,
                       "kernel": 1.,
                       "allow_autapses": False}
        tp.ConnectLayers(in_cells, ex_cells, conndict_IE)

        nest.SetKernelStatus({'rng_seeds': range(self.params["gen"]["base_seed"] + N_vp + 1,
                                                 self.params["gen"]["base_seed"] + 2 * N_vp + 1)})
        nest.SetKernelStatus({'grng_seed': self.params["gen"]["base_seed"] + N_vp})
        print "Seed for connectivity (should stay the same): ", nest.GetKernelStatus()["rng_seeds"]

        # noise on membrane parameters
        if self.params["gen"][mpr.EL_NOISE] > 0. or not self.params["gen"][mpr.EL_NOISE_ARRAY] == None:
            print "Using nonzero leak noise:"
            if not self.params["gen"][mpr.EL_NOISE_ARRAY] == None:
                leak_variance = self.params["gen"][mpr.EL_NOISE_ARRAY]
                print "Setting variable leak currents with given noise array, mean: %.2f, std: %.2f" % (
                np.mean(leak_variance), np.std(leak_variance))
                assert len(leak_variance) == n_e, "noise array has wrong size"
            else:
                np.random.seed(self.params["gen"]["base_seed"] + N_vp)
                self.params["gen"][mpr.EL_SEED] = self.params["gen"]["base_seed"] + N_vp
                leak_variance = np.random.normal(self.handler.p_e[npr.VL], self.params["gen"][mpr.EL_NOISE], n_e)
                self.params["gen"][mpr.EL_NOISE_ARRAY] = leak_variance
                print "Setting variable leak currents: mean %.1f, std %.1f (seed=%i)" % (
                np.mean(leak_variance), np.std(leak_variance), self.params["gen"][mpr.EL_SEED])
            print leak_variance[:10], "..."
            nest.SetStatus(nest.GetNodes(ex_cells)[0], "E_L", leak_variance)

        # EE CONNECTIONS ##################################

        # TO SAVE CONNECTIONS WE DO AMPA DIST DEPENDENT CONNS WITH THE SAME SEED
        if self.params["gen"]["return_w_mat"]:
            print "Using fake synapses for weight matrix output"
            nest.CopyModel('static_synapse', 'NMDA_EE', {'receptor_type': ex_receptors['AMPA']})
        else:
            # STP CONNECTIONS
            if self.handler.network.e.nmda_synapse.is_stp:
                print "Using E->E stp (%s) with params: U=%f tau_r=%f rau_f=%f" % (
                tsodyks_model, self.handler.network.e.nmda_synapse.U, self.handler.network.e.nmda_synapse.tau_r,
                self.handler.network.e.nmda_synapse.tau_f)
                nest.CopyModel(tsodyks_model, 'NMDA_EE', {
                    'U': self.handler.network.e.nmda_synapse.U,
                    'tau_fac': self.handler.network.e.nmda_synapse.tau_f,
                    'tau_rec': self.handler.network.e.nmda_synapse.tau_r,
                    'receptor_type': ex_receptors['NMDA']
                })
            # REGULAR CONNECTIONS
            else:
                print "Using static synapses"
                nest.CopyModel('static_synapse', 'NMDA_EE', {'receptor_type': ex_receptors['NMDA']})

        if self.params["gen"][mpr.P_EE] < 1.:
            print "Using connecticity p_ee: %f" % self.params["gen"][mpr.P_EE]

        if self.params["gen"][mpr.W_NOISE] > 0.:
            print "Using w noise with level (multiplied after by w1:%f): %f" % (
            self.params["gen"][mpr.W_NOISE], self.handler.p_e[npr.W_1])

        print "Using w parameters - w_j: %.2f --> w_0: %.2f, w_1: %.2f, w_sigma: %.2f" % (
        self.handler.p_mf[mpr.W_J], self.handler.p_e[npr.W_0], self.handler.p_e[npr.W_1], self.handler.p_e[npr.W_SIGMA])
        conndict_EE = {
            "connection_type": "divergent",
            "synapse_model": "NMDA_EE",
            "mask": {"grid": {"rows": 1, "columns": n_e}},
            "weights": {"gaussian_noisy": {"c": self.handler.p_e[npr.W_0], "p_center": self.handler.p_e[npr.W_1],
                                           "sigma": self.handler.p_e[npr.W_SIGMA],
                                           "sigma_noise": self.params["gen"][mpr.W_NOISE] * self.handler.p_e[npr.W_1]}},
            "kernel": self.params["gen"][mpr.P_EE],
            "allow_autapses": False
        }
        if self.params["gen"][mpr.P_EE] < 1.:
            conndict_EE["kernel"] = self.params["gen"][mpr.P_EE]

        nest.SetKernelStatus({'rng_seeds': range(self.params["gen"]["base_seed"] + N_vp + 1,
                                                 self.params["gen"]["base_seed"] + 2 * N_vp + 1)})
        nest.SetKernelStatus({'grng_seed': self.params["gen"]["base_seed"] + N_vp})

        tp.ConnectLayers(ex_cells, ex_cells, conndict_EE)

        if self.params["gen"]["return_w_mat"]:

            w_mat = np.zeros((n_e, n_e))
            cids = nest.GetNodes(ex_cells)[0]
            conns = nest.FindConnections(cids, synapse_type="NMDA_EE")
            stats = nest.GetStatus(conns)
            id_min = min(cids)
            id_max = max(cids)
            assert id_max - id_min == n_e - 1, "something is wrong with weight estimation"
            for stat in stats:
                pre = stat['source'] - id_min
                post = stat['target'] - id_min
                weight = stat['weight']
                w_mat[post, pre] = weight

            return w_mat

        if self.params["gen"]["show_weights"]:

            target = []
            weight = []

            a = nest.GetNodes(ex_cells)[0][0]
            b = nest.FindConnections([a])
            for con in b:

                aha = nest.GetStatus([con])[0]
                if aha['target'] < n_e:
                    target.append(aha['target'])
                    weight.append(aha['weight'])

            print "Total weight is: %f" % (np.sum(weight) / float(n_e))
            pl.figure()
            pl.scatter((np.array(target) / float(n_e) * 360.), weight)
            tmp = np.zeros(n_e)
            for i in range(n_e):
                tmp[i] = self.handler.p_e[npr.W_0] + self.handler.p_e[npr.W_1] * np.exp(
                    -((i - n_e - 2) / (1. * n_e / (2 * np.pi))) ** 2 / 2. / (self.handler.p_e[npr.W_SIGMA]) ** 2)

            pl.plot((np.arange(n_e) / float(n_e) * 360.), tmp, 'r')
            pl.show()
            print "Targets (should stay constant)"
            print target

        if self.params["gen"]["base_seed_run"] is not None:
            print "Using preset seed for run: %i" % self.params["gen"]["base_seed_run"]
            msd = self.params["gen"]["base_seed_run"]
        else:
            print "Using random seed for run: seeding by time."
            np.random.seed(int(time.time()))
            msd = np.random.randint(100000000000)

        self.params["gen"]["base_seed_run"] = msd
        nest.SetKernelStatus({'rng_seeds': range(msd + N_vp + 1, msd + 2 * N_vp + 1)})
        nest.SetKernelStatus({'grng_seed': msd + N_vp})

        # READOUTS ----------------------- #
        # spike detectors
        ex_spikes = nest.Create("spike_detector")
        nest.SetStatus(ex_spikes, [{"label": "ex", "withtime": True, "withgid": True}])
        nest.ConvergentConnect(nest.GetNodes(ex_cells)[0], ex_spikes, model="static_synapse")

        in_spikes = nest.Create("spike_detector")
        nest.SetStatus(in_spikes, [{"label": "in", "withtime": True, "withgid": True}])
        nest.ConvergentConnect(nest.GetNodes(in_cells)[0], in_spikes, model="static_synapse")

        spikes = nest.Create("spike_detector")
        nest.SetStatus(spikes, [{"withtime": True, "withgid": True}])
        nest.ConvergentConnect(nest.GetNodes(ex_cells)[0], spikes, model="static_synapse")
        nest.ConvergentConnect(nest.GetNodes(in_cells)[0], spikes, model="static_synapse")

        # write cell ids to object
        self.ex_cells = np.array(nest.GetNodes(ex_cells)[0])
        self.in_cells = np.array(nest.GetNodes(in_cells)[0])

        if do_only_setup:
            return True

        # SIMULATE ----------------------- #
        print "Seeds for run: ", nest.GetKernelStatus()["rng_seeds"]
        print "Running %.1f ms" % self.params["gen"]["tmax"]
        starttime = time.time()

        # COOLDOWN FIRST
        ###########
        if self.params["gen"]["do_bistab_pert"]:
            ex_signal_init = nest.Create('poisson_generator', 1, params={'rate': 1500., 'start': 1000., 'stop': 1500.})
            sig_range = []
            d = 4
            for i in range(d):
                sig_range += get_cue_neurons(.025, i / float(d), n_e)
            sig_len = len(sig_range)
            nest.OneToOneConnect(ex_signal_init * sig_len, [ex_neurons[i] for i in sig_range],
                                 params={'weight': 2., 'receptor_type': ex_receptors['AMPA']}, model='static_synapse')

        # FULL SIM
        ###########

        nest.Simulate(self.params["gen"]["tmax"])
        print "Done (%.2fs): " % ((time.time() - starttime)),

        # Output --------------------------- #
        ex_ev = nest.GetStatus(ex_spikes, "events")[0]
        self.spikes["e"] = process_spikes(ex_ev, self.params["gen"]["tmax"])
        in_ev = nest.GetStatus(in_spikes, "events")[0]
        self.spikes["i"] = process_spikes(in_ev, self.params["gen"]["tmax"])
        all_ev = nest.GetStatus(spikes, "events")[0]
        self.spikes["all"] = process_spikes(all_ev, self.params["gen"]["tmax"])

        window = 100.
        ex_spiketrains = self.spikes["e"]["spike_trains"].spiketrains
        kernel, norm, m_idx = analysis.make_kernel('exp', window, 1)
        act_len = len(self.ex_cells)
        act = np.zeros((act_len, int(self.params["gen"]["tmax"])))

        for l, spktr in enumerate(ex_spiketrains):
            nrnid = np.where(self.ex_cells == spktr)[0]
            if len(nrnid) == 0:
                continue
            taxis, raxis = ex_spiketrains[spktr].instantaneous_rate(1, kernel, norm, m_idx, trim=False)
            act[nrnid, :] = raxis

        # end of bump activity only

        bump_t_end = self.params["gen"]["tmax"] - self.params["gen"]["bump_shape_enddist"]
        bump_t_start = self.params["gen"]["sig_start"] + self.params["gen"][
            "sig_len"] + 1000.  # 1s of buffer between start and avg
        if bump_t_end <= bump_t_start:
            bump_t_start_old = bump_t_start
            bump_t_start = self.params["gen"]["sig_start"] + self.params["gen"]["sig_len"]
            print "Run is not long enough (t_max: %.1f, bump_t_start: %.1f, bump_t_end: %.1f). Setting bump_t_start = %.1f" % (
            self.params["gen"]["tmax"], bump_t_start_old, bump_t_end, bump_t_start)

        dir_vec = np.exp(2. * 1.j * np.pi * np.arange(act_len, dtype=float) / float(act_len))
        dirs = (np.angle(np.dot(np.transpose(act), dir_vec))) / np.pi * act_len / 2. % act_len

        act_bump = act[:, int(bump_t_start):int(bump_t_end)]
        dirs_bump = dirs[int(bump_t_start):int(bump_t_end)]

        # rectify bump to center and fit gauss curve
        shift_n = (len(act_bump) / 2. - dirs_bump).round().astype(int)
        out_means = []
        points = np.arange(0, act_bump.shape[1], 20)
        for j in points:
            out_means.append(np.roll(act_bump[:, j], shift_n[j], 0))
        out_means = np.array(out_means).T

        if self.params["gen"]["show_results"] or self.params["gen"]["show_spikes"]:
            nest.raster_plot.from_device(ex_spikes, hist=True)
            pl.title("Exctitatory Population")
            pl.savefig('lastrun_spikes_e.pdf')
            nest.raster_plot.from_device(in_spikes, hist=True)
            pl.title("Inhibitory Population")
            pl.savefig('lastrun_spikes_i.pdf')

        return {
            "shape_mean": np.mean(out_means, 1),
            "shape_std": np.std(out_means, 1),
            "shape": out_means,
            "dirs": dirs,
            "bump_rate": act,
            "spikes": self.spikes,
            "pop_rate": self.get_rates(bump_t_start, bump_t_end)
        }

    def set_paramset(self, todo):
        """Set parameterset to this instance

        :param todo: "bump", "base", or "drift"
        :return:
        """
        if todo == "bump":
            self.params["gen"] = gen_params_bump.copy()
        elif todo == "base":
            self.params["gen"] = gen_params_baserate.copy()
        elif todo == "drift":
            self.params["gen"] = gen_params_drift.copy()
        else:
            raise ValueError("Unknown todo")

    def set_params(self, dic, **kwargs):
        """Directly set a parameter"""
        self.params[dic].update(kwargs)

    def get_rates(self, t_start, t_stop):
        """Calculates and returns mean rates for e and i populations.

        This uses the internally saved spikes which are available after calling run()

        :param t_start: start time
        :param t_stop: end time
        :return: dictionary of two rates
        """
        print "Returning rates from %f to %f" % (t_start, t_stop)
        if t_start >= t_stop:
            raise ValueError("t_stop should be bigger than t_start.")
        ex_spiketrains_base = signals.SpikeList(self.spikes["e"]["spike_list"], self.spikes["e"]["ids"],
                                                t_start=t_start, t_stop=t_stop)
        ex_base = ex_spiketrains_base.mean_rate()
        in_spiketrains_base = signals.SpikeList(self.spikes["i"]["spike_list"], self.spikes["i"]["ids"],
                                                t_start=t_start, t_stop=t_stop)
        in_base = in_spiketrains_base.mean_rate()
        return {"e": ex_base, "i": in_base}

    def get_baserates(self):
        """Returns base rates of time period before external input stimulus was applied.

        Evaluating this only makes sense if you simulated am extended period before the network
        received a stimulus.

        :return: dictionary of two rates
        """
        base_end = self.params["gen"]["tmax"]

        # should we be simulating a bump, get rates before the bump
        if self.params["gen"]["sig_len"] > 0:
            base_end = self.params["gen"]["sig_start"]

        base_start = self.params["gen"]["tstart_spikes_base"]
        return self.get_rates(base_start, base_end)
