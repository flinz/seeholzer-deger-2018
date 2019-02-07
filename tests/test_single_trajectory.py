import matplotlib
matplotlib.use('pdf')
import numpy as np
from simulation.nest_simulation import NestSimulator
import meanfield.error_functions as erfs
import tools.db
from classes.static import MeanfieldParameters as mpr


def test_network_sim():
    """
    Test single trajectory without using simtypes
    """

    session, net = tools.db.get_network(45)
    handler = erfs.MeanfieldHandler(net)

    # somewhat legacy parameters for this test run
    std_pars = {
        "show_results": False,
        "id": net.id,
        "sig_start": 0.,

        # noise parameters
        mpr.W_NOISE: 0.,  # sigma_w
        mpr.EL_NOISE: 0.,  # = sigma_L
        mpr.P_EE: 0.5,   # recurrent connectivity p

        "cores": 4,  # cores to use for sim
        "sig_center": 0.3,
        "sig_len": 1000.,
        "sig_weight": .4,
        "sig_width": .2,
        "sig_rate": 2000.,
        "base_seed_run": 25757712619,
        "base_seed": 29529492,

        # runtime
        "tmax": 2000.,
    }

    handler.set_states({ mpr.W_J: net.e.paramset.parameters['w_0'] + net.e.paramset.parameters['w_1']})
    sim = NestSimulator(handler, cores=4)
    sim.set_paramset("bump")
    sim.set_params("gen", **std_pars)

    sim_ret = sim.run()

    print "E-rate after offset: ", sim_ret["pop_rate"]["e"]
    print "I-rate after offset: ", sim_ret["pop_rate"]["i"]
    print "Averaged maximal E-rate after offset: ", max(sim_ret["shape_mean"])

    np.testing.assert_almost_equal(sim_ret["pop_rate"]["e"], 6.07852487681, decimal=10)
    np.testing.assert_almost_equal(sim_ret["pop_rate"]["i"], 4.81666666667, decimal=10)
    np.testing.assert_almost_equal(max(sim_ret["shape_mean"]), 42.245566508, decimal=10)
