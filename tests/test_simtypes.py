import numpy as np
import classes.base as cb
import meanfield.error_functions as erfs
import tools.db
from classes.static import MeanfieldParameters as mpr


def test_simtype_bump():
    """Test the SimBump simtype"""

    session, net = tools.db.get_network(45)

    # get meanfield handler (abstracts network parameter handling) and set connectivity parameter
    handler = erfs.MeanfieldHandler(net)

    # get simulation type
    btype = session.query(cb.SimType).filter(cb.SimType.name == "bump").one()

    # get simulation wrapper
    # set to_memory = False to save simulation to hdf5 file
    wrapper = btype.get_wrapper(to_memory=True, session=session)

    tmax = 2000.

    # parameters used for this example
    gen_params = {
        # noise parameters
        mpr.EL_NOISE: 0.,
        mpr.P_EE: 1.,  # 1mV of noise on the leak reversal potentials
        mpr.W_NOISE: 0.,

        # run parameters
        "cores": 4,  # cores to use for sim
        "tmax": tmax,  # runtime
        "base_seed": 1,  # seed for connectivity noise
        "base_seed_run": 1,  # seed for run noise

        "sig_center": 0.3,  # center of signal
    }

    # run the simulation
    wrapper.run(handler, **gen_params)

    assert wrapper.data_file['dirs'].shape == (int(tmax),)
    np.testing.assert_almost_equal(
        wrapper.data_file['dirs'][-1],
        239.6555633544922)


def test_simtype_drift():
    """Test the SimDrift simtype"""

    session, net = tools.db.get_network(45)

    # get meanfield handler (abstracts network parameter handling) and set connectivity parameter
    handler = erfs.MeanfieldHandler(net)

    # get simulation type
    btype = session.query(cb.SimType).filter(cb.SimType.name == "drift").one()

    # get simulation wrapper
    # set to_memory = False to save simulation to hdf5 file
    wrapper = btype.get_wrapper(to_memory=True, session=session)

    reps = 1
    initials = 2
    tmax = 2000.

    # parameters used for this example
    gen_params = {
        # noise parameters
        mpr.EL_NOISE: 0.,
        mpr.P_EE: 1.,  # 1mV of noise on the leak reversal potentials
        mpr.W_NOISE: 0.,

        # run parameters
        "cores": 4,  # cores to use for sim
        "tmax": tmax,  # runtime
        "base_seed": 1,  # seed for connectivity noise
        "base_seed_run": 1,  # seed for run noise

        "reps": reps,  # do only 1 repetition
        "initials": initials
    }

    # run the simulation
    wrapper.run(handler, **gen_params)
    assert wrapper.data_file['dirs'].shape == (initials, reps, int(tmax))
    np.testing.assert_array_almost_equal(
        wrapper.data_file['dirs'][:, :, -1].flatten(),
        [6.210974216461182, 396.7076110839844])
