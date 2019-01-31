import matplotlib
matplotlib.use('pdf')
import tools.db
import meanfield.error_functions as erfs
from classes.static import MeanfieldParameters as mpr
import classes.base as cb
from matplotlib import pyplot as pl

session, net = tools.db.get_network(45)

# get meanfield handler (abstracts network parameter handling) and set connectivity parameter
handler = erfs.MeanfieldHandler(net)

# get simulation type
btype = session.query(cb.SimType).filter(cb.SimType.name == "drift").one()

# get simulation wrapper
# set to_memory = False to save simulation to hdf5 file
wrapper = btype.get_wrapper(to_memory=True, session=session)

# parameters for simulations as used in the publication
gen_params = {
    "reps": 20,  # 20 repetitions per initial position
    "tmax": 8000.,  # 8s of simulation time
    "initials": 10,  # 10 initial positions
}

# parameters used for this example
gen_params.update({
    # noise parameters
    mpr.EL_NOISE: 1.,
    mpr.P_EE: 1.,  # 1mV of noise on the leak reversal potentials
    mpr.W_NOISE: 0.,

    # run parameters
    "cores": 4,  # cores to use for sim
    "tmax": 8000.,  # runtime
    "base_seed": 1,  # seed for connectivity noise
    "base_seed_run": 1,  # seed for run noise

    "reps": 5  # do only 5 repetitions
})

# run the simulation
wrapper.run(handler, **gen_params)

wrapper.plot_drift()
pl.savefig('out_drift.png', dpi=300)
