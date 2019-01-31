import matplotlib.pyplot as pl
import classes.base as cb
import meanfield.error_functions as erfs
import tools.db
from classes.static import MeanfieldParameters as mpr

# to show all network ids: tools.db.print_networks()
session, net = tools.db.get_network(48)

# get meanfield handler (abstracts network parameter handling) and set connectivity parameter
handler = erfs.MeanfieldHandler(net)

# get simulation type
btype = session.query(cb.SimType).filter(cb.SimType.name == "bump").one()

# get simulation wrapper
# set to_memory = False to save simulation to hdf5 file
wrapper = btype.get_wrapper(to_memory=True, session=session)

# parameters for simulation
sim_params = {
    "cores": 4,  # cores to use for sim
    "sig_center": 0.3,  # center of signal
    "tmax": 8000.,  # runtime
    "base_seed": 1,  # seed for connectivity noise
    "base_seed_run": 1,  # seed for run noise

    # network noise parameters
    mpr.W_NOISE: 0.,  # sigma_w
    mpr.EL_NOISE: 0.,  # = sigma_L
    mpr.P_EE: 0.5,  # recurrent connectivity p
}

wrapper.run(handler, **sim_params)

# plot mean
wrapper.plot_mean()
pl.savefig('out_mean.pdf')

# plot rates
wrapper.plot_rates()
pl.savefig('out_rates.pdf')