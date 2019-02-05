import matplotlib

matplotlib.use('pdf')
import tools.db
import meanfield.error_functions as erfs
import classes.base as cb
from classes.static import MeanfieldParameters as mpr
import numpy as np

def test_network_sim():
    networks = tools.db.print_networks([4,5,6])
    np.testing.assert_equal(len(networks), 34)