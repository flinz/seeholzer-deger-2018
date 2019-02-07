import matplotlib
matplotlib.use('pdf')
import tools.db
import numpy as np


def test_network_sim():
    """Test DB setup by getting and printing all networks"""
    networks = tools.db.print_networks([4,5,6])
    np.testing.assert_equal(len(networks), 34)