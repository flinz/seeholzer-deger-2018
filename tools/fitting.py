import numpy as np
import pylab as pl
from scipy.optimize import minimize

from classes.static import MeanfieldParameters as mpr


def fit_gauss_curve(ydata, plot=False, ax=None, dots=False, label=None, **kwargs):
    """Fit a generalized gaussian to given uniformly spaced datapoints on [-pi,pi)

    :param ydata: data points of measurements uniformly spaced on [-pi,pi)
    :param plot: do a plot of result
    :param ax: matplotlib axis to plot to
    :param dots: scatterplot the original data
    :param label: plot label
    :param kwargs: arguments to pass to plotting
    :return: dict of coefficients and the fitted function
    """
    def gauss(x, *p):
        g_0, g_1, center, sigma, r = p
        if sigma <= 0:
            return g_0 + g_1 * np.ones(x.shape)
        x_new = x - center
        return g_0 + g_1 * np.exp(-(abs(x_new) / sigma) ** r)

    xlen = len(ydata)
    xdata = np.arange(0, xlen) * (2. * np.pi) / float(xlen) - np.pi

    p0 = [10., 0., 0., 18. / 360. * 2. * np.pi, 2.]

    err = lambda p: np.mean((ydata - gauss(xdata, *p)) ** 2)

    try:
        for i in range(1):
            min_err = 1e10
            min_sol = None
            p0_use = p0
            p0_use[1] = np.random.rand(1) * 10.

            p_opt = minimize(
                err,  # minimize wrt to the noisy data
                p0_use,
                bounds=[(0, None), (0, None), (-np.pi, np.pi), (.2, np.pi), (1.5, 30.)],  # set the bounds
                method="L-BFGS-B"  # this method supports bounds
            )
            if p_opt.success and p_opt.fun < min_err:
                min_err = p_opt.fun
                min_sol = p_opt.x

        coeff = min_sol
    except Exception, e:
        print e
        print "not converged"
        coeff = p0

    # get rid of fitting fragments, only influences the base (flat) state
    if coeff[1] / coeff[0] < 1.5:
        coeff[0] = np.mean(ydata)
        coeff[1] = 0.
    fittedcurve = gauss(xdata, *coeff)

    if plot:
        if ax is None:
            pl.figure()
            ax = pl.subplot(111)

        if dots:
            ax.scatter(xdata, ydata, c="#cccccc", s=1., alpha=.5)
        ax.plot(xdata, fittedcurve, label=label, **kwargs)

    return {mpr.G_0: coeff[0], mpr.G_1: coeff[1], mpr.G_SIGMA: coeff[3], mpr.G_R: coeff[4],
            "func": lambda x: gauss(x, *coeff), "x": xdata}
