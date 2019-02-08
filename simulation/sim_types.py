from classes.static import MeanfieldParameters as mpr
from classes.static import NeuronParameters as npr
import h5py
import time
import dill as pickle
import pickle as realpickle
import numpy as np
from classes.base import SimType, SimFile
from sqlalchemy.orm.exc import NoResultFound
from tools.db import similar_or_new_net
import matplotlib.pyplot as pl


class SimWrapper(object):
    """Wraps h5py HDF5 DataFiles and handles their creation and opening.

    Wraps either classes.base.SimFiles (which contain h5py DataFiles) or h5py DataFiles directly
    [provided by either their data_filename (path) or a DataFile instance].
    """

    def __init__(self, sim_file=None, sim_type=None, data_filename=None, data_file=None, session=None, to_memory=False):
        """

        :param sim_file: classes.base.SimFile instance
        :param sim_type: classes.base.SimType instance
        :param data_filename: filename (path) to the datafile
        :param data_file: DataFile object
        :param session: SqlAlchemy session
        :param to_memory: Save to memory instead of physical hdf5 file
        """
        self.data_filename = data_filename

        self.sim_file = sim_file
        self.data_file = data_file

        self.to_memory = to_memory

        self.sim = None
        self.session = session

        if sim_file:
            self.sim_type = sim_file.sim_type
        elif sim_type:
            self.sim_type = sim_type
        else:
            ValueError, "Can't find sim type"

    """Getters and setters for data_file"""

    @property
    def data_file(self):
        return self.get_data_file()

    @data_file.setter
    def data_file(self, value):
        if value and not value.__class__ == h5py._hl.files.File:
            raise ValueError, "%s is not a h5py file." % value
        self.data_file_ = value
        if value:
            self.data_filename = value.filename

    """Getters and setters for sim_file"""

    @property
    def sim_file(self):
        return self.get_sim_file()

    @sim_file.setter
    def sim_file(self, value):
        if value and not value.__class__ == SimFile:
            raise ValueError, "%s is not a SimFile." % value
        self.sim_file_ = value

    def get_data_file(self):
        """Returns (in this order)
        1) existing DataFile in this instance
        2) h5py file from data_filename
        3) h5py from obtained from sim_file

        Return value is set to self.data_file.

        :return: h5py DataFile instance
        """
        if self.data_file_:
            return self.data_file_

        if self.data_filename:
            self.data_file_ = h5py.File(self.data_filename, 'r+')
            return self.data_file_

        if self.sim_file_:
            self.data_file_ = h5py.File(self.sim_file.file_name, 'r+')
            return self.data_file_
        else:
            return None

    def get_data_file_write(self):
        """Deprecated"""
        if self.data_filename:
            return h5py.File(self.data_filename, 'r+')

        if self.sim_file_:
            return h5py.File(self.sim_file_.file_name, 'r+')
        return None

    def get_sim_file(self):
        """Gets a SimFile instance, either from directly attachd instance or creating one for the DataFile."""
        if self.sim_file_:
            return self.sim_file_
        if self.data_file or self.data_filename:
            sim_file = self.get_sim_file_from_data_file()
            self.sim_file = sim_file
            return sim_file

        else:
            return None

    def get_sim_file_from_data_file(self):
        """Creates a SimFile instance using DataFile"""
        if not self.data_filename:
            raise ValueError, "Needs data filename, please provide in constructor by kwarg 'data_filename'."

        if not self.session:
            raise ValueError, "Needs session, please provide in constructor by kwarg 'session'."

        session = self.session

        # return an existing SimFile for the given path
        try:
            ret = session.query(SimFile).filter(SimFile.file_name == self.data_filename).one()
            if ret:
                print "File with name %s already imported" % self.data_filename
                return ret
        except NoResultFound:
            pass

        # load handler, simtype and parameters from datafile
        record_file = self.data_file
        handler = realpickle.loads(record_file.attrs["handler"])
        sim_type = record_file.attrs["wrapper"]
        gen_parameters = self.get_gen_params()

        # get simtype object from DB
        try:
            type_db = session.query(SimType).filter(SimType.name == sim_type).one()
        except NoResultFound:
            raise ValueError, "Simulation type %s not known" % sim_type

        # get the network from the handler
        net = similar_or_new_net(session, handler.network)

        # commented to prevent database changes
        # session.add(net)

        # create a simfile, set paramaters from the datafile
        simfile = SimFile()
        if not self.to_memory:
            simfile.file_name = self.data_filename
        simfile.network = net
        simfile.sim_type = type_db
        simfile.parameters = {}
        simfile.parameters["gen_params"] = gen_parameters

        return simfile

    def get_data_fields(self):
        """hdf5 fields from DataFile"""
        return self.data_file.keys()

    def get_data_set(self):
        """Returns dict of hdf5 fields

        :return: dict of fields
        """
        keys = self.get_data_fields()
        out = {}
        for key in keys:
            print key
            out[key] = self.data_file[key]
        return out

    def prep_file(self, handler, suffix=''):
        """Prepare new DataFiles.

        :param handler: meanfield.error_functions.MeanfieldHandler instance
        :param suffix: suffix to add to the data_filename
        :return:
        """

        # Error if we already have a simfile
        if self.sim_file:
            raise RuntimeError, "The wrapped SimType has been run already (file exists)"

        # if to memory, we use a in-memory hdf5
        if self.to_memory:
            fn = 'memory_file_' + str(int(time.time()))
            print "[%s] creating memory only file" % self.__class__.__name__
            record_file = h5py.File(fn, 'a', driver="core", backing_store=False)
        else:
            # create a unique filename
            fn = "data/" + self.sim_type.name + "/" + self.sim_type.name + "_"
            if suffix:
                fn += suffix + "_"
            fn += str(int(time.time())) + ".hdf5"
            print "[%s] creating file: %s" % (self.__class__.__name__, fn)
            record_file = h5py.File(fn, 'a')

        # standard attachments
        record_file.attrs["handler"] = np.void(pickle.dumps(handler))
        record_file.attrs["wrapper"] = self.sim_type.name

        self.data_file = record_file
        print self.data_file

    def set_gen_params(self, gen_params):
        """sets gen_params - should probably be a property instead"""
        self.data_file.attrs["gen_params"] = np.void(pickle.dumps(gen_params))

    def get_gen_params(self):
        """gets gen_params - should probably be a property instead"""
        try:
            return pickle.loads(self.data_file.attrs["gen_params"])
        except AttributeError, e:
            try:
                return realpickle.loads(self.data_file.attrs["gen_params"])
            except AttributeError, e:
                raise e

    def close(self):
        """Close the wrapped datafile"""
        if self.data_file_:
            self.data_file_.close()
            return True
        return False

    def run(self, handler):
        """Run this simulation wrapper"""
        raise NotImplementedError


class SimBump(SimWrapper):
    """Runs a single bump trajectory storing in the same datafile:
    - recorded rates ("bump_rate")
    - the resulting centers ("dirs")
    - averaged shape ("shape_mean" and "shape_std")
    """

    def __init__(self, **kwargs):
        SimWrapper.__init__(self, **kwargs)

    def run(self, handler, cores=2, tmax=1000., do_only_setup=False, **kwargs):
        """Implements running and storage of a single bump.

        :param handler: meanfield.error_functions.MeanfieldHandler instance
        :param cores: cores to use for simulation
        :param tmax: maximal time for simulation
        :param do_only_setup:  only setup simulation but do not run
        :param kwargs: passed on to simulation.nest_simulation.NestSimulator
        :return: DataFile
        """
        self.prep_file(handler)

        from simulation.nest_simulation import NestSimulator

        sim = NestSimulator(handler, cores=cores)
        sim.set_paramset("bump")
        sim.set_params("gen", tmax=tmax, **kwargs)

        sim_ret = sim.run(do_only_setup)

        self.data_file.attrs["gen_params"] = np.void(pickle.dumps(sim.params["gen"]))

        act_set = self.data_file.create_dataset("bump_rate", sim_ret["bump_rate"].shape)
        act_set[:, :] = sim_ret["bump_rate"]

        act_set = self.data_file.create_dataset("shape_mean", sim_ret["shape_mean"].shape)
        act_set[:] = sim_ret["shape_mean"]
        act_set = self.data_file.create_dataset("shape_std", sim_ret["shape_std"].shape)
        act_set[:] = sim_ret["shape_std"]

        act_set = self.data_file.create_dataset("dirs", sim_ret["dirs"].shape)
        act_set[:] = sim_ret["dirs"]

        self.sim = sim

        return self.data_file

    def plot_mean(self, ax=None):
        """Plots the averaged rectified (zero-centered) shape of the excitatory population activity.

        :param ax: matplotlib axis for output.
        :return:
        """
        if ax is None:
            pl.figure()
            ax = pl.subplot(111)

        shape_mean = self.data_file["shape_mean"].value
        shape_std = self.data_file["shape_std"].value

        n_act = shape_mean.shape[0]
        xv_data = np.arange(n_act) / float(n_act) * 2. * np.pi - np.pi

        ax.plot(xv_data, shape_mean, lw=1., alpha=1.)
        ax.fill_between(
            xv_data, shape_mean - shape_std, shape_mean + shape_std,
            alpha=0.2, edgecolor="#222222", facecolor="#dddddd", antialiased=True)
        ax.set_title('Mean rates (centered)')
        ax.set_xlabel("Relative bump position")
        ax.set_ylabel("Mean firing rate [Hz]")
        ax.set_xticks([-np.pi, 0., np.pi], [r'$-\pi$', r'$0$', r'$\pi$'])
        ax.set_xlim([-np.pi, np.pi])
        ax.set_ylim([0, ax.get_ylim()[1]])

    def plot_rates(self, ax=None):
        """Plots excitatory population firing rates.

        :param ax: matplotlib axis for output.
        :return:
        """
        if ax is None:
            pl.figure(figsize=(10, 5))
            ax = pl.subplot(111)
        dirs = self.data_file["dirs"].value
        bump_rate = self.data_file["bump_rate"].value
        mappable = ax.imshow(bump_rate, aspect='auto', origin='lower')
        cbar = pl.colorbar(mappable, ax=ax)
        cbar.set_label('Firing rate [Hz]', rotation=270)
        dirs[:500] = np.nan
        ax.plot(dirs, lw=1.5, c='#FF0000', label="center")
        ax.set_title('Excitatory firing rates (filtered spikes)')
        ax.legend()
        ax.set_xlim([0, len(dirs)])
        ax.set_ylim([0, bump_rate.shape[0]])
        ax.set_xlabel("Time [ms]")
        ax.set_ylabel("Excitatory neurons")


def get_suffix(dic):
    """Create filename suffix for a given run, including on noise parameters in filename

    :param dic: dictionary with simulation parameters
    :return:
    """
    suffix = []
    keys = [mpr.P_EE, mpr.EL_NOISE, mpr.W_NOISE, "id", "fn_suffix"]
    fs = [lambda x: x < 1, lambda x: x > 0., lambda x: x > 0., lambda x: int, lambda x: str]

    for i, k in enumerate(keys):
        if dic.has_key(k) and fs[i](dic[k]) == int:
            suffix.append(k + "%i" % dic[k])
        elif dic.has_key(k) and fs[i](dic[k]) == str:
            suffix.append("%s" % dic[k])
        elif dic.has_key(k) and fs[i](dic[k]):
            suffix.append(k + "%.2f" % dic[k])

    return '_'.join(suffix)


class SimDrift(SimWrapper):
    """Runs several bump trajectories after another, storing in the same datafile:
    - the resulting centers ("dirs")
    - averaged shape ("shape_mean" and "shape_std")
    - population rates ("pop_rate_e" and "pop_rate_i")
    """

    def __init__(self, **kwargs):
        SimWrapper.__init__(self, **kwargs)

    def run(self, handler, initials=2, reps=1, cores=2, tmax=1000., base_seed_run=None, do_only_setup=False, **kwargs):
        """
        Implements running and storage of a several bumps on different initial positions with several repetitions.

        The base
        :param handler: meanfield.error_functions.MeanfieldHandler instance
        :param initials: how many initial positions (equally spaced)
        :param reps: how many repetitions
        :param cores: cores to use for simulation
        :param tmax: maximal time for simulation
        :param base_seed_run: base seed used for the running of each repetition (not connectivity, but used to control the Poisson noise)
        :param do_only_setup:  only setup simulation but do not run
        :param kwargs: passed on to simulation.nest_simulation.NestSimulator
        :return: DataFile
        """
        self.prep_file(handler, suffix=get_suffix(kwargs))

        data_file = self.data_file

        from simulation.nest_simulation import NestSimulator

        init_positions = np.arange(0., 1., 1. / float(initials))

        data_file.attrs["init_positions"] = init_positions
        data_file.attrs["reps"] = reps

        dir_set = data_file.create_dataset("dirs", (initials, reps, int(tmax)))

        shape_set_m = data_file.create_dataset("shape_mean", (initials, reps, handler.p_e[npr.C_E]))
        shape_set_std = data_file.create_dataset("shape_std", (initials, reps, handler.p_e[npr.C_E]))

        pop_rate_e = data_file.create_dataset("pop_rate_e", (initials, reps,))
        pop_rate_i = data_file.create_dataset("pop_rate_i", (initials, reps,))

        for i, init_pos in enumerate(init_positions):

            for r in range(reps):

                print "\n" + "#" * 15 + \
                      "\nRunning repetition %i/%i of center position %g" % (r + 1, reps, init_pos) + \
                      "\n" + "#" * 15 + "\n"

                sim = NestSimulator(handler, cores=cores)
                sim.set_paramset("bump")
                gen_params = {
                    "sig_center": init_pos,
                }
                gen_params.update(kwargs)

                # set base seed for run to something reproducible
                # but changing over each repetition and position
                if base_seed_run is not None:
                    base_seed_run_rep = base_seed_run + r * cores + i * reps * cores
                    gen_params.update({"base_seed_run": base_seed_run_rep})

                sim.set_params("gen", tmax=tmax, cores=cores, **gen_params)

                sim_ret = sim.run(do_only_setup)

                if not do_only_setup:
                    dir_set[i, r, :] = sim_ret["dirs"]
                    shape_set_m[i, r, :] = sim_ret["shape_mean"]
                    shape_set_std[i, r, :] = sim_ret["shape_std"]
                    pop_rate_e[i, r] = sim_ret["pop_rate"]['e']
                    pop_rate_i[i, r] = sim_ret["pop_rate"]['i']

        self.sim = sim
        self.set_gen_params(sim.params["gen"])
        return self.data_file

    def get_w_mat(self, handler, cores=2, **kwargs):
        """Returns the weight matrix generated by NEST for the used parameters.

        This is done in NEST to ensure exactly the same connectivity is generated by the seed used.
        Since this also depends on the number of cores, changing the number of cores will also change
        the resulting connectivity - therefore, to get the same weights, use the same number of cores.

        :param handler: meanfield.error_functions.MeanfieldHandler
        :param cores: number of cores (see note above)
        :param kwargs: passed to simulation.nest_simulation.NestSimulator
        :return: dictionary with weight matrix and a (deprecated) dump of the data used for plots at some point
        """
        from simulation.nest_simulation import NestSimulator

        print "################"
        print "Using %i cores! Made sure the number is the same as the original run?" % cores
        print "################"
        sim = NestSimulator(handler, cores=cores)
        sim.set_paramset("bump")
        gen_params = {
            "sig_start": 100.,
            "sig_len": 200.,
            "sig_weight": .4,
            "sig_center": 1,
            "return_w_mat": True
        }
        gen_params.update(kwargs)
        sim.set_params("gen", tmax=1., **gen_params)

        w_mat = sim.run()

        offset = 2
        test_dump = []
        n_e = w_mat.shape[1]
        for j in range(n_e):
            for i in range(n_e):
                if w_mat[i, j] > 0.:
                    tmp = (i - j) / float(n_e) * 2. * np.pi
                    dist = tmp - round(tmp / 2. / np.pi) * 2. * np.pi
                    tmp = "%i %i %.2f 1 %.6f" % (j + offset, i + offset, w_mat[i, j], dist)
                    test_dump.append(tmp)

        return {"w_mat": w_mat, "dump": test_dump}

    def plot_drift(
            self, ax=None, selectors=None, flip_y=0.,
            tmax_s=None, tmin_s=0., every_nth=1, every_nth_init=1,
            rasterized=True, color_lr=None, lw=.8, **kwargs):
        """Plots a drift profile

        :param ax: matplotlib axis
        :param selectors: 2 dim array with [inits,run_numbers] that should be plotted. None to plot all.
        :param flip_y: flip the axis orientation
        :param tmax_s: maximum to plot in seconds
        :param tmin_s: minimum to plot in seconds
        :param every_nth: plot only every nth trajectory (useful if you have loads)
        :param every_nth_init: plot only every nth initial position (useful if you have loads)
        :param rasterized: rasterize plot, significantly reduces filesize for many trajectories saved in pdf
        :param color_lr: colors by left/right of middle
        :param lw: line width
        :param kwargs: passed to matplotlib
        :return:
        """
        if ax is None:
            pl.figure(figsize=(10, 5))
            ax = pl.subplot(111)

        init_positions = self.data_file.attrs["init_positions"]
        gen_params = self.get_gen_params()

        dirs = self.data_file["dirs"].value

        inits = dirs.shape[0]
        reps = dirs.shape[1]
        times = dirs.shape[2]
        colval = lambda j: pl.cm.rainbow(j / float(inits))
        if color_lr:
            colval = lambda j: color_lr[int(j >= inits / 2)]

        tmax = times / 1e3
        if tmax_s is not None:
            tmax = tmax_s
        tmin = tmin_s
        print(tmax)
        if selectors:
            assert len(selectors) == inits

        for i in range(inits)[::every_nth_init]:
            for r in range(reps):

                if selectors and r not in selectors[i]:
                    continue

                if r % every_nth == 0:

                    xdim = np.arange(times, dtype=np.float) / 1e3
                    dirs_ = dirs[i, r, :].copy()

                    xdim = xdim[::10]
                    dirs_ = dirs_[::10]

                    pos = np.where(np.abs(np.diff(dirs_)) >= 400.)[0]

                    if len(pos) > 0:
                        dirs_[pos] = np.nan
                        xdim[pos] = np.nan

                    idxmax = np.where(xdim > tmax)
                    if len(idxmax[0]) > 0:
                        idxmax = idxmax[0][0]
                    else:
                        idxmax = -1
                    idxmin = np.where(xdim > tmin)[0][0]

                    if flip_y:
                        xdim, dirs_ = dirs_, xdim

                    ax.plot(xdim[idxmin:idxmax], dirs_[idxmin:idxmax], c=colval(i), rasterized=rasterized, lw=lw)

        xlim = (tmin, tmax)
        ylim = (0, 799)
        xlabel = r"Time [s]"
        ylabel = r"Center position $\varphi$"
        xticks = [0, 400, 800]
        yticks = np.arange(xlim[0], np.floor(xlim[1]), 2.)

        if flip_y:
            xlim, ylim = ylim, xlim
            ylim = (ylim[1], ylim[0])
            xlabel, ylabel = ylabel, xlabel
            xticks, yticks = yticks, yticks

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_yticks(xticks)
        ax.set_xticks(yticks)

        if flip_y:
            ax.set_xticklabels([r'$-\pi$', '$0$', r'$\pi$'])
        else:
            ax.set_yticklabels([r'$-\pi$', '$0$', r'$\pi$'])

    def get_initial_final(self, t0=2000, t1=-1):
        """Calculate initial and final positions in degrees

        :param t0: initial time [ms]
        :param t1: final time [ms]
        :return: dictionary of initials, finals and directories in degrees
        """
        convraddeg = 360. / (2. * np.pi)  # 1./800.*360 # convert from rad/ms to degree/ms
        convidxdeg = 360. / 800.  # 1./800.*360 # convert from idx/ms to degree/ms

        dirs = self.data_file["dirs"].value

        initials = dirs[:, :, t0].flatten()
        finals = dirs[:, :, t1].flatten()

        return {
            "initial": initials * convidxdeg / convraddeg,
            "final": finals * convidxdeg / convraddeg,
            "dirs": dirs * convidxdeg / convraddeg,
        }


def get_wrapper_class(sim_type):
    """Register simulation types here - not beautiful, but it works."""
    if sim_type.name == "bump":
        return SimBump
    elif sim_type.name == "drift":
        return SimDrift
    else:
        raise NameError, "Unknown simulation type (%s)." % sim_type.name


def get_wrapper_from_file(sim_file, **kwargs):
    """Helper to get wrapper from a classes.base.SimFile - should probably be a method there"""
    return get_wrapper_class(sim_file.sim_type)(sim_file=sim_file, **kwargs)


def get_wrapper_from_type(sim_type, **kwargs):
    """Helper to get wrapper from a classes.base.SimType - should probably be a method there"""
    return get_wrapper_class(sim_type)(sim_type=sim_type, **kwargs)
