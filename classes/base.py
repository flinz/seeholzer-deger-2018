import cPickle
import json

import dill
import numpy as np
import pylab as pl
from sqlalchemy import Table, Column, Float, ForeignKey, Integer, String, Boolean, PickleType
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, backref

from classes.static import NeuronParameters as npr

Base = declarative_base()
prefix = "EIBrunel"


class NameRepr(object):
    def __repr__(self):
        out = self.__class__.__name__
        if self.name:
            out += " (%s)" % self.name
        return out


class DictMixin(object):
    """Adds obj[key] access to a mapped class.

    This class proxies dictionary access to an attribute
    called ``_proxied``.  The class which inherits this class
    should have an attribute called ``_proxied`` which points to a dictionary.

    """

    def __len__(self):
        return len(self.parameters)

    def __iter__(self):
        return iter(self.parameters)

    def __getitem__(self, key):
        return self.parameters[key]

    def __contains__(self, key):
        return key in self.parameters

    def __setitem__(self, key, value):
        self.parameters[key] = value

    def __delitem__(self, key):
        del self.parameters[key]


class Paramset(Base, NameRepr, DictMixin):
    __tablename__ = prefix + '_paramset'
    __table_args__ = {'mysql_engine': 'InnoDB', 'sqlite_autoincrement': True}

    id = Column(Integer, primary_key=True)
    name = Column(String(512), nullable=True)

    parameters = Column(PickleType(pickler=json), nullable=False)

    def __init__(self, **kwargs):
        self.parameters = kwargs

    def __eq__(self, other):
        if not other or (other.__class__ != self.__class__):
            return False
        for attr in self.parameters:
            if attr in [npr.W_0, npr.W_1, npr.W_SIGMA]:
                continue
            if attr in [npr.G_GABA, npr.G_AMPA, npr.G_NMDA, npr.NU_EXT]:
                if np.abs(self.parameters[attr] - other.parameters[attr]) > 1e-4:
                    return False
            if not self.parameters[attr] == other.parameters[attr]:
                return False
        return True


class SynapseFitTypes(object):

    def fit_linear_zero(x, a):
        return a * x

    def fit_1(x, a, b, c, d, e):
        return e + b / (c + a * np.exp(-x * d))

    def fit_2(x, a, b, c, d):
        return a * x ** b / (c + x ** b) + d

    types = {
        0: fit_linear_zero,
        1: fit_1,
        2: fit_2,
    }


def stp_dgl_u(U, tau_f, tau_r, rate_ms):
    return float(U) * (1. + rate_ms * tau_f) * (1. + float(U) * rate_ms * tau_f) ** (-1.)


def stp_dgl_r(U, tau_f, tau_r, rate_ms):
    return (1. + float(U) * rate_ms * tau_f) * (
                1. + float(U) * rate_ms * (tau_f + tau_r + rate_ms * tau_f * tau_r)) ** (-1.)


def stp_dgl_u_nobase(U, tau_f, tau_r, rate_ms):
    return float(U) * (rate_ms * tau_f) / (1. + float(U) * rate_ms * tau_f)


def stp_dgl_r_nobase(U, tau_f, tau_r, rate_ms):
    return 1. / (1. + float(U) * rate_ms ** 2 * tau_f * tau_r / (1. + float(U) * tau_f * rate_ms))


class Synapse(Base):
    __tablename__ = prefix + '_synapse'
    __table_args__ = {'mysql_engine': 'InnoDB', 'sqlite_autoincrement': True}

    id = Column(Integer, primary_key=True)
    name = Column(String(512), nullable=True)

    U = Column(Float, default=1.)
    tau_r = Column(Float, default=0.)
    tau_f = Column(Float, default=0.)
    tau_post = Column(Float, default=100.)
    cv = Column(Float, default=1.)
    is_nmda = Column(Boolean, default=False)

    def __eq__(self, other):
        if not other or (other.__class__ != self.__class__):
            return False
        for attr in ['U', 'tau_r', 'tau_f', 'tau_post', 'cv', 'is_nmda']:
            if not self.__getattribute__(attr) == other.__getattribute__(attr):
                return False
        return True

    @property
    def is_stp(self):
        return self.U < 1.0 or self.tau_r > 0.

    fit_tau_nmda = Column(Float, default=100.)
    fit_p_opt = Column(PickleType(pickler=json), nullable=True)
    fit_data_x = Column(PickleType(pickler=cPickle), nullable=True)
    fit_data_y = Column(PickleType(pickler=cPickle), nullable=True)
    fit_type = Column(Integer, default=2)

    fun = None

    def __repr__(self):
        return "Synapse (U: %.3f, tau_f: %.1f, tau_r: %.1f, is_nmda: %s)" % (
        self.U, self.tau_f, self.tau_r, self.is_nmda)

    def __init__(self, U=1., tau_f=0., tau_r=0., tau_post=100., is_nmda=False, fit_type=2):
        self.U = U
        self.tau_f = tau_f
        self.tau_r = tau_r
        self.tau_post = tau_post
        self.is_nmda = is_nmda
        self.fit_type = fit_type

    def do_fit(self, fit_type):
        # fit the synapse
        # change fit type to comply
        # change fit_nmda, fit_data_x, fit_data_y
        raise NotImplementedError

    def __call__(self, x):
        if self.fun == None:
            self.fun = self.make_fun()
        return self.fun(x)

    def plot(self, **param):
        try:
            param['label']
        except KeyError:
            param['label'] = self.__str__()

        x = np.arange(0., 181., .1) * 1e-3
        if self.fit_data_x:
            x = self.fit_data_x
            if self.fit_data_y:
                pl.scatter(x * 1e3, self.fit_data_y, alpha=.2, **param)
        pl.plot(x * 1e3, [self(xv) for xv in x], **param)
        pl.xlim((0, 180))
        pl.legend(loc="best")

        pl.xlabel("Input rate [Hz]")
        pl.ylabel("Channel activation [1]")

    def stp_u(self, x):
        return stp_dgl_u(self.U, self.tau_f, self.tau_r, x)

    def stp_r(self, x):
        return stp_dgl_r(self.U, self.tau_f, self.tau_r, x)

    def stp_du_dx(self, x):
        u = self.stp_u(x)
        return -(
                (self.tau_f * (-1. + u) * u)
                / (1. + x * self.tau_f)
        )

    def stp_dr_dx(self, x):
        r = self.stp_r(x)
        return -(
                (self.U * r *
                 (
                         -self.tau_f + (self.tau_f + self.tau_r + 2. * x * self.tau_f * self.tau_r) * r
                 )
                 ) / (1. + self.U * x * self.tau_f)
        )

    def stp_dur_dx(self, x):
        return self.stp_u(x) * self.stp_dr_dx(x) + self.stp_r(x) * self.stp_du_dx(x)

    def stp_ur(self, x):
        return self.stp_r(x) * self.stp_u(x)

    def make_fun(self):
        # output the analytical trace
        if self.fit_p_opt == None:
            if self.is_nmda:
                raise Exception("NMDA+STP synapse that has not been fitted yet.")
            # calculate the synaptic activation as a linear synapse scaled by STP
            # tau in ms, x in spikes/ms
            return lambda x: self.tau_post * x * self.stp_ur(x)
        # output the fitted trace
        else:
            return lambda x: max(0, SynapseFitTypes.types[self.fit_type](x, *self.fit_p_opt))

        raise Exception("Should not have gotten here.")

    def __eq__(self, other):
        if not other or (other.__class__.__name__ != self.__class__.__name__):
            return False
        for attr in ['U', 'tau_r', 'tau_f', 'tau_post', 'fit_type', 'is_nmda']:
            if not self.__getattribute__(attr) == other.__getattribute__(attr):
                print attr, "is different"
                return False
        return True


class NeuronModel(Base, NameRepr):
    __table_args__ = {'mysql_engine': 'InnoDB', 'sqlite_autoincrement': True}
    __tablename__ = prefix + '_neuron_model'

    id = Column(Integer, primary_key=True)
    name = Column(String(512), nullable=False)

    # extend this later to others
    simulator = "nest"

    def __eq__(self, other):
        if not other or (other.__class__ != self.__class__):
            return False
        return self.name == other.name

    def __init__(self, name):
        self.name = name


class Population(Base, NameRepr):
    __table_args__ = {'mysql_engine': 'InnoDB', 'sqlite_autoincrement': True}
    __tablename__ = prefix + '_population'

    id = Column(Integer, primary_key=True)
    name = Column(String(512), nullable=True)

    nmda_synapse_id = Column(ForeignKey(prefix + '_synapse.id'))
    nmda_synapse = relationship("Synapse", backref=backref("populations", order_by=id))

    paramset_id = Column(ForeignKey(prefix + '_paramset.id'))
    paramset = relationship("Paramset", backref=backref("populations", order_by=id))

    neuron_model_id = Column(Integer, ForeignKey(prefix + '_neuron_model.id'))
    neuron_model = relationship("NeuronModel", foreign_keys=neuron_model_id)

    def __init__(self, neuron_model):
        self.neuron_model = neuron_model

    def __eq__(self, other):
        if not other or (other.__class__ != self.__class__):
            return False
        for attr in ['nmda_synapse', 'paramset', 'neuron_model']:
            if not self.__getattribute__(attr) == other.__getattribute__(attr):
                return False
        return True


class Network(Base, NameRepr):
    __table_args__ = {'mysql_engine': 'InnoDB', 'sqlite_autoincrement': True}
    __tablename__ = prefix + '_network'

    id = Column(Integer, primary_key=True)
    name = Column(String(512), nullable=True)

    e_id = Column(ForeignKey(prefix + '_population.id'))
    e = relationship(Population, foreign_keys=e_id)

    i_id = Column(ForeignKey(prefix + '_population.id'))
    i = relationship(Population, foreign_keys=i_id)

    @property
    def p_e(self):
        return self.e.paramset

    @property
    def p_i(self):
        return self.i.paramset

    def p_(self, pop_key, key):
        return self.__getattribute__(pop_key).paramset[key]

    def __eq__(self, other):
        if not other or (other.__class__ != self.__class__):
            return False
        for attr in ['e', 'i']:
            if not self.__getattribute__(attr) == other.__getattribute__(attr):
                return False
        return True


post_keywords = Table(prefix + '_experiment_sets_networks', Base.metadata,
                      Column('exp_id', Integer, ForeignKey(prefix + '_experiment_set.id'), primary_key=True),
                      Column('net_id', Integer, ForeignKey(prefix + '_network.id'), primary_key=True)
                      )

post_keywords = Table(prefix + '_experiments_simfiles', Base.metadata,
                      Column('exp_id', Integer, ForeignKey(prefix + '_experiment.id'), primary_key=True),
                      Column('simfile_id', Integer, ForeignKey(prefix + '_sim_file.id'), primary_key=True)
                      )


class ExperimentSet(Base, NameRepr):
    __table_args__ = {'mysql_engine': 'InnoDB', 'sqlite_autoincrement': True}
    __tablename__ = prefix + '_experiment_set'

    id = Column(Integer, primary_key=True)
    name = Column(String(255))
    description = Column(String(1024))

    networks = relationship('Network', secondary=prefix + '_experiment_sets_networks', backref='experiment_sets')


class Experiment(Base, NameRepr):
    __table_args__ = {'mysql_engine': 'InnoDB', 'sqlite_autoincrement': True}
    __tablename__ = prefix + '_experiment'

    id = Column(Integer, primary_key=True)
    name = Column(String(255))
    description = Column(String(1024))

    experiment_set_id = Column(Integer, ForeignKey(prefix + '_experiment_set.id'))
    experiment_set = relationship(ExperimentSet, foreign_keys=experiment_set_id,
                                  backref=backref("experiments", order_by=id))

    simfiles = relationship('SimFile', secondary=prefix + '_experiments_simfiles', backref='experiments')


class SimType(Base, NameRepr):
    __table_args__ = {'mysql_engine': 'InnoDB', 'sqlite_autoincrement': True}
    __tablename__ = prefix + '_sim_type'

    id = Column(Integer, primary_key=True)
    name = Column(String(255))  # , primary_key=True)
    description = Column(String(1024))

    def get_wrapper(self, **kwargs):
        import simulation.sim_types
        return simulation.sim_types.get_wrapper_from_type(self, **kwargs)


class SimFile(Base):
    __table_args__ = {'mysql_engine': 'InnoDB', 'sqlite_autoincrement': True}
    __tablename__ = prefix + '_sim_file'

    id = Column(Integer, primary_key=True)

    file_name = Column(String(255), unique=True, nullable=False)

    network_id = Column(ForeignKey(prefix + '_network.id'))
    network = relationship(Network, foreign_keys=network_id, backref=backref("files", order_by=id))

    sim_type_id = Column(ForeignKey(prefix + '_sim_type.id'))
    sim_type = relationship(SimType, foreign_keys=sim_type_id, backref=backref("files", order_by=id))

    parameters = Column(PickleType(pickler=dill), nullable=True)
    data = Column(PickleType(pickler=cPickle), nullable=True)
    choice = Column(Integer, default=1)

    def get_wrapper(self, **kwargs):
        import simulation.sim_types
        return simulation.sim_types.get_wrapper_from_file(self, **kwargs)

    @property
    def is_imported(self):
        return len(self.runs) > 0


class SimRun(Base):
    __tablename__ = prefix + '_sim_run'
    __table_args__ = {'mysql_engine': 'InnoDB'}

    id = Column(Integer, primary_key=True)

    parameters = Column(PickleType(pickler=dill), nullable=True)
    data = Column(PickleType(pickler=cPickle), nullable=True)

    sim_file_id = Column(ForeignKey(prefix + '_sim_file.id'))
    sim_file = relationship(SimFile, backref=backref("runs", order_by=id))

    def get_wrapper(self, **kwargs):
        return self.sim_file.get_wrapper(**kwargs)
