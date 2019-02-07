import time
from sqlalchemy import create_engine
from sqlalchemy.orm import joinedload, sessionmaker
import classes.base as cb
from classes.base import (Base, Network, NeuronModel, Paramset, Population,
                          Synapse)
from classes.static import NeuronParameters as npr

metadata = Base.metadata

dbs = {
    'mysql': {
        'schema': "network_db",
        'url': "mysql://root:root@localhost/"
    },
    'sqlite': {
        'schema': "db/database_published.sqlite3",
        'url': "sqlite:///"
    },
}


def db_connection(db_name='sqlite'):
    """
    Creates SqlAlchemy db connection
    :param db_name: name of database, supported: "sqlite" and "mysql"
    :return: SqlAlchemy engine, SqlAlchemy sessionmaker, SqlAlchemy session
    """
    try:
        engine = create_engine(dbs[db_name]['url'] + dbs[db_name]['schema'], echo=False)
    except ImportError:
        print "Failed to connect to MYSQL database. Loading from SQLite3."
        db_name = 'sqlite'
        engine = create_engine(dbs[db_name]['url'] + dbs[db_name]['schema'], echo=False)
    print "Using database '%s': " % db_name, engine
    Session = sessionmaker(bind=engine)
    session = Session()
    return (engine, Session, session)


def get_network(id, session=None):
    """
    Gets network object with specified id
    :param id: id of network
    :param session: SqlAlchemy session
    :return: SqlAlchemy session, network
    """
    if session == None:
        _, _, session = db_connection()
    network = session.query(Network).filter(Network.id == id).options(joinedload('*')).one()
    return session, network


def print_latex_tables_synpars(exp_ids=[4,5,6]):
    """
    Print latex table of parameters
    :param exp_ids: Experiment ids
    """
    variables = ["U", "tau_f", "tau_r", npr.G_NMDA, npr.G_GABA, "w_0", "w_1", "w_sigma"]

    out = """\\begin{tabular}{lcc}\n"""
    out += r"""\\toprule\n"""

    _, _, session = db_connection()

    for i, exp_id in enumerate(exp_ids):

        experiment = session.query(cb.ExperimentSet).filter(
            cb.ExperimentSet.id == exp_id).one()

        for n, net in enumerate(experiment.networks):

            outkeys = {x: "\_".join(x.split("_")) for x in variables}
            outkeys["tau_f"] = "\tau_u"
            outkeys["tau_r"] = "\tau_x"

            print net.p_e["w_0"] + net.p_e["w_1"]
            for key in variables:

                try:
                    e_val = net.p_e.parameters[key]
                    if key not in ['w_0', 'w_1', 'w_sigma']:
                        i_val = net.p_i.parameters[key]
                except KeyError:
                    e_val = getattr(net.e.nmda_synapse, key)
                    if key not in ['w_0', 'w_1', 'w_sigma']:
                        i_val = getattr(net.i.nmda_synapse, key)

                if int(e_val) == e_val:
                    e_val = int(e_val)
                if key not in ['w_0', 'w_1', 'w_sigma']:
                    if int(i_val) == i_val:
                        i_val = int(i_val)

                out += " & $" + "%.4g" % e_val + "$ "
                if key in [npr.G_GABA, npr.G_NMDA]:
                    out += " & $" + "%.4g" % i_val + "$ "

            out += "\\\\ \n"
    out += """\\bottomrule\n\end{tabular}"""
    print out


def print_networks(exp_ids=[4,5,6]):
    """
    Print networks
    :param exp_ids: Experiment ids
    :return: list of networks that were printed
    """
    variables = ["U", "tau_f", "tau_r"]
    output = ["U", "tau_f", "tau_x"]

    networks = []
    _, _, session = db_connection()

    for i, exp_id in enumerate(exp_ids):

        experiment = session.query(cb.ExperimentSet).filter(
            cb.ExperimentSet.id == exp_id).one()
        print "Set: %s" % (experiment.name)
        for n, net in enumerate(experiment.networks):
            networks.append(net)
            print "id: %i -" % (net.id),
            for i, v in enumerate(variables):
                val = getattr(net.e.nmda_synapse, v)
                print "%s: %g," % (output[i], val),
            print ""
        print ""
    return networks


def similar_or_new_net(session, net):
    """
    For a given network, compare if we have a similar one in the database, and
    return this if found, otherwise return new network.

    :param session: SqlAlchemy session
    :param net: network
    :return: network
    """
    # try by id
    ret = session.query(Network).filter(Network.id == net.id).one()
    if ret:
        return ret

    # try by full net
    ret = session.query(Network).all()
    for net_db in ret:
        if net_db == net:
            return net_db

    # otherwise, new network
    net_new = Network()
    net_new.name = net.name + "_" + str(int(time.time()))

    # try the populations
    pops = session.query(Population).all()
    net_pop = {"e": None, "i": None}
    for pop in pops:
        if pop == net.e:
            net_pop["e"] = pop
        if pop == net.i:
            net_pop["i"] = pop

    for pop_key in net_pop.keys():

        if net_pop[pop_key] == None:

            pop_old = net.__getattribute__(pop_key)
            # new population
            pop_new = Population(None)
            if pop_old.name:
                pop_new.name = pop_old.name + "_" + str(int(time.time()))

            syn_old = pop_old.nmda_synapse
            print "---------"
            print syn_old
            print "---------"
            # synapse model from existing ones
            syns = session.query(Synapse).all()
            for syn_ in syns:
                print syn_, "?"
                if syn_ == syn_old:
                    pop_new.nmda_synapse = syn_
                    print "found existing synapse"
                    break

            # if still not set make new one
            if not pop_new.nmda_synapse:
                syn_new = Synapse(U=syn_old.U, tau_f=syn_old.tau_f, tau_r=syn_old.tau_r, tau_post=syn_old.tau_post,
                                  is_nmda=syn_old.is_nmda, fit_type=syn_old.fit_type)
                assert syn_new == syn_old, "Something went wrong while making new synapse"
                pop_new.nmda_synapse = syn_new

            # neuron model
            nmod = session.query(NeuronModel).filter(NeuronModel.id == pop_old.neuron_model.id).one()
            if nmod and nmod == pop_old.neuron_model:
                pop_new.neuron_model = nmod
            else:
                raise ValueError, "Neuronmodel not known, not supported yet."

            # paramset
            psets = session.query(Paramset).all()
            for pset in psets:
                if pop_old.paramset == pset:
                    pop_new.paramset = pset
                    break
            if not pop_new.paramset:
                pop_new.paramset = Paramset()
                if pop_old.paramset.name:
                    pop_new.paramset.name = pop_old.paramset.name + "_" + str(int(time.time()))
                pop_new.paramset.parameters = pop_old.paramset.parameters

            net_pop[pop_key] = pop_new

        net_new.__setattr__(pop_key, net_pop[pop_key])

    return net_new
