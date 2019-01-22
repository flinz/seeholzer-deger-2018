import sys

import numpy as np
from NeuroTools import signals


def process_spikes(ev, tmax, tmin=0.):
    senders = ev['senders']
    times = ev['times']
    idxs = np.where(times >= tmin)

    ev['senders'] = senders[idxs]
    ev['times'] = times[idxs]
    senders = senders[idxs]
    times = times[idxs]

    spikelist = [(senders[i], ev['times'][i]) for i in range(len(senders))]

    if len(senders) > 0:
        ids = np.arange(senders.min(), senders.max() + 1)
    else:
        ids = np.array([])
    spike_trains = signals.SpikeList(spikelist, ids, t_start=tmin, t_stop=tmax)

    return {"ids": ids, "spike_list": spikelist, "spike_trains": spike_trains, 'ev': ev}


def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """
    valid = {"yes": True, "y": True, "ye": True, "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' (or 'y' or 'n').\n")


def query_value(question, values=['u', 'l', 'r'], default="r"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """

    if default not in values:
        raise ValueError("invalid default answer: '%s'" % default)

    prompt = ' [' + ', '.join(values) + ']'

    while True:
        sys.stdout.write("%s: %s" % (question, prompt))
        choice = raw_input().lower()
        if default is not None and choice == '':
            return default
        elif choice in values:
            return choice
        else:
            sys.stdout.write("Please respond with one of the values.\n")
