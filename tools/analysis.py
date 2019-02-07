import numpy as np
from NeuroTools import signals


def process_spikes(ev, tmax, tmin=0.):
    """
    Converts NEST events to dictionary of numpy arrays

    :param ev: NEST events
    :param tmax: max time until when to convert
    :param tmin: min time from which to convert
    :return: dict of parsed spikes
    """
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
