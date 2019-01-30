[![Build Status](https://travis-ci.com/flinz/seeholzer-deger-2018.svg?branch=master)](https://travis-ci.com/flinz/seeholzer-deger-2018)

# seeholzer-deger-2018

This repository contains code accompanying the publication:

A. Seeholzer, M. Deger, and W. Gerstner, “[Stability of working memory in continuous attractor networks under the control of short-term plasticity][1],” bioRxiv, p. 424515, Sep. 2018.


## Installation

### Minimum dependencies
* gcc-5 (tested) or higher
* libgsl
* cython 0.23.4 (tested) or higher
* python 2.7

### Build via docker image (the easy way)
```
git clone --recursive https://github.com/flinz/seeholzer-deger-2018.git
docker build -t seeholzer/2019 .
docker run -it --rm seeholzer/2019 /bin/bash
```

### Build locally (the harder way)
Makes Nest & installs python dependencies.
```
git clone --recursive https://github.com/flinz/seeholzer-deger-2018.git
make
``` 

## Networks

### Available networks

All networks used in [the publication][1] are included as presets in this repository. To list the available networks, call

```python
import tools.db
tools.db.print_networks()
```  

which will output
<sub>
```
Using database 'sqlite':  Engine(sqlite:///db/database_published.sqlite3)
Set: tau650 stablewidth
id: 45 - U: 0.6, tau_f: 650, tau_x: 150,
id: 46 - U: 0.4, tau_f: 650, tau_x: 150,
id: 47 - U: 0.2, tau_f: 650, tau_x: 150,
id: 48 - U: 0.1, tau_f: 650, tau_x: 150,
id: 49 - U: 0.08, tau_f: 650, tau_x: 150,
id: 50 - U: 0.06, tau_f: 650, tau_x: 150,
id: 51 - U: 0.04, tau_f: 650, tau_x: 150,
id: 52 - U: 0.02, tau_f: 650, tau_x: 150,
id: 53 - U: 0.8, tau_f: 650, tau_x: 150,
id: 66 - U: 1, tau_f: 650, tau_x: 150,

Set: tau650 depressions
id: 54 - U: 0.1, tau_f: 650, tau_x: 140,
id: 55 - U: 0.1, tau_f: 650, tau_x: 120,
id: 56 - U: 0.1, tau_f: 650, tau_x: 180,
id: 57 - U: 0.1, tau_f: 650, tau_x: 200,
id: 58 - U: 0.4, tau_f: 650, tau_x: 140,
id: 59 - U: 0.4, tau_f: 650, tau_x: 120,
id: 60 - U: 0.4, tau_f: 650, tau_x: 180,
id: 61 - U: 0.4, tau_f: 650, tau_x: 200,
id: 62 - U: 0.8, tau_f: 650, tau_x: 200,
id: 63 - U: 0.8, tau_f: 650, tau_x: 180,
id: 64 - U: 0.8, tau_f: 650, tau_x: 140,
id: 65 - U: 0.8, tau_f: 650, tau_x: 120,
id: 67 - U: 0.1, tau_f: 650, tau_x: 160,
id: 68 - U: 0.4, tau_f: 650, tau_x: 160,
id: 69 - U: 0.8, tau_f: 650, tau_x: 160,

Set: tau1000 stablewidth
id: 66 - U: 1, tau_f: 650, tau_x: 150,
id: 71 - U: 0.8, tau_f: 1000, tau_x: 150,
id: 72 - U: 0.6, tau_f: 1000, tau_x: 150,
id: 73 - U: 0.4, tau_f: 1000, tau_x: 150,
id: 74 - U: 0.2, tau_f: 1000, tau_x: 150,
id: 75 - U: 0.1, tau_f: 1000, tau_x: 150,
id: 76 - U: 0.08, tau_f: 1000, tau_x: 150,
id: 77 - U: 0.06, tau_f: 1000, tau_x: 150,
id: 78 - U: 0.04, tau_f: 1000, tau_x: 150,
```
</sub>

To obtain a single network (e.g. `id=48`) from this list

```python
net_id = 48
session, net = tools.db.get_network(net_id)
```

which returns

```
Network (Scan3_U=0.1_sig=0.5)
```

### Neuron parameters

To access (and change) single neuron parameters use the dictionaries

```python
net.e.paramset.parameters  # for excitatory neurons
net.i.paramset.parameters  # for inhibitory neurons

```

Note: the excitatory parameter additionally contains the connectivity parameters `w_0, w_1, w_sigma`  


## Running simulations

Simulations use
- a handler class `meanfield.error_functions.MeanfieldHandler`: access for simulator to network parameters
- a simulation type from `simulation.sim_types`: implements different simulation protocols and provides a simulation wrapper
- a simulation wrapper from `simulation.sim_types`: handles access to data files and reading/writing of simulation data  

See below for examples of how to run a single bump trajectory, or several bump trajectories at once. 

### Single bump trajectory

The simulation type `SimBump` from `simulation.sim_types` implements running of a single bump trajectory, 
and saves all simulation data (including spike times). This can be used to investigate 
the shape of bump firing-rate distributions, for example.

See [scripts/run_bump.py](tree/master/scripts/run_bump.py) for an example of how to use the interface:

```python
# get meanfield handler (abstracts network parameter handling) and set connectivity parameter
handler = erfs.MeanfieldHandler(net)

# get simulation type
btype = session.query(cb.SimType).filter(cb.SimType.name == "bump").one()

# get simulation wrapper
wrapper = btype.get_wrapper(to_memory=True, session=session)

# parameters for simulation
sim_params = {
    "show_results": True,  # plot results
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
```

#### Plots

The wrapper provides the following plots after running

##### Excitatory firing rates (filtered spikes)

```python
wrapper.plot_rates()
pl.savefig('out_rates.png')
```
![Excitatory firing rates (filtered spikes)](docs/img/run_bump_rates.png)


##### Mean firing rates (rectified)

```python
wrapper.plot_mean()
pl.savefig('out_mean.pdf')
```

![Mean firing rates (rectified)](docs/img/run_bump_mean.png)

# Repeated bump trajectories (drift profiles)

The simulation type `SimDrift` from `simulation.sim_types` implements 
running several bump trajectories from equally spaced initial positions with repetitions. 
To save space, this simulation type saves only bump center positions together with aggregate data, but no spike times.

[1]: https://www.biorxiv.org/content/10.1101/424515v2
