# seeholzer-deger-2018

This repository contains code accompanying the publication by Alexander Seeholzer, Moritz Deger & Wulfram Gerstner: [Stability of working memory in continuous attractor networks under the control of short-term plasticity](https://www.biorxiv.org/content/early/2018/09/23/424515.1).

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
