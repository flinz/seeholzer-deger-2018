#!/bin/bash

CURDIR=$(pwd);

source nest/nest_build/bin/nest_vars.sh;
export DYLD_LIBRARY_PATH=$CURDIR/nest/nest_build/lib/nest;
export LD_LIBRARY_PATH=$CURDIR/nest/nest_build/lib/nest;