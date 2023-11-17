#!/bin/bash

# TURBOMOLE paths
export TURBODIR=~/COSMOlogic/TmoleX16/TURBOMOLE
source $TURBODIR/Config_turbo_env
export TMP_DIR=.tmp

# parallel
export PARA_ARCH=SMP
export PATH=$TURBODIR/bin/‘sysname‘:$PATH
export PARNODES=24

# setup timing script
cp timing.sh ~/bin
chmod +x ~/bin/timing.sh
echo "export PATH=:~/bin:$PATH" >> ~/.bashrc