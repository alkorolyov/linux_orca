# TURBOMOLE paths
TURBODIR=~/COSMOlogic/TmoleX16/TURBOMOLE
source $TURBODIR/Config_turbo_env


# parallel
export PARA_ARCH=SMP
export PATH=$TURBODIR/bin/‘sysname‘:$PATH
export PARNODES=24
