#!/bin/bash

if [ "$EUID" -ne 0 ]
  then echo "Please run as root"
  exit
fi

echo "download orca-5.0.3_openmpi-4.1.1 from https://orcaforum.kofo.mpg.de/app.php/dlext/?view=detail&df_id=180"
echo "copy it to /tmp"
read -p "Press any key to continue"

# openmpi compile and install
cd /usr/src
apt update
apt install build-essential gfortran -y
wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.1.tar.gz -q --show-progress
tar xvf openmpi-4.1.1.tar.gz ; rm openmpi-4.1.1.tar.gz
cd openmpi-4.1.1
mkdir -p /opt/openmpi-4.1.1
./configure --prefix=/opt/openmpi-4.1.1
make install

# check openmpi
# mpiexec --version

# extract orca
cd /tmp
tar xvf orca_5_0_3*.tar.xz
mv orca_5_0_3_linux_x86-64_shared_openmpi411 /opt/orca-5.0.3

# add to path
source orcainit

# test run
orca test.inp

echo "Add orcainit to .bashrc"