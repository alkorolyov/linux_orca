#!/bin/bash
#####################################
#### Prepare working environment ####
#####################################
env_name=linux_qm

# miniforge
if [ ! -f Miniforge3-Linux-x86_64.sh ]; then
  wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh  -q --show-progress
fi
bash Miniforge3-Linux-x86_64.sh -b
source ~/miniforge3/etc/profile.d/conda.sh
source ~/miniforge3/etc/profile.d/mamba.sh
rm Miniforge3-Linux-x86_64.sh
conda config --set auto_activate_base false
conda config --set ssl_verify false

# create env
mamba env create --name $env_name --file environment.yml
conda activate $env_name

pip config set global.trusted-host "pypi.org files.pythonhosted.org pypi.python.org"
pip install cclib

# molli
# wget https://github.com/SEDenmarkLab/molli_firstgen/archive/refs/heads/main.zip -O molli.zip -q --show-progress
# pip install molli.zip #molli package
# rm molli.zip

# install as a local copy for dev
pip install -e .

# add ~/bin to PATH
# and setup scripts
# mkdir ~/bin
# cp scripts/timing ~/bin
# chmod +x ~/bin/timing
# export PATH=:~/bin:$PATH
# echo "export PATH=:$PATH:~/bin:" >> ~/.bashrc
