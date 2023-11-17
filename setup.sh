env_name=linux_qm

# miniforge
if [[! -e Miniforge3-Linux-x86_64.sh]]
  then
  wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
fi
bash Miniforge3-Linux-x86_64.sh -b
echo "source ~/miniforge3/etc/profile.d/conda.sh" >> ~/.bashrc
echo "source ~/miniforge3/etc/profile.d/mamba.sh" >> ~/.bashrc
rm Miniforge3-Linux-x86_64.sh
conda config --set auto_activate_base false

# create env
mamba env create --name $env_name --file environment.yml
conda activate $env_name

# molli
wget https://github.com/SEDenmarkLab/molli_firstgen/archive/refs/heads/main.zip -O molli.zip
~/miniforge3/envs/$env_name/bin/pip install molli.zip #molli package
rm molli.zip