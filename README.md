# linux_orca
Installation of orca in Linux

# ORCA install

# download https://orcaforum.kofo.mpg.de/app.php/dlext/?view=detail&df_id=180
# copy it to /tmp

# openmpi install in install.sh
cd /tmp
apt update
apt install build-essential mc -y
wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.1.tar.gz
tar xvf openmpi-4.1.1.tar.gz
cd openmpi-4.1.1
./configure --prefix=/opt/orca-5.0.3/openmpi-4.1.1
mkdir -p /opt/orca-5.0.3
make install

cd /opt
tar xvf /tmp/orca_5_0_3*.tar.xz
mv orca_5_0_3* orca-5.0.3/orca
cd /orca-5.0.3
wget https://raw.githubusercontent.com/lijikun/orca-linux-deployment/master/orcainit5

# execute it in you user shell and restart it
echo "source /opt/orca-5.0.3/orcainit5 >> ~/.bashrc"

export ORCA_PATH=/opt/orca-5.0.3/orca
export PATH="/opt/orca-5.0.3/orca:$PATH"; export LD_LIBRARY_PATH="/opt/orca-5.0.3/orca:$LD_LIBRARY_PATH"
export PATH="/opt/orca-5.0.3/openmpi-4.1.1/bin:$PATH"; export LD_LIBRARY_PATH="/opt/orca-5.0.3/openmpi-4.1.1/lib:$LD_LIBRARY_PATH"
