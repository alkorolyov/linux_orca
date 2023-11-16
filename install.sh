if [ "$EUID" -ne 0 ]
  then echo "Please run as root"
  exit
fi

echo "download https://orcaforum.kofo.mpg.de/app.php/dlext/?view=detail&df_id=180"
echo "copy it to /tmp"
read -p "Press any key to continue"

# openmpi install
cd /tmp
apt update
apt install build-essential -y
wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.1.tar.gz
tar xvf openmpi-4.1.1.tar.gz
cd openmpi-4.1.1
./configure --prefix=/opt/orca-5.0.3/openmpi-4.1.1
mkdir -p /opt/orca-5.0.3
make install

cd /opt
tar xvf /tmp/orca_5_0_3*.tar.xz
mv orca_5_0_3* orca-5.0.3/orca

ORCA_PATH="/opt/orca-5.0.3"
export PATH="$ORCA_PATH/orca:$PATH"; export LD_LIBRARY_PATH="$ORCA_PATH/orca:$LD_LIBRARY_PATH"
export PATH="$ORCA_PATH/openmpi-4.1.1/bin:$PATH"; export LD_LIBRARY_PATH="$ORCA_PATH/openmpi-4.1.1/lib:$LD_LIBRARY_PATH"

