#!/bin/bash

# Check if Java is installed
if type -p java >/dev/null 2>&1; then
    echo "Java is already installed."
else
    echo "Java is not installed. Installing..."
    sudo apt update  # Use 'sudo yum update' for CentOS/RHEL systems
    sudo apt install default-jre  # Use 'sudo yum install java' for CentOS/RHEL systems
fi

wget https://sourceforge.net/projects/janpa/files/latest/download -q --show-progress -O ~/janpa_binaries.zip
unzip ~/janpa_binaries.zip
rm ~/janpa_binaries.zip
