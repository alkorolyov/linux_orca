#!/bin/bash

# setup timing script
mkdir ~/bin
cp timing ~/bin
chmod +x ~/bin/timing
export PATH=:~/bin:$PATH
echo "export PATH=:$PATH:~/bin:" >> ~/.bashrc