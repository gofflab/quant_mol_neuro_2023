#!/bin/bash
set -e
curl -LSsO https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
chmod +x Miniforge3-Linux-x86_64.sh
./Miniforge3-Linux-x86_64.sh -b
miniforge3/condabin/conda init
rm Miniforge3-Linux-x86_64.sh

source ~/.bashrc
curl -LSsO https://update.code.visualstudio.com/1.75.1/linux-x64/stable
tar xvzf stable
rm stable
