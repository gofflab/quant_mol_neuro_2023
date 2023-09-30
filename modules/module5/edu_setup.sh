#!/bin/bash
set -e
echo "Hi $USER! Let's get you set up for the course."
curl -LSsO https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
chmod +x Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b
miniforge3/condabin/conda init
rm Miniforge3-Linux-x86_64.sh

source ~/.bashrc
mamba install notebook matplotlib plotnine pandas
pip install pyfastx

curl -LSsO https://update.code.visualstudio.com/1.75.1/linux-x64/stable
tar xvzf stable
rm stable
ln -s $PWD/VSCode-linux-x64/bin/code $PWD/code
