#!/bin/bash
set -e

# Loop over all .ipynb files in the current directory and its subdirectories
find ./modules -type f -name '*.ipynb' | while read file; do
    abs_path=$(realpath $file)
    dir=$(dirname $file)
    abs_dir=$(realpath $dir)

    # Change the working directory to the file's directory
    (cd $abs_dir && jupyter nbconvert --to notebook --execute $abs_path)
done
