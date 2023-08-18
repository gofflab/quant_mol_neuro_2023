# quant_mol_neuro_2023

ME400.825 FA 23 Quantitative Molecular Neuroscience - JHU

## Setup environment

1. Install VSCode from `https://code.visualstudio.com/download`
2. Install conda
    <details>
        <summary>

            **Windows**

        </summary>

        1. Use all the defaults settings. Download <https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Windows-x86_64.exe> and install.
        2. Then, open PowerShell and run

        ```
        %HOMEDRIVE%%HOMEPATH%\mambaforge\condabin\conda init pwsh
        ```

    </details>

    <details>
        <summary>

            **Mac**

        </summary>

        Copy the following and paste in your terminal

        ```sh
        if [ $(arch) = "arm64" ]; then
            file="Mambaforge-MacOSX-arm64.sh"
        else
            file="Mambaforge-MacOSX-x86_64.sh"
        fi

        curl -LO "https://github.com/conda-forge/miniforge/releases/latest/download/${file}"
        chmod +x $file
        ./$file -b
        rm $file
        $HOME/mambaforge/condabin/conda init zsh
        $HOME/mambaforge/condabin/mamba init zsh
        source $HOME/.zshrc
        ```

    </details>
