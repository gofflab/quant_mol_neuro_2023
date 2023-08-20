# quant_mol_neuro_2023

ME400.825 FA 23 Quantitative Molecular Neuroscience - JHU

## Setup environment


<details>
    <summary>Windows</summary>
    
    1. Install VSCode https://code.visualstudio.com/sha/download?build=stable&os=win32-x64-user
    
    2. Run PowerShell as admin and run `wsl --install` without the backticks.

    3. Use all the defaults settings. Download <https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Windows-x86_64.exe> and install.
    
    4. Then, open PowerShell and run

    %HOMEDRIVE%%HOMEPATH%\mambaforge\condabin\conda init pwsh
    
    $extensions = @(
    "ms-python.python",
    "charliermarsh.ruff",
    "ms-python.isort",
    "ms-toolsai.jupyter",
    "ms-python.black-formatter",
    "kevinrose.vsc-python-indent",
    "reageyao.biosyntax",
    "usernamehw.errorlens",
    "yellpika.latex-input",
    "christian-kohler.path-intellisense",
    "mechatroner.rainbow-csv",
    "ms-vscode-remote.remote-ssh",
    "ms-vscode-remote.remote-wsl",
    "tetradresearch.vscode-h2o",
    "foxundermoon.shell-format",
    "tomoki1207.pdf",
    "redhat.vscode-yaml"
    )
    
    foreach ($ext in $extensions) {
        code --install-extension $ext
    }
        

    3. Close and reopen PowerShell.

</details>

<details>
    <summary>Mac: copy the following and paste in your terminal</summary>

    xcode-select --install
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

    if [ $(arch) = "arm64" ]; then
        file="Mambaforge-MacOSX-arm64.sh"
    else
        file="Mambaforge-MacOSX-x86_64.sh"
    fi

    brew install --cask visual-studio-code

    for ext in ms-python.python charliermarsh.ruff ms-python.isort ms-toolsai.jupyter ms-python.black-formatter kevinrose.vsc-python-indent reageyao.biosyntax usernamehw.errorlens yellpika.latex-input christian-kohler.path-intellisense mechatroner.rainbow-csv ms-vscode-remote.remote-ssh ms-vscode-remote.remote-wsl tetradresearch.vscode-h2o foxundermoon.shell-format tomoki1207.pdf redhat.vscode-yaml; do
      code --install-extension $ext
    done
    
    curl -LO "https://github.com/conda-forge/miniforge/releases/latest/download/${file}"
    chmod +x $file
    ./$file -b
    rm $file
    $HOME/mambaforge/condabin/conda init zsh
    $HOME/mambaforge/condabin/mamba init zsh
    source $HOME/.zshrc

</details>
