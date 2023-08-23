# Mac

1. Open the terminal and copy the following and paste in your terminal

    ```bash
    xcode-select --install
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

    if [ $(arch) = "arm64" ]; then
        file="Mambaforge-MacOSX-arm64.sh"
    else
        file="Mambaforge-MacOSX-x86_64.sh"
    fi

    brew install --cask visual-studio-code

    for ext in ms-python.python \
           charliermarsh.ruff \
           ms-python.isort \
           ms-toolsai.jupyter \
           ms-python.black-formatter \
           kevinrose.vsc-python-indent \
           reageyao.biosyntax \
           usernamehw.errorlens \
           yellpika.latex-input \
           christian-kohler.path-intellisense \
           mechatroner.rainbow-csv \
           ms-vscode-remote.remote-ssh \
           ms-vscode-remote.remote-wsl \
           tetradresearch.vscode-h2o \
           foxundermoon.shell-format \
           tomoki1207.pdf \
           redhat.vscode-yaml; do

        code --install-extension $ext
    done

    mkdir -p "$HOME/Library/Application Support/Code/User/"
    if test -f "$HOME/Library/Application Support/Code/User/settings.json"; then
        cp "$HOME/Library/Application Support/Code/User/settings.json" "$HOME/Library/Application Support/Code/User/settings.json.old"
    fi
    curl https://raw.githubusercontent.com/gofflab/quant_mol_neuro_2023/main/setup/settings.json > "$HOME/Library/Application Support/Code/User/settings.json"

    curl -LO "https://github.com/conda-forge/miniforge/releases/latest/download/${file}"
    chmod +x $file
    ./$file -b
    rm $file
    $HOME/mambaforge/condabin/conda init zsh
    $HOME/mambaforge/condabin/mamba init zsh
    source $HOME/.zshrc

    ```
