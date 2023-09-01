# Windows Setup

1. Install VSCode <https://code.visualstudio.com/sha/download?build=stable&os=win32-x64-user>

2. Run PowerShell as admin and run `wsl --install`.

3. Use all the defaults settings. Download <https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Windows-x86_64.exe> and install.

4. Then, open PowerShell and run

    ```pwsh
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

    New-Item -ItemType Directory -Force -Path "$env:HOME/AppData/Roaming/Code/User"
    if (Test-Path "$env:HOME/AppData/Roaming/Code/User/settings.json") {
        Copy-Item "$env:HOME/AppData/Roaming/Code/User/settings.json" "$env:HOME/AppData/Roaming/Code/User/settings.json.old"

    }

    Invoke-WebRequest -Uri https://raw.githubusercontent.com/gofflab/quant_mol_neuro_2023/main/setup/settings.json -OutFile "$env:HOME/AppData/Roaming/Code/User/settings.json"
    ```

5. Close and reopen PowerShell.
