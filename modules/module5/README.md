## How to setup edu Rockfish

1. Go to this folder in bash
2. Run `ssh me440_<**REPLACE THIS WITH YOUR JHED ID IN SMALL CAPS**>@edulogin.arch.jhu.edu 'srun --time=00:10:00 --mem-per-cpu=3GB --cpus-per-task=1 bash -s' < edu_setup.sh`
    > For example, `ssh me440_lgoff137@edulogin.arch.jhu.edu 'srun --time=00:10:00 --mem-per-cpu=3GB --cpus-per-task=1 bash -s' < edu_setup.sh`
4. Enter your password
5. Wait for the script to finish

---

---

## How to use VSCode on Rockfish
5. Run `ssh me440_<REPLACE THIS WITH YOUR JHED ID>@edulogin.arch.jhu.edu "srun -n1 --time=2:00:00 --mem-per-cpu=4GB --cpus-per-task=2 VSCode-linux-x64/bin/code tunnel --accept-server-license-terms"`
    > For example, `ssh me440_lgoff137@edulogin.arch.jhu.edu "srun -n1 --time=2:00:00 --mem-per-cpu=4GB --cpus-per-task=2 VSCode-linux-x64/bin/code tunnel --accept-server-license-terms"`
7. It will tell you to go to <https://github.com/login/device> and enter the code it gives you. You will need a GitHub account.
8. Follow the guide at <https://marketplace.visualstudio.com/items?itemName=ms-vscode.remote-server> in the section `Connect to a tunnel`

### Click on the image to play

![stretched](https://github.com/gofflab/quant_mol_neuro_2023/assets/34997334/e83b3259-34ab-4a93-81b5-add9a847886c)

---

### Read below as well (first time only)

8. Now that you are connected, you'll need to transfer all your local extensions up. Press <Cmd+Shift+P> and type the first few words of `Remote: Install Local Extensions in '[SERVER-NAME]'`, press Enter, select everything, and click ok.
