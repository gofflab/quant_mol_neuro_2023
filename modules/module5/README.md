# How to setup edu Rockfish

1. Go to this folder in bash
2. Run `ssh [USERNAME]@edulogin.arch.jhu.edu 'srun --time=00:10:00 --mem-per-cpu=3GB --cpus-per-task=1 bash -s' < edu_setup.sh`
3. Enter your password
4. Wait for the script to finish

---

**Click on the image to play**

5. Run `ssh me440_<REPLACE THIS WITH YOUR JHED ID>@edulogin.arch.jhu.edu srun --time=2:00:00 --mem-per-cpu=4GB --cpus-per-task=2 ./code tunnel --accept-server-license-terms`
6. It will tell you to go to <https://github.com/login/device> and enter the code it gives you. You will need a GitHub account.
7. Follow the guide at <https://marketplace.visualstudio.com/items?itemName=ms-vscode.remote-server> in the section `Connect to a tunnel`

> ![stretched](https://github.com/gofflab/quant_mol_neuro_2023/assets/34997334/95c7c9a4-f7ac-40e5-ab57-95a43809fec1)

---

### Read below as well

8. Now that you are connected, you'll need to transfer all your local extensions up. Press <Cmd+Shift+P> and type the first few words of `Remote: Install Local Extensions in '[SERVER-NAME]'`, press Enter, select everything, and click ok.
