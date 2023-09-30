# How to setup edu Rockfish

1. Go to this folder in bash
2. Run `ssh [USERNAME]@edulogin.arch.jhu.edu 'srun --time=00:10:00 --mem-per-cpu=3GB --cpus-per-task=1 bash -s' < edu_setup.sh`
3. Enter your password
4. Wait for the script to finish
5. Run `ssh [USERNAME]@edulogin.arch.jhu.edu srun --time=2:00:00 --mem-per-cpu=4GB --cpus-per-task=2 ./code tunnel`
6. It will tell you to go to <https://github.com/login/device> and enter the code it gives you. You will need a GitHub account.
7. Follow the guide at <https://marketplace.visualstudio.com/items?itemName=ms-vscode.remote-server> in the section `Connect to a tunnel`
8. Now that you are connected, you'll need to transfer all your local extensions up. Press <Cmd+Shift+P> and type the first few words of `Remote: Install Local Extensions in '[SERVER-NAME]'`, press Enter, select everything, and click ok.
