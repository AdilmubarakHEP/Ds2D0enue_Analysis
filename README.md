# Ds2D0enue_Analysis



## Rare Decay
## D_s+ -> [D0 -> K- pi+] e+ nu_e


## Folder Description

- [ ] 00-Setup: 
- [ ] 01-Event_Generation: The steering and .dec file used to generate signal events
- [ ] 02-Reconstruction_Scripts: The background folder contain steering scripts that is used to simulate the peaking background in the analysis. The last file in the folder is a steering script to reconstruct my decay mode.
- [ ] 03-Grid: 
- [ ] 04-ML: The basf2 folder is where I set up and train a BDT using the basf2 software. This folder contains the BDT for fake $D^{0}$ suppression and background suppression. The local folder is where I attempt to use a python machine learning package separate from the basf2 software.
- [ ] 05-ROOT: All the fitting in the analysis is done here
- [ ] 06-Python: The machine learning folder is where I check for correlations and check all the variables between signal and background to see which has the most discriminating power. The "Plots_and_Tables" folder is where I keep all my jupyter scripts for my analysis that create plots. The name of the files describe what is mostly plotted
- [ ] 07-Python_Functions: Self made python functions for tasks that I often repeat. This is so I don't keeo copy and pasting long pieces of code
- [ ] 08-Images:
- [ ] 09-Sources:
- [ ] 10-Save: Exercises
- [ ] Unused: Codes that aren't being used, but may be usefull in the future