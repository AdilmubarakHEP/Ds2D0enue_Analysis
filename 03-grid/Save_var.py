import ROOT

# Open the ROOT file
file = ROOT.TFile("/home/belle2/amubarak/C01-Simulated_Events/Signal/output_test_kmpippi0_eff10_May2020.root")
tree = file.Get("DstreeCh2")

# Check if the tree was loaded properly
if not tree:
    raise RuntimeError("Could not find 'Dstree' in the file.")

# Extract and write all branch names to a text file
output_path = "/home/belle2/amubarak/Ds2D0enue_Analysis/03-grid/Save_var_Mode2.txt"
with open(output_path, "w") as f:
    for branch in tree.GetListOfBranches():
        name = branch.GetName()
        f.write(name + "\n")

print(f"Saved {tree.GetListOfBranches().GetEntries()} branch names to {output_path}")