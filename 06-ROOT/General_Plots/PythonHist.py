import ROOT

# Disable multi-threading if it causes issues
ROOT.ROOT.DisableImplicitMT()

# Define settings
Bins = 50
var = "e_electronID"  # Variable to plot
Range = (0.0, 1.0)
BS = 2  # Electron ID cut
factor = 0.7

def create_hist(df, var_name, bins, range_min, range_max, color, title):
    hist = ROOT.TH1F(title, title, bins, range_min, range_max)
    filtered_df = df.Filter(f"{var_name} <= {BS}")
    
    # Fill histogram manually
    for row in filtered_df.AsNumpy([var_name])[var_name]:
        hist.Fill(row)
    
    hist.SetLineColor(color)
    hist.SetLineWidth(2)
    return hist

# Load ROOT files into RDataFrame
signal_file = "/home/belle2/amubarak/C01-Simulated_Events/Ds2D0enu-Signal.root"
all_file = "/home/belle2/amubarak/C03-Grid/Completed/Ds2D0e-Generic_Ds_031325_0_All.root"
tree_name = "Dstree"

df_signal = ROOT.RDataFrame(tree_name, signal_file)
df_all = ROOT.RDataFrame(tree_name, all_file)

# Create histograms for different categories
hist_signal = create_hist(df_signal, var, Bins, Range[0], Range[1], ROOT.kBlack, "Signal")

hist1 = create_hist(df_all.Filter("abs(Ds_D0_Dstarplus) == 1 && abs(D0_isSignal) == 1"), var, Bins, Range[0], Range[1], ROOT.kGreen+2, "D*+ -> D0 pi+")
hist2 = create_hist(df_all.Filter("abs(Ds_D0_Dstar0) == 1 && abs(D0_isSignal) == 1"), var, Bins, Range[0], Range[1], ROOT.kRed, "D*0 -> D0 gamma/pi0")
hist3 = create_hist(df_all.Filter("abs(Ds_D0_NoDstarplusDstar0) == 1 && abs(D0_isSignal) == 1"), var, Bins, Range[0], Range[1], ROOT.kBlue, "D0")
hist4 = create_hist(df_all.Filter("abs(Ds_D0_Other) == 1 || (abs(D0_mcPDG) == 421 && abs(D0_isSignal) == 0)"), var, Bins, Range[0], Range[1], ROOT.kMagenta, "Other")

# Scale signal histogram
hist_signal.Scale(factor)

# Stack histograms
hs = ROOT.THStack("hs", "")
hs.Add(hist1)
hs.Add(hist2)
hs.Add(hist3)
hs.Add(hist4)

# Draw histograms
canvas = ROOT.TCanvas("canvas", "Histogram", 800, 600)
hs.Draw("hist")
hist_signal.Draw("hist same")

# Customize plot
hs.GetXaxis().SetTitle("e^{-} Electron ID")
hs.GetYaxis().SetTitle("Entries / Bin")
ROOT.gStyle.SetOptStat(0)

# Legend
legend = ROOT.TLegend(0.6, 0.7, 0.9, 0.9)
legend.AddEntry(hist_signal, "Signal", "l")
legend.AddEntry(hist1, "D*+ -> D0 pi+", "l")
legend.AddEntry(hist2, "D*0 -> D0 gamma/pi0", "l")
legend.AddEntry(hist3, "D0", "l")
legend.AddEntry(hist4, "Other", "l")
legend.Draw()

# Save and show plot
canvas.SaveAs("histogram.png")
canvas.Draw()