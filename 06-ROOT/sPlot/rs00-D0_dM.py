import ROOT
from ROOT import RooFit
import uproot
from sweights import SWeight
from sweights import convert_rf_pdf

#=========================================================================================================================================================================
# Import Data into Dataframe:
#----------------------------------
# In this notebook we only process the main signal and the generic events,
# for illustration purposes.
# You can add other backgrounds after if you wish.
samples = ["Signal", "ccbar"]
# samples = ["Signal","All"]

DataFrames = {}  # define empty dictionary to hold dataframes
Date = "1112"
Attempt = "0"

# Signal:
DataFrames[samples[0]] =  uproot.concatenate("/home/belle2/amubarak/C01-Simulated_Events/Ds2D0enu-Signal.root:Dstree",library='pd')
# Background
for s in samples[1:]: # loop over samples
    DataFrames[s] =  uproot.concatenate("/home/belle2/amubarak/C03-Grid/Completed/Ds2D0e-Generic_Ds_" + Date +"24_"+ Attempt +"_"+ s +".root:Dstree",library='pd')
#=========================================================================================================================================================================
# S e t u p   t h e   F r a m e   f o r   t h e   P l o t
#------------------------------------------------------------
c = ROOT.TCanvas("fit", "fit", 650, 700)
pad1 = ROOT.TPad("pad1","pad1",0,0.25,1,1)
pad2 = ROOT.TPad("pad2","pad2",0,0,1,0.25)
pad1.SetTopMargin(0.05)
pad1.SetBottomMargin(0.05)
pad1.SetLeftMargin(0.15)
pad1.SetRightMargin(0.05)
pad2.SetLeftMargin(0.15)
pad2.SetRightMargin(0.05)
pad2.SetTopMargin(0.01)
pad2.SetBottomMargin(0.4)
# pad2.SetBorderMode(0);
pad1.Draw()
pad2.Draw()
pad1.cd()
ROOT.gStyle.SetOptStat(0)
#=========================================================================================================================================================================
# I n i t i a l i z a t i o n:
# --------------------------------------------------------------
# The lines below initializes all my values
bin = 50
xL = -0.05
xR = 0.05
#=========================================================================================================================================================================
# I m p o r t s   D a t a   f r o m   t r e e   b r a n c h:
# --------------------------------------------------------------
f = ROOT.TFile('/home/belle2/amubarak/C03-Grid/Completed/Ds2D0e-Generic_Ds_111224_0_ccbar.root')
tree = f.Get('Dstree')
#=========================================================================================================================================================================
# D a t a:
# --------------------------------------------------------------
dM = ROOT.RooRealVar("D0_dM", "Mass Difference", xL, xR)
Veto = ROOT.RooRealVar("Ds_gammaveto_M_Correction","Veto",0.1,1000)
BCS = ROOT.RooRealVar("Ds_chiProb_rank","BCS",-1,1.25)

data = ROOT.RooDataSet("data", "dataset with mass", ROOT.RooArgSet(dM,Veto,BCS), RooFit.Import(tree))
print(data.Print())
#=========================================================================================================================================================================
# B u i l d   M o d e l
# -------------------------------------------------------------- 
# Signal Peak:
#---------------
# # First Attempt:
# sigmean1 = ROOT.RooRealVar("sigmean1", "mean", -0.005, 0.005)
# sigwidth1 = ROOT.RooRealVar("sigwidth1", "width", 0, 0.01)
# gaussm1 = ROOT.RooGaussModel("gaussm1", "Signal PDF", dM, sigmean1, sigwidth1)

# sigmean2 = ROOT.RooRealVar("sigmean2", "mean", -0.005, 0.005)
# sigwidth2 = ROOT.RooRealVar("sigwidth2", "width", 0, 0.01)
# gaussm2 = ROOT.RooGaussModel("gaussm2", "Signal PDF", dM, sigmean2, sigwidth2)

# # Add the components
# g1frac = ROOT.RooRealVar("g1frac","fraction of gauss1", 0., 1.)
# g2frac = ROOT.RooRealVar("g2frac","fraction of gauss2", 0., 1.)
# Signal = ROOT.RooAddPdf("Signal","g1+g2",ROOT.RooArgSet(gaussm1,gaussm2),ROOT.RooArgSet(g1frac,g2frac))

# Second Attempt
sig_mean = ROOT.RooRealVar("sig_mean", "mean of crystal", 2.0582e-04, -0.005, 0.005)
sig_sigmaL = ROOT.RooRealVar("sig_sigmaL", "sig_sigmaL", 4.4691e-03, 1e-4, 0.02)
sig_sigmaR =  ROOT.RooRealVar("sig_sigmaR", "sig_sigmaR", 4.0743e-03, 1e-4, 0.02)
sig_alphaL = ROOT.RooRealVar("sig_alphaL", "sig_alphaL", 1.6108e+00, 1e-4, 4.)
sig_alphaR = ROOT.RooRealVar("sig_alphaR", "sig_alphaR", 1.3859e+00, 1e-4, 4.)
sig_nL = ROOT.RooRealVar("sig_nL", "sig_nL", 1.7096e+00, 1e-4, 5)
sig_nR = ROOT.RooRealVar("sig_nR", "sig_nR", 5.1600e+00, 1e-4, 5)
sig_mean.setConstant()
sig_sigmaL.setConstant()
sig_sigmaR.setConstant()
sig_alphaL.setConstant()
sig_alphaR.setConstant()
sig_nL.setConstant()
sig_nR.setConstant()
Signal = ROOT.RooCrystalBall("Signal", "Signal", dM, sig_mean,sig_sigmaL,sig_sigmaR,sig_alphaL,sig_nL,sig_alphaR,sig_nR)

# Background PDF:
#---------------------
# Constant Background:
c0 = ROOT.RooRealVar("c0", "coefficient #0", 1)
Background = ROOT.RooPolynomial("Background", "background p.d.f.", dM, c0)

# C o n s t r u c t   a   S i g na l   a n d   B a c k g r o u n d   P D F:
# --------------------------------------------------------------
nsig = ROOT.RooRealVar("nsig", "#Signal events", 2.9470e+05, 0., 600000)
nbkg = ROOT.RooRealVar("nbkg", "#background events", 7.3477e+04, 0., 600000)
nsig.setConstant()
nbkg.setConstant()
model = ROOT.RooAddPdf("model", "g+a", ROOT.RooArgSet(Signal, Background), ROOT.RooArgSet(nsig, nbkg))
#=========================================================================================================================================================================
# F i t
#-------------
# Perform fit and save result
# r = model.fitTo(data, RooFit.Save(True), RooFit.PrintLevel(-1))
# If the fit is too broad and does not seem to fit a peak use the code below. Use it
# once to determine where the peak should be and then adjust the mean above.
# model.fitTo(data, Minos(kTRUE), PrintLevel(-1), Range(2,2.1));
#=========================================================================================================================================================================
# P l o t   M o d e l   a n d   P r e p a r e   P a r a m e t e r   V a l u e s
# ---------------------------------------------------------------------------
# Create a RooPlot to draw on. We don't manage the memory of the returned
# pointer. Instead we let it leak such that the plot still exists at the end of
# the macro and we can take a look at it.
framedM = dM.frame()
data.plotOn( framedM, RooFit.Binning(bin))
model.plotOn(framedM, RooFit.Components("Background"), RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(46), RooFit.Name("Background"))
model.plotOn(framedM, RooFit.Components("Signal"), RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(28), RooFit.Name("Signal"))
model.plotOn(framedM, RooFit.LineColor(28), RooFit.Name("model"))
#=========================================================================================================================================================================
# C a l c u l a t e   c h i ^ 2
#------------------------------
# Show the chi^2 of the curve w.r.t. the histogram
# If multiple curves or datasets live in the frame you can specify
# the name of the relevant curve and/or dataset in chiSquare()
npar = model.getParameters(data).selectByAttrib("Constant",ROOT.kFALSE).getSize() #select floating parameters and count their number

print("NDF = ",npar)
print("chi^2/NDF = ",framedM.chiSquare(npar))

chi2ndf = framedM.chiSquare(npar)
#=========================================================================================================================================================================
# P l o t    L a b e l s
# ---------------------------------------------------------------------------
# Add text to frame
latex = ROOT.TLatex()
lum = "#splitline{Generic Events}{#scale[0.5]{#int}#it{L} dt=200 fb^{-1}}"
latex.SetNDC()
latex.SetTextSize(0.03)

PerBin = ((xR-xL)/(bin))*1000
y_axis = r'$Entries/(\; {width:.2f}\;MeV/c^2)$'.format(width = PerBin)

framedM.SetTitle("")
framedM.GetYaxis().SetTitle(y_axis)
# framedM.GetYaxis().SetRangeUser(0,22000)
framedM.GetXaxis().SetTitle("")
# framedM.GetXaxis().CenterTitle(ROOT.ktrue)
framedM.Draw()

# The line below prints the luminosity onto the screen.
latex.DrawLatex(0.18,0.89, lum)
#=========================================================================================================================================================================
# P a r a m e t e r   B o x
# ---------------------------------------------------------------------------
xIF1 = pad1.GetUxmin() + pad1.GetLeftMargin()
xIF2 = pad1.GetUxmax() - pad1.GetRightMargin()
yIF1 = pad1.GetUymin() + pad1.GetBottomMargin()
yIF2 = pad1.GetUymax() - pad1.GetTopMargin()
dx = xIF2 - xIF1
dy = yIF2 - yIF1
x1 = xIF1 + 0.4*dx
x2 = xIF2
y1 = yIF1 + 0.75*dy
y2 = 0.93
legend = ROOT.TLegend(x1, y1, x2, y2)

chi2 = r'$\chi^2 /ndf = {chi:.2f}$'.format(chi = chi2ndf) 
Nsig_para = r'{sig:.2f} $\pm$ {error:.2f}'.format(sig = nsig.getValV(), error = nsig.getError())
Nbkg_para = r'{bkg:.2f} $\pm$ {errorb:.2f}'.format(bkg = nbkg.getValV(), errorb = nbkg.getError())

legend.AddEntry(framedM.FindObject("data"),"data","pe")
legend.AddEntry("", chi2, "")
legend.AddEntry("Signal", "D^{0} #rightarrow K^{-} #pi^{+}", "pl" )
legend.AddEntry("", Nsig_para, "")
legend.AddEntry("Background", "Background", "pl" )
legend.AddEntry("", Nbkg_para, "")
legend.AddEntry("model", "Total Fit", "pl" )
legend.AddEntry("", " ", "")

legend.SetBorderSize(0)
legend.SetTextFont(132)
legend.SetTextSize(0.03)
legend.SetNColumns(2)
legend.SetMargin(0.1)
legend.Draw()
#=========================================================================================================================================================================
# S h o w   P a r a m e t e r/ F i t    R e s u l t s
# -------------------------------------------------------

# Verbose printing: Basic info, values of constant parameters, initial and
# final values of floating parameters, global correlations
# r.Print("v")

# # Access basic information
# cout << "EDM = " << r->edm() << endl;
# cout << "-log(L) at minimum = " << r->minNll() << endl;

# # Extract covariance and correlation matrix as TMatrixDSym
# const TMatrixDSym &cor = r->correlationMatrix();
# const TMatrixDSym &cov = r->covarianceMatrix();

# # Print correlation, covariance matrix
# cout << "correlation matrix" << endl;
# cor.Print();
# cout << "covariance matrix" << endl;
# cov.Print();
#=========================================================================================================================================================================
# S h o w   r e s i d u a l   d i s t s
# -------------------------------------------------------
pad2.cd()
# Construct a histogram with the residuals of the data w.r.t. the curve
hresid = framedM.residHist()
for iPoint in range(hresid.GetN()):
    hresid.SetPointEYlow(iPoint, 0.0)
    hresid.SetPointEYhigh(iPoint, 0.0)

hresid.SetFillColor(38)

# Create a new frame to draw the residual distribution and add the distribution to the frame
frame1 = dM.frame(RooFit.Title("Residual Distribution"))
frame1.addPlotable(hresid, "B")
frame1.Draw()

# Plot Titles
#---------------
frame1.SetTitle("")

# Y-axis
frame1.SetYTitle("Residual")
# frame1.GetYaxis().CenterTitle(true)
frame1.GetYaxis().SetTitleSize(0.10)
frame1.GetYaxis().SetLabelSize(0.10)

# X-axis
x_axis = r"$\Delta m  [GeV/c^{2}]$"
frame1.SetXTitle(x_axis)
# frame1.GetXaxis().CenterTitle(true)
frame1.GetXaxis().SetTitleSize(0.13)
frame1.GetXaxis().SetLabelSize(0.10)

frame1.SetTitleOffset(0.5, "Y")

# Draw Horizontal Line On Residual Plot
#----------------------------------------
# double xmin1 = xIF1; 
# double xmax1 = 1.0;
# TLine *line11 = new TLine(xIF1,0.0,xIF2,0.0);
# line11->SetLineColor(kRed); line11->SetLineWidth(3);
# line11->Draw("SAME");

c.cd()
#=========================================================================================================================================================================
# S a v e   F i l e s
#-------------------------------------------------------------
c.SaveAs("/home/belle2/amubarak/Ds2D0enue_Analysis/10-Images/D0_dM_SignalRegion.png")
#=========================================================================================================================================================================
# SWeights:
#---------------
# pick the fitting variable and convert it into a numpy array
df_sw = DataFrames[samples[1]]['D0_dM'].to_numpy()

# convert the PDFs from RooFit to python callables
SigPDF = convert_rf_pdf(Signal, dM)
BkgPDF = convert_rf_pdf(Background, dM)

# get the yields    
SigYield = Signal.expectedEvents(dM)
BkgYield = Background.expectedEvents(dM)

# calculate the sweights
mrange = (xL,xR)
sweighter = SWeight(df_sw, [SigPDF, BkgPDF], [nsig.getValV(), nbkg.getValV()],
            (mrange,), method="summation", verbose=True, checks=True)

# create columns for the signal and background weights
DataFrames[samples[1]]["D0_SigWeight"] = sweighter.get_weight(0, df_sw)
DataFrames[samples[1]]["D0_BkgWeight"] = sweighter.get_weight(1, df_sw)

# store the dataframe with the sweight columns in a new file
with uproot.recreate('/home/belle2/amubarak/output.root') as f:
    f['Dstree'] = DataFrames[samples[1]]
#=========================================================================================================================================================================