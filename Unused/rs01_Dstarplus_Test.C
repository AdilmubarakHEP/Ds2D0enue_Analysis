// #include "RooRealVar.h"
// #include "RooDataSet.h"
// #include "RooDataHist.h"
// #include "RooGaussian.h"
// #include "RooChebychev.h"
// #include "RooPolynomial.h"
// #include "TCanvas.h"
// #include "RooPlot.h"
// #include "TTree.h"
// #include "TH1D.h"
using namespace RooFit;

void rs01_Dstarplus_Test()
{
   // S e t u p   t h e   F r a m e   f o r   t h e   P l o t
   //------------------------------------------------------------
   TCanvas *c = new TCanvas("fit", "fit", 650, 700);
   TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
   TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.25);
   pad1->SetTopMargin(0.05);
   pad1->SetBottomMargin(0.05);
   pad1->SetLeftMargin(0.15);
   pad1->SetRightMargin(0.05);
   pad2->SetLeftMargin(0.15);
   pad2->SetRightMargin(0.05);
   pad2->SetTopMargin(0.01);
   pad2->SetBottomMargin(0.4);
   // pad2->SetBorderMode(0);
   pad1->Draw();
   pad2->Draw();
   pad1->cd();
   gStyle->SetOptStat(0);
   //=========================================================================================================================================================================//
   // I n i t i a l i z a t i o n:
   // --------------------------------------------------------------
   // The lines below initializes all my values
   int bin = 50;
   double xL = 0.142;
   double xR = 0.149;
   //=========================================================================================================================================================================//
   // I m p o r t s   D a t a   f r o m   t r e e   b r a n c h:
   // --------------------------------------------------------------
   TFile *file = new TFile("/home/belle2/amubarak/C03-Grid/TopoAna.root");
   TTree *tree = (TTree*)file->Get("Dstarplus");
   //=========================================================================================================================================================================//
   // D a t a:
   // --------------------------------------------------------------
   RooRealVar DeltaM("Ds_diff_D0pi", "Mass Difference", xL, xR);
   // RooRealVar Veto("Ds_gammaveto_M_Correction","Veto",-1,100);
   // RooRealVar isSignal("Ds_ifNANgiveX_isSignal_5","BCS",0,10);
   RooRealVar BS("Ds_extraInfo_BkgBDT", "Background Suppression", 0, 1);
   RooDataSet data("data", "dataset with mass", tree, RooArgSet(DeltaM,BS));
   // RooDataSet data("data", "dataset with mass", tree, RooArgSet(DeltaM,isSignal,BS));

   // Apply multiple selection cuts: 
   // - BS >= 0.531
   // - isSignal==1
   // - DeltaM within [0.1, 0.55] GeV/c^2
   RooDataSet *data_Dstarplus = (RooDataSet*)data.reduce("Ds_extraInfo_BkgBDT >= 0.531 && Ds_diff_D0pi >= 0.1 && Ds_diff_D0pi <= 0.55");
   // RooDataSet *data_Dstarplus = (RooDataSet*)data.reduce("Ds_extraInfo_BkgBDT >= 0.531 && Ds_ifNANgiveX_isSignal_5==1 && Ds_diff_D0pi >= 0.1 && Ds_diff_D0pi <= 0.55");
   //=========================================================================================================================================================================//
   // B u i l d   M o d e l
   // -------------------------------------------------------------- 
   // // First Attempt
   // //---------------------
   // // Gaussian Resolution model 1
	// RooRealVar sigmean1_PBkg("sigmean1_PBkg", "mean", 1.4543e-01, 0.135, 0.16);
	// RooRealVar sigwidth1_PBkg("sigwidth1_PBkg", "width", 2.9107e-04, 0, 0.01);
   //  sigmean1_PBkg.setConstant();
   //  sigwidth1_PBkg.setConstant();
	// RooGaussModel gaussm1_PBkg("gaussm1_PBkg", "Signal PDF", DeltaM, sigmean1_PBkg, sigwidth1_PBkg);
	// // Gaussian Resolution model 2
	// RooRealVar sigmean2_PBkg("sigmean2_PBkg", "mean", 1.4545e-01, 0.135, 0.16);
	// RooRealVar sigwidth2_PBkg("sigwidth2_PBkg", "width", 7.0138e-04, 0, 0.01);
   //  sigmean2_PBkg.setConstant();
   //  sigwidth2_PBkg.setConstant();
	// RooGaussModel gaussm2_PBkg("gaussm2_PBkg", "Signal PDF", DeltaM, sigmean2_PBkg, sigwidth2_PBkg);
	// // Gaussian Resolution model 2
	// RooRealVar sigmean3_PBkg("sigmean3_PBkg", "mean", 1.4604e-01, 0.135, 0.16);
	// RooRealVar sigwidth3_PBkg("sigwidth3_PBkg", "width", 2.6258e-03, 0, 0.01);
   //  sigmean3_PBkg.setConstant();
   //  sigwidth3_PBkg.setConstant();
	// RooGaussModel gaussm3_PBkg("gaussm3_PBkg", "Signal PDF", DeltaM, sigmean3_PBkg, sigwidth3_PBkg);

	// // Add the components
	// RooRealVar g1frac_PBkg("g1frac_PBkg","fraction of gauss1", 3.9035e-01, 0., 1.);
	// RooRealVar g2frac_PBkg("g2frac_PBkg","fraction of gauss2", 2.9710e-01, 0., 1.);
   // RooRealVar g3frac_PBkg("g3frac_PBkg","fraction of gauss2", 1.3614e-01, 0., 1.);
   //  g1frac_PBkg.setConstant();
   //  g2frac_PBkg.setConstant();
   //  g3frac_PBkg.setConstant();
	// RooAddPdf PeakingBkgModel("PeakingBkgModel","g1+g2",RooArgList(gaussm1_PBkg,gaussm2_PBkg,gaussm3_PBkg), RooArgList(g1frac_PBkg,g2frac_PBkg,g3frac_PBkg));

   // Second Model
   //---------------------
   // RooRealVar sig_mean("sig_mean", "mean of crystal", 1.4540e-01, 0.144, 0.147);
   // RooRealVar sig_sigmaL("sig_sigmaL", "sig_sigmaL", 3.0553e-04, 1e-4, 0.2);
   // RooRealVar sig_sigmaR("sig_sigmaR", "sig_sigmaR", 3.3288e-04, 1e-4, 0.2);
   // RooRealVar sig_alphaL("sig_alphaL", "sig_alphaL", 1.0972e+00, 1e-4, 5.);
   // RooRealVar sig_alphaR("sig_alphaR", "sig_alphaR", 1.0905e+00, 1e-4, 5.);
   // RooRealVar sig_nL("sig_nL", "sig_nL", 6.6200e+00, 1e-4,50);
   // RooRealVar sig_nR("sig_nR", "sig_nR", 7.1221e+00, 1e-4,50);
   //  sig_mean.setConstant();
   //  sig_sigmaL.setConstant();
   //  sig_sigmaR.setConstant();
   //  sig_alphaL.setConstant();
   //  sig_alphaR.setConstant();
   //  sig_nL.setConstant();
   //  sig_nR.setConstant();
   // RooCrystalBall PeakingBkgModel("PeakingBkgModel", "sig4", DeltaM, sig_mean,sig_sigmaL,sig_sigmaR,sig_alphaL,sig_nL,sig_alphaR,sig_nR);

   // Peaking Component
   //---------------------
   RooRealVar sig_mean("sig_mean", "mean of crystal", 1.4542e-01, 1.45e-01, 1.46e-01);
   RooRealVar sig_sigmaL("sig_sigmaL", "sig_sigmaL", 3.3295e-04, 2.0e-4, 6.0e-04);
   RooRealVar sig_sigmaR("sig_sigmaR", "sig_sigmaR", 3.4044e-04, 2.0e-4, 6.0e-04);
   RooRealVar sig_alphaL("sig_alphaL", "sig_alphaL", 1.3745e+00, 1e-4, 3.0e+00);
   RooRealVar sig_alphaR("sig_alphaR", "sig_alphaR", 1.3085e+00, 1e-4, 3.0e+00);
   RooRealVar sig_nL("sig_nL", "sig_nL", 2.8184e+00, 1e-4, 5.0e+00);
   RooRealVar sig_nR("sig_nR", "sig_nR", 2.8435e+00, 1e-4, 5.0e+00);
    sig_mean.setConstant();
    sig_sigmaL.setConstant();
    sig_sigmaR.setConstant();
    sig_alphaL.setConstant();
    sig_alphaR.setConstant();
    sig_nL.setConstant();
    sig_nR.setConstant();
   RooCrystalBall DstarPModel("DstarPModel", "sig4", DeltaM, sig_mean,sig_sigmaL,sig_sigmaR,sig_alphaL,sig_nL,sig_alphaR,sig_nR);

   // // Combinatorial Part
   // //----------------------
   // // First Gamma Function
   // RooRealVar x00_DstarPBkg("x00_DstarPBkg","x00_DstarPBkg", 1.4200e-01, 0.13, 0.15);
   // RooRealVar a0_DstarPBkg("a0_DstarPBkg","a0_DstarPBkg", 1.2268e-04, 0.0, 5.);
   // RooRealVar b0_DstarPBkg("b0_DstarPBkg","b0_DstarPBkg", 6.0771e-01, 0.0, 5.);    
   // x00_DstarPBkg.setConstant();
   // a0_DstarPBkg.setConstant();
   // b0_DstarPBkg.setConstant();
   // RooFormulaVar x00prime_DstarPBkg("#x00prime_DstarPBkg", "x[0] < x[1] ? x[0] : x[1]", RooArgList(DeltaM, x00_DstarPBkg));     
   // RooFormulaVar gamma0_DstarPBkg("#gamma0_DstarPBkg", "a0_DstarPBkg + 1", RooArgList(a0_DstarPBkg));    
   // RooFormulaVar beta0_for_gamma_DstarPBkg("beta0_for_gamma_DstarPBkg", "1./b0_DstarPBkg", RooArgList(b0_DstarPBkg));    
   // RooGamma GammaModel0_DstarPBkg("GammaModel1_DstarPBkg", "Gamma pdf", DeltaM, gamma0_DstarPBkg, beta0_for_gamma_DstarPBkg, x00prime_DstarPBkg);

   // // Second Gamma Function
   // RooRealVar x01_DstarPBkg("x01_DstarPBkg","x01_DstarPBkg", 1.3944e-01, 0.13, 0.15);
   // RooRealVar a1_DstarPBkg("a1_DstarPBkg","a1_DstarPBkg", 1.2337e+00, 0.0, 5.);
   // RooRealVar b1_DstarPBkg("b1_DstarPBkg","b1_DstarPBkg", 2.4856e+02, 100.0, 500.);     
   // x01_DstarPBkg.setConstant();
   // a1_DstarPBkg.setConstant();
   // b1_DstarPBkg.setConstant();
   // RooFormulaVar x01prime_DstarPBkg("#x01prime_S", "x[0] < x[1] ? x[0] : x[1]", RooArgList(DeltaM, x01_DstarPBkg));     
   // RooFormulaVar gamma1_DstarPBkg("#gamma1_DstarPBkg", "a1_DstarPBkg + 1", RooArgList(a1_DstarPBkg));    
   // RooFormulaVar beta1_for_gamma_DstarPBkg("beta1_for_gamma_DstarPBkg", "1./b1_DstarPBkg", RooArgList(b1_DstarPBkg));    
   // RooGamma GammaModel1_DstarPBkg("GammaModel2_DstarPBkg", "Gamma pdf", DeltaM, gamma1_DstarPBkg, beta1_for_gamma_DstarPBkg, x01prime_DstarPBkg);

   // // Third Gamma Function
   // RooRealVar x02_DstarPBkg("x02_DstarPBkg","x02_DstarPBkg", 1.4396e-01, 0.13, 0.15);
   // RooRealVar a2_DstarPBkg("a2_DstarPBkg","a2_DstarPBkg", 2.2879e+00, 0.0, 6.);
   // RooRealVar b2_DstarPBkg("b2_DstarPBkg","b2_DstarPBkg", 1.5561e+01, 0.0, 20.);     
   // x02_DstarPBkg.setConstant();
   // a2_DstarPBkg.setConstant();
   // b2_DstarPBkg.setConstant();
   // RooFormulaVar x02prime_DstarPBkg("#x02prime_S", "x[0] < x[1] ? x[0] : x[1]", RooArgList(DeltaM, x02_DstarPBkg));     
   // RooFormulaVar gamma2_DstarPBkg("#gamma2_DstarPBkg", "a2_DstarPBkg + 1", RooArgList(a2_DstarPBkg));    
   // RooFormulaVar beta2_for_gamma_DstarPBkg("beta2_for_gamma_DstarPBkg", "1./b2_DstarPBkg", RooArgList(b2_DstarPBkg));    
   // RooGamma GammaModel2_DstarPBkg("GammaModel3_DstarPBkg", "Gamma pdf", DeltaM, gamma2_DstarPBkg, beta2_for_gamma_DstarPBkg, x02prime_DstarPBkg);

   // // Add the components
   // RooRealVar f0_DstarPBkg("f0_DstarPBkg", "f0", 1.1043e+00, 0, 2);
   // RooRealVar f1_DstarPBkg("f1_DstarPBkg", "f1", 1.1058e-01, 0, 2);
   // RooRealVar f2_DstarPBkg("f2_DstarPBkg", "f2", 9.5479e-01, 0, 2);
   // f0_DstarPBkg.setConstant();
   // f1_DstarPBkg.setConstant();
   // f2_DstarPBkg.setConstant();
   // RooAddPdf CombBkgModel("CombBkgModel","g1+g2",RooArgList(GammaModel0_DstarPBkg,GammaModel1_DstarPBkg,GammaModel2_DstarPBkg), RooArgList(f0_DstarPBkg,f1_DstarPBkg,f2_DstarPBkg));

   // Total Model
   // Add the components
   // RooRealVar f_DstarPBkg_Peak("f_DstarPBkg_Peak", "f0", 1);
   // // RooRealVar f_DstarPBkg_Comb("f_DstarPBkg_Comb", "f1", 1);
   // f_DstarPBkg_Peak.setConstant();
   // // f_DstarPBkg_Comb.setConstant();
   // RooAddPdf DstarPModel("DstarPModel","g1+g2",RooArgList(PeakingBkgModel,CombBkgModel), RooArgList(f_DstarPBkg_Peak,f_DstarPBkg_Comb));

   // C o n s t r u c t   a   S i g na l   a n d   B a c k g r o u n d   P D F:
   // --------------------------------------------------------------
   RooRealVar Nsig_PBkg("nsig_PBkg", "#Signal events", 2.0573e+05, 0.5e+05, 4.0e+05);
   // RooAddPdf model("model", "g+a", RooArgList(PeakingBkgModel), RooArgList(Nsig_PBkg));
   RooAddPdf model("model", "g+a", RooArgList(DstarPModel), RooArgList(Nsig_PBkg));
   //=========================================================================================================================================================================//
   // S e t   t h e   C o n s t a n t   f o r   t h e    V a r i a b l e s
   //--------------------------------------------------------------------------
   // The code below fixes the value for a certain value.
   // B0_PBkg.setVal();
   // P0_PBkg.setVal();
   // g1frac_PBkg.setVal();
   // g2frac_PBkg.setVal();
   // g3frac_PBkg.setVal();
   // sigmean1_PBkg.setVal();
   // sigmean2_PBkg.setVal();
   // sigmean3_PBkg.setVal();
   // sigwidth1_PBkg.setVal();
   // sigwidth2_PBkg.setVal();
   // sigwidth3_PBkg.setVal();

   // The code below removes the range and fixes the value.
   // B0_PBkg.setConstant(true);
   // P0_PBkg.setConstant(true);
   // g1frac_PBkg.setConstant(true);
   // g2frac_PBkg.setConstant(true);
   // g3frac_PBkg.setConstant(true);
   // sigmean1_PBkg.setConstant(true);
   // sigmean2_PBkg.setConstant(true);
   // sigmean3_PBkg.setConstant(true);
   // sigwidth1_PBkg.setConstant(true);
   // sigwidth2_PBkg.setConstant(true);
   // sigwidth3_PBkg.setConstant(true);

   // DeltaM.setVal(0.325);
   // DeltaM.setConstant(true);
   //=========================================================================================================================================================================//
   // C o n s t r u c t   a   S i g na l   a n d   B a c k g r o u n d   P D F:
   // --------------------------------------------------------------
   // RooRealVar N("N", "PBackground events", 0., 500000);
   // RooAddPdf Signal("Signal", "g+a", RooArgList(SignalModel), RooArgList(N));
   //=========================================================================================================================================================================//
   // F i t
   //-------------
   // Perform extended ML fit of data:
   std::unique_ptr<RooFitResult> r{model.fitTo(*data_Dstarplus, Save(), PrintLevel(-1))};
   // If the fit is too broad and does not seem to fit a peak use the code below. Use it
   // once to determine where the peak should be and then adjust the mean above.
   // model.fitTo(data, Minos(kTRUE), PrintLevel(-1), Range(2,2.1));
   //=========================================================================================================================================================================//
   // P l o t   M o d e l   a n d   P r e p a r e   P a r a m e t e r   V a l u e s
   // ---------------------------------------------------------------------------
   // Create a RooPlot to draw on. We don't manage the memory of the returned
   // pointer. Instead we let it leak such that the plot still exists at the end of
   // the macro and we can take a look at it.
   // Create binning object with range (0.1,0.55)
   RooBinning tbins(xL, xR);
 
   // Add 60 bins with uniform spacing in range (0.1,0.16)
   tbins.addUniform(50, xL, 0.16);
 
   // Add 15 bins with uniform spacing in range (0.15,0.55)
   tbins.addUniform(20, 0.16, xR);

   RooPlot *frame = DeltaM.frame();
   *data_Dstarplus->plotOn(frame, Binning(50));
   // data_Dstarplus.plotOn(frame, Binning(tbins));
   model.plotOn(frame, LineColor(kGreen+1), RooFit::Name("model"));
   //=========================================================================================================================================================================//
   // C a l c u l a t e   c h i ^ 2
   // ------------------------------
   // Show the chi^2 of the curve w.r.t. the histogram
   // If multiple curves or datasets live in the frame you can specify
   // the name of the relevant curve and/or dataset in chiSquare()
   Int_t npar = model.getParameters(data_Dstarplus)->selectByAttrib("Constant",kFALSE)->getSize(); //select floating parameters and count their number

   cout << "NDF = " << npar << endl;
   cout << "chi^2/NDF = " << frame->chiSquare(npar) << endl;

   Double_t chi2ndf = frame->chiSquare(npar);
   //=========================================================================================================================================================================//
   // P l o t    L a b e l s
   // ---------------------------------------------------------------------------
   // // Add text to frame
   // TLatex *latex = new TLatex();
   // TString lum ("");
   // lum.Form("#splitline{Simulated Events}{2M Events}");
   // latex->SetNDC();
   // latex->SetTextSize(0.035);
   // Add text to frame
   TLatex *latex = new TLatex();
   TString lum ("");
   lum.Form("#splitline{Generic c#bar{c} Events}{#scale[0.5]{#int}#it{L} dt=1.44 ab^{-1}}");
   // lum.Form("#splitline{Generic c#bar{c} Events}{#scale[0.5]{#int}#it{L} dt=200 fb^{-1}}");
   latex->SetNDC();
   latex->SetTextSize(0.03); 
   
   // Add "chi-squared fit" below "2M Events"
   TString Likelihood ("");
   Likelihood.Form("#splitline{}{Likelihood Fit}");

   auto PerBin = ((xR-xL)/(bin))*1000;
   // TString title ("");
   TString y_axis ("");
   // TString x_axis ("");

   // title.Form("Generic d#bar{d} Events");
   y_axis.Form("Entries/(%.2f MeV/c^{2})", PerBin);
   // x_axis.Form("#Delta m  [GeV/c^{2}]");

   frame->SetTitle("");
   frame->GetYaxis()->SetTitle(y_axis);
   frame->GetYaxis()->SetRangeUser(0,28000);
   frame->GetXaxis()->SetTitle("");
   // frame->GetXaxis()->CenterTitle(true);
   // pad1->SetLogx();
   // pad1->SetLogy();
   frame->Draw();

   // The line below prints the luminosity onto the screen.
   // latex->DrawLatex(0.6,0.2, lum);
   latex->DrawLatex(0.18,0.89, lum);
   latex->DrawLatex(0.18,0.85, Likelihood);
   //=========================================================================================================================================================================//
   // DeltaM.setRange("window",0.15,1) ;
   // RooAbsReal* fracSigRange = Signal.createIntegral(DeltaM,DeltaM,"window") ; 
   // Double_t NWindow = N.getVal() * fracSigRange->getVal() ;
   // Double_t NWindow_error = N.getError() * fracSigRange->getVal() ;
   //=========================================================================================================================================================================//
   // P a r a m e t e r s
   // ---------------------------------------------------------------------------
   // The code below grabs the parameters and then places them in the legend
   // Float_t N_value = Nsig_PBkg.getValV();
   // Float_t b_value = b.getValV();
   // Float_t P_value = P.getValV();
   //=========================================================================================================================================================================//
   // P a r a m e t e r   B o x
   // ---------------------------------------------------------------------------
   auto xIF1 = pad1 -> GetUxmin() + pad1 -> GetLeftMargin();
   auto xIF2 = pad1 -> GetUxmax() - pad1 -> GetRightMargin();
   auto yIF1 = pad1 -> GetUymin() + pad1 -> GetBottomMargin();
   auto yIF2 = pad1 -> GetUymax() - pad1 -> GetTopMargin();
   auto dx = xIF2 - xIF1;
   auto dy = yIF2 - yIF1;
   auto x1 = xIF1 + 0.35*dx;
   auto x2 = xIF2;
   auto y1 = yIF1 + 0.75*dy;
   auto y2 = 0.93;
   auto legend = new TLegend(x1, y1, x2, y2);

   TString chi2("");
   TString Nsig_para("");
   TString Nbkg_para("");

   chi2.Form("#chi^{2}/ndf = %.2f", chi2ndf);
   Nsig_para.Form("N_{D^{*+}} = %.2f #pm %.2f", Nsig_PBkg.getValV(), Nsig_PBkg.getError());

   legend->AddEntry(frame->FindObject("data"),"data","pe");
   legend->AddEntry((TObject *)0, chi2, "");
   legend->AddEntry("model", "D^{*+} #rightarrow [D^{0} #rightarrow K^{-} #pi^{+}] #pi^{+}", "pl" );
   legend->AddEntry((TObject *)0, Nsig_para, "");

   legend -> SetBorderSize(0);
   legend -> SetTextFont(132);
   legend -> SetTextSize(0.03);
   legend -> SetNColumns(2);
   legend -> SetMargin(0.1);
   legend->Draw();
   //=========================================================================================================================================================================//
   // S h o w   P a r a m e t e r/ F i t    R e s u l t s
   // -------------------------------------------------------
   // When I use rooworkspace this first part will make sure that the parameters of the 
   // model will stay fixed.
   // DeltaM.setConstant(true);
   // First Try:
   // B0_PBkg.setConstant(true);
   // B0_PBkg.removeError();
   // P0_PBkg.setConstant(true);
   // P0_PBkg.removeError();
   // g1frac_PBkg.setConstant(true);
   // g1frac_PBkg.removeError();
   // g2frac_PBkg.setConstant(true);
   // g2frac_PBkg.removeError();
   // g3frac_PBkg.setConstant(true);
   // g3frac_PBkg.removeError();
   // sigmean1_PBkg.setConstant(true);
   // sigmean1_PBkg.removeError();
   // sigmean2_PBkg.setConstant(true);
   // sigmean2_PBkg.removeError();
   // sigmean3_PBkg.setConstant(true);
   // sigmean3_PBkg.removeError();
   // sigwidth1_PBkg.setConstant(true);
   // sigwidth1_PBkg.removeError();
   // sigwidth2_PBkg.setConstant(true);
   // sigwidth2_PBkg.removeError();
   // sigwidth3_PBkg.setConstant(true);
   // sigwidth3_PBkg.removeError();
   // Second Try:
   // sig_mean.setConstant(true);
   // sig_mean.removeError();
   // sig_sigmaL.setConstant(true);
   // sig_sigmaL.removeError();
   // sig_sigmaR.setConstant(true);
   // sig_sigmaR.removeError();
   // sig_alphaL.setConstant(true);
   // sig_alphaL.removeError();
   // sig_alphaR.setConstant(true);
   // sig_alphaR.removeError();
   // sig_nL.setConstant(true);
   // sig_nL.removeError();
   // sig_nR.setConstant(true);
   // sig_nR.removeError();

   // Verbose printing: Basic info, values of constant parameters, initial and
   // final values of floating parameters, global correlations
   r->Print("v");

   // Access basic information
   cout << "EDM = " << r->edm() << endl;
   cout << "-log(L) at minimum = " << r->minNll() << endl;

   // Extract covariance and correlation matrix as TMatrixDSym
   const TMatrixDSym &cor = r->correlationMatrix();
   const TMatrixDSym &cov = r->covarianceMatrix();
 
   // Print correlation, covariance matrix
   cout << "correlation matrix" << endl;
   cor.Print();
   cout << "covariance matrix" << endl;
   cov.Print();
   //=========================================================================================================================================================================//
   // S h o w   r e s i d u a l   d i s t s
   // -------------------------------------------------------
   pad2->cd();
   // Construct a histogram with the residuals of the data w.r.t. the curve
   RooHist *hresid = frame->residHist();
   for(int iPoint = 0; iPoint < hresid->GetN(); ++iPoint) {
      hresid->SetPointEYlow(iPoint, 0.0);
      hresid->SetPointEYhigh(iPoint, 0.0);
   }
   hresid->SetFillColor(38);
   // hresid->SetLineWidth(-802);
   // hresid->SetFillStyle(3002);
   // hresid->SetFillColor(2);

   // Create a new frame to draw the residual distribution and add the distribution to the frame
   RooPlot *frame1 = DeltaM.frame(Title("Residual Distribution"));
   frame1->addPlotable(hresid, "B");

   // frame1->GetYaxis()->SetRangeUser(-200,200);

   frame1->Draw();
   
   // Plot Titles
   //---------------
   frame1->SetTitle("");

   // Y-axis
   frame1->SetYTitle("Residual");
   frame1->GetYaxis()->CenterTitle(true);
   frame1->GetYaxis()->SetTitleSize(0.10);
   frame1->GetYaxis()->SetLabelSize(0.10);
   frame1->SetTitleOffset(0.5, "Y");

   // X-axis
   TString x_axis ("");
   x_axis.Form("#Delta m  [GeV/c^{2}]");
   frame1->SetXTitle(x_axis);
   frame1->GetXaxis()->CenterTitle(true);
   frame1->GetXaxis()->SetTitleSize(0.13);
   frame1->GetXaxis()->SetLabelSize(0.10);

   c->cd();
   //=========================================================================================================================================================================//
   // C r e a t e   w o r k s p a c e ,   i m p o r t   d a t a   a n d   m o d e l
   // -----------------------------------------------------------------------------
   // Create a new empty workspace
   RooWorkspace *w = new RooWorkspace("w_PBkg", "workspace");
   // Import model and all its components into the workspace
   w->import(DstarPModel);

   // E n c o d e   d e f i n i t i o n   o f   p a r a m e t e r s   i n   w o r k s p a c e
   // ---------------------------------------------------------------------------------------
   // Define named sets "parameters" and "observables", which list which variables should be considered
   // parameters and observables by the users convention
   //
   // Variables appearing in sets _must_ live in the workspace already, or the autoImport flag
   // of defineSet must be set to import them on the fly. Named sets contain only references
   // to the original variables, therefore the value of observables in named sets already
   // reflect their 'current' value
   RooArgSet *params = DstarPModel.getVariables();
   w->defineSet("parameters", *params);
   w->defineSet("observables", DeltaM);

   // E n c o d e   r e f e r e n c e   v a l u e   f o r   p a r a m e t e r s   i n   w o r k s p a c e
   // ---------------------------------------------------------------------------------------------------
 
   // Define a parameter 'snapshot' in the pdf
   // Unlike a named set, a parameter snapshot stores an independent set of values for
   // a given set of variables in the workspace. The values can be stored and reloaded
   // into the workspace variable objects using the loadSnapshot() and saveSnapshot()
   // methods. A snapshot saves the value of each variable, any errors that are stored
   // with it as well as the 'Constant' flag that is used in fits to determine if a
   // parameter is kept fixed or not.
 
   // The true flag imports the values of the objects in (*params) into the workspace
   // If not set, the present values of the workspace parameters objects are stored
   w->saveSnapshot("reference_fit_PeakingBkg", *params, true);

   // Print workspace contents
   w->Print();

   // S a v e   w o r k s p a c e   i n   f i l e
   // -------------------------------------------
   // Save the workspace into a ROOT file
   // w->writeToFile("/home/belle2/amubarak/Ds2D0enue_Analysis/05-ROOT/Roofit_Suggestion_Likehood/Fit_Parameter/rs01_DstarP.root");
   //=========================================================================================================================================================================//
   c->SaveAs("/home/belle2/amubarak/Ds2D0enue_Analysis/05-ROOT/Roofit_Suggestion_Likehood/Fit_Image/rs01_DeltaM_DstarP_Test.png");
   //=========================================================================================================================================================================//
   // Final Fit Status
   //---------------------------
   int fitStatus = r->status();  // Get the fit status
   int covQual = r->covQual();   // Get covariance quality

   std::cout << "Fit Status: " << fitStatus << std::endl;
   std::cout << "Covariance Matrix Quality: " << covQual << std::endl << endl;
}