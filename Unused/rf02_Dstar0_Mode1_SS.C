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

void rf02_Dstar0_Mode1_SS()
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
   double xL = 0.1;
   double xR = 0.55;
   //=========================================================================================================================================================================//
   // I m p o r t s   D a t a   f r o m   t r e e   b r a n c h:
   // --------------------------------------------------------------
   // TFile *file = new TFile("C01-Simulated_Events/Ds2D0enu-Background.root");
   // TTree *tree = (TTree*)file->Get("Dtree");

   // I m p o r t s   D a t a   f r o m   t r e e   b r a n c h:
   // --------------------------------------------------------------
   TFile *file = new TFile("/home/belle2/amubarak/C01-Simulated_Events/Dstar0/Dstar0-Background_Mode_1.root");
   TTree *tree = (TTree*)file->Get("Dstree");
   //=========================================================================================================================================================================//
   // D a t a:
   // --------------------------------------------------------------
   RooRealVar DeltaM("Ds_diff_D0pi", "Mass Difference", xL, xR);
   // RooRealVar Veto("Ds_gammaveto_M_Correction","Veto",-1,100);
   RooRealVar isSignal("Ds_ifNANgiveX_isSignal_5","BCS",0,10);
   RooRealVar BS("Ds_extraInfo_BkgBDT", "Background Suppression", 0, 1);
   // RooDataSet data("data", "dataset with mass", tree, RooArgSet(DeltaM,BS));
   RooDataSet data("data", "dataset with mass", tree, RooArgSet(DeltaM,isSignal,BS));

   // Apply multiple selection cuts: 
   // - BS >= 0.531
   // - isSignal==1
   // - DeltaM within [0.1, 0.55] GeV/c^2
   // RooDataSet *data_Mode1 = (RooDataSet*)data.reduce("Ds_extraInfo_BkgBDT >= 0.531 && Ds_diff_D0pi >= 0.1 && Ds_diff_D0pi <= 0.55");
   RooDataSet *data_Mode1 = (RooDataSet*)data.reduce("Ds_extraInfo_BkgBDT >= 0.531 && Ds_ifNANgiveX_isSignal_5==1 && Ds_diff_D0pi >= 0.1 && Ds_diff_D0pi <= 0.55");
   //=========================================================================================================================================================================//
   // B u i l d   M o d e l
   // -------------------------------------------------------------- 
   // Peaking 
   // -------------------------------------------------------------- 
   // First Gamma Function
   RooRealVar x00_Mode1("x00_Mode1","x00_Mode1", 1.4686e-01, xL, 0.2);
   RooRealVar a0_Mode1("a0_Mode1","a0_Mode1", 2.1711e+00, 1.0e+00, 5.0e+00);
   RooRealVar b0_Mode1("b0_Mode1","b0_Mode1", 6.3342e+01, 2.0e+01, 9.0e+01);    
   x00_Mode1.setConstant(true);
   a0_Mode1.setConstant(true);
   b0_Mode1.setConstant(true);
   RooFormulaVar x00prime_Mode1("#x00prime_Mode1", "x[0] < x[1] ? x[0] : x[1]", RooArgList(DeltaM, x00_Mode1));     
   RooFormulaVar gamma0_Mode1("#gamma0_Mode1", "a0_Mode1 + 1", RooArgList(a0_Mode1));    
   RooFormulaVar beta0_for_gamma_Mode1("beta0_for_gamma_Mode1", "1./b0_Mode1", RooArgList(b0_Mode1));    
   RooGamma GammaModel0_Mode1("GammaModel0_Mode1", "Gamma pdf", DeltaM, gamma0_Mode1, beta0_for_gamma_Mode1, x00prime_Mode1);

   // Second Gamma Function
   RooRealVar x01_Mode1("x01_Mode1","x01_Mode1",  1.3958e-01, xL, 0.2);
   RooRealVar a1_Mode1("a1_Mode1","a1_Mode1", 6.0323e+00, 2.0e+00, 9.0e+00);
   RooRealVar b1_Mode1("b1_Mode1","b1_Mode1", 8.0368e+01, 5.0e+01, 10.0e+01);     
   x01_Mode1.setConstant(true);
   a1_Mode1.setConstant(true);
   b1_Mode1.setConstant(true);
   RooFormulaVar x01prime_Mode1("#x01prime_Mode1", "x[0] < x[1] ? x[0] : x[1]", RooArgList(DeltaM, x01_Mode1));     
   RooFormulaVar gamma1_Mode1("#gamma1_Mode1", "a1_Mode1 + 1", RooArgList(a1_Mode1));    
   RooFormulaVar beta1_for_gamma_Mode1("beta1_for_gamma_Mode1", "1./b1_Mode1", RooArgList(b1_Mode1));    
   RooGamma GammaModel1_Mode1("GammaModel1_Mode1", "Gamma pdf", DeltaM, gamma1_Mode1, beta1_for_gamma_Mode1, x01prime_Mode1);

   // Second Gamma Function
   RooRealVar x02_Mode1("x02_Mode1","x02_Mode1",  1.1740e-01, xL, 0.2);
   RooRealVar a2_Mode1("a2_Mode1","a2_Mode1", 5.5652e-01, 2.0e-01, 7.0e-01);
   RooRealVar b2_Mode1("b2_Mode1","b2_Mode1", 4.7553e+01, 1.0e+01, 8.0e+01);     
   x02_Mode1.setConstant(true);
   a2_Mode1.setConstant(true);
   b2_Mode1.setConstant(true);
   RooFormulaVar x02prime_Mode1("#x01prime_Mode1", "x[0] < x[1] ? x[0] : x[1]", RooArgList(DeltaM, x02_Mode1));     
   RooFormulaVar gamma2_Mode1("#gamma2_Mode1", "a2_Mode1 + 1", RooArgList(a2_Mode1));    
   RooFormulaVar beta2_for_gamma_Mode1("beta2_for_gamma_Mode1", "1./b2_Mode1", RooArgList(b2_Mode1));    
   RooGamma GammaModel2_Mode1("GammaModel2_Mode1", "Gamma pdf", DeltaM, gamma2_Mode1, beta2_for_gamma_Mode1, x02prime_Mode1);

   // Add the components
   RooRealVar f0_Mode1("f0_Mode1", "f0", 2.0709e-01, 0.0e-01, 5.0e-01);
   RooRealVar f1_Mode1("f1_Mode1", "f1", 1.7216e-01, 0.0e-01, 5.0e-01);
   RooRealVar f2_Mode1("f2_Mode1", "f1", 2.4599e-01, 0.0e-01, 5.0e-01);
   f0_Mode1.setConstant(true);
   f1_Mode1.setConstant(true);
   f2_Mode1.setConstant(true);
	RooAddPdf Dstar0Mode1Model("Dstar0Mode1Model","g1+g2",RooArgList(GammaModel0_Mode1,GammaModel1_Mode1,GammaModel2_Mode1), RooArgList(f0_Mode1,f1_Mode1,f2_Mode1));

   // // -------------------------------------------------------------- 
   // // Combinatorial Part
   // // -------------------------------------------------------------- 
   // // First Attempt:
   // // RooRealVar B0_Mode1("B0_Mode1", "B0_Mode1", -40., 20.);
   // // RooRealVar P0_Mode1("P0_Mode1", "P0_Mode1", 0., 5.);
   // // RooGenericPdf comb_Mode1("comb_Mode1","Generic PDF","(Ds_diff_D0pi - 0.13953 > 0 ? ( ((Ds_diff_D0pi-0.13953)^P0_Mode1)*exp( B0_Mode1*(Ds_diff_D0pi-0.13953))) : 0.0)", RooArgSet( DeltaM, B0_Mode1, P0_Mode1));
   // // Second Attempt:
   // // RooRealVar c0("c0", "coefficient #0", -1, 1);
   // // RooRealVar c1("c1", "coefficient #1", -1, 1);
   // // RooRealVar c2("c2", "coefficient #2", -1, 1);
   // // RooPolynomial comb_Mode1("comb_Mode1", "background p.d.f.", DeltaM, RooArgList(c0,c1,c2));
   // // Third Attempt:
   // // First Gamma Function
   // RooRealVar x00_comb("x00_comb","x00_comb", 1.4043e-01, xL, 0.2);
   // RooRealVar a0_comb("a0_comb","a0_comb", 5.6201e-01, 2.0e-01, 7.0e-01);
   // RooRealVar b0_comb("b0_comb","b0_comb", 2.7753e+00, 1.0+00, 5.0+00);    
   // // x00_comb.setConstant(true);
   // // a0_comb.setConstant(true);
   // // b0_comb.setConstant(true);
   // RooFormulaVar x00prime_comb("#x00prime_comb", "x[0] < x[1] ? x[0] : x[1]", RooArgList(DeltaM, x00_comb));     
   // RooFormulaVar gamma0_comb("#gamma0_comb", "a0_comb + 1", RooArgList(a0_comb));    
   // RooFormulaVar beta0_for_gamma_comb("beta0_for_gamma_comb", "1./b0_comb", RooArgList(b0_comb));    
   // RooGamma GammaModel0_comb("GammaModel0_comb", "Gamma pdf", DeltaM, gamma0_comb, beta0_for_gamma_comb, x00prime_comb);

   // // Second Gamma Function
   // RooRealVar x01_comb("x01_comb","x01_comb", 1.4260e-01, xL, 0.2);
   // RooRealVar a1_comb("a1_comb","a1_comb", 3.1636e-01, 1.0e-01, 6.0e-01);
   // RooRealVar b1_comb("b1_comb","b1_comb", 4.1849e+01, 1.0e+01, 6.0e+01);     
   // // x01_comb.setConstant(true);
   // // a1_comb.setConstant(true);
   // // b1_comb.setConstant(true);
   // RooFormulaVar x01prime_comb("#x01prime_comb", "x[0] < x[1] ? x[0] : x[1]", RooArgList(DeltaM, x01_comb));     
   // RooFormulaVar gamma1_comb("#gamma1_comb", "a1_comb + 1", RooArgList(a1_comb));    
   // RooFormulaVar beta1_for_gamma_comb("beta1_for_gamma_comb", "1./b1_comb", RooArgList(b1_comb));    
   // RooGamma GammaModel1_comb("GammaModel1_comb", "Gamma pdf", DeltaM, gamma1_comb, beta1_for_gamma_comb, x01prime_comb);

   // // Second Gamma Function
   // RooRealVar x02_comb("x02_comb","x02_comb", 1.5000e-01, xL, 0.2);
   // RooRealVar a2_comb("a2_comb","a2_comb", 9.8836e+00, 5.0e+00, 12.0e+00);
   // RooRealVar b2_comb("b2_comb","b2_comb", 9.6892e+01, 5.0e+01, 12.0e+01);     
   // // x02_comb.setConstant(true);
   // // a2_comb.setConstant(true);
   // // b2_comb.setConstant(true);
   // RooFormulaVar x02prime_comb("#x01prime_comb", "x[0] < x[1] ? x[0] : x[1]", RooArgList(DeltaM, x02_comb));     
   // RooFormulaVar gamma2_comb("#gamma2_comb", "a2_comb + 1", RooArgList(a2_comb));    
   // RooFormulaVar beta2_for_gamma_comb("beta2_for_gamma_comb", "1./b2_comb", RooArgList(b2_comb));    
   // RooGamma GammaModel2_comb("GammaModel2_comb", "Gamma pdf", DeltaM, gamma2_comb, beta2_for_gamma_comb, x02prime_comb);

   // // Add the components
   // RooRealVar f0_comb("f0_comb", "f0", 9.4008e-01, 5.0e-01, 12.0e-01);
   // RooRealVar f1_comb("f1_comb", "f1", 4.0968e-06, 2.0e-06,  6.0e-06);
   // RooRealVar f2_comb("f2_comb", "f2", 8.2837e-02, 4.0e-02, 10.0e-02);
   // // f0_comb.setConstant(true);
   // // f1_comb.setConstant(true);
   // // f2_comb.setConstant(true);
   // RooAddPdf comb_Mode1("comb_Mode1","g0+g1+g2",RooArgList(GammaModel0_comb,GammaModel1_comb,GammaModel2_comb), RooArgList(f0_comb,f1_comb,f2_comb));

   // C o n s t r u c t   a   S i g na l   a n d   B a c k g r o u n d   P D F:
   // --------------------------------------------------------------
   RooRealVar Nsig_Mode1("nsig_Mode1", "#Signal events", 4.1247e+03, 1.0e+03, 8.0e+03);
   // RooRealVar Nbkg_Mode1("nbkg_Mode1", "#background events", 2.1807e+04, 1.0e+04, 5.0e+04);
   // RooAddPdf model("model", "g+a", RooArgList(PeakingBkgModel), RooArgList(Nsig_Mode1));
   RooAddPdf model("model", "g+a", RooArgList(Dstar0Mode1Model), RooArgList(Nsig_Mode1));
   // RooAddPdf model("model", "g+a", RooArgList(comb_Mode1), RooArgList(Nbkg_Mode1));
   // RooAddPdf model("model", "g+a", RooArgList(Dstar0Mode1Model, comb_Mode1), RooArgList(Nsig_Mode1, Nbkg_Mode1));
   //=========================================================================================================================================================================//
   // S e t   t h e   C o n s t a n t   f o r   t h e    V a r i a b l e s
   //--------------------------------------------------------------------------
   // The code below fixes the value for a certain value.
   // B0_Mode1.setVal();
   // P0_Mode1.setVal();
   // g1frac_Mode1.setVal();
   // g2frac_Mode1.setVal();
   // g3frac_Mode1.setVal();
   // sigmean1_Mode1.setVal();
   // sigmean2_Mode1.setVal();
   // sigmean3_Mode1.setVal();
   // sigwidth1_Mode1.setVal();
   // sigwidth2_Mode1.setVal();
   // sigwidth3_Mode1.setVal();

   // The code below removes the range and fixes the value.
   // B0_Mode1.setConstant(true);
   // P0_Mode1.setConstant(true);
   // g1frac_Mode1.setConstant(true);
   // g2frac_Mode1.setConstant(true);
   // g3frac_Mode1.setConstant(true);
   // sigmean1_Mode1.setConstant(true);
   // sigmean2_Mode1.setConstant(true);
   // sigmean3_Mode1.setConstant(true);
   // sigwidth1_Mode1.setConstant(true);
   // sigwidth2_Mode1.setConstant(true);
   // sigwidth3_Mode1.setConstant(true);

   // DeltaM.setVal(0.143);
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
   std::unique_ptr<RooFitResult> r{model.fitTo(*data_Mode1, Save(), PrintLevel(-1))};
   // If the fit is too broad and does not seem to fit a peak use the code below. Use it
   // once to determine where the peak should be and then adjust the mean above.
   // model.fitTo(data, Minos(kTRUE), PrintLevel(-1), Range(2,2.1));
   //=========================================================================================================================================================================//
   // P l o t   M o d e l   a n d   P r e p a r e   P a r a m e t e r   V a l u e s
   // ---------------------------------------------------------------------------
   // Create a RooPlot to draw on. We don't manage the memory of the returned
   // pointer. Instead we let it leak such that the plot still exists at the end of
   // the macro and we can take a look at it.
   RooPlot *frame = DeltaM.frame();
   *data_Mode1->plotOn(frame, Binning(bin));
   // model.plotOn(frame, Components(comb_Mode1),       LineStyle(kDashed), LineColor(kRed+1), RooFit::Name("Background"));  // Background in deep red
   model.plotOn(frame, Components(Dstar0Mode1Model), LineStyle(kDashed), LineColor(kBlue+1), RooFit::Name("Signal"));  // Signal in deep blue
   model.plotOn(frame, LineColor(kGreen+3), LineWidth(2), RooFit::Name("model"));  // Total fit in dark green
   //=========================================================================================================================================================================//
   // C a l c u l a t e   c h i ^ 2
   // ------------------------------
   // Show the chi^2 of the curve w.r.t. the histogram
   // If multiple curves or datasets live in the frame you can specify
   // the name of the relevant curve and/or dataset in chiSquare()
   Int_t npar = model.getParameters(*data_Mode1)->selectByAttrib("Constant",kFALSE)->getSize(); //select floating parameters and count their number

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
   lum.Form("#splitline{Simulation}{200k Events}");
   // lum.Form("#splitline{Generic c#bar{c} Events}{#scale[0.5]{#int}#it{L} dt=200 fb^{-1}}");
   latex->SetNDC();
   latex->SetTextSize(0.03); 
   
   auto PerBin = ((xR-xL)/(bin))*1000;
   // TString title ("");
   TString y_axis ("");
   // TString x_axis ("");

   // title.Form("Generic d#bar{d} Events");
   y_axis.Form("Entries/(%.2f MeV/c^{2})", PerBin);
   // x_axis.Form("#Delta m  [GeV/c^{2}]");

   frame->SetTitle("");
   frame->GetYaxis()->SetTitle(y_axis);
   frame->GetYaxis()->SetRangeUser(0,800);
   frame->GetXaxis()->SetTitle("");
   // frame->GetXaxis()->CenterTitle(true);
   frame->Draw();

   // The line below prints the luminosity onto the screen.
   // latex->DrawLatex(0.6,0.2, lum);
   latex->DrawLatex(0.18,0.89, lum);
   //=========================================================================================================================================================================//
   // DeltaM.setRange("window",0.15,1) ;
   // RooAbsReal* fracSigRange = Signal.createIntegral(DeltaM,DeltaM,"window") ; 
   // Double_t NWindow = N.getVal() * fracSigRange->getVal() ;
   // Double_t NWindow_error = N.getError() * fracSigRange->getVal() ;
   //=========================================================================================================================================================================//
   // P a r a m e t e r s
   // ---------------------------------------------------------------------------
   // The code below grabs the parameters and then places them in the legend
   // Float_t N_value = Nsig_Mode1.getValV();
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
   auto x1 = xIF1 + 0.4*dx;
   auto x2 = xIF2;
   auto y1 = yIF1 + 0.75*dy;
   auto y2 = 0.93;
   auto legend = new TLegend(x1, y1, x2, y2);

   TString chi2("");
   TString Nsig_para("");
   TString Nbkg_para("");

   chi2.Form("#chi^{2}/ndf = %.2f", chi2ndf);
   Nsig_para.Form("N_{sig} = %.2f #pm %.2f", Nsig_Mode1.getValV(), Nsig_Mode1.getError());
   // Nbkg_para.Form("N_{bkg} = %.2f #pm %.2f", Nbkg_Mode1.getValV(), Nbkg_Mode1.getError());

   legend->AddEntry(frame->FindObject("data"),"data","pe");
   legend->AddEntry((TObject *)0, chi2, "");
   legend->AddEntry("Signal", "D^{*+} #rightarrow [D^{0} #rightarrow K^{-} #pi^{+}] #pi^{+}", "pl" );
   legend->AddEntry((TObject *)0, Nsig_para, "");
   legend->AddEntry("Background", "Background", "pl" );
   legend->AddEntry((TObject *)0, Nbkg_para, "");
   legend->AddEntry("model", "Total Fit", "pl" );
   legend->AddEntry((TObject *)0, " ", "");

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
   // B0_Mode1.setConstant(true);
   // B0_Mode1.removeError();
   // P0_Mode1.setConstant(true);
   // P0_Mode1.removeError();
   // g1frac_Mode1.setConstant(true);
   // g1frac_Mode1.removeError();
   // g2frac_Mode1.setConstant(true);
   // g2frac_Mode1.removeError();
   // g3frac_Mode1.setConstant(true);
   // g3frac_Mode1.removeError();
   // sigmean1_Mode1.setConstant(true);
   // sigmean1_Mode1.removeError();
   // sigmean2_Mode1.setConstant(true);
   // sigmean2_Mode1.removeError();
   // sigmean3_Mode1.setConstant(true);
   // sigmean3_Mode1.removeError();
   // sigwidth1_Mode1.setConstant(true);
   // sigwidth1_Mode1.removeError();
   // sigwidth2_Mode1.setConstant(true);
   // sigwidth2_Mode1.removeError();
   // sigwidth3_Mode1.setConstant(true);
   // sigwidth3_Mode1.removeError();
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
   RooWorkspace *w = new RooWorkspace("w_Mode1", "workspace");
   // Import model and all its components into the workspace
   w->import(Dstar0Mode1Model);

   // E n c o d e   d e f i n i t i o n   o f   p a r a m e t e r s   i n   w o r k s p a c e
   // ---------------------------------------------------------------------------------------
   // Define named sets "parameters" and "observables", which list which variables should be considered
   // parameters and observables by the users convention
   //
   // Variables appearing in sets _must_ live in the workspace already, or the autoImport flag
   // of defineSet must be set to import them on the fly. Named sets contain only references
   // to the original variables, therefore the value of observables in named sets already
   // reflect their 'current' value
   RooArgSet *params = Dstar0Mode1Model.getVariables();
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
   w->saveSnapshot("reference_fit_Mode1", *params, true);

   // Print workspace contents
   w->Print();

   // S a v e   w o r k s p a c e   i n   f i l e
   // -------------------------------------------
   // Save the workspace into a ROOT file
   w->writeToFile("/home/belle2/amubarak/Ds2D0enue_Analysis/05-ROOT/Roofit_New/Fit_Parameter/rf02-Dstar0_Mode1.root");
   //=========================================================================================================================================================================//
   c->SaveAs("/home/belle2/amubarak/Ds2D0enue_Analysis/05-ROOT/Roofit_New/Fit_Image/rf02-Dstar0_Mode1.png");

   // cout << "Total Number After Peak = " << NWindow << endl;
   //=========================================================================================================================================================================//
   // Final Fit Status
   //---------------------------
   int fitStatus = r->status();  // Get the fit status
   int covQual = r->covQual();   // Get covariance quality

   std::cout << "Fit Status: " << fitStatus << std::endl;
   std::cout << "Covariance Matrix Quality: " << covQual << std::endl << endl;
}