#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
using namespace RooFit;
using namespace TMath;

void rf06_D0_dM()
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
   // I n i t i a l i z a t i o n
   // --------------------------------------------------------------
   // The lines below initializes all my values
   int bin = 150;
   double xL = -0.05;
   double xR = 0.05;
   //=========================================================================================================================================================================//
   // I m p o r t s   D a t a   f r o m   t r e e   b r a n c h:
   // --------------------------------------------------------------
   TFile *file = new TFile("C01-Simulated_Events/Ds2D0enu-Signal.root");
   TTree *tree = (TTree*)file->Get("D02kmpiptree");
   //=========================================================================================================================================================================//
   RooRealVar M("D0_kmpip_dM", "m(D0)", xL, xR);

   RooDataSet data("data", "dataset with mass", tree, M);
   // RooDataSet data("data", "dataset with mass", tree, RooArgSet(M,z));
   //=========================================================================================================================================================================//
	// C r e a t e   m o d e l   f o r   s i g n a l   s a m p l e
   // -------------------------------------------------------------- 
   RooRealVar sig_mean("sig_mean", "mean of crystal", -0.005, 0.005);
   RooRealVar sig_sigmaL("sig_sigmaL", "sig_sigmaL", 1e-4, 0.02);
   RooRealVar sig_sigmaR("sig_sigmaR", "sig_sigmaR", 1e-4, 0.02);
   RooRealVar sig_alphaL("sig_alphaL", "sig_alphaL", 1e-4, 10.);
   RooRealVar sig_alphaR("sig_alphaR", "sig_alphaR", 1e-4, 10.);
   RooRealVar sig_nL("sig_nL", "sig_nL", 1e-4, 30);
   RooRealVar sig_nR("sig_nR", "sig_nR", 1e-4, 30);
   //  sig_mean.setConstant();
   //  sig_sigmaL.setConstant();
   //  sig_sigmaR.setConstant();
   //  sig_alphaL.setConstant();
   //  sig_alphaR.setConstant();
   //  sig_nL.setConstant();
   //  sig_nR.setConstant();
   RooCrystalBall dcb("dcb", "sig4", M, sig_mean,sig_sigmaL,sig_sigmaR,sig_alphaL,sig_nL,sig_alphaR,sig_nR);

   // Gaussian Resolution model
	RooRealVar sigmean1("sigmean1", "mean", -0.005, 0.005);
	RooRealVar sigwidth1("sigwidth1", "width", 0, 0.01);
   // sigmean1.setConstant(true);
   // sigwidth1.setConstant(true);
	RooGaussModel gaussm1("gaussm1", "Signal PDF", M, sig_mean, sigwidth1);

	// Gaussian Resolution model
	RooRealVar sigmean2("sigmean2", "mean", -0.005, 0.005);
	RooRealVar sigwidth2("sigwidth2", "width", 0, 0.01);
   // sigmean2.setConstant(true);
   // sigwidth2.setConstant(true);
	RooGaussModel gaussm2("gaussm2", "Signal PDF", M, sig_mean, sigwidth2);

   RooRealVar c0("c0", "coefficient #0", -100, 100);
   RooRealVar c1("c1", "coefficient #1", -100, 100);
   RooRealVar c2("c2", "coefficient #1", -100, 100);
   RooRealVar c3("c3", "coefficient #1", -100, 100);
   // RooRealVar c2("c2", "coefficient #2", -0.99, 1);
   // RooRealVar c3("c3", "coefficient #3", -0.99, 1);
   // RooChebychev comb("comb", "background p.d.f.", x, RooArgList(c0));
   RooPolynomial comb("comb", "background p.d.f.", M, RooArgList(c0,c1,c2,c3));

	// Add the components
	RooRealVar g1frac("g1frac","fraction of gauss1", 0., 1.);
	RooRealVar g2frac("g2frac","fraction of gauss2", 0., 1.);
   RooRealVar g3frac("g3frac","fraction of gauss3", 0., 1.);
   // g1frac.setConstant(true);
   // g2frac.setConstant(true);
	RooAddPdf Signal("Signal","g1+g2",RooArgList(gaussm1,gaussm2,dcb), RooArgList(g1frac,g2frac,g3frac));
   //=========================================================================================================================================================================//
   // C r e a t e   m o d e l   f o r   b a c k g r o u n d   s a m p l e
   // -------------------------------------------------------------- 
   // Constant Background:
   RooRealVar c5("c5", "coefficient #5", 0., 1);
   // c5.setConstant(true);
   // RooPolynomial Polynomial("Polynomial", "background p.d.f.", x, c5);
   RooPolynomial Background("Background", "background p.d.f.", M, c5);
   //=========================================================================================================================================================================//
   // C o n s t r u c t   a   S i g na l   a n d   B a c k g r o u n d   P D F:
   // --------------------------------------------------------------
   RooRealVar nsig("nsig", "#Signal events", 0., 1e+08);
   RooRealVar nbkg("nbkg", "#background events", 0., 1e+08);
   RooAddPdf model("model", "g+a", RooArgList(Signal, Background), RooArgList(nsig, nbkg));
   //=========================================================================================================================================================================//
   // F i t
   //-------------
   // Perform extended ML fit of data:
   std::unique_ptr<RooFitResult> r{model.fitTo(data, Save(), PrintLevel(-1), Range(xL,xR))};
   // If the fit is too broad and does not seem to fit a peak use the code below. Use it
   // once to determine where the peak should be and then adjust the mean above.
   // model.fitTo(data, Minos(kTRUE), PrintLevel(-1), Range(2,2.1));
   //=========================================================================================================================================================================//
   // P l o t   M o d e l   a n d   P r e p a r e   P a r a m e t e r   V a l u e s
   // ---------------------------------------------------------------------------
   // Create a RooPlot to draw on. We don't manage the memory of the returned
   // pointer. Instead we let it leak such that the plot still exists at the end of
   // the macro and we can take a look at it.
   RooPlot *frame = M.frame();
   data.plotOn(frame, Binning(50));

   // model.plotOn(frame, Components(*SignalModel),     LineStyle(kDashed), LineColor(kAzure-2),    RooFit::Name("Signal"));
   model.plotOn(frame, Components(Signal),   LineStyle(kDashed), LineColor(kRed+2),    RooFit::Name("Peak"));
   model.plotOn(frame, Components(Background),  LineStyle(kDashed), LineColor(kGreen+3),   RooFit::Name("Comb"));

   model.plotOn(frame, LineColor(kBlue+2), RooFit::Name("model"));
   //=========================================================================================================================================================================//
   // C a l c u l a t e   c h i ^ 2
   // ------------------------------
   // Show the chi^2 of the curve w.r.t. the histogram
   // If multiple curves or datasets live in the frame you can specify
   // the name of the relevant curve and/or dataset in chiSquare()
   Int_t npar = model.getParameters(data)->selectByAttrib("Constant",kFALSE)->getSize(); //select floating parameters and count their number

   cout << "NDF = " << npar << endl;
   cout << "chi^2/NDF = " << frame->chiSquare(npar) << endl;

   Double_t chi2ndf = frame->chiSquare(npar);
   //=========================================================================================================================================================================//
   // P l o t    L a b e l s
   // ---------------------------------------------------------------------------
   // Add text to frame
   TLatex *latex = new TLatex();
   TString lum ("");
   lum.Form("#splitline{Signal Events}{2M Events}");
   latex->SetNDC();
   latex->SetTextSize(0.035); 
   
   // // Add "Likelihood fit" below "2M Events"
   // TString Likelihood ("");
   // Likelihood.Form("#splitline{}{Likelihood Fit}");

   auto PerBin = ((xR-xL)/(bin))*1000;
   // TString title ("");
   TString y_axis ("");
   // TString x_axis ("");

   // title.Form("Generic d#bar{d} Events");
   y_axis.Form("Entries/(%.2f MeV/c^{2})", PerBin);
   // x_axis.Form("#Delta m  [GeV/c^{2}]");

   frame->SetTitle("");
   frame->GetYaxis()->SetTitle(y_axis);
   frame->GetYaxis()->SetRangeUser(0,95000);
   frame->GetXaxis()->SetTitle("");
   // frame->GetXaxis()->CenterTitle(true);
   frame->Draw();

   // The line below prints the luminosity onto the screen.
   // latex->DrawLatex(0.6,0.2, lum);
   latex->DrawLatex(0.18,0.89, lum);
   // latex->DrawLatex(0.18,0.85, Likelihood);
   //=========================================================================================================================================================================//
   // DeltaM.setRange("window",0.15,1) ;
   // RooAbsReal* fracSigRange = Signal.createIntegral(DeltaM,DeltaM,"window") ; 
   // Double_t NWindow = N.getVal() * fracSigRange->getVal() ;
   // Double_t NWindow_error = N.getError() * fracSigRange->getVal() ;
   //=========================================================================================================================================================================//
   // S t a n d a r d   D e v i a t i o n
   // ---------------------------------------------------------------------------
   // Extract uncertainties
   double df1 = g1frac.getError();
   double df2 = g2frac.getError();
   double ds1 = sigwidth1.getError();
   double ds2 = sigwidth2.getError();
   double dsL = sig_sigmaL.getError();
   double dsR = sig_sigmaR.getError();

   // Extract values
   double f1 = g1frac.getVal();                     // gauss1 fraction
   double f2 = g2frac.getVal();                     // gauss2 fraction
   double f3 = 1.0 - f1 - f2;                        // dcb fraction

   double sigma1 = sigwidth1.getVal();              // gauss1 sigma
   double sigma2 = sigwidth2.getVal();              // gauss2 sigma
   double sigma_cb_L = sig_sigmaL.getVal();         // left width of DCB
   double sigma_cb_R = sig_sigmaR.getVal();         // right width of DCB

   // Estimate effective sigma for DCB (weighted by symmetry)
   double sigma_cb = 0.5 * (sigma_cb_L + sigma_cb_R);
   double ds_cb = 0.5 * sqrt(dsL*dsL + dsR*dsR);  // error on average

   // Variance formula
   double dvar_df1 = sigma1 * sigma1 - sigma_cb * sigma_cb;
   double dvar_df2 = sigma2 * sigma2 - sigma_cb * sigma_cb;
   double dvar_ds1 = 2 * f1 * sigma1;
   double dvar_ds2 = 2 * f2 * sigma2;
   double dvar_ds_cb = 2 * f3 * sigma_cb;

   // Error propagation on variance
   double var_error_sq =
      (dvar_df1 * df1) * (dvar_df1 * df1) +
      (dvar_df2 * df2) * (dvar_df2 * df2) +
      (dvar_ds1 * ds1) * (dvar_ds1 * ds1) +
      (dvar_ds2 * ds2) * (dvar_ds2 * ds2) +
      (dvar_ds_cb * ds_cb) * (dvar_ds_cb * ds_cb);

   // Compute weighted variance
   double var_total = f1 * (sigma1 * sigma1) + f2 * (sigma2 * sigma2) + f3 * (sigma_cb * sigma_cb);
   double sigma_total = sqrt(var_total);
   double sigma_error = 0.5 * var_error_sq / sigma_total;

   std::cout << "-------------------------------------------" << std::endl;
   std::cout << "Combined Signal σ (from fractions): " << sigma_total << std::endl;
   std::cout << "-------------------------------------------" << std::endl;
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
   TString mu("");
   TString sigma("");

   chi2.Form("#chi^{2}/ndf = %.2f", chi2ndf);
   mu.Form("#mu = %.2f #pm %.2f", sig_mean.getValV(), sig_mean.getError());
   sigma.Form("#sigma_{eff} = %.2f #pm %.2f", sigma_total, sigma_error);
   // b_para.Form("b = %.6f #pm %.6f", b, b.getError());
   // P_para.Form("P = %.6f #pm %.6f", P, P.getError());

   legend->AddEntry(frame->FindObject("*data_sig"),"data","pe");
   legend->AddEntry((TObject *)0, chi2, "");
   legend->AddEntry("Peak", "Signal", "pl" );
   legend->AddEntry((TObject *)0, mu, "");
   legend->AddEntry("Comb", "Background", "pl" );
   legend->AddEntry((TObject *)0, sigma, "");
   // legend->AddEntry((TObject *)0, " ", "");
   // legend->AddEntry((TObject *)0, b_para, "");
   // legend->AddEntry((TObject *)0, " ", "");
   // legend->AddEntry((TObject *)0, P_para, "");

   legend -> SetBorderSize(0);
   legend -> SetTextFont(132);
   legend -> SetTextSize(0.04);
   legend -> SetNColumns(2);
   legend -> SetMargin(0.1);
   legend->Draw();
   //=========================================================================================================================================================================//
   // S h o w   P a r a m e t e r/ F i t    R e s u l t s
   // -------------------------------------------------------

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
   // S h o w   P u l l   d i s t s
   // -------------------------------------------------------
   pad2->cd();
   // Construct a histogram with the residuals of the data w.r.t. the curve
   RooHist *hpull = frame->pullHist();
   for(int iPoint = 0; iPoint < hpull->GetN(); ++iPoint) {
      hpull->SetPointEYlow(iPoint, 0.0);
      hpull->SetPointEYhigh(iPoint, 0.0);
   }
   hpull->SetFillColor(38);
   // hpull->SetLineWidth(-802);
   // hpull->SetFillStyle(3002);
   // hpull->SetFillColor(2);

   // Create a new frame to draw the residual distribution and add the distribution to the frame
   RooPlot *frame1 = M.frame(Title("Pull Distribution"));
   frame1->addPlotable(hpull, "B");
   frame1->Draw();
   
   // Plot Titles
   //---------------
   frame1->SetTitle("");

   // Y-axis
   frame1->SetYTitle("Pull");
   frame1->GetYaxis()->CenterTitle(true);
   frame1->GetYaxis()->SetTitleSize(0.10);
   frame1->GetYaxis()->SetLabelSize(0.12);

   // X-axis
   TString x_axis ("");
   x_axis.Form("m(D^{0})-m_{PDG}(D^{0})  [GeV/c^{2}]");
   frame1->SetXTitle(x_axis);
   frame1->GetXaxis()->CenterTitle(true);
   frame1->GetXaxis()->SetTitleSize(0.13);
   frame1->GetXaxis()->SetLabelSize(0.10);

   frame1->SetTitleOffset(0.5, "Y");

   // Draw Horizontal Line On Residual Plot
   //----------------------------------------
   // double xmin1 = xIF1; 
   // double xmax1 = 1.0;
   // TLine *line11 = new TLine(xIF1,0.0,xIF2,0.0);
   // line11->SetLineColor(kRed); line11->SetLineWidth(3);
   // line11->Draw("SAME");

   c->cd();
   //=========================================================================================================================================================================//
   // S a v e   F i l e s
   // -------------------------------------------------------------
   c->SaveAs("/home/belle2/amubarak/Ds2D0enue_Analysis/05-ROOT/Roofit_Original/Fit_Image/D0_dM_SignalRegion.png");
   //=========================================================================================================================================================================//
   // Final Fit Status
   //---------------------------
   int fitStatus = r->status();  // Get the fit status
   int covQual = r->covQual();   // Get covariance quality

   std::cout << "Fit Status: " << fitStatus << std::endl;
   std::cout << "Covariance Matrix Quality: " << covQual << std::endl << endl;
}