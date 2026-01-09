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

void rf10_Signal()
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
   double xL = 0.;
   double xR = 0.25;
   //=========================================================================================================================================================================//
   // I m p o r t s   D a t a   f r o m   t r e e   b r a n c h:
   // --------------------------------------------------------------
   TFile *file = new TFile("C01-Simulated_Events/Ds2D0enu-Signal.root");
   TTree *tree = (TTree*)file->Get("Dstree");
   //=========================================================================================================================================================================//
   // D a t a:
   // --------------------------------------------------------------
   RooRealVar DeltaM("Ds_massDifference_0", "Mass Difference", xL, xR);
   RooRealVar Veto("Ds_gammaveto_M_Correction","Veto",0.1,100);
   RooRealVar BCS("Ds_Ds_chiProb_rank","BCS",-1,1.25);
   RooRealVar dM("D0_dM", "Mass Difference", -0.02, 0.02);

   RooDataSet data_sig("data_sig", "dataset with mass", tree, RooArgSet(DeltaM,Veto,BCS,dM));
   // RooDataSet data("data", "dataset with mass", tree, RooArgSet(M,z));
   //=========================================================================================================================================================================//
   // B u i l d   M o d e l
   // -------------------------------------------------------------- 
   // Gaussian Resolution model
	RooRealVar sigmean1("sigmean1", "mean", 0, 0.1);
	RooRealVar sigwidth1("sigwidth1", "width", 1.0687e-02, 0.0, 0.1);
	RooGaussModel gaussm1("gaussm1", "Signal PDF", DeltaM, sigmean1, sigwidth1);
	// Gaussian Resolution model
	RooRealVar sigmean2("sigmean2", "mean", 0, 0.1);
	RooRealVar sigwidth2("sigwidth2", "width", 2.6965e-02, 0.0, 0.1);
	RooGaussModel gaussm2("gaussm2", "Signal PDF", DeltaM, sigmean2, sigwidth2);
   // Gaussian Resolution model
	RooRealVar sigmean3("sigmean3", "mean", 0, 0.1);
	RooRealVar sigwidth3("sigwidth3", "width", 1.7712e-02, 0.0, 0.1);
	RooGaussModel gaussm3("gaussm3", "Signal PDF", DeltaM, sigmean3, sigwidth3);

	// Add the components
	RooRealVar g1frac("g1frac","fraction of gauss1", 0., 1.);
	RooRealVar g2frac("g2frac","fraction of gauss2", 0., 1.);
   RooRealVar g3frac("g3frac","fraction of gauss2", 0., 1.);
	RooAddPdf SignalModel("SignalModel","g1+g2",RooArgList(gaussm1,gaussm2,gaussm3), RooArgList(g1frac,g2frac,g3frac));

   // Comb.
   RooRealVar x00_S("x00_S","x00_S", xL, xL, xR);
   RooRealVar a0_S("a0_S","a0_S",  0.0, 50.);
   RooRealVar b0_S("b0_S","b0_S", 0.0, 100.);    
   x00_S.setConstant(true);
   // a0_S.setConstant(true);
   // b0_S.setConstant(true);
   RooFormulaVar x00prime_S("#x00prime_S", "x[0] < x[1] ? x[0] : x[1]", RooArgList(DeltaM, x00_S));     
   RooFormulaVar gamma0_S("#gamma0_S", "a0_S + 1", RooArgList(a0_S));    
   RooFormulaVar beta0_for_gamma_S("beta0_for_gamma_S", "1./b0_S", RooArgList(b0_S));    
   RooGamma BkgModel("BkgModel", "Gamma pdf", DeltaM, gamma0_S, beta0_for_gamma_S, x00prime_S);

   // F i n a l    S i g n a l   C o u n t
   //-----------------------------------------
   RooRealVar Nsig_S("Nsig", "#Signal events", 2.4233e+05, 0., 500000);
   RooRealVar Nbkg_S("Nbkg", "#Signal events", 2.4233e+05, 0., 500000);
   RooAddPdf model("model", "g+a", RooArgList(SignalModel,BkgModel), RooArgList(Nsig_S,Nbkg_S));
   //=========================================================================================================================================================================//
   // S e t   t h e   C o n s t a n t   f o r   t h e    V a r i a b l e s
   //--------------------------------------------------------------------------
   // // The code below fixes the value for a certain value.
   // B0_S.setVal(-3.9039e+01);
   // B1_S.setVal(-4.7587e+01);
   // P0_S.setVal(4.8075e-01);
   // P1_S.setVal(2.7543e+00);
   // f0_S.setVal(1.0394e+01);
   // f1_S.setVal(8.1720e+00);

   // // // The code below removes the range and fixes the value.
   // B0_S.setConstant(true);
   // B1_S.setConstant(true);
   // P0_S.setConstant(true);
   // P1_S.setConstant(true);
   // f0_S.setConstant(true);
   // f1_S.setConstant(true);

   // DeltaM.setVal(0.0065);
   // DeltaM.setConstant(true);
   //=========================================================================================================================================================================//
   // F i t
   //-------------
   // Perform extended ML fit of data:
   std::unique_ptr<RooFitResult> r{model.fitTo(data_sig, Save(), PrintLevel(-1))};
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
   data_sig.plotOn(frame, Binning(bin));
   model.plotOn(frame, Components(SignalModel),       LineStyle(kDashed), LineColor(46), RooFit::Name("Background"));
   model.plotOn(frame, Components(BkgModel), LineStyle(kDashed), LineColor(28), RooFit::Name("Signal"));
   model.plotOn(frame, LineColor(28), RooFit::Name("model"));
   //=========================================================================================================================================================================//
   // C a l c u l a t e   c h i ^ 2
   // ------------------------------
   // Show the chi^2 of the curve w.r.t. the histogram
   // If multiple curves or datasets live in the frame you can specify
   // the name of the relevant curve and/or dataset in chiSquare()
   Int_t npar = model.getParameters(data_sig)->selectByAttrib("Constant",kFALSE)->getSize(); //select floating parameters and count their number

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
   
   auto PerBin = ((xR-xL)/(bin))*1000;
   // TString title ("");
   TString y_axis ("");
   // TString x_axis ("");

   // title.Form("Generic d#bar{d} Events");
   y_axis.Form("Entries/(%.2f MeV/c^{2})", PerBin);
   // x_axis.Form("#Delta m  [GeV/c^{2}]");

   frame->SetTitle("");
   frame->GetYaxis()->SetTitle(y_axis);
   frame->GetYaxis()->SetRangeUser(0,22000);
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
   // Float_t N_value = Nsig_S.getValV();
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
   TString Ntot_para("");
   TString N_para("");
   TString b_para("");
   TString P_para("");

   chi2.Form("#chi^{2}/ndf = %.2f", chi2ndf);
   Ntot_para.Form("N_{sig} = %.2f #pm %.2f", Nsig_S.getValV(), Nsig_S.getError());
   // N_para.Form("N^{*} = %.2f #pm %.2f", NWindow, NWindow_error);
   // b_para.Form("b = %.6f #pm %.6f", b_value, b.getError());
   // P_para.Form("P = %.6f #pm %.6f", P_value, P.getError());

   legend->AddEntry(frame->FindObject("data_sig"),"data","pe");
   legend->AddEntry((TObject *)0, chi2, "");
   legend->AddEntry("Signal", "Total Fit", "pl" );
   legend->AddEntry((TObject *)0, Ntot_para, "");
   // legend->AddEntry((TObject *)0, " ", "");
   // legend->AddEntry((TObject *)0, b_para, "");
   // legend->AddEntry((TObject *)0, " ", "");
   // legend->AddEntry((TObject *)0, P_para, "");

   legend -> SetBorderSize(0);
   legend -> SetTextFont(132);
   legend -> SetTextSize(0.03);
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

   // X-axis
   TString x_axis ("");
   x_axis.Form("#Delta m  [GeV/c^{2}]");
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
   // C r e a t e   w o r k s p a c e ,   i m p o r t   d a t a   a n d   m o d e l
   // -----------------------------------------------------------------------------
   // Create a new empty workspace
   RooWorkspace *w = new RooWorkspace("w_S", "workspace");
   // Import model and all its components into the workspace
   w->import(SignalModel);

   // E n c o d e   d e f i n i t i o n   o f   p a r a m e t e r s   i n   w o r k s p a c e
   // ---------------------------------------------------------------------------------------
   // Define named sets "parameters" and "observables", which list which variables should be considered
   // parameters and observables by the users convention
   //
   // Variables appearing in sets _must_ live in the workspace already, or the autoImport flag
   // of defineSet must be set to import them on the fly. Named sets contain only references
   // to the original variables, therefore the value of observables in named sets already
   // reflect their 'current' value
   RooArgSet *params = SignalModel.getVariables();
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
   w->saveSnapshot("reference_fit_Signal", *params, true);

   // Print workspace contents
   w->Print();

   // S a v e   w o r k s p a c e   i n   f i l e
   // -------------------------------------------
   // Save the workspace into a ROOT file
   w->writeToFile("Ds2D0enue_Analysis/05-ROOT/Fit_Parameter/rf00_SignalPara.root");
   //=========================================================================================================================================================================//
   c->SaveAs("Ds2D0enue_Analysis/10-Images/rf10_DeltaM_Signal.png");

   cout << endl << "Efficiency = " << Nsig_S.getValV()/(2000000) << " +/- " << Nsig_S.getError()/(200000) << endl;
}