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

void rs03_D0()
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
   TFile* file = new TFile("/home/belle2/amubarak/C03-Grid/TopoAna_Binned.root");
   TH1F* hh = (TH1F*)file->Get("hist_D0");
   //=========================================================================================================================================================================//
   // D a t a:
   // --------------------------------------------------------------
   RooRealVar DeltaM("Ds_diff_D0pi", "Mass Difference", xL, xR);
   RooDataHist data_D0("data_D0", "dataset with mass", DeltaM, hh);
   //=========================================================================================================================================================================//
   // B u i l d   M o d e l
   // -------------------------------------------------------------- 
   // // First Try
   // //--------------------------------
   // // Generic PDF
   // RooRealVar B0_S("B0_S", "B0_S", -100., 10.);
   // RooRealVar P0_S("P0_S", "P0_S", 0., 5.);
   // RooGenericPdf Generic0_S("Generic0_S","Generic PDF","(Ds_diff_D0pi - 0.13953 > 0 ? ( ((Ds_diff_D0pi-0.13953)^P0_S)*exp( B0_S*(Ds_diff_D0pi-0.13953))) : 0.0)", RooArgSet( DeltaM, B0_S, P0_S));

   // // Generic PDF
   // RooRealVar B1_S("B1_S", "B1_S", -100., 10.);
   // RooRealVar P1_S("P1_S", "P1_S", 0., 5.);
   // RooGenericPdf Generic1_S("Generic1_S","Generic PDF","(Ds_diff_D0pi - 0.13953 > 0 ? ( ((Ds_diff_D0pi-0.13953)^P1_S)*exp( B1_S*(Ds_diff_D0pi-0.13953))) : 0.0)", RooArgSet( DeltaM, B1_S, P1_S));

   // // Combined Background
   // RooRealVar f0_S("f0_S", "f0", 1.0394e+01, 0, 20);
   // RooRealVar f1_S("f1_S", "f1", 8.1720e+00, 0, 20);
   // f0_S.setConstant(true);
   // f1_S.setConstant(true);
   // RooAddPdf SignalModel("SignalModel", "b+c", RooArgList(Generic0_S,Generic1_S), RooArgList(f0_S,f1_S));

   // // Second Try
   // //--------------------------------
   // RooRealVar mean("mean", "mean", -1., 1.);
   // RooRealVar sigma("sigma", "sigma", 0., 0.5);
   // RooLandau SignalModel("SignalModel","Signal", DeltaM, mean, sigma);

   // // Third Try
   // //--------------------------------
   // RooRealVar m0_S("m0_S", "m0", 1e-4, 1.);
   // RooRealVar k_S("k_S", "k", 1e-4, 2.);
   // RooLognormal SignalModel("SignalModel","Signal", DeltaM, m0_S, k_S);

   // Second Attempt
   //---------------------
   // First Gamma Function
   RooRealVar x00_D0Bkg("x00_D0Bkg","x00_D0Bkg", 1.3980e-01, 0.13, 0.2);
   RooRealVar a0_D0Bkg("a0_D0Bkg","a0_D0Bkg", 3.7829e-01, 2.0e-01, 6.0e-01);
   RooRealVar b0_D0Bkg("b0_D0Bkg","b0_D0Bkg", 1.7939e+00, 1.0e+00, 4.0e+00);    
   x00_D0Bkg.setConstant();
   a0_D0Bkg.setConstant();
   b0_D0Bkg.setConstant();
   RooFormulaVar x00prime_D0Bkg("#x00prime_D0Bkg", "x[0] < x[1] ? x[0] : x[1]", RooArgList(DeltaM, x00_D0Bkg));     
   RooFormulaVar gamma0_D0Bkg("#gamma0_D0Bkg", "a0_D0Bkg + 1", RooArgList(a0_D0Bkg));    
   RooFormulaVar beta0_for_gamma_D0Bkg("beta0_for_gamma_D0Bkg", "1./b0_D0Bkg", RooArgList(b0_D0Bkg));    
   RooGamma GammaModel0_D0Bkg("GammaModel1_D0Bkg", "Gamma pdf", DeltaM, gamma0_D0Bkg, beta0_for_gamma_D0Bkg, x00prime_D0Bkg);

   // Second Gamma Function
   RooRealVar x01_D0Bkg("x01_D0Bkg","x01_D0Bkg", 2.3937e-01, 0.13, 0.25);
   RooRealVar a1_D0Bkg("a1_D0Bkg","a1_D0Bkg", 1.0241e+01, 0.0e+01, 3.0e+01);
   RooRealVar b1_D0Bkg("b1_D0Bkg","b1_D0Bkg", 6.0650e+00, 4.0e+00, 8.0e+00);     
   x01_D0Bkg.setConstant();
   a1_D0Bkg.setConstant();
   b1_D0Bkg.setConstant();
   RooFormulaVar x01prime_D0Bkg("#x01prime_S", "x[0] < x[1] ? x[0] : x[1]", RooArgList(DeltaM, x01_D0Bkg));     
   RooFormulaVar gamma1_D0Bkg("#gamma1_D0Bkg", "a1_D0Bkg + 1", RooArgList(a1_D0Bkg));    
   RooFormulaVar beta1_for_gamma_D0Bkg("beta1_for_gamma_D0Bkg", "1./b1_D0Bkg", RooArgList(b1_D0Bkg));    
   RooGamma GammaModel1_D0Bkg("GammaModel2_D0Bkg", "Gamma pdf", DeltaM, gamma1_D0Bkg, beta1_for_gamma_D0Bkg, x01prime_D0Bkg);

   // Third Gamma Function
   RooRealVar x02_D0Bkg("x02_D0Bkg","x02_D0Bkg", 1.4278e-01, 0.13, 0.2);
   RooRealVar a2_D0Bkg("a2_D0Bkg","a2_D0Bkg", 1.5341e+00, 1.0e+00, 3.0e+00);
   RooRealVar b2_D0Bkg("b2_D0Bkg","b2_D0Bkg", 1.0385e+01, 0.0e+01, 3.0e+01);     
   x02_D0Bkg.setConstant();
   a2_D0Bkg.setConstant();
   b2_D0Bkg.setConstant();
   RooFormulaVar x02prime_D0Bkg("#x02prime_S", "x[0] < x[1] ? x[0] : x[1]", RooArgList(DeltaM, x02_D0Bkg));     
   RooFormulaVar gamma2_D0Bkg("#gamma2_D0Bkg", "a2_D0Bkg + 1", RooArgList(a2_D0Bkg));    
   RooFormulaVar beta2_for_gamma_D0Bkg("beta2_for_gamma_D0Bkg", "1./b2_D0Bkg", RooArgList(b2_D0Bkg));    
   RooGamma GammaModel2_D0Bkg("GammaModel3_D0Bkg", "Gamma pdf", DeltaM, gamma2_D0Bkg, beta2_for_gamma_D0Bkg, x02prime_D0Bkg);

   // Add the components
   RooRealVar f0_D0Bkg("f0_D0Bkg", "f0", 7.1325e-01, 5.0e-01, 9.0e-01);
   RooRealVar f1_D0Bkg("f1_D0Bkg", "f1", 2.0819e-04, 0.0e-04, 4.0e-04);
   RooRealVar f2_D0Bkg("f2_D0Bkg", "f2", 9.6634e-01, 7.0e-01, 11.0e-01);
   f0_D0Bkg.setConstant();
   f1_D0Bkg.setConstant();
   f2_D0Bkg.setConstant();
   RooAddPdf D0BkgModel("D0BkgModel","g1+g2",RooArgList(GammaModel0_D0Bkg,GammaModel1_D0Bkg,GammaModel2_D0Bkg), RooArgList(f0_D0Bkg,f1_D0Bkg,f2_D0Bkg));

   // Final Background Count
   //-----------------------------------------
   RooRealVar Nsig_D0Bkg("Nsig_D0Bkg", "#Signal events", 7.1232e+04, 5.0e+04, 9.0e+04);
   RooAddPdf model("model", "g+a", RooArgList(D0BkgModel), RooArgList(Nsig_D0Bkg));
   //=========================================================================================================================================================================//
   // C h i - s q u a r e d   F i t
   //--------------------------------------------------------------------------
   // Perform chi-squared fit
   RooChi2Var chi2("chi2", "chi2", model, data_D0, DataError(RooAbsData::SumW2));
   //=========================================================================================================================================================================//
   // F i t
   //-------------
   // Perform extended ML fit of data:
   // std::unique_ptr<RooFitResult> r{model.fitTo(*data_D0, Save(), PrintLevel(-1))};
   std::unique_ptr<RooFitResult> r{model.chi2FitTo(data_D0, Save(), PrintLevel(-1))};
   //=========================================================================================================================================================================//
   // P l o t   M o d e l   a n d   P r e p a r e   P a r a m e t e r   V a l u e s
   // ---------------------------------------------------------------------------
   // Create a RooPlot to draw on. We don't manage the memory of the returned
   // pointer. Instead we let it leak such that the plot still exists at the end of
   // the macro and we can take a look at it.
   RooPlot *frame = DeltaM.frame();
   data_D0.plotOn(frame, Binning(bin));
   model.plotOn(frame, LineColor(kViolet+3), RooFit::Name("Signal"));
   //=========================================================================================================================================================================//
   // C a l c u l a t e   c h i ^ 2
   // ------------------------------
   // Show the chi^2 of the curve w.r.t. the histogram
   // If multiple curves or datasets live in the frame you can specify
   // the name of the relevant curve and/or dataset in chiSquare()
   Int_t npar = model.getParameters(data_D0)->selectByAttrib("Constant",kFALSE)->getSize(); //select floating parameters and count their number

   cout << "NDF = " << npar << endl;
   cout << "chi^2/NDF = " << frame->chiSquare(npar) << endl;

   Double_t chi2ndf = frame->chiSquare(npar);
   //=========================================================================================================================================================================//
   // P l o t    L a b e l s
   // ---------------------------------------------------------------------------
   // Add text to frame
   // TLatex *latex = new TLatex();
   // TString lum ("");
   // lum.Form("#splitline{Signal Events}{2M Events}");
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
   TString chi_label ("");
   chi_label.Form("#splitline{}{#chi^{2} Fit}");

   auto PerBin = ((xR-xL)/(bin))*1000;
   // TString title ("");
   TString y_axis ("");
   // TString x_axis ("");

   // title.Form("Generic d#bar{d} Events");
   y_axis.Form("Entries/(%.2f MeV/c^{2})", PerBin);
   // x_axis.Form("#Delta m  [GeV/c^{2}]");

   frame->SetTitle("");
   frame->GetYaxis()->SetTitle(y_axis);
   frame->GetYaxis()->SetRangeUser(0,4000);
   frame->GetXaxis()->SetTitle("");
   // frame->GetXaxis()->CenterTitle(true);
   frame->Draw();

   // The line below prints the luminosity onto the screen.
   // latex->DrawLatex(0.6,0.2, lum);
   latex->DrawLatex(0.18,0.89, lum);
   latex->DrawLatex(0.18,0.85, chi_label);
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

   TString chi2_ndf("");
   TString Ntot_para("");
   TString N_para("");
   TString b_para("");
   TString P_para("");

   chi2_ndf.Form("#chi^{2}/ndf = %.2f", chi2ndf);
   Ntot_para.Form("N_{D^{0}} = %.2f #pm %.2f", Nsig_D0Bkg.getValV(), Nsig_D0Bkg.getError());
   // N_para.Form("N^{*} = %.2f #pm %.2f", NWindow, NWindow_error);
   // b_para.Form("b = %.6f #pm %.6f", b_value, b.getError());
   // P_para.Form("P = %.6f #pm %.6f", P_value, P.getError());

   legend->AddEntry(frame->FindObject("*data_D0"),"data","pe");
   legend->AddEntry((TObject *)0, chi2_ndf, "");
   legend->AddEntry("Signal", "D^{0} #rightarrow K^{-} #pi^{+}", "pl" );
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
   // When I use rooworkspace this first part will make sure that the parameters of the 
   // model will stay fixed.
   DeltaM.setConstant(true);
   // B0_S.removeError();

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

   // Verbose printing: Basic info, values of constant parameters, initial and
   // final values of floating parameters, global correlations
   r->Print("v");
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
   RooWorkspace *w = new RooWorkspace("w_D0", "workspace");
   // Import model and all its components into the workspace
   w->import(D0BkgModel);

   // E n c o d e   d e f i n i t i o n   o f   p a r a m e t e r s   i n   w o r k s p a c e
   // ---------------------------------------------------------------------------------------
   // Define named sets "parameters" and "observables", which list which variables should be considered
   // parameters and observables by the users convention
   //
   // Variables appearing in sets _must_ live in the workspace already, or the autoImport flag
   // of defineSet must be set to import them on the fly. Named sets contain only references
   // to the original variables, therefore the value of observables in named sets already
   // reflect their 'current' value
   RooArgSet *params = D0BkgModel.getVariables();
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

   // Retrieve the final fitted parameters
   RooArgList fitParams = r->floatParsFinal();
   TIterator* iter = fitParams.createIterator();
   RooRealVar* param;

   // Loop over parameters and set them as constants
   while ((param = (RooRealVar*)iter->Next())) {
      param->setVal(param->getVal());  // Update value to the fitted result
      param->setConstant(kTRUE);       // Set as constant
   }

   delete iter; // Clean up memory

   w->import(fitParams);  // Import the fitted and now constant parameters
   w->saveSnapshot("reference_fit_D0", fitParams, true);

   // Print workspace contents
   w->Print();

   // S a v e   w o r k s p a c e   i n   f i l e
   // -------------------------------------------
   // Save the workspace into a ROOT file
   w->writeToFile("/home/belle2/amubarak/Ds2D0enue_Analysis/10-Suggestion/Roofit_Suggestion_Chi/Fit_Parameter/rs03_D0.root");
   //=========================================================================================================================================================================//
   c->SaveAs("/home/belle2/amubarak/Ds2D0enue_Analysis/10-Suggestion/Roofit_Suggestion_Chi/Fit_Image/rs03_DeltaM_D0.png");
   //=========================================================================================================================================================================//
   // Final Fit Status
   //---------------------------
   int fitStatus = r->status();  // Get the fit status
   int covQual = r->covQual();   // Get covariance quality

   std::cout << "Fit Status: " << fitStatus << std::endl;
   std::cout << "Covariance Matrix Quality: " << covQual << std::endl << endl;
   //=========================================================================================================================================================================//
   double chi2_value = chi2.getVal();
   int ndof = model.getParameters(data_D0)->selectByAttrib("Constant", kFALSE)->getSize(); 
   cout << endl << "Chi2 = " << chi2_value << ", NDF = " << ndof << ", Chi2/NDF = " << chi2_value / ndof << endl; 
}