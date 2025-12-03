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

void rf02_Combinatorial()
{
   // This code is used for the case where the electron was give a pion mass. Unless more work is done
   // it would be difficult to extract a yield from this code.
   
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
   double xR = 1;
   //=========================================================================================================================================================================//
   // I m p o r t s   D a t a   f r o m   t r e e   b r a n c h:
   // --------------------------------------------------------------
   TFile *file = new TFile("C03-Grid/Completed/Ds2D0e-Generic_Ds_110324_0_All.root");
   TTree *tree = (TTree*)file->Get("Dstree");
   //=========================================================================================================================================================================//
   // D a t a:
   // --------------------------------------------------------------
   RooRealVar DeltaM("Ds_diff_D0pi", "Mass Difference", xL, xR);
   RooRealVar Veto("Ds_gammaveto_M_Correction","Veto",0.1,100);
   RooRealVar BCS("Ds_chiProb_rank","BCS",-1,1.25);
   RooRealVar dM("D0_dM", "Mass Difference", -0.02, 0.02);

   RooRealVar CombFailed("Ds_CombFailed","Ds_CombFailed",0.5,1.1);

   // RooDataSet data("data", "dataset with mass", tree, RooArgSet(M,x,y));
   RooDataSet data_CF("data_CF", "dataset with mass", tree, RooArgSet(DeltaM,Veto,BCS,dM,CombFailed));
   //=========================================================================================================================================================================//
   // C r e a t e   C o m b i n e d   M o d e l
   // ------------------------------------------------------------------------------
   // // First Attempt
   // //---------------------
   // // Generic PDF
   // RooRealVar B0_CBkg("B0_CBkg", "B0_CBkg", -8.3111e+00, -40., 20.);
   // RooRealVar P0_CBkg("P0_CBkg", "P0_CBkg", 8.1845e-01, 0., 5.);
   // B0_CBkg.setConstant(true);
   // P0_CBkg.setConstant(true);
   // RooGenericPdf Generic0_CBkg("Generic0_CBkg","Generic PDF","(Ds_diff_D0pi - 0.13953 > 0 ? ( ((Ds_diff_D0pi-0.13953)^P0_CBkg)*exp( B0_CBkg*(Ds_diff_D0pi-0.13953))) : 0.0)", RooArgSet( DeltaM, B0_CBkg, P0_CBkg));

   // // Generic PDF
   // RooRealVar B1_CBkg("B1_CBkg", "B1_CBkg", -3.2604e+00, -40., 20.);
   // RooRealVar P1_CBkg("P1_CBkg", "P1_CBkg", 1.6060e+00, 0., 5.);
   // B1_CBkg.setConstant(true);
   // P1_CBkg.setConstant(true);
   // RooGenericPdf Generic1_CBkg("Generic1_CBkg","Generic PDF","(Ds_diff_D0pi - 0.13953 > 0 ? ( ((Ds_diff_D0pi-0.13953)^P1_CBkg)*exp( B1_CBkg*(Ds_diff_D0pi-0.13953))) : 0.0)", RooArgSet( DeltaM, B1_CBkg, P1_CBkg));

   // // Combined Background
   // RooRealVar f0_CBkg("f0_CBkg", "f0", 5.3990e-01, 0, 1);
   // RooRealVar f1_CBkg("f1_CBkg", "f1", 4.2536e-01, 0, 1);
   // f0_CBkg.setConstant(true);
   // f1_CBkg.setConstant(true);
   // RooAddPdf CFBkgModel("CFBkgModel", "b+c", RooArgList(Generic0_CBkg,Generic1_CBkg), RooArgList(f0_CBkg,f1_CBkg));

   // Second Attempt
   //---------------------
   // First Gamma Function
   RooRealVar x00_CBkg("x00_CBkg","x00_CBkg", 0.13957, xL, 0.2);
   RooRealVar a0_CBkg("a0_CBkg","a0_CBkg", 2.1254e-01, 0.0, 1.);
   RooRealVar b0_CBkg("b0_CBkg","b0_CBkg", 2.6059e+00, 1.0, 5.);    
   x00_CBkg.setConstant();
   a0_CBkg.setConstant();
   b0_CBkg.setConstant();
   RooFormulaVar x00prime_CBkg("#x00prime_CBkg", "x[0] < x[1] ? x[0] : x[1]", RooArgList(DeltaM, x00_CBkg));     
   RooFormulaVar gamma0_CBkg("#gamma0_CBkg", "a0_CBkg + 1", RooArgList(a0_CBkg));    
   RooFormulaVar beta0_for_gamma_CBkg("beta0_for_gamma_CBkg", "1./b0_CBkg", RooArgList(b0_CBkg));    
   RooGamma GammaModel1_CBkg("GammaModel1_CBkg", "Gamma pdf", DeltaM, gamma0_CBkg, beta0_for_gamma_CBkg, x00prime_CBkg);

   // Second Gamma Function
   RooRealVar x01_CBkg("x01_CBkg","x01_CBkg", 0.13957, xL, 0.2);
   RooRealVar a1_CBkg("a1_CBkg","a1_CBkg", 4.6803e+00, 1.0, 5.);
   RooRealVar b1_CBkg("b1_CBkg","b1_CBkg", 6.9240e+00, 2.0, 10.);     
   x01_CBkg.setConstant();
   a1_CBkg.setConstant();
   b1_CBkg.setConstant();
   RooFormulaVar x01prime_CBkg("#x01prime_S", "x[0] < x[1] ? x[0] : x[1]", RooArgList(DeltaM, x01_CBkg));     
   RooFormulaVar gamma1_CBkg("#gamma1_CBkg", "a1_CBkg + 1", RooArgList(a1_CBkg));    
   RooFormulaVar beta1_for_gamma_CBkg("beta1_for_gamma_CBkg", "1./b1_CBkg", RooArgList(b1_CBkg));    
   RooGamma GammaModel2_CBkg("GammaModel2_CBkg", "Gamma pdf", DeltaM, gamma1_CBkg, beta1_for_gamma_CBkg, x01prime_CBkg);

   // Third Gamma Function
   RooRealVar x02_CBkg("x02_CBkg","x02_CBkg", 0.13957, xL, 0.2);
   RooRealVar a2_CBkg("a2_CBkg","a2_CBkg", 2.0633e+00, 1.0, 5.);
   RooRealVar b2_CBkg("b2_CBkg","b2_CBkg", 1.2933e+01, 10., 15.);     
   x02_CBkg.setConstant();
   a2_CBkg.setConstant();
   b2_CBkg.setConstant();
   RooFormulaVar x02prime_CBkg("#x02prime_S", "x[0] < x[1] ? x[0] : x[1]", RooArgList(DeltaM, x02_CBkg));     
   RooFormulaVar gamma2_CBkg("#gamma2_CBkg", "a2_CBkg + 1", RooArgList(a2_CBkg));    
   RooFormulaVar beta2_for_gamma_CBkg("beta2_for_gamma_CBkg", "1./b2_CBkg", RooArgList(b2_CBkg));    
   RooGamma GammaModel3_CBkg("GammaModel3_CBkg", "Gamma pdf", DeltaM, gamma2_CBkg, beta2_for_gamma_CBkg, x02prime_CBkg);

   // Add the components
   RooRealVar f0_CBkg("f0_CBkg", "f0", 9.9573e-01, 0, 1);
   RooRealVar f1_CBkg("f1_CBkg", "f1", 2.7963e-01, 0, 1);
   RooRealVar f2_CBkg("f2_CBkg", "f2", 4.9076e-01, 0, 1);
   f0_CBkg.setConstant();
   f1_CBkg.setConstant();
   f2_CBkg.setConstant();
   RooAddPdf CFBkgModel("CFBkgModel","g1+g2",RooArgList(GammaModel1_CBkg,GammaModel2_CBkg,GammaModel3_CBkg), RooArgList(f0_CBkg,f1_CBkg,f2_CBkg));

   // Final Background Count
   //-----------------------------------------
   RooRealVar Nsig_CBkg("Nsig_CBkg", "#Signal events", 3.2983e+05, 100000., 500000);
   RooAddPdf model("model", "g+a", RooArgList(CFBkgModel), RooArgList(Nsig_CBkg));
   //=========================================================================================================================================================================//
   // S e t   t h e   C o n s t a n t   f o r   t h e    V a r i a b l e s
   //--------------------------------------------------------------------------
   // // The code below fixes the value for a certain value.
   // B0_CBkg.setVal(-8.3580e+00);
   // B1_CBkg.setVal(-3.2331e-01);
   // P0_CBkg.setVal(1.1092e+00);
   // P1_CBkg.setVal(7.1018e-02);
   // f0_CBkg.setVal(2.0123e-04);
   // f1_CBkg.setVal(2.1654e-04);

   // // The code below removes the range and fixes the value.
   // B0_CBkg.setConstant(true);
   // B1_CBkg.setConstant(true);
   // P0_CBkg.setConstant(true);
   // P1_CBkg.setConstant(true);
   // f0_CBkg.setConstant(true);
   // f1_CBkg.setConstant(true);

   DeltaM.setVal(0.0065);
   DeltaM.setConstant(true);
   //=========================================================================================================================================================================//
   // C o n s t r u c t   a   S i g na l   a n d   B a c k g r o u n d   P D F:
   // --------------------------------------------------------------
   // RooRealVar N("N", "PBackground events", 0., 500000);
   // RooAddPdf Signal("Signal", "g+a", RooArgList(SignalModel), RooArgList(N));
   //=========================================================================================================================================================================//
   // F i t
   //-------------
   // Perform extended ML fit of data:
   std::unique_ptr<RooFitResult> r{model.fitTo(data_CF, Save(), PrintLevel(-1))};
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
   data_CF.plotOn(frame, Binning(bin));
   // model.plotOn(frame, Components(Generic0_CBkg),       LineStyle(kDashed), LineColor(46), RooFit::Name("Background"));
   // model.plotOn(frame, Components(Generic1_CBkg), LineStyle(kDashed), LineColor(38), RooFit::Name("Signal"));
   model.plotOn(frame, LineColor(kRed), RooFit::Name("model"));
   //=========================================================================================================================================================================//
   // C a l c u l a t e   c h i ^ 2
   // ------------------------------
   // Show the chi^2 of the curve w.r.t. the histogram
   // If multiple curves or datasets live in the frame you can specify
   // the name of the relevant curve and/or dataset in chiSquare()
   Int_t npar = model.getParameters(data_CF)->selectByAttrib("Constant",kFALSE)->getSize(); //select floating parameters and count their number

   cout << "NDF = " << npar << endl;
   cout << "chi^2/NDF = " << frame->chiSquare(npar) << endl;

   Double_t chi2ndf = frame->chiSquare(npar);
   //=========================================================================================================================================================================//
   // P l o t    L a b e l s
   // ---------------------------------------------------------------------------
      // Add text to frame
   TLatex *latex = new TLatex();
   TString lum ("");
   lum.Form("#splitline{Generic Events}{#scale[0.5]{#int}#it{L} dt=200 fb^{-1}}");
   // lum.Form("#splitline{Generic c#bar{c} Events}{#scale[0.5]{#int}#it{L} dt=200 fb^{-1}}");
   latex->SetNDC();
   latex->SetTextSize(0.03); 
   // // Add text to frame
   // TLatex *latex = new TLatex();
   // TString lum ("");
   // lum.Form("#splitline{Simulated Events}{100k Events}");
   // latex->SetNDC();
   // latex->SetTextSize(0.035); 
   
   auto PerBin = ((xR-xL)/(bin))*1000;
   // TString title ("");
   TString y_axis ("");
   // TString x_axis ("");

   // title.Form("Generic d#bar{d} Events");
   y_axis.Form("Entries/(%.2f MeV/c^{2})", PerBin);
   // x_axis.Form("#Delta m  [GeV/c^{2}]");

   frame->SetTitle("");
   frame->GetYaxis()->SetTitle(y_axis);
   frame->GetYaxis()->SetRangeUser(0,15000);
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
   auto x1 = xIF1 + 0.4*dx;
   auto x2 = xIF2;
   auto y1 = yIF1 + 0.75*dy;
   auto y2 = 0.93;
   auto legend = new TLegend(x1, y1, x2, y2);

   TString chi2("");
   TString Nsig_para("");
   TString Nbkg_para("");

   chi2.Form("#chi^{2}/ndf = %.2f", chi2ndf);
   Nsig_para.Form("N_{sig} = %.2f #pm %.2f", Nsig_CBkg.getValV(), Nsig_CBkg.getError());

   legend->AddEntry(frame->FindObject("data"),"data","pe");
   legend->AddEntry((TObject *)0, chi2, "");
   legend->AddEntry("model", "combinatorial+failed", "pl" );
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
   // B0_CBkg.setConstant(true);
   // B0_CBkg.removeError();
   // B1_CBkg.setConstant(true);
   // B1_CBkg.removeError();
   // P0_CBkg.setConstant(true);
   // P0_CBkg.removeError();
   // P1_CBkg.setConstant(true);
   // P1_CBkg.removeError();
   // f0_CBkg.setConstant(true);
   // f0_CBkg.removeError();
   // f1_CBkg.setConstant(true);
   // f1_CBkg.removeError();

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
   RooWorkspace *w = new RooWorkspace("w_CBkg", "workspace");
   // Import model and all its components into the workspace
   w->import(CFBkgModel);

   // E n c o d e   d e f i n i t i o n   o f   p a r a m e t e r s   i n   w o r k s p a c e
   // ---------------------------------------------------------------------------------------
   // Define named sets "parameters" and "observables", which list which variables should be considered
   // parameters and observables by the users convention
   //
   // Variables appearing in sets _must_ live in the workspace already, or the autoImport flag
   // of defineSet must be set to import them on the fly. Named sets contain only references
   // to the original variables, therefore the value of observables in named sets already
   // reflect their 'current' value
   RooArgSet *params = CFBkgModel.getVariables();
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
   w->saveSnapshot("reference_fit_Combinatorial_Bkg", *params, true);

   // Print workspace contents
   w->Print();

   // S a v e   w o r k s p a c e   i n   f i l e
   // -------------------------------------------
   // Save the workspace into a ROOT file
   w->writeToFile("Ds2D0enue_Analysis/05-ROOT/Fit_Parameter/rf02_CombinatorialBkgPara.root");
   //=========================================================================================================================================================================//
   c->SaveAs("Ds2D0enue_Analysis/09-Images/rf02_DeltaM_Combinatorial.png");

   // cout << "Total Number After Peak = " << NWindow << endl;
}