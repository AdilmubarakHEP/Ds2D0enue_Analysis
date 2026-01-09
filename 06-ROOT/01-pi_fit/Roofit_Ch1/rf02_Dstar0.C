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

void rf02_Dstar0()
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
   double xL = 0.15;
   double xR = 0.38;
   //=========================================================================================================================================================================//
   // I m p o r t s   D a t a   f r o m   t r e e   b r a n c h:
   // --------------------------------------------------------------
   TFile *file = new TFile("/group/belle/users/amubarak/03-ML/BkgBDT/Ds2D0e-Generic_Ds_053025_0_All_withBkgBDT.root");
   TTree *tree = (TTree*)file->Get("Dstree");
   //=========================================================================================================================================================================//
   // D a t a:
   // --------------------------------------------------------------
   RooRealVar x("Ds_diff_D0pi", "Mass Difference", xL, xR);
   RooRealVar PDG("Ds_mcPDG","PDG Value", -9999, 9999);
   RooRealVar BS("Ds_BkgBDT", "Background Suppression", 0, 1);
   RooDataSet data("data", "dataset with mass", tree, RooArgSet(x,PDG,BS));

   // Apply multiple selection cuts: 
   RooDataSet *filteredData = (RooDataSet*)data.reduce("Ds_BkgBDT >= 0.7 && abs(Ds_mcPDG)==423");
   //=========================================================================================================================================================================//
   // B u i l d   M o d e l
   // --------------------------------------------------------------

   // Method 1:
   //============================
   // Base gamma near threshold
   RooRealVar Dstar0_x00 ("Dstar0_x00","threshold", 1.3956e-01, 0.135, 0.145);
   RooRealVar Dstar0_a0  ("Dstar0_a0", "a0",        3.6092e-01, 2.6092e-01, 4.6092e-01);
   RooRealVar Dstar0_b0  ("Dstar0_b0", "b0",        4.1477e+01, 3.1477e+01, 5.1477e+01);
   Dstar0_x00.setConstant();
   Dstar0_a0.setConstant();
   Dstar0_b0.setConstant();
   RooFormulaVar Dstar0_x00p ("Dstar0_x00p", "@0<@1?@0:@1", RooArgList(x, Dstar0_x00));
   RooFormulaVar Dstar0_gam0 ("Dstar0_gam0","@0+1.0",      RooArgList(Dstar0_a0));
   RooFormulaVar Dstar0_bet0 ("Dstar0_bet0","1.0/@0",      RooArgList(Dstar0_b0));
   RooGamma     Dstar0_G0    ("Dstar0_G0", "Gamma0", x, Dstar0_gam0, Dstar0_bet0, Dstar0_x00p);

   // Sharpening gamma defined *relative* to the base
   RooRealVar   Dstar0_dA ("Dstar0_dA","shape increment", 2.2718e+00, 1.2718e+00, 3.2718e+00); // a1 = a0 + dA
   RooRealVar   Dstar0_sB ("Dstar0_sB","scale ratio",     1.2558e+00, 0.2558e+00, 2.2558e+00); // b1 = sB * b0
   Dstar0_dA.setConstant();
   Dstar0_sB.setConstant();
   RooFormulaVar Dstar0_a1F ("Dstar0_a1F","@0+@1",  RooArgList(Dstar0_a0, Dstar0_dA));
   RooFormulaVar Dstar0_b1F ("Dstar0_b1F","@0*@1",  RooArgList(Dstar0_b0, Dstar0_sB));
   RooFormulaVar Dstar0_gam1("Dstar0_gam1","@0+1.0",RooArgList(Dstar0_a1F));
   RooFormulaVar Dstar0_bet1("Dstar0_bet1","1.0/@0",RooArgList(Dstar0_b1F));
   RooGamma      Dstar0_G1  ("Dstar0_G1","Gamma1", x, Dstar0_gam1, Dstar0_bet1, Dstar0_x00p);

   // Mixture
   RooRealVar Dstar0_f1("Dstar0_f1","frac(G1)", 5.7185e-01, 4.7185e-01, 6.7185e-01);
   Dstar0_f1.setConstant();
   RooAddPdf  Dstar0Model("Dstar0Model","G1+G0", RooArgList(Dstar0_G1, Dstar0_G0),
                        RooArgList(Dstar0_f1));

   // // Method 2:
   // //============================
   // RooRealVar Dstar0_mu  ("Dstar0_mu","mean",   0.1454,  0.1452, 0.1457);
   // RooRealVar Dstar0_sigL("Dstar0_sigL","sigL", 3.6e-4,  1.5e-4, 8.0e-4);
   // RooRealVar Dstar0_sigR("Dstar0_sigR","sigR", 3.4e-4,  1.5e-4, 8.0e-4);
   // RooRealVar Dstar0_aL  ("Dstar0_aL","alphaL", 1.1,     0.001,    5.0);
   // RooRealVar Dstar0_aR  ("Dstar0_aR","alphaR", 1.1,     0.001,    5.0);
   // RooRealVar Dstar0_nL  ("Dstar0_nL","nL",     12.0,    2.0,    60.0);
   // RooRealVar Dstar0_nR  ("Dstar0_nR","nR",     12.0,    2.0,    60.0);

   // RooCrystalBall Dstar0Model("Dstar0Model","DSCB",
   //                            x, Dstar0_mu, Dstar0_sigL, Dstar0_sigR,
   //                            Dstar0_aL, Dstar0_nL, Dstar0_aR, Dstar0_nR);

   // Background 
   //======================
   // Fermi function using RooGenericPdf
   RooRealVar x0("x0", "Turn-on", 0.135, 0.145);
   RooRealVar k("k", "Slope", 0.000001, 0.05);
   RooGenericPdf bkg("bkg", "1.0 / (1.0 + exp(-(Ds_diff_D0pi - x0)/k))", RooArgSet(x, x0, k));

   // Total model
   //----------------------------
   RooRealVar nsig("nsig", "# signal", 0, 1e6);
   RooRealVar nbkg("nbkg", "# background", 0, 1e6);
   RooAddPdf model("model", "sig + bkg", RooArgList(Dstar0Model, bkg), RooArgList(nsig, nbkg));
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
   std::unique_ptr<RooFitResult> r{model.fitTo(*filteredData, Save(), PrintLevel(-1))};
   // If the fit is too broad and does not seem to fit a peak use the code below. Use it
   // once to determine where the peak should be and then adjust the mean above.
   // model.fitTo(data, Minos(kTRUE), PrintLevel(-1), Range(2,2.1));
   //=========================================================================================================================================================================//
   // P l o t   M o d e l   a n d   P r e p a r e   P a r a m e t e r   V a l u e s
   // ---------------------------------------------------------------------------
   // Create a RooPlot to draw on. We don't manage the memory of the returned
   // pointer. Instead we let it leak such that the plot still exists at the end of
   // the macro and we can take a look at it.
   RooPlot *frame = x.frame();
   filteredData->plotOn(frame, Binning(50));

   // Peak
   model.plotOn(frame, Components(Dstar0Model),   LineStyle(kDashed), LineColor(28),    RooFit::Name("Peak"));
   model.plotOn(frame, Components(bkg),           LineStyle(kDashed), LineColor(kGreen+3),   RooFit::Name("Comb"));
   model.plotOn(frame, LineColor(kBlue+2), RooFit::Name("model"));
   //=========================================================================================================================================================================//
   // C a l c u l a t e   c h i ^ 2
   // ------------------------------
   // Show the chi^2 of the curve w.r.t. the histogram
   // If multiple curves or datasets live in the frame you can specify
   // the name of the relevant curve and/or dataset in chiSquare()
   Int_t npar = model.getParameters(*filteredData)->selectByAttrib("Constant",kFALSE)->getSize(); //select floating parameters and count their number

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
   // lum.Form("#splitline{Simulation}{200k Events}");
   lum.Form("#splitline{Generic Events}{#scale[0.5]{#int}#it{L} dt=1.4 ab^{-1}}");
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
   // frame->GetYaxis()->SetRangeUser(0,1500);
   double maxY = frame->GetMaximum();
   // frame->GetYaxis()->SetRangeUser(0, 5000);
   frame->GetYaxis()->SetRangeUser(0, 1.4 * maxY);
   // frame->GetYaxis()->SetTitleSize(0.07);

   frame->GetXaxis()->SetTitle("");
   // frame->GetXaxis()->CenterTitle(true);

   frame->Draw();

   // The line below prints the luminosity onto the screen.
   // latex->DrawLatex(0.6,0.2, lum);
   latex->DrawLatex(0.18,0.87, lum);
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
   auto x1 = xIF1 + 0.29*dx;
   auto x2 = xIF2;
   auto y1 = yIF1 + 0.77*dy;
   auto y2 = 0.93;
   auto legend = new TLegend(x1, y1, x2, y2);

   TString chi2("");
   TString Nsig_para("");
   TString Nbkg_para("");

   chi2.Form("#chi^{2}/ndf = %.2f", chi2ndf);
   Nsig_para.Form("N_{sig} = %.2f #pm %.2f", nsig.getValV(), nsig.getError());
   Nbkg_para.Form("N_{bkg} = %.2f #pm %.2f", nbkg.getValV(), nbkg.getError());

   legend->AddEntry(frame->FindObject("*filteredData"),"data","pe");
   legend->AddEntry((TObject *)0, chi2, "");
   legend->AddEntry("model", "Total Fit", "pl" );
   legend->AddEntry((TObject *)0, Nsig_para, "");
   legend->AddEntry("Peak", "D^{*0} #rightarrow D^{0} #pi^{0}/#gamma", "pl");
   legend->AddEntry((TObject *)0, Nbkg_para, "");
   legend->AddEntry("Comb", "Background", "pl" );
   legend->AddEntry((TObject *)0, " ", "");
   // legend->AddEntry("Comb", "Background", "pl" );

   // legend->AddEntry("Peak", "D^{*0} #rightarrow D^{0} #pi^{0}/#gamma", "pl" );

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
   RooPlot *frame1 = x.frame(Title("Pull Distribution"));
   frame1->addPlotable(hpull, "B");
   frame1->Draw();
   
   // Plot Titles
   //---------------
   frame1->SetTitle("");

   // Y-axis
   frame1->SetYTitle("Pull");
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
   RooWorkspace *w = new RooWorkspace("w_Dstar0", "workspace");
   // Import model and all its components into the workspace
   w->import(Dstar0Model);

   // E n c o d e   d e f i n i t i o n   o f   p a r a m e t e r s   i n   w o r k s p a c e
   // ---------------------------------------------------------------------------------------
   // Define named sets "parameters" and "observables", which list which variables should be considered
   // parameters and observables by the users convention
   //
   // Variables appearing in sets _must_ live in the workspace already, or the autoImport flag
   // of defineSet must be set to import them on the fly. Named sets contain only references
   // to the original variables, therefore the value of observables in named sets already
   // reflect their 'current' value
   RooArgSet *params = Dstar0Model.getVariables();
   w->defineSet("parameters", *params);
   w->defineSet("observables", x);

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

   RooArgList fitParams = r->floatParsFinal();

   // Loop over each parameter using ROOT's modern iteration style
   for (auto obj : fitParams) {
      RooRealVar* param = dynamic_cast<RooRealVar*>(obj);
      if (param) {
         param->setVal(param->getVal());  // Update value to fitted result
         param->setConstant(kTRUE);       // Set as constant
      }
   }

   w->import(fitParams);  // Import the fitted and now constant parameters
   w->saveSnapshot("reference_fit_Dstar0", fitParams, true);

   // Print workspace contents
   w->Print();

   // S a v e   w o r k s p a c e   i n   f i l e
   // -------------------------------------------
   // Save the workspace into a ROOT file
   w->writeToFile("/home/belle2/amubarak/Ds2D0enue_Analysis/05-ROOT/01-pi_fit/Roofit_Ch1/Fit_Parameter/rf02_Dstar0.root");
   //=========================================================================================================================================================================//
   c->SaveAs("/home/belle2/amubarak/Ds2D0enue_Analysis/05-ROOT/01-pi_fit/Roofit_Ch1/Fit_Image/rf02_Dstar0.png");
   //=========================================================================================================================================================================//
   // Final Fit Status
   //---------------------------
   int fitStatus = r->status();  // Get the fit status
   int covQual = r->covQual();   // Get covariance quality

   std::cout << "Fit Status: " << fitStatus << std::endl;
   std::cout << "Covariance Matrix Quality: " << covQual << std::endl << endl;
}