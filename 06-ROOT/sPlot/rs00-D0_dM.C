#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooAbsPdf.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
 
#include "TCanvas.h"
#include "TLegend.h"
 
#include <iomanip>

// use this order for safety on library loading
using namespace RooFit;
using namespace RooStats;
using namespace TMath;

void rs00_D0_dM()
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
   double xL = -0.05;
   double xR = 0.05;
   //=========================================================================================================================================================================//
   // I m p o r t s   D a t a   f r o m   t r e e   b r a n c h:
   // --------------------------------------------------------------
   TFile *file = new TFile("C03-Grid/Completed/Ds2D0e-Generic_Ds_111424_0_All.root");
   TTree *tree = (TTree*)file->Get("Dstree");
   //=========================================================================================================================================================================//
   // D a t a:
   // --------------------------------------------------------------
   RooRealVar dM("D0_dM", "Mass Difference", xL, xR);
   RooRealVar DeltaM("Ds_diff_D0pi", "Mass Difference", 0.1, 1);
   RooRealVar Veto("Ds_gammaveto_M_Correction","Veto",0.1,100);
   RooRealVar BCS("Ds_chiProb_rank","BCS",-1,1.25);
   // RooRealVar Category("Ds_D0_sideband","Category",0.5,1);

   RooDataSet data("data", "dataset with mass", tree, RooArgSet(dM,DeltaM,Veto,BCS));
   //=========================================================================================================================================================================//
	// C r e a t e   m o d e l
   // -------------------------------------------------------------- 
   // Signal Peak:
   //-----------------
   // // First Attempt:
   // // Gaussian Resolution model
	// RooRealVar sigmean1("sigmean1", "mean", -0.005, 0.005);
	// RooRealVar sigwidth1("sigwidth1", "width", 0, 0.01);
	// RooGaussModel gaussm1("gaussm1", "Signal PDF", dM, sigmean1, sigwidth1);
	// // Gaussian Resolution model
	// RooRealVar sigmean2("sigmean2", "mean", -0.005, 0.005);
	// RooRealVar sigwidth2("sigwidth2", "width", 0, 0.01);
	// RooGaussModel gaussm2("gaussm2", "Signal PDF", dM, sigmean2, sigwidth2);

	// // Add the components
	// RooRealVar g1frac("g1frac","fraction of gauss1", 0., 1.);
	// RooRealVar g2frac("g2frac","fraction of gauss2", 0., 1.);
	// RooAddPdf Signal("Signal","g1+g2",RooArgList(gaussm1,gaussm2), RooArgList(g1frac,g2frac));

   // Second Attempt:
   RooRealVar sig_mean("sig_mean", "mean of crystal", 1.2035e-04, -0.005, 0.005);
   RooRealVar sig_sigmaL("sig_sigmaL", "sig_sigmaL", 4.4147e-03, 1e-4, 0.02);
   RooRealVar sig_sigmaR("sig_sigmaR", "sig_sigmaR", 4.2152e-03, 1e-4, 0.02);
   RooRealVar sig_alphaL("sig_alphaL", "sig_alphaL", 1.4540e+00, 1e-4, 4.);
   RooRealVar sig_alphaR("sig_alphaR", "sig_alphaR", 1.3880e+00, 1e-4, 4.);
   RooRealVar sig_nL("sig_nL", "sig_nL", 1.0961e+00, 1e-4, 10);
   RooRealVar sig_nR("sig_nR", "sig_nR", 1.9424e+00, 1e-4, 10);
    sig_mean.setConstant();
    sig_sigmaL.setConstant();
    sig_sigmaR.setConstant();
    sig_alphaL.setConstant();
    sig_alphaR.setConstant();
   //  sig_nL.setConstant();
   //  sig_nR.setConstant();
   RooCrystalBall Signal("Signal", "Signal", dM, sig_mean,sig_sigmaL,sig_sigmaR,sig_alphaL,sig_nL,sig_alphaR,sig_nR);

   // Background PDF:
   //---------------------
   // Constant Background:
   // RooRealVar c0("c0", "coefficient #0", -0.99, 0.99);
   // RooRealVar c1("c1", "coefficient #1", -0.99, 0.99);
   // RooRealVar c2("c2", "coefficient #2", -0.99, 0.99);
   // RooRealVar c3("c3", "coefficient #3", -0.99, 0.99);
   // RooRealVar c4("c4", "coefficient #4", -4.3222e-02, -1, 1);
   // c0.setConstant(true);
   // c1.setConstant(true);
   // c2.setConstant(true);
   // c2.setConstant(true);
   // c3.setConstant(true);
   // c4.setConstant(true);
   // RooChebychev Background("Background", "Other", DeltaM, RooArgList(c0,c1,c2,c3));
   RooRealVar c0("c0", "coefficient #0", 1);
   RooPolynomial Background("Background", "background p.d.f.", dM, c0);

   // C o n s t r u c t   a   S i g na l   a n d   B a c k g r o u n d   P D F:
   // --------------------------------------------------------------
   RooRealVar nsig("nsig", "#Signal events", 3.2225e+05, 0., 600000);
   RooRealVar nbkg("nbkg", "#background events", 2.4334e+05, 0., 600000);
   RooAddPdf model("model", "g+a", RooArgList(Signal, Background), RooArgList(nsig, nbkg));
   //=========================================================================================================================================================================//
   // // S e t   t h e   C o n s t a n t   f o r   t h e    V a r i a b l e s
   // //--------------------------------------------------------------------------
   // // The code below fixes the value for a certain value.
   // sigmean1.setVal(-0.00119229);
   // sigmean2.setVal(0.00018854);
   // sigwidth1.setVal(0.00814808);
   // sigwidth2.setVal(0.00375633);
   // g1frac.setVal(0.487768);
   // g2frac.setVal(0.995866);
   // c5.setVal(1.58085e-05);

   // // // The code below removes the range and fixes the value.
   // sigmean1.setConstant(true);
   // // sigmean1.removeError();
   // sigmean2.setConstant(true);
   // // sigmean2.removeError();
   // sigwidth1.setConstant(true);
   // // sigwidth1.removeError();
   // sigwidth2.setConstant(true);
   // // sigwidth2.removeError();
   // g1frac.setConstant(true);
   // // g1frac.removeError();
   // g2frac.setConstant(true);
   // // g2frac.removeError();
   // c5.setConstant(true);
   // // c5.removeError();
   //=========================================================================================================================================================================//
   // Perform extended ML fit of data: The first line shows more information
   std::unique_ptr<RooFitResult> r{model.fitTo(data, Save(), PrintLevel(-1))};
   // If the fit is too broad and does not seem to fit a peak use the code below. Use it 
   // once to determine where the peak should be and then adjust the mean above.
   // model.fitTo(data, Minos(kTRUE), PrintLevel(-1), Range(2,2.1));
   //=========================================================================================================================================================================//
   // P l o t   M o d e l   a n d   P r e p a r e   P a r a m e t e r   V a l u e s
   // -------------------------------------------------------------------------------
   // Create a RooPlot to draw on. We don't manage the memory of the returned
   // pointer. Instead we let it leak such that the plot still exists at the end of
   // the macro and we can take a look at it.
   RooPlot* frame = dM.frame();
   data.plotOn(frame, Binning(bin));
   model.plotOn(frame, Components(Background), LineStyle(kDashed), LineColor(46), RooFit::Name("Background"));
   model.plotOn(frame, Components(Signal),     LineStyle(kDashed), LineColor(28), RooFit::Name("Signal"));
   model.plotOn(frame, LineColor(28), RooFit::Name("model"));
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
   // // Add text to frame
   // TLatex *latex = new TLatex();
   // TString lum ("");
   // lum.Form("#splitline{Simulated Events}{2M Events}");
   // latex->SetNDC();
   // latex->SetTextSize(0.035);
   // Add text to frame
   TLatex *latex = new TLatex();
   TString lum ("");
   lum.Form("#splitline{Generic Events}{#scale[0.5]{#int}#it{L} dt=200 fb^{-1}}");
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
   // frame->GetYaxis()->SetRangeUser(0,90000);
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
   Nsig_para.Form("N_{sig} = %.2f #pm %.2f", nsig.getValV(), nsig.getError());
   Nbkg_para.Form("N_{bkg} = %.2f #pm %.2f", nbkg.getValV(), nbkg.getError());

   legend->AddEntry(frame->FindObject("data"),"data","pe");
   legend->AddEntry((TObject *)0, chi2, "");
   legend->AddEntry("Signal", "D^{0} #rightarrow K^{-} #pi^{+}", "pl" );
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
   RooPlot *frame1 = dM.frame(Title("Residual Distribution"));
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
   x_axis.Form("dM(K^{-}#pi^{+})  [GeV/c^{2}]");
   frame1->SetXTitle(x_axis);
   // frame1->GetXaxis()->CenterTitle(true);
   frame1->GetXaxis()->SetTitleSize(0.13);
   frame1->GetXaxis()->SetLabelSize(0.10);

   c->cd();
   //=========================================================================================================================================================================//
   // S a v e   F i l e s
   // -------------------------------------------------------------
   c->SaveAs("Ds2D0enue_Analysis/05-ROOT/sPlot_results/rs00_D0_dM_SignalRegion.png");
   //=========================================================================================================================================================================//
   // SPlot
   // -------------------------------------------------------------
   //......after a normal unbinned ML fit on `RooDataSet *Data` with `RooAddPdf sumpdf` (sigPDF, bkgPDF) with yields Nsig and Nbkg ......
   //using sPlot technique
   
   //---------------------------------------------------------------------------------------------------------------------
   // Cross-check Section:
   //----------------------
   std::cout << endl << "Calculate sWeights" << std::endl;

   std::cout << "\n\n------------------------------------------\nThe dataset before creating sWeights:\n";
   data.Print();

   // Now we use the SPlot class to add SWeights for the isolation variable to
   // our data set based on fitting the yields to the invariant mass variable
   RooStats::SPlot splot("splot", "splot", data, &model, RooArgList(nsig,nbkg));
 
   std::cout << "\n\nThe dataset after creating sWeights:\n";
   data.Print();

   // Check that our weights have the desired properties
 
   std::cout << "\n\n------------------------------------------\n\nCheck SWeights:" << std::endl;
 
   std::cout << std::endl
             << "Yield of Signal is\t" << nsig.getVal() << ".  From sWeights it is "
             << splot.GetYieldFromSWeight("nsig") << std::endl;
 
   std::cout << "Yield of Background is\t" << nbkg.getVal() << ".  From sWeights it is "
             << splot.GetYieldFromSWeight("nbkg") << std::endl
             << std::endl;

   for (Int_t i = 0; i < 10; i++) {
      std::cout << "Signal Weight for event " << i << std::right << std::setw(12) << splot.GetSWeight(i, "nsig") << "  Background Weight"
                << std::setw(12) << splot.GetSWeight(i, "nbkg") << "  Total Weight" << std::setw(12) << splot.GetSumOfEventSWeight(i)
                << std::endl;
   }
 
   std::cout << std::endl;
   //---------------------------------------------------------------------------------------------------------------------

   RooDataSet data_signal(data.GetName(), data.GetTitle(), &data, *data.get(), nullptr, "nsig_sw");
   RooDataSet data_background(data.GetName(), data.GetTitle(), &data, *data.get(), nullptr, "nbkg_sw");

   // Here we make plots of the discriminating variable (invMass) after the fit
   // and of the control variable (isolation) after unfolding with sPlot.
 
   // make our canvas
   TCanvas *cdata = new TCanvas("sPlot", "sPlot demo", 600, 450);

   TH1F* hsig_dM = new TH1F("hsig_dM", "", 30, 0, 4.5);
   TH1F* hbkg_dM = new TH1F("hbkg_dM", "", 30, 0, 4.5);
   data_signal.fillHistogram( hsig_dM, RooArgList(dM));
   data_background.fillHistogram( hbkg_dM, RooArgList(dM));

   cout << "Making splots of the signal and background" <<endl;
   RooPlot* f_splot = dM.frame(30) ;
   f_splot->SetTitle(" ");
   data_signal.plotOn(f_splot, RooFit::DataError(RooAbsData::SumW2), LineColor(kRed), MarkerColor(kRed),LineWidth(2),MarkerStyle(20), RooFit::Name("sig")) ;
   data_background.plotOn(f_splot, RooFit::DataError(RooAbsData::SumW2), LineColor(kBlue), MarkerColor(kBlue),LineWidth(2),MarkerStyle(22), RooFit::Name("bkg")) ;
   cdata->cd();
   f_splot->Draw("hist");
   TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
   leg->AddEntry(f_splot->findObject("sig"),"signal","EP");
   leg->AddEntry(f_splot->findObject("bkg"),"background","EP");
   leg->Draw();  
   cdata->SaveAs("Ds2D0enue_Analysis/05-ROOT/sPlot_results/sPlot.png");

   // TFile *sfile = new TFile("Ds2D0enue_Analysis/05-ROOT/SPlot_results/sPlot_SigBkg.root", "recreate");
   // cdata->Write();  
   // hsig_dM->Write(); 
   // hbkg_dM->Write();
   //=========================================================================================================================================================================//
   // C r e a t e   w o r k s p a c e ,   i m p o r t   d a t a   a n d   m o d e l
   // -----------------------------------------------------------------------------
   // Create a new empty workspace
   RooWorkspace *w = new RooWorkspace("w_splot", "workspace");
   
   // import this new dataset with sWeights
   std::cout << "import new dataset with sWeights" << std::endl;
   w->import(data, Rename("dataWithSWeights"));
   
   w->Print();

   // S a v e   w o r k s p a c e   i n   f i l e
   // -------------------------------------------
   // Save the workspace into a ROOT file
   w->writeToFile("Ds2D0enue_Analysis/05-ROOT/sPlot_results/rs00-RealD0.root");
   //=========================================================================================================================================================================//
   // ROOT File
   //---------------
   // TFile f("rf402_datahandling.root", "RECREATE");
   // data.Write(0,TObject::kOverwrite);

   TFile outputFile("filename.root", "RECREATE");
   RooAbsData::setDefaultStorageType(RooAbsData::Tree);
   data.convertToTreeStore();
   data.SetName("Dstree");
   data.Write();
   // data.Write(0,TObject::kOverwrite);
   //=========================================================================================================================================================================//
}