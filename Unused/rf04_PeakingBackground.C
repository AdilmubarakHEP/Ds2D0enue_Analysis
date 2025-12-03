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

void rf04_PeakingBackground()
{
   // S e t u p   t h e   F r a m e   f o r   t h e   P l o t
   //------------------------------------------------------------
   TCanvas *c = new TCanvas("fit", "fit", 650, 700);
   TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
   TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.25);
   pad1->SetBottomMargin(0.08);
   pad1->SetLeftMargin(0.15);
   pad2->SetLeftMargin(0.15);
   pad2->SetTopMargin(0.01);
   pad2->SetBottomMargin(0.25);
   // pad2->SetBorderMode(0);
   pad1->Draw();
   pad2->Draw();
   pad1->cd();
   gStyle->SetOptStat(0);
   //=========================================================================================================================================================================//
   // I n i t i a l i z a t i o n
   // --------------------------------------------------------------
   // The lines below initializes all my values
   int bin = 50;
   //=========================================================================================================================================================================//
   // I m p o r t s   D a t a   f r o m   t r e e   b r a n c h:
   // --------------------------------------------------------------
   TFile *file = new TFile("C:/Users/Adilm/Desktop/High_Energy_Physics/Analysis_Code/KEKCC_File/Grid/Ds2D0e-Generic_Ds_07252024_ccbar.root");
   TTree *tree = (TTree*)file->Get("dstree");
   //=========================================================================================================================================================================//
   // RooRealVar M("dM_D0pi", "m (GeV)", 0, 0.25);
   // RooRealVar x("e_electronID", "electronID", 0.8, 1);
   // RooRealVar y("D0_useCMSFrame_p", "p*(D0)", 2.5, 5);
   RooRealVar z("Ds_dM_D0pi", "m(D0)", 0.135, 0.16);

   RooDataSet data("data", "dataset with mass", tree, z);
   // RooDataSet data("data", "dataset with mass", tree, RooArgSet(M,z));
   //=========================================================================================================================================================================//
	// C r e a t e   m o d e l   f o r   s i g n a l   s a m p l e
   // -------------------------------------------------------------- 
   double xL = 0.135;
   double xR = 0.16;
   // Gaussian Resolution model
	RooRealVar sigmean1("sigmean1", "mean", xL, xR);
	RooRealVar sigwidth1("sigwidth1", "width", 0, 0.01);
	RooGaussModel gaussm1("gaussm1", "Signal PDF", z, sigmean1, sigwidth1);
	// Gaussian Resolution model
	RooRealVar sigmean2("sigmean2", "mean", xL, xR);
	RooRealVar sigwidth2("sigwidth2", "width", 0, 0.01);
	RooGaussModel gaussm2("gaussm2", "Signal PDF", z, sigmean2, sigwidth2);

	// Add the components
	RooRealVar g1frac("g1frac","fraction of gauss1", 0., 1.);
	RooRealVar g2frac("g2frac","fraction of gauss2", 0., 1.);
	RooAddPdf Signal("Signal","g1+g2",RooArgList(gaussm1,gaussm2), RooArgList(g1frac,g2frac));
   //=========================================================================================================================================================================//
   // C r e a t e   m o d e l   f o r   b a c k g r o u n d   s a m p l e
   // -------------------------------------------------------------- 
//    // Constant Background:
//    RooRealVar c5("c5", "coefficient #5", 0., 1);
//    // RooPolynomial Polynomial("Polynomial", "background p.d.f.", x, c5);
//    RooPolynomial Background("Background", "background p.d.f.", z, c5);
   //-----------
   RooRealVar b("b", "b", -40., 20.);
   RooRealVar N("N", "N", 1000., 100000.);
   RooRealVar P("P", "P", 0., 1.);
   // RooGenericPdf Signal("Signal","Signal","N*((dM_D0pi-0.13957)^P)*exp( b*(dM_D0pi-0.13957))", RooArgSet(M,N,b,P));
   RooGenericPdf Background("Background", "Background", "(Ds_dM_D0pi - 0.13957 > 0 ? ( N*((Ds_dM_D0pi-0.13957)^P)*exp( b*(Ds_dM_D0pi-0.13957))) : 0.0)", RooArgSet(z, N, b, P));
   //=========================================================================================================================================================================//
   // C o n s t r u c t   a   S i g na l   a n d   B a c k g r o u n d   P D F:
   // --------------------------------------------------------------
   RooRealVar nsig("nsig", "#Signal events", 0., 100000);
   RooRealVar nbkg("nbkg", "#background events", 0., 100000);
   RooAddPdf model("model", "g+a", RooArgList(Signal, Background), RooArgList(nsig, nbkg));
   //=========================================================================================================================================================================//
   // Perform extended ML fit of data: The first line shows more information
   model.fitTo(data,PrintLevel(0));
   // If the fit is too broad and does not seem to fit a peak use the code below. Use it 
   // once to determine where the peak should be and then adjust the mean above.
   // model.fitTo(data, Minos(kTRUE), PrintLevel(-1), Range(2,2.1));
   //=========================================================================================================================================================================//
   // P l o t   M o d e l   a n d   P r e p a r e   P a r a m e t e r   V a l u e s
   // ---------------------------------------------------------------------------
   // Create a RooPlot to draw on. We don't manage the memory of the returned
   // pointer. Instead we let it leak such that the plot still exists at the end of
   // the macro and we can take a look at it.
   RooPlot* frame = z.frame();
   data.plotOn(frame, Binning(bin));
   model.plotOn(frame, Components(Background), LineStyle(kDashed), LineColor(kRed));
   model.plotOn(frame, LineColor(kBlue),RooFit::Name("model"));
   //=========================================================================================================================================================================//
   // P l o t    L a b e l s
   // ---------------------------------------------------------------------------
   // Add text to frame
   TLatex *latex = new TLatex();
   TString lum ("");
   lum.Form("#splitline{Generic c#bar{c} Events}{#scale[0.5]{#int}#it{L} dt=100.12 fb^{-1}}");
   // lum.Form("#splitline{Generic c#bar{c} Events}{#scale[0.5]{#int}#it{L} dt=107.33 fb^{-1}}");
   latex->SetNDC();
   latex->SetTextSize(0.03); 
   
   auto PerBin = ((xR-xL)/(bin))*1000;
   // TString title ("");
   TString y_axis ("");
   TString x_axis ("");

   // title.Form("Generic d#bar{d} Events");
   y_axis.Form("Entries/(%.2f MeV/c^{2})", PerBin);
   x_axis.Form("#Deltam  [GeV/c^{2}]");

   frame->SetTitle("");
   frame->GetYaxis()->SetTitle(y_axis);
   frame->GetYaxis()->SetRangeUser(0,14000);
   frame->GetXaxis()->SetTitle(x_axis);
   // frame->GetXaxis()->CenterTitle(true);
   frame->Draw();

   // The line below prints the luminosity onto the screen.
   // latex->DrawLatex(0.6,0.2, lum);
   latex->DrawLatex(0.18,0.85, lum);
   //=========================================================================================================================================================================//
   // P a r a m e t e r s
   // ---------------------------------------------------------------------------
   // The code below grabs the parameters and then places them in the legend
   Float_t N_value = data.sumEntries();
   Float_t nbkg_value = nbkg.getValV();
   Float_t nbkg_error = nbkg.getError();
   Float_t nsig_value = nsig.getValV();
   Float_t nsig_error = nsig.getError();
   // Float_t sigmean1_value = sigmean1.getValV();
   // Float_t sigmean2_value = sigmean2.getValV();
   // Float_t sigwidth1_value = sigwidth1.getValV();
   // Float_t sigwidth2_value = sigwidth2.getValV();
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
   auto y1 = yIF1 + 0.85*dy;
   auto y2 = 0.895;
   auto legend = new TLegend(x1, y1, x2, y2);

   // Parameter String
   TString nbkg_para ("");
   TString nsig_para ("");
   // TString sigmean1_para ("");
   // TString sigmean2_para ("");
   // TString sigwidth1_para ("");
   // TString sigwidth2_para ("");
   // N_para.Form("N = %.2f", N_value);
//    nbkg_para.Form("nbkg = %.2f #pm %.2f", nbkg.getValV(), nbkg.getError());
//    nsig_para.Form("nsig = %.2f #pm %.2f", nsig.getValV(), nsig.getError());
   // sigmean1_para.Form("mean_{1} = %.8f", sigmean1_value);
   // sigmean2_para.Form("mean_{2} = %.8f", sigmean2_value);
   // sigwidth1_para.Form("#sigma_{1} = %.8f", sigwidth1_value);
   // sigwidth2_para.Form("#sigma_{2} = %.8f", sigwidth2_value);

   legend->AddEntry(frame->FindObject("data"),"data","pe");
   legend->AddEntry((TObject*)nullptr, "", "");
//    legend->AddEntry((TObject*)nullptr, nbkg_para, "");
   legend->AddEntry("model", "Total Fit", "pl" );
   legend->AddEntry((TObject*)nullptr, "", "");
//    legend->AddEntry((TObject*)nullptr, nsig_para, "");

   legend -> SetBorderSize(0);
   legend -> SetTextFont(132);
   legend -> SetTextSize(0.030);
   legend -> SetNColumns(2);
   legend -> SetMargin(0.3);
   legend->Draw();
   //=========================================================================================================================================================================//
   // S h o w   P a r a m e t e r   R e s u l t s
   // -------------------------------------------------------
   // nbkg.setConstant(true);
   // nbkg.removeError();
   RooArgSet* params = model.getVariables() ; 
   params->Print("v") ;

   data.Print();
   //=========================================================================================================================================================================//
   // S h o w   r e s i d u a l   d i s t s
   // -------------------------------------------------------
   pad2->cd();
   // Construct a histogram with the residuals of the data w.r.t. the curve
   RooHist *hresid = frame->residHist();
   // Create a new frame to draw the residual distribution and add the distribution to the frame
   RooPlot *frame1 = z.frame(Title("Residual Distribution"));
   frame1->addPlotable(hresid, "P");
   frame1->Draw();
   frame1->SetTitle("");
   frame1->SetYTitle("Residual");
   frame1->GetYaxis()->CenterTitle(true);
   frame1->SetXTitle("");
   // frame1->GetXaxis()->SetTitleSize(.25);
   frame1->GetYaxis()->SetTitleSize(0.10);
   frame1->GetYaxis()->SetLabelSize(0.10);
   frame1->GetXaxis()->SetLabelSize(0.10);
   // frame1->SetTitleOffset(1);
   frame1->SetTitleOffset(0.35, "Y");

   c->cd();
   //=========================================================================================================================================================================//
   // S a v e   F i l e s
   // -------------------------------------------------------------
   c->SaveAs("C:/Users/Adilm/Desktop/Update/Next_Meeting/rf04_PeakingBackground.png");
   //=========================================================================================================================================================================//
}