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

void rs01_Projection()
{
    // S e t u p   t h e   F r a m e   f o r   t h e   P l o t
    //------------------------------------------------------------
    // make our canvas
    TCanvas *cdata = new TCanvas("sPlot", "sPlot demo", 400, 600);
    TPad *pad1 = new TPad("pad1","pad1",0,0.5,1,1);
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.5);
    pad1->SetTopMargin(0.05);
    pad2->SetTopMargin(0.05);
    // pad1->SetBottomMargin(0.05);
    pad1->SetLeftMargin(0.15);
    pad1->SetRightMargin(0.05);
    pad2->SetLeftMargin(0.15);
    pad2->SetRightMargin(0.05);
    // pad2->SetTopMargin(0.01);
    pad2->SetBottomMargin(0.12);
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
    // sPlot:
    //------------
    TFile *f = new TFile("Ds2D0enue_Analysis/05-ROOT/sPlot_results/rs00-RealD0.root");
    RooWorkspace *w = (RooWorkspace *)f->Get("w_splot");
    RooRealVar *dM = w->var("D0_dM");
    RooRealVar *DeltaM = w->var("Ds_diff_D0pi");

    // note, we get the dataset with sWeights
    auto& data = static_cast<RooDataSet&>(*w->data("dataWithSWeights"));

    RooDataSet data_signal(data.GetName(), data.GetTitle(), &data, *data.get(), nullptr, "nsig_sw");
    RooDataSet data_background(data.GetName(), data.GetTitle(), &data, *data.get(), nullptr, "nbkg_sw");
    //=========================================================================================================================================================================//
    // I m p o r t s   D a t a   f r o m   t r e e   b r a n c h:
    // --------------------------------------------------------------
    TFile *file = new TFile("C03-Grid/Completed/Ds2D0e-Generic_Ds_111424_0_All.root");
    TTree *tree = (TTree*)file->Get("Dstree");
    //=========================================================================================================================================================================//
    // Original Data:
    // --------------------------------------------------------------
    RooRealVar dM_Orig("D0_dM", "Mass Difference", -0.05, 0.05);
    RooRealVar DeltaM_Orig("Ds_diff_D0pi", "Mass Difference", 0.1, 1);
    RooRealVar Veto_Orig("Ds_gammaveto_M_Correction","Veto",0.1,100);
    RooRealVar BCS_Orig("Ds_chiProb_rank","BCS",-1,1.25);
    RooRealVar Signal("D0_isSignal","BCS",0.5,1.25);
    RooRealVar NoSignal("D0_isSignal","BCS",-1,0.5);

    RooDataSet data_sig("data_sig", "dataset with mass", tree, RooArgSet(dM_Orig,DeltaM_Orig,Veto_Orig,BCS_Orig,Signal));
    RooDataSet data_bkg("data_bkg", "dataset with mass", tree, RooArgSet(dM_Orig,DeltaM_Orig,Veto_Orig,BCS_Orig,NoSignal));
    //=========================================================================================================================================================================//
    // Plot Model:
    // --------------------
    auto PerBin = ((xR-xL)/(bin))*1000;

    TString x_axis ("");
    TString y_axis ("");
    x_axis.Form("#Delta m  [GeV/c^{2}]");
    y_axis.Form("Entries/(%.2f MeV/c^{2})", PerBin);

    // Here we make plots of the discriminating variable (invMass) after the fit
    // and of the control variable (isolation) after unfolding with sPlot.
    std::cout << "make plots" << std::endl;

    // Now use the sWeights to show isolation distribution for Z and QCD.
    // The SPlot class can make this easier, but here we demonstrate in more
    // detail how the sWeights are used.  The SPlot class should make this
    // very easy and needs some more development.

    // Plot isolation for Z component.
    // Do this by plotting all events weighted by the sWeight for the Z component.
    // The SPlot class adds a new variable that has the name of the corresponding
    // yield + "_sw".
    // cdata->cd(1);

    RooPlot *frame1 = DeltaM->frame(Title("Isolation distribution with s weights to project out Signal"));
    // Since the data are weighted, we use SumW2 to compute the errors.
    data_signal.plotOn(frame1, DataError(RooAbsData::SumW2), MarkerColor(kRed), Binning(bin), RooFit::Name("sweightsS"));
    data_sig.plotOn(frame1, DataError(RooAbsData::SumW2), MarkerColor(kBlue), Binning(bin), RooFit::Name("Signal"));

    // pad1->SetLogy();
    frame1->Draw();
    frame1->SetTitle("");
    frame1->GetYaxis()->SetTitle(y_axis);
    frame1->GetXaxis()->SetTitle(x_axis);
    frame1->GetXaxis()->SetTitleSize(0.05);
    // frame1->GetXaxis()->SetLabelSize(0.10);

    // Legend
    //---------------------------------------------------------------------------------------------------------------------------
    auto xIF11 = pad1 -> GetUxmin() + pad1 -> GetLeftMargin();
    auto xIF21 = pad1 -> GetUxmax() - pad1 -> GetRightMargin();
    auto yIF11 = pad1 -> GetUymin() + pad1 -> GetBottomMargin();
    auto yIF21 = pad1 -> GetUymax() - pad1 -> GetTopMargin();
    auto dx1 = xIF21 - xIF11;
    auto dy1 = yIF21 - yIF11;
    auto x11 = xIF11 + 0.5*dx1;
    auto x21 = xIF21;
    auto y11 = yIF11 + 0.75*dy1;
    auto y21 = 0.93;
    auto legend1 = new TLegend(x11, y11, x21, y21);

    legend1->AddEntry((TObject *)0, " ", "");
    legend1->AddEntry("sweightsS","isSignal=1","p");
    legend1->AddEntry((TObject *)0, " ", "");
    legend1->AddEntry("Signal","sWeights","p");

    legend1 -> SetBorderSize(0);
    legend1 -> SetTextFont(132);
    legend1 -> SetTextSize(0.04);
    legend1 -> SetNColumns(2);
    legend1 -> SetMargin(0.3);
    legend1 -> Draw();
    //---------------------------------------------------------------------------------------------------------------------------

    // Plot isolation for QCD component.
    // Eg. plot all events weighted by the sWeight for the QCD component.
    // The SPlot class adds a new variable that has the name of the corresponding
    // yield + "_sw".
    pad2->cd();
    RooPlot *frame2 = DeltaM->frame(Title("Isolation distribution with s weights to project out Background"));
    data_background.plotOn(frame2, DataError(RooAbsData::SumW2), MarkerColor(kRed), Binning(bin), RooFit::Name("sweightsB"));
    data_bkg.plotOn(frame2, DataError(RooAbsData::SumW2), MarkerColor(kBlue), Binning(bin), RooFit::Name("Background"));

    frame2->Draw();
    frame2->SetTitle("");
    frame2->GetYaxis()->SetTitle(y_axis);
    frame2->GetXaxis()->SetTitle(x_axis);
    frame2->GetXaxis()->SetTitleSize(0.05);
    // frame2->GetXaxis()->SetLabelSize(0.10);

    // Legend
    //---------------------------------------------------------------------------------------------------------------------------
    auto xIF12 = pad2 -> GetUxmin() + pad2 -> GetLeftMargin();
    auto xIF22 = pad2 -> GetUxmax() - pad2 -> GetRightMargin();
    auto yIF12 = pad2 -> GetUymin() + pad2 -> GetBottomMargin();
    auto yIF22 = pad2 -> GetUymax() - pad2 -> GetTopMargin();
    auto dx2 = xIF22 - xIF12;
    auto dy2 = yIF22 - yIF12;
    auto x12 = xIF12 + 0.5*dx2;
    auto x22 = xIF22;
    auto y12 = yIF12 + 0.75*dy2;
    auto y22 = 0.93;
    auto legend2 = new TLegend(x12, y12, x22, y22);

    legend2->AddEntry((TObject *)0, " ", "");
    legend2->AddEntry("sweightsB","isSignal=0","p");
    legend2->AddEntry((TObject *)0, " ", "");
    legend2->AddEntry("Background","sWeights","p");

    legend2 -> SetBorderSize(0);
    legend2 -> SetTextFont(132);
    legend2 -> SetTextSize(0.04);
    legend2 -> SetNColumns(2);
    legend2 -> SetMargin(0.3);
    legend2 -> Draw();
    //---------------------------------------------------------------------------------------------------------------------------

    cdata->cd();

    cdata->SaveAs("/home/belle2/amubarak/Ds2D0enue_Analysis/05-ROOT/sPlot_results/rs01-Projection.png");
}