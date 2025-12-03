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
#include <RooLinearVar.h>
#include <RooConstVar.h>
using namespace RooFit;

void rs02_Dstar0()
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
   double xL = 0.16;
   double xR = 0.45;
   //=========================================================================================================================================================================//
   // I m p o r t s   D a t a   f r o m   t r e e   b r a n c h:
   // --------------------------------------------------------------
   TFile *file = new TFile("/home/belle2/amubarak/C03-Grid/TopoAna.root");
   TTree *tree = (TTree*)file->Get("Dstar0");
   //=========================================================================================================================================================================//
   // D a t a:
   // --------------------------------------------------------------
   RooRealVar DeltaM("Ds_diff_D0pi", "Mass Difference", xL, xR);
   // RooRealVar Veto("Ds_gammaveto_M_Correction","Veto",-1,100);
   // RooRealVar Ds_Dstar0("Ds_Dstar0","Ds_Dstar0",0.0,1.0);
   RooRealVar BS("Ds_extraInfo_BkgBDT", "Background Suppression", 0.0, 1.0);
   RooDataSet data("data", "dataset with mass", tree, RooArgSet(DeltaM,BS));
   // RooDataSet data("data", "dataset with mass", tree, RooArgSet(DeltaM,isSignal,BS));

   // Apply multiple selection cuts: 
   // - BS >= 0.531
   // - isSignal==1
   // - DeltaM within [0.1, 0.55] GeV/c^2
   RooDataSet *data_Dstar0 = (RooDataSet*)data.reduce("Ds_extraInfo_BkgBDT >= 0.531 && Ds_diff_D0pi >= 0.16 && Ds_diff_D0pi <= 0.45");
   // RooDataSet **data_Dstar0 = (RooDataSet*)data.reduce("Ds_extraInfo_BkgBDT >= 0.531 && Ds_ifNANgiveX_isSignal_5==1 && Ds_diff_D0pi >= 0.1 && Ds_diff_D0pi <= 0.55");
   //=========================================================================================================================================================================//
   // C r e a t e   C o m b i n e d   M o d e l
   // ------------------------------------------------------------------------------
   // // First Attempt
   // //---------------------
   // // Generic PDF
   // RooRealVar B0_Dstar0Bkg("B0_Dstar0Bkg", "B0_Dstar0Bkg", -8.3111e+00, -40., 20.);
   // RooRealVar P0_Dstar0Bkg("P0_Dstar0Bkg", "P0_Dstar0Bkg", 8.1845e-01, 0., 5.);
   // B0_Dstar0Bkg.setConstant(true);
   // P0_Dstar0Bkg.setConstant(true);
   // RooGenericPdf Generic0_Dstar0Bkg_Dstar0Bkg("Generic0_Dstar0Bkg_Dstar0Bkg","Generic PDF","(Ds_diff_D0pi - 0.13953 > 0 ? ( ((Ds_diff_D0pi-0.13953)^P0_Dstar0Bkg)*exp( B0_Dstar0Bkg*(Ds_diff_D0pi-0.13953))) : 0.0)", RooArgSet( DeltaM, B0_Dstar0Bkg, P0_Dstar0Bkg));

   // // Generic PDF
   // RooRealVar B1_Dstar0Bkg("B1_Dstar0Bkg", "B1_Dstar0Bkg", -3.2604e+00, -40., 20.);
   // RooRealVar P1_Dstar0Bkg("P1_Dstar0Bkg", "P1_Dstar0Bkg", 1.6060e+00, 0., 5.);
   // B1_Dstar0Bkg.setConstant(true);
   // P1_Dstar0Bkg.setConstant(true);
   // RooGenericPdf Generic1_Dstar0Bkg_Dstar0Bkg("Generic1_Dstar0Bkg_Dstar0Bkg","Generic PDF","(Ds_diff_D0pi - 0.13953 > 0 ? ( ((Ds_diff_D0pi-0.13953)^P1_Dstar0Bkg)*exp( B1_Dstar0Bkg*(Ds_diff_D0pi-0.13953))) : 0.0)", RooArgSet( DeltaM, B1_Dstar0Bkg, P1_Dstar0Bkg));

   // // Combined Background
   // RooRealVar f0_Dstar0Bkg("f0_Dstar0Bkg", "f0", 5.3990e-01, 0, 1);
   // RooRealVar f1_Dstar0Bkg("f1_Dstar0Bkg", "f1", 4.2536e-01, 0, 1);
   // f0_Dstar0Bkg.setConstant(true);
   // f1_Dstar0Bkg.setConstant(true);
   // RooAddPdf CFBkgModel("CFBkgModel", "b+c", RooArgList(Generic0_Dstar0Bkg_Dstar0Bkg,Generic1_Dstar0Bkg_Dstar0Bkg), RooArgList(f0_Dstar0Bkg,f1_Dstar0Bkg));

   // // Second Attempt
   // //---------------------
   // // First Gamma Function
   // RooRealVar x00_Dstar0Bkg("x00_Dstar0Bkg","x00_Dstar0Bkg", 1.3124e-01, 0.12, 0.2);
   // RooRealVar a0_Dstar0Bkg("a0_Dstar0Bkg","a0_Dstar0Bkg", 1.1902e+00, 0.5e+00, 2.0e+00);
   // RooRealVar b0_Dstar0Bkg("b0_Dstar0Bkg","b0_Dstar0Bkg", 6.1763e+01, 5.0e+01, 7.0e+01);    
   // x00_Dstar0Bkg.setConstant();
   // a0_Dstar0Bkg.setConstant();
   // b0_Dstar0Bkg.setConstant();
   // RooFormulaVar x00prime_Dstar0Bkg("#x00prime_Dstar0Bkg", "x[0] < x[1] ? x[0] : x[1]", RooArgList(DeltaM, x00_Dstar0Bkg));     
   // RooFormulaVar gamma0_Dstar0Bkg("#gamma0_Dstar0Bkg", "a0_Dstar0Bkg + 1", RooArgList(a0_Dstar0Bkg));    
   // RooFormulaVar beta0_for_gamma_Dstar0Bkg("beta0_for_gamma_Dstar0Bkg", "1./b0_Dstar0Bkg", RooArgList(b0_Dstar0Bkg));    
   // RooGamma GammaModel0_Dstar0Bkg("GammaModel0_Dstar0Bkg", "Gamma pdf", DeltaM, gamma0_Dstar0Bkg, beta0_for_gamma_Dstar0Bkg, x00prime_Dstar0Bkg);

   // // Second Gamma Function
   // RooRealVar x01_Dstar0Bkg("x01_Dstar0Bkg","x01_Dstar0Bkg", 1.5263e-01, 0.13, 0.2);
   // RooRealVar a1_Dstar0Bkg("a1_Dstar0Bkg","a1_Dstar0Bkg", 6.7404e-01, 6.0e-01, 8.0e-01);
   // RooRealVar b1_Dstar0Bkg("b1_Dstar0Bkg","b1_Dstar0Bkg", 4.9962e+00, 4.0e+00, 6.0e+00);     
   // x01_Dstar0Bkg.setConstant();
   // a1_Dstar0Bkg.setConstant();
   // b1_Dstar0Bkg.setConstant();
   // RooFormulaVar x01prime_Dstar0Bkg("#x01prime_S", "x[0] < x[1] ? x[0] : x[1]", RooArgList(DeltaM, x01_Dstar0Bkg));     
   // RooFormulaVar gamma1_Dstar0Bkg("#gamma1_Dstar0Bkg", "a1_Dstar0Bkg + 1", RooArgList(a1_Dstar0Bkg));    
   // RooFormulaVar beta1_for_gamma_Dstar0Bkg("beta1_for_gamma_Dstar0Bkg", "1./b1_Dstar0Bkg", RooArgList(b1_Dstar0Bkg));    
   // RooGamma GammaModel1_Dstar0Bkg("GammaModel1_Dstar0Bkg", "Gamma pdf", DeltaM, gamma1_Dstar0Bkg, beta1_for_gamma_Dstar0Bkg, x01prime_Dstar0Bkg);

   // // Third Gamma Function
   // RooRealVar x02_Dstar0Bkg("x02_Dstar0Bkg","x02_Dstar0Bkg", 1.6356e-01, 0.13, 0.2);
   // RooRealVar a2_Dstar0Bkg("a2_Dstar0Bkg","a2_Dstar0Bkg", 1.4099e+00, 0.5e+00, 2.0e+00);
   // RooRealVar b2_Dstar0Bkg("b2_Dstar0Bkg","b2_Dstar0Bkg", 2.9852e+01, 2.0e+01, 4.0e+01);     
   // x02_Dstar0Bkg.setConstant();
   // a2_Dstar0Bkg.setConstant();
   // b2_Dstar0Bkg.setConstant();
   // RooFormulaVar x02prime_Dstar0Bkg("#x02prime_S", "x[0] < x[1] ? x[0] : x[1]", RooArgList(DeltaM, x02_Dstar0Bkg));     
   // RooFormulaVar gamma2_Dstar0Bkg("#gamma2_Dstar0Bkg", "a2_Dstar0Bkg + 1", RooArgList(a2_Dstar0Bkg));    
   // RooFormulaVar beta2_for_gamma_Dstar0Bkg("beta2_for_gamma_Dstar0Bkg", "1./b2_Dstar0Bkg", RooArgList(b2_Dstar0Bkg));    
   // RooGamma GammaModel2_Dstar0Bkg("GammaModel2_Dstar0Bkg", "Gamma pdf", DeltaM, gamma2_Dstar0Bkg, beta2_for_gamma_Dstar0Bkg, x02prime_Dstar0Bkg);

   // // Fourth Gamma Function
   // RooRealVar x03_Dstar0Bkg("x03_Dstar0Bkg","x03_Dstar0Bkg", 1.3018e-01, 0.12, 0.15);
   // RooRealVar a3_Dstar0Bkg("a3_Dstar0Bkg","a3_Dstar0Bkg", 1.0441e+01, 0.5e+01, 2.0e+01);
   // RooRealVar b3_Dstar0Bkg("b3_Dstar0Bkg","b3_Dstar0Bkg", 4.9884e+01, 4.0e+01, 6.0e+01);     
   // x03_Dstar0Bkg.setConstant();
   // a3_Dstar0Bkg.setConstant();
   // b3_Dstar0Bkg.setConstant();
   // RooFormulaVar x03prime_Dstar0Bkg("#x03prime_S", "x[0] < x[1] ? x[0] : x[1]", RooArgList(DeltaM, x03_Dstar0Bkg));     
   // RooFormulaVar gamma3_Dstar0Bkg("#gamma3_Dstar0Bkg", "a3_Dstar0Bkg + 1", RooArgList(a3_Dstar0Bkg));    
   // RooFormulaVar beta3_for_gamma_Dstar0Bkg("beta3_for_gamma_Dstar0Bkg", "1./b3_Dstar0Bkg", RooArgList(b3_Dstar0Bkg));    
   // RooGamma GammaModel3_Dstar0Bkg("GammaModel3_Dstar0Bkg", "Gamma pdf", DeltaM, gamma3_Dstar0Bkg, beta3_for_gamma_Dstar0Bkg, x03prime_Dstar0Bkg);

   // // Add the components
   // RooRealVar f0_Dstar0Bkg("f0_Dstar0Bkg", "f0", 3.1274e-02, 0, 1);
   // RooRealVar f1_Dstar0Bkg("f1_Dstar0Bkg", "f1", 6.4932e-01, 0, 1);
   // RooRealVar f2_Dstar0Bkg("f2_Dstar0Bkg", "f2", 9.2640e-02, 0, 1);
   // RooRealVar f3_Dstar0Bkg("f3_Dstar0Bkg", "f3", 3.5501e-02, 0, 1);
   // f0_Dstar0Bkg.setConstant();
   // f1_Dstar0Bkg.setConstant();
   // f2_Dstar0Bkg.setConstant();
   // f3_Dstar0Bkg.setConstant();
   // // RooAddPdf Dstar0BkgModel("Dstar0BkgModel","g1+g2",RooArgList(GammaModel0_Dstar0Bkg,GammaModel1_Dstar0Bkg,GammaModel2_Dstar0Bkg), RooArgList(f0_Dstar0Bkg,f1_Dstar0Bkg,f2_Dstar0Bkg));
   // RooAddPdf Dstar0BkgModel("Dstar0BkgModel","g1+g2",RooArgList(GammaModel0_Dstar0Bkg,GammaModel1_Dstar0Bkg,GammaModel2_Dstar0Bkg,GammaModel3_Dstar0Bkg), RooArgList(f0_Dstar0Bkg,f1_Dstar0Bkg,f2_Dstar0Bkg,f3_Dstar0Bkg));

   // Third Attempt
   //----------------
   RooRealVar c0_Dstar0Bkg("c0_Dstar0Bkg", "coefficient #0", -1.1996e-01, -1.0, 1.0);
   RooRealVar c1_Dstar0Bkg("c1_Dstar0Bkg", "coefficient #1", -1.9521e-01, -1.0, 1.0);
   RooRealVar c2_Dstar0Bkg("c2_Dstar0Bkg", "coefficient #2",  1.0892e-01, -1.0, 1.0);
   RooRealVar c3_Dstar0Bkg("c3_Dstar0Bkg", "coefficient #3", -5.8254e-02, -1.0, 1.0);
   RooRealVar c4_Dstar0Bkg("c4_Dstar0Bkg", "coefficient #4",  2.7045e-02, -1.0, 1.0);
   RooRealVar c5_Dstar0Bkg("c5_Dstar0Bkg", "coefficient #5", -3.4500e-03, -1.0, 1.0);
   RooRealVar c6_Dstar0Bkg("c6_Dstar0Bkg", "coefficient #6",  -7.2964e-03, -1.0, 1.0);
   RooRealVar c7_Dstar0Bkg("c7_Dstar0Bkg", "coefficient #7",  2.9229e-03, -1.0, 1.0);
   RooRealVar c8_Dstar0Bkg("c8_Dstar0Bkg", "coefficient #8", -2.7331e-03, -1.0, 1.0);
   c0_Dstar0Bkg.setConstant(true);
   c1_Dstar0Bkg.setConstant(true);
   c2_Dstar0Bkg.setConstant(true);
   c3_Dstar0Bkg.setConstant(true);
   c4_Dstar0Bkg.setConstant(true);
   // c5_Dstar0Bkg.setConstant(true);
   // c6_Dstar0Bkg.setConstant(true);
   // c7_Dstar0Bkg.setConstant(true);
   // c8_Dstar0Bkg.setConstant(true);
   RooChebychev Dstar0BkgModel("Dstar0BkgModel", "Dstar0BkgModel", DeltaM, RooArgList(c0_Dstar0Bkg,c1_Dstar0Bkg,c2_Dstar0Bkg,c3_Dstar0Bkg,c4_Dstar0Bkg));

   // Fourth Way
   //------------------
   // // Reverse Argus transformation using RooLinearVar (instead of RooFormulaVar)
   // RooLinearVar DeltaM_reversed("DeltaM_reversed", "Reversed Mass", DeltaM, RooConst(-1.0), RooConst(xR));

   // // Define Reverse Argus function parameters
   // RooRealVar c0_Dstar0Bkg("c0_Dstar0Bkg", "Cutoff", xR, 4);
   // RooRealVar m0("m0", "m0", 1, 100);
   // RooRealVar p0("p0", "p0", 0.0, 10.0);
   // // RooFormulaVar x_flip0("x_flip0", "c0_Dstar0Bkg - Ds_diff_D0pi", RooArgList(c0_Dstar0Bkg, DeltaM));
   // RooArgusBG Argus_0("Argus_0", "1st Argus Function", DeltaM, c0_Dstar0Bkg, m0, p0);

   // // Combine Reverse Argus Functions
   // RooRealVar f0("f0", "Fraction #1", 1.0);
   // RooAddPdf Dstar0BkgModel("Dstar0BkgModel", "Combined Reverse Argus", RooArgList(Argus_0), RooArgList(f0));

   // // Fifth Way
   // //------------------------
   // RooRealVar var0("var0", "var0", 0.13, 0.15);
   // RooRealVar var1("var1", "var1", -10, 10);
   // RooRealVar var2("var2", "var2", -10, 10);
   // RooRealVar var3("var3", "var3", -10, 10);
   // RooDstD0BG Dstar0BkgModel("Dstar0BkgModel","Dstar0BkgModel", DeltaM, var0, var1, var2, var3);

   // // Sixth Attempt
   // //---------------------
   // // First Gamma Function
   // RooRealVar x00_Dstar0Bkg("x00_Dstar0Bkg","x00_Dstar0Bkg", 1.3955e-01, 0.13, 0.15);
   // RooRealVar a0_Dstar0Bkg("a0_Dstar0Bkg","a0_Dstar0Bkg", 5.3309e-01, 0, 200);
   // RooRealVar b0_Dstar0Bkg("b0_Dstar0Bkg","b0_Dstar0Bkg", 5.2582e+00, 0, 200);    
   // x00_Dstar0Bkg.setConstant();
   // a0_Dstar0Bkg.setConstant();
   // b0_Dstar0Bkg.setConstant();
   // RooFormulaVar x00prime_Dstar0Bkg("#x00prime_Dstar0Bkg", "x[0] < x[1] ? x[0] : x[1]", RooArgList(DeltaM, x00_Dstar0Bkg));     
   // RooFormulaVar gamma0_Dstar0Bkg("#gamma0_Dstar0Bkg", "a0_Dstar0Bkg + 1", RooArgList(a0_Dstar0Bkg));    
   // RooFormulaVar beta0_for_gamma_Dstar0Bkg("beta0_for_gamma_Dstar0Bkg", "1./b0_Dstar0Bkg", RooArgList(b0_Dstar0Bkg));    
   // RooGamma GammaModel0_Dstar0Bkg("GammaModel0_Dstar0Bkg", "Gamma pdf", DeltaM, gamma0_Dstar0Bkg, beta0_for_gamma_Dstar0Bkg, x00prime_Dstar0Bkg);

   // RooRealVar c0_Dstar0Bkg("c0_Dstar0Bkg", "coefficient #0", -7.8719e+01, -8.0e+01, -3.0e+01);
   // RooRealVar c1_Dstar0Bkg("c1_Dstar0Bkg", "coefficient #1", 1.3065e+03, 0.0e+03, 3.0e+03);
   // RooRealVar c2_Dstar0Bkg("c2_Dstar0Bkg", "coefficient #2", 1.1131e-02, 0.0e+01, 5.0e+01);
   // RooRealVar c3_Dstar0Bkg("c3_Dstar0Bkg", "coefficient #3", -7.1648e+03 -9.0e+03, -5.0e+03);
   // RooRealVar c4_Dstar0Bkg("c4_Dstar0Bkg", "coefficient #4", 9.5703e+03, 6.0e+03, 10.0e+03);
   // RooRealVar c5_Dstar0Bkg("c5_Dstar0Bkg", "coefficient #5", 9.9666e+03, -10000, 10000);
   // // RooRealVar c6_Dstar0Bkg("c6_Dstar0Bkg", "coefficient #6", -10000, 10000);
   // c0_Dstar0Bkg.setConstant(true);
   // c1_Dstar0Bkg.setConstant(true);
   // c2_Dstar0Bkg.setConstant(true);
   // c2_Dstar0Bkg.setConstant(true);
   // c3_Dstar0Bkg.setConstant(true);
   // c4_Dstar0Bkg.setConstant(true);
   // c5_Dstar0Bkg.setConstant(true);
   // RooPolynomial PolykgModel("PolykgModel", "Other", DeltaM, RooArgList(c0_Dstar0Bkg,c1_Dstar0Bkg,c2_Dstar0Bkg,c3_Dstar0Bkg,c4_Dstar0Bkg,c5_Dstar0Bkg));

   // // Add the components
   // RooRealVar f0_Dstar0Bkg("f0_Dstar0Bkg", "f0", 7.9662e-01, 0, 1);
   // RooRealVar f1_Dstar0Bkg("f1_Dstar0Bkg", "f1", 1.4062e-08, 0, 1);
   // f0_Dstar0Bkg.setConstant();
   // f1_Dstar0Bkg.setConstant();
   // // f2_Dstar0Bkg.setConstant();
   // // f3_Dstar0Bkg.setConstant();
   // // RooAddPdf Dstar0BkgModel("Dstar0BkgModel","g1+g2",RooArgList(GammaModel0_Dstar0Bkg,GammaModel1_Dstar0Bkg,GammaModel2_Dstar0Bkg), RooArgList(f0_Dstar0Bkg,f1_Dstar0Bkg,f2_Dstar0Bkg));
   // RooAddPdf Dstar0BkgModel("Dstar0BkgModel","g1+g2",RooArgList(GammaModel0_Dstar0Bkg,PolykgModel), RooArgList(f0_Dstar0Bkg,f1_Dstar0Bkg));

   // Final Background Count
   //-----------------------------------------
   RooRealVar Nsig_Dstar0Bkg("Nsig_Dstar0Bkg", "#Signal events", 1.4364e+05, 1.0e+04, 3.0e+05);
   RooAddPdf model("model", "g+a", RooArgList(Dstar0BkgModel), RooArgList(Nsig_Dstar0Bkg));
   //=========================================================================================================================================================================//
   // S e t   t h e   C o n s t a n t   f o r   t h e    V a r i a b l e s
   //--------------------------------------------------------------------------
   // // The code below fixes the value for a certain value.
   // B0_Dstar0Bkg.setVal(-8.3580e+00);
   // B1_Dstar0Bkg.setVal(-3.2331e-01);
   // P0_Dstar0Bkg.setVal(1.1092e+00);
   // P1_Dstar0Bkg.setVal(7.1018e-02);
   // f0_Dstar0Bkg.setVal(2.0123e-04);
   // f1_Dstar0Bkg.setVal(2.1654e-04);

   // // The code below removes the range and fixes the value.
   // B0_Dstar0Bkg.setConstant(true);
   // B1_Dstar0Bkg.setConstant(true);
   // P0_Dstar0Bkg.setConstant(true);
   // P1_Dstar0Bkg.setConstant(true);
   // f0_Dstar0Bkg.setConstant(true);
   // f1_Dstar0Bkg.setConstant(true);

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
   std::unique_ptr<RooFitResult> r{model.fitTo(*data_Dstar0, Save(), PrintLevel(-1), Range(xL,xR))};
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
   *data_Dstar0->plotOn(frame, Binning(bin));
   // model.plotOn(frame, Components(Generic0_Dstar0Bkg_Dstar0Bkg),       LineStyle(kDashed), LineColor(46), RooFit::Name("Background"));
   // model.plotOn(frame, Components(Generic1_Dstar0Bkg_Dstar0Bkg), LineStyle(kDashed), LineColor(38), RooFit::Name("Signal"));
   model.plotOn(frame, LineColor(kOrange+8), RooFit::Name("model"));
   //=========================================================================================================================================================================//
   // C a l c u l a t e   c h i ^ 2
   // ------------------------------
   // Show the chi^2 of the curve w.r.t. the histogram
   // If multiple curves or datasets live in the frame you can specify
   // the name of the relevant curve and/or dataset in chiSquare()
   Int_t npar = model.getParameters(*data_Dstar0)->selectByAttrib("Constant",kFALSE)->getSize(); //select floating parameters and count their number

   cout << "NDF = " << npar << endl;
   cout << "chi^2/NDF = " << frame->chiSquare(npar) << endl;

   Double_t chi2ndf = frame->chiSquare(npar);
   //=========================================================================================================================================================================//
   // P l o t    L a b e l s
   // ---------------------------------------------------------------------------
      // Add text to frame
   TLatex *latex = new TLatex();
   TString lum ("");
   lum.Form("#splitline{Generic Events}{#scale[0.5]{#int}#it{L} dt=1.44 ab^{-1}}");
   // lum.Form("#splitline{Generic c#bar{c} Events}{#scale[0.5]{#int}#it{L} dt=200 fb^{-1}}");
   latex->SetNDC();
   latex->SetTextSize(0.03); 
   // // Add text to frame
   // TLatex *latex = new TLatex();
   // TString lum ("");
   // lum.Form("#splitline{Simulated Events}{100k Events}");
   // latex->SetNDC();
   // latex->SetTextSize(0.035); 
   
   // // Add "chi-squared fit" below "2M Events"
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
   frame->GetYaxis()->SetRangeUser(0,6000);
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
   Nsig_para.Form("N_{D^{*0}} = %.2f #pm %.2f", Nsig_Dstar0Bkg.getValV(), Nsig_Dstar0Bkg.getError());

   legend->AddEntry(frame->FindObject("data"),"data","pe");
   legend->AddEntry((TObject *)0, chi2, "");
   legend->AddEntry("model", "D^{*0} #rightarrow D^{0} #pi^{0}/#gamma", "pl" );
   // legend->AddEntry((TObject *)0, Nsig_para, "");

   legend -> SetBorderSize(0);
   legend -> SetTextFont(132);
   legend -> SetTextSize(0.032);
   legend -> SetNColumns(2);
   legend -> SetMargin(0.1);
   legend->Draw();
   //=========================================================================================================================================================================//
   // S h o w   P a r a m e t e r/ F i t    R e s u l t s
   // -------------------------------------------------------
   // When I use rooworkspace this first part will make sure that the parameters of the 
   // model will stay fixed.
   // DeltaM.setConstant(true);
   // B0_Dstar0Bkg.setConstant(true);
   // B0_Dstar0Bkg.removeError();
   // B1_Dstar0Bkg.setConstant(true);
   // B1_Dstar0Bkg.removeError();
   // P0_Dstar0Bkg.setConstant(true);
   // P0_Dstar0Bkg.removeError();
   // P1_Dstar0Bkg.setConstant(true);
   // P1_Dstar0Bkg.removeError();
   // f0_Dstar0Bkg.setConstant(true);
   // f0_Dstar0Bkg.removeError();
   // f1_Dstar0Bkg.setConstant(true);
   // f1_Dstar0Bkg.removeError();

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
   RooPlot *frame1 = DeltaM.frame(Title("Pull Distribution"));
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
   w->import(Dstar0BkgModel);

   // E n c o d e   d e f i n i t i o n   o f   p a r a m e t e r s   i n   w o r k s p a c e
   // ---------------------------------------------------------------------------------------
   // Define named sets "parameters" and "observables", which list which variables should be considered
   // parameters and observables by the users convention
   //
   // Variables appearing in sets _must_ live in the workspace already, or the autoImport flag
   // of defineSet must be set to import them on the fly. Named sets contain only references
   // to the original variables, therefore the value of observables in named sets already
   // reflect their 'current' value
   RooArgSet *params = Dstar0BkgModel.getVariables();
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
   w->saveSnapshot("reference_fit_Dstar0_Bkg", *params, true);

   // Print workspace contents
   w->Print();

   // S a v e   w o r k s p a c e   i n   f i l e
   // -------------------------------------------
   // Save the workspace into a ROOT file
   w->writeToFile("/home/belle2/amubarak/Ds2D0enue_Analysis/10-Save/Roofit_Suggestion_Likehood_Window/Fit_Parameter/rs02_Dstar0.root");
   //=========================================================================================================================================================================//
   c->SaveAs("/home/belle2/amubarak/Ds2D0enue_Analysis/10-Save/Roofit_Suggestion_Likehood_Window/Fit_Image/rs02_DeltaM_Dstar0.png");
   //=========================================================================================================================================================================//
   // Final Fit Status
   //---------------------------
   int fitStatus = r->status();  // Get the fit status
   int covQual = r->covQual();   // Get covariance quality

   std::cout << "Fit Status: " << fitStatus << std::endl;
   std::cout << "Covariance Matrix Quality: " << covQual << std::endl << endl;
}