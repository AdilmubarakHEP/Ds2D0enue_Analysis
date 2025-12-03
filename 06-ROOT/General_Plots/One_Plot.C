void One_Plot()
{
    TCanvas *c = new TCanvas();
    gStyle->SetOptStat(0);

    int bin = 50;
    double xL = 1.8;
    double xR = 2.2;
    // double xL = 0.;
    // double xR = 0.25;
    //=========================================================================================================================================================================//
    // G r a b   N e e d e d   F i l e   a n d   T r e e
    // ---------------------------------------------------------------------------
    TFile *file = new TFile("C:/Users/Adilm/Desktop/High_Energy_Physics/Analysis_Code/KEKCC_File/Ds2D0e-Generic_Ds_100124_0_ccbar_D.root");
    TTree *dst = (TTree *)file->Get("dstree");

    double m1;
    // dst->SetBranchAddress("dM_D0pi", &m1);
    dst->SetBranchAddress("Ds_M", &m1);
    // dst->SetBranchAddress("D0_dM", &m2);
    // dst->SetBranchAddress("e_electronID", &m3);
    // dst->SetBranchAddress("D0_useCMSFrame_p", &m4);
    //=========================================================================================================================================================================//
    // F i l l   H i s t o g r a m
    // ---------------------------------------------------------------------------
    TH1F *h = new TH1F("h", "h", bin, xL, xR);

    double N = 0;
    for (int iEntry = 0; dst->LoadTree(iEntry) >= 0; ++iEntry)
    {
      dst->GetEntry(iEntry);
      // if (m2 >= -0.25 && m2 <= -0.02 && m3 >= 0.8 && m4 >= 2.5) {
      //   draw->Fill(m2);
      // } else if (m2 >= 0.02 && m2 <= 0.25 && m3 >= 0.8 && m4 >= 2.5){
      //   draw->Fill(m2);
      // }
      // if (m3 >= 0.8 && m4 >= 2.5) {
      //   draw->Fill(m2);
      // }
      // if (m2 >= 0.04 && m2 <= 0.25 && m3 >= 0.8 && m4 >= 2.5){
      //   h->Fill(m1);
      //   N += 1;
      // }
      // if (m2 >= -0.02 && m2 <= 0.02 && m3 >= 0.8 && m4 >= 2.5)
      // {
        h->Fill(m1);
        N += 1;
      // }
      // if (m3 >= 0.8 && m4 >= 2.5){
      //   h->Fill(m2);
      //   N += 1;
      // }
    }

    // h->SetFillColor(38);
    //=========================================================================================================================================================================//
    // A d d   t e x t   t o   f r a m e
    // ---------------------------------------------------------------------------
    TLatex *latex = new TLatex();
    TString lum("");
    // lum.Form("#int #it{L} dt=107.33 fb^{-1}");
    lum.Form("#splitline{Generic c#bar{c} Events}{#scale[0.5]{#int}#it{L} dt=107.33 fb^{-1}}");
    latex->SetNDC();
    latex->SetTextSize(0.03);
    latex->DrawLatex(0.4, 0.4, lum);

    auto PerBin = ((xR - xL) / (bin)) * 100.;
    TString str("");
    str.Form("Entries/(%.2f MeV/c^{2})", PerBin);
    h->SetTitle("");
    h->GetXaxis()->SetTitle("#Deltam [GeV/c^{2}]");
    h->GetXaxis()->CenterTitle(true);
    h->GetYaxis()->SetTitle(str);
    // h->GetYaxis()->SetRangeUser(0,1400);

    h->Draw();
    latex->DrawLatex(0.7,0.55, lum);
    //=========================================================================================================================================================================//
    // L e g e n d
    // ---------------------------------------------------------------------------
    // The code below grabs the parameters and then places them in the legend
    double Area = h->Integral();

    // auto legend = new TLegend(0.1,0.7,0.48,0.9);
    auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->SetHeader("Parameters", "C"); // option "C" allows to center the header
    TString N_para("");
    TString Area_para("");
    N_para.Form("N = %.2f", N);
    Area_para.Form("Area = %.2f", Area);
    legend->AddEntry((TObject *)0, N_para, "");
    legend->AddEntry((TObject *)0, Area_para, "");
    legend->Draw();
    gStyle->SetOptStat("e");
    //=========================================================================================================================================================================//
    c->SaveAs("C:/Users/Adilm/Desktop/Update/Next_Meeting/DeltaM_WithoutD0_ccbar_Better.png");
}