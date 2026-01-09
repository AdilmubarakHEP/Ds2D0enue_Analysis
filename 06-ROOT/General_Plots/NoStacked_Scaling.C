double Scaling(int bin, double xL, double xR)
{
    //=========================================================================================================================================================================//
    // G r a b   N e e d e d   F i l e   a n d   T r e e
    // ---------------------------------------------------------------------------
    TFile *file = new TFile("C:/Users/Adilm/Desktop/High_Energy_Physics/Analysis_Code/KEKCC_File/Grid/Ds2D0e-Reconstruction_Ds_052024_ccbar.root");
    TTree *dst = (TTree *)file->Get("dstree");

    double m1, m2, m3, m4;
    dst->SetBranchAddress("dM_D0pi", &m1);
    // dst->SetBranchAddress("Ds_daughterMotherDiffOf_0_M", &m1);
    dst->SetBranchAddress("D0_dM", &m2);
    dst->SetBranchAddress("e_electronID", &m3);
    dst->SetBranchAddress("D0_useCMSFrame_p", &m4);
    //=========================================================================================================================================================================//
    // F i l l   H i s t o g r a m
    // ---------------------------------------------------------------------------
    // int bin = 50;
    // double xL = 0.135;
    // double xR = 0.25;
    TH1F *h1 = new TH1F("h1", "With D0", bin, xL, xR);
    TH1F *h2 = new TH1F("h2", "Without D0", bin, xL, xR);

    // double N1 = 0;
    // double N2 = 0;
    for (int iEntry = 0; dst->LoadTree(iEntry) >= 0; ++iEntry)
    {
      dst->GetEntry(iEntry);
      if (m2 >= -0.02 && m2 <= 0.02 && m3 >= 0.8 && m4 >= 2.5)
      {
        h1->Fill(m1);
        // N1 += 1;
      }
      if (m2 >= 0.04 && m2 <= 0.25 && m3 >= 0.8 && m4 >= 2.5){
        h2->Fill(m1);
        // N2 += 1;
      }
    }
    
    double Area1 = h1->Integral();
    double Area2 = h2->Integral();
    
    double scale_factor = Area1/Area2;

    return scale_factor;
}

void NoStacked_Scaling()
{
    //=========================================================================================================================================================================//
    // S e t u p
    // ---------------------------------------------------------------------------
    TCanvas *c = new TCanvas();
    gPad->SetLeftMargin(0.13);
    auto hs = new THStack("hs","");
    hs->SetMinimum(0.);
    hs->SetMaximum(2000.);

    int bin = 50;
    double xL = 0.135;
    double xR = 0.25;
    // double xL = 0.;
    // double xR = 0.25;
    //=========================================================================================================================================================================//
    // G r a b   N e e d e d   F i l e   a n d   T r e e
    // ---------------------------------------------------------------------------
    TFile *file = new TFile("C:/Users/Adilm/Desktop/High_Energy_Physics/Analysis_Code/KEKCC_File/Grid/Ds2D0e-Reconstruction_Ds_052024_ccbar.root");
    TTree *dst = (TTree *)file->Get("dstree");

    double m1, m2, m3, m4, m5, m6;
    dst->SetBranchAddress("dM_D0pi", &m1);
    // dst->SetBranchAddress("Ds_daughterMotherDiffOf_0_M", &m1);
    dst->SetBranchAddress("D0_dM", &m2);
    dst->SetBranchAddress("e_electronID", &m3);
    dst->SetBranchAddress("D0_useCMSFrame_p", &m4);
    //=========================================================================================================================================================================//
    // S e t u p   H i s t o g r a m
    // ---------------------------------------------------------------------------
    // TH1F *h = new TH1F("h", "Number of Candidates per Event", bin, xL, xR);
    auto h0st = new TH1F("Stacked", "With D0", bin, xL, xR);
    auto h1st = new TH1F("Stacked", "Without D0", bin, xL, xR);

    for (int iEntry = 0; dst->LoadTree(iEntry) >= 0; ++iEntry)
    {
      dst->GetEntry(iEntry);
      if (m2 >= -0.02 && m2 <= 0.02 && m3 >= 0.8 && m4 >= 2.5){
        h0st->Fill(m1);
      } else if (m2 >= 0.04 && m2 <= 0.25 && m3 >= 0.8 && m4 >= 2.5){
        h1st->Fill(m1);
      }
    }
  
    h0st->SetFillColor(38);
    h0st->SetFillStyle(1001);
    // h0st->GetYaxis()->SetRangeUser(0,1000);

    h1st->Scale(0.27787*Scaling(bin,xL,xR));
    h1st->SetFillColor(46);
    h1st->SetFillStyle(3001);
    // h1st->GetYaxis()->SetRangeUser(0,1000);

    hs->Add(h0st,"");
    hs->Add(h1st,"");

    hs->Draw("nostack,HIST");
    // gPad->SetLogy();
    //=========================================================================================================================================================================//
    // A d d   t e x t   t o   f r a m e
    // ---------------------------------------------------------------------------
    TLatex *latex = new TLatex();
    TString lum("");
    lum.Form("#splitline{Generic c#bar{c} Events}{#scale[0.5]{#int}#it{L} dt=107.33 fb^{-1}}");
    latex->SetNDC();
    latex->SetTextSize(0.03);
    latex->DrawLatex(0.7,0.55, lum);

    auto PerBin = ((xR - xL) / (bin)) * 100.;
    TString str("");
    str.Form("Entries/(%.2f MeV/c^{2})", PerBin);
    hs->SetTitle("");
    hs->GetXaxis()->SetTitle("#Deltam [GeV/c^{2}]");
    hs->GetXaxis()->CenterTitle(true);
    hs->GetYaxis()->SetTitle(str);

    c->Modified();
    //=========================================================================================================================================================================//
    // L e g e n d
    // ---------------------------------------------------------------------------
    // gPad->BuildLegend();
    TLegend *leg = new TLegend(0.78,0.695,0.980,0.935);
    leg->AddEntry(h0st,"Inside D^{0} Signal Region","f");
    leg->AddEntry(h1st,"Outside D^{0} Signal Region","f");
    leg->Draw();
    //=========================================================================================================================================================================//
    // A r e a s
    // ---------------------------------------------------------------------------
    double Area1 = h0st->Integral();
    double Area2 = h1st->Integral();

    cout << "With D0: " << "Area1 = " << Area1 << endl;
    cout << "Without D0: " << "Area2 = " << Area2 << endl;
    //=========================================================================================================================================================================//
    // S a v e   P i c t u r e   o f  P l o t
    c->SaveAs("C:/Users/Adilm/Desktop/Update/Next_Meeting/D0_dM_Scaled_Comparison_ccbar.png");
}