void CompareSB()
{
    TCanvas *c = new TCanvas("fit", "fit", 650, 700);
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.25);
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
    auto hs = new THStack("hs", "");
    //=========================================================================================================================================================================//
    // G r a b   N e e d e d   F i l e   a n d   T r e e
    // ---------------------------------------------------------------------------
    // Simulated Data
    TFile *file1 = new TFile("C:/Users/Adilm/Desktop/High_Energy_Physics/Analysis_Code/KEKCC_File/Grid/Ds2D0e-Generic_Ds_052324_ccbar.root");
    TTree *sim = (TTree *)file1->Get("dstree");
    // Generic Events
    TFile *file2 = new TFile("C:/Users/Adilm/Desktop/High_Energy_Physics/Analysis_Code/KEKCC_File/Simulated/Ds2D0enu-Signal.root");
    TTree *gen = (TTree *)file2->Get("dstree");

    double m1, m2, m3, m4;
    sim->SetBranchAddress("Ds_daughterMotherDiffOf_0_M", &m1);
    sim->SetBranchAddress("dM_D0pi", &m2);
    gen->SetBranchAddress("Ds_daughterMotherDiffOf_0_M", &m3);
    gen->SetBranchAddress("dM_D0pi", &m4);
    //=========================================================================================================================================================================//
    // S e t u p   H i s t o g r a m
    // ---------------------------------------------------------------------------
    int bin = 50;
    double xL = 0.135;
    double xR = 0.25;
    auto h0st = new TH1F("Stacked", "", bin, xL, xR);
    auto h1st = new TH1F("Stacked", "", bin, xL, xR);

    for (int iEntry = 0; sim->LoadTree(iEntry) >= 0; ++iEntry)
    {
        sim->GetEntry(iEntry);
        h0st->Fill(m1);
    }

    for (int iEntry = 0; gen->LoadTree(iEntry) >= 0; ++iEntry)
    {
        gen->GetEntry(iEntry);
        h1st->Fill(m4);
    }

    h0st->SetFillColor(38);
    h0st->SetFillStyle(1001);

    h1st->SetFillColor(46);
    h1st->SetFillStyle(3001);

    hs->Add(h0st, "");
    hs->Add(h1st, "");

    hs->Draw("nostack");
}