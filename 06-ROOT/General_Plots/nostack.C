// The line below unstacks the two histograms.
void nostack()
{
    // The everything before the breaker sets everything up.
    // gStyle->SetTimeOffset(0); here
    TCanvas *c1 = new TCanvas();

    auto hs = new THStack("hs","");

	// The line below initializes everything.
    int bin = 50;
    double xL = 0;
    double xR = 0.26;

    TFile *file = new TFile("C:/Users/Adilm/Desktop/Particle/Analysis Code/KEKCC_File/Grid/Ds2D0e-Reconstruction_Ds_052024_ccbar.root");
    TTree *dst1 = (TTree *)file->Get("dstree");
    double m1, m2, m3, m4, m5, m6, m7, m8, m9;
    dst1->SetBranchAddress("Ds_daughterMotherDiffOf_0_M", &m1);
    dst1->SetBranchAddress("D0_mcPDG", &m2);
    dst1->SetBranchAddress("e_mcPDG", &m3);
    dst1->SetBranchAddress("pi_mcPDG", &m4);
    dst1->SetBranchAddress("K_mcPDG", &m5);

    dst1->SetBranchAddress("D0_genMotherPDG", &m6);
    dst1->SetBranchAddress("e_genMotherPDG", &m7);
    dst1->SetBranchAddress("pi_genMotherPDG", &m8);
    dst1->SetBranchAddress("K_genMotherPDG", &m9);

    gStyle->SetOptStat(0);
//--------------------------------------------------------------------------------------------------//
    auto h0st = new TH1F("Stats", "generic ccbar events", bin, xL, xR);
    auto h1st = new TH1F("Stats", "D^{*+} #rightarrow [D^{0} #rightarrow K^{-} #pi^{+}] #pi^{+} ", bin, xL, xR);

    for (int iEntry = 0; dst1->LoadTree(iEntry) >= 0; ++iEntry)
    {
        dst1->GetEntry(iEntry);
        h0st->Fill(m1);

        if (abs(m2) == 421 && abs(m3) == 211 && abs(m4) == 211 && abs(m5) == 321 && abs(m6) == 413 && abs(m7) == 413 && abs(m8) == 421 && abs(m9) == 421) {
            h1st->Fill(m1);
        }

        // if (abs(m2) != 421 && abs(m3) != 211 && abs(m4) != 211 && abs(m5) != 321 && abs(m6) != 413 && abs(m7) != 413 && abs(m8) != 421 && abs(m9) != 421) {
        //     h1st->Fill(m1);
        // }

    }
    h0st->SetFillColor(38);
    h0st->SetFillStyle(1001);
    h0st->SetMarkerStyle(21);

    h1st->SetFillColor(46);
    h1st->SetFillStyle(3001);
    h1st->SetMarkerStyle(21);

    hs->Add(h0st,"");
    hs->Add(h1st,"");
//--------------------------------------------------------------------------------------------------//
    hs->Draw("nostack");
    
    TLatex *latex = new TLatex();
    TString lum("");
    lum.Form("#splitline{Generic c#bar{c} Events}{#scale[0.5]{#int}#it{L} dt=107.33 fb^{-1}}");
    latex->SetNDC();
    latex->SetTextSize(0.03);
    latex->DrawLatex(0.7,0.55, lum);

    auto PerBin = ((xR-xL)/(bin))*100.;
    TString str ("");
    str.Form("Entries/(%.2f MeV/c^{2})", PerBin);
    hs->SetTitle("");
    hs->GetXaxis()->SetTitle("m(D_{s}^{+} - D^{0}) [GeV/c^{2}]");
    hs->GetXaxis()->CenterTitle(true);
    hs->GetYaxis()->SetTitle(str);
    hs->GetYaxis()->SetRangeUser(0,1400);
 
    c1->Modified();
    // c1.Update();
    
    // gPad->BuildLegend();
    TLegend *leg = new TLegend(0.78,0.695,0.980,0.935);
    leg->AddEntry(h0st,"With D^{0}","l");
    leg->AddEntry(h1st,"Without D^{0}","l");
    leg->Draw();

}


// // The code below combines files from the KEKCC server and from the grid. It can be utilized to combine 
// // two different files.
// void nostack()
// {
//     // The everything before the breaker sets everything up.
//     // gStyle->SetTimeOffset(0); here
//     TCanvas *c1 = new TCanvas();

//     auto hs = new THStack("hs","Mass Distribution of D_s+");

// 	// The line below initializes everything.
//     int bin = 50;
//     double xL = 0;
//     double xR = 0.26;

//     TTree* dst0 = _file0->Get<TTree>("dstree");
//     TTree* dst1 = _file1->Get<TTree>("dstree");
//     TTree* dst2 = _file2->Get<TTree>("dstree");

//     double m0, m1, m2;
//     dst0->SetBranchAddress("Ds_daughterMotherDiffOf_0_M", &m0);
//     dst1->SetBranchAddress("Ds_daughterMotherDiffOf_0_M", &m1);
//     dst2->SetBranchAddress("Ds_daughterMotherDiffOf_0_M", &m2);

//     gStyle->SetOptStat(0);
// //--------------------------------------------------------------------------------------------------//
//     auto h0st = new TH1F("Stats", "Generic c#bar{c} Events (All)", bin, xL, xR);
//     for (int iEntry = 0; dst0->LoadTree(iEntry) >= 0; ++iEntry)
//     {
//         dst0->GetEntry(iEntry);
//         h0st->Fill(m0);
//     }
//     h0st->SetFillColor(46);
//     h0st->SetFillStyle(1001);
//     h0st->GetYaxis()->SetRangeUser(0,1400);

//     hs->Add(h0st,"");
// //--------------------------------------------------------------------------------------------------//
//     auto h1st = new TH1F("Stats", "Generic c#bar{c} Events !(D^{*+} -> [D^{0}->K^{-}#pi^{+}] #pi^{+} )  ", bin, xL, xR);
//     for (int iEntry = 0; dst1->LoadTree(iEntry) >= 0; ++iEntry)
//     {
//         dst1->GetEntry(iEntry);
//         h1st->Fill(m1);
//     }
//     h1st->SetFillColor(38);
//     h1st->SetFillStyle(3002);
//     h1st->GetYaxis()->SetRangeUser(0,1400);

//     hs->Add(h1st,"");
// //--------------------------------------------------------------------------------------------------//
//     auto h2st = new TH1F("Stats", "Generic c#bar{c} Events !(D^{*+})  ", bin, xL, xR);
//     for (int iEntry = 0; dst2->LoadTree(iEntry) >= 0; ++iEntry)
//     {
//         dst2->GetEntry(iEntry);
//         h2st->Fill(m2);
//     }
//     h2st->SetFillColor(29);
//     h2st->SetFillStyle(1001);
//     h2st->SetFillStyle(3001);
//     h2st->GetYaxis()->SetRangeUser(0,1400);

//     hs->Add(h2st,"");
// //--------------------------------------------------------------------------------------------------//
//     hs->Draw("nostack");
    
//     auto PerBin = ((xR-xL)/(bin))*100.;
//     TString str ("");
//     str.Form("Entries/(%.2f MeV/c^{2})", PerBin);
//     hs->SetTitle("");
//     hs->GetXaxis()->SetTitle("m(D_{s}^{+} - D^{0}) [GeV/c^{2}]");
//     hs->GetXaxis()->CenterTitle(true);
//     hs->GetYaxis()->SetTitle(str);
 
//     c1->Modified();
//     hs->GetYaxis()->SetRangeUser(0,1400);

//     gPad->BuildLegend(0.78,0.695,0.980,0.935,"","f");
// }