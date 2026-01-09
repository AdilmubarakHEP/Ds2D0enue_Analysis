// The line below is used to split the histogram between isSignal=1 and isSignal=0.
void hstack()
{
    // The everything before the breaker sets everything up.
    // gStyle->SetTimeOffset(0); here
    TCanvas *c1 = new TCanvas();

    auto hs = new THStack("hs","D_s+ Mass Distribution");

	// The line below initializes everything.
    int bin = 50;
    double xL = 0.5;
    double xR = 3;

    TTree* dst1 = _file0->Get<TTree>("dstree");
    double m1, m2, m3, m4, m5;
    dst1->SetBranchAddress("Ds_M", &m1);
    dst1->SetBranchAddress("D0_mcPDG", &m2);
    dst1->SetBranchAddress("e_mcPDG", &m3);
    dst1->SetBranchAddress("pi_mcPDG", &m4);
    dst1->SetBranchAddress("K_mcPDG", &m5);

    gStyle->SetOptStat(0);
//--------------------------------------------------------------------------------------------------//
    auto h1st = new TH1F("Stats", "isSignal=0", bin, xL, xR);
    auto h2st = new TH1F("Stats", "isSignal=1", bin, xL, xR);

    for (int iEntry = 0; dst1->LoadTree(iEntry) >= 0; ++iEntry)
    {
        dst1->GetEntry(iEntry);

        if (m2 == 0) {
            h1st->Fill(m1);
        } else if (m2 == 1) {
            h2st->Fill(m1);
        } else{
            h3st->Fill(m1);
        }
    }
    // h1st->SetFillColor(38);
    // h1st->SetMarkerStyle(21);
    // h1st->SetMarkerColor(kBlue);

    h1st->SetFillColor(38);
    h1st->SetFillStyle(1001);

    h2st->SetFillColor(46);
    h2st->SetFillStyle(1001);

    h3st->SetFillColor(14);
    h3st->SetFillStyle(1001);

    hs->Add(h1st,"");
    hs->Add(h2st,"");
    // hs->Add(h3st,"");
//--------------------------------------------------------------------------------------------------//
    // This lower section sets up the plot and gives it its labels.
    hs->Draw();
    gPad->BuildLegend(0.78,0.695,0.980,0.935,"","f");
 
    hs->GetXaxis()->SetTitle("mass [GeV]");
    hs->GetYaxis()->SetTitle("Entry");
 
    //c1->Modified();
    //gPad->BuildLegend(0.78,0.695,0.980,0.935,"","f");
}


// The code below combines files from the KEKCC server and from the grid. It can be utilized to combine 
// two different files and stack the plots.
// void hstack()
// {
//     // The everything before the breaker sets everything up.
//     TCanvas *c1 = new TCanvas();

//     auto hs = new THStack("hs","Mass Distribution of D_s+");

// 	// The line below initializes everything.
// 	double xL = 1.8;
// 	double xR = 2.05;
//     int bin = 50;

//     TTree* dst0 = _file0->Get<TTree>("dstree");
//     TTree* dst1 = _file1->Get<TTree>("dstree");

//     double m0, m1;
//     dst0->SetBranchAddress("Ds_M", &m0);
//     dst1->SetBranchAddress("Ds_M", &m1);

//     gStyle->SetOptStat(0);
// //--------------------------------------------------------------------------------------------------//
//     auto h0st = new TH1F("Stats", "generic events", bin, xL, xR);
//     for (int iEntry = 0; dst0->LoadTree(iEntry) >= 0; ++iEntry)
//     {
//         dst0->GetEntry(iEntry);
//         h0st->Fill(m0);
//     }
//     // h2st->SetFillColor(46);
//     // h2st->SetMarkerStyle(21);
//     // h2st->SetMarkerColor(kRed);

//     h0st->SetFillColor(46);
//     h0st->SetFillStyle(1001);

//     hs->Add(h0st,"");
// //--------------------------------------------------------------------------------------------------//
//     auto h1st = new TH1F("Stats", "signal events", bin, xL, xR);
//     for (int iEntry = 0; dst1->LoadTree(iEntry) >= 0; ++iEntry)
//     {
//         dst1->GetEntry(iEntry);
//         h1st->Fill(m1);
//     }
//     // h1st->SetFillColor(38);
//     // h1st->SetMarkerStyle(21);
//     // h1st->SetMarkerColor(kBlue);

//     h1st->SetFillColor(38);
//     h1st->SetFillStyle(3002);

//     hs->Add(h1st,"");
// //--------------------------------------------------------------------------------------------------//
//     hs->Draw();
    
//     hs->GetXaxis()->SetTitle("mass [GeV]");
//     hs->GetYaxis()->SetTitle("Entry");
 
//     c1->Modified();

//     gPad->BuildLegend(0.78,0.695,0.980,0.935,"","f");
// }