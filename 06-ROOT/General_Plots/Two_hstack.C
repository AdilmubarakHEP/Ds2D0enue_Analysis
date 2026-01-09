void Two_hstack()
{
    // gStyle->SetTimeOffset(0); here
    TCanvas *c1 = new TCanvas();
    c1->Divide(2, 1);

    int bin = 50;

    auto hs1 = new THStack("hv","D0 Mass Distribution");
    auto hs2 = new THStack("hc","Number of Candidates");

	TTree* dst1 = _file0->Get<TTree>("d0tree");
    double m, s;
    int c;
    dst1->SetBranchAddress("M", &m);
    dst1->SetBranchAddress("__ncandidates__", &c);
    
    dst1->SetBranchAddress("isSignal", &s);

    gStyle->SetOptStat(0);
//--------------------------------------------------------------------------------------------------//
    auto h11st = new TH1F("Stats", "isSignal=0", bin, 1.4, 2.2);
    auto h12st = new TH1F("Stats", "isSignal=1", bin, 1.4, 2.2);

    auto h21st = new TH1F("Stats", "isSignal=0", 10, 0.5, 5);
    auto h22st = new TH1F("Stats", "isSignal=1", 10, 0.5, 5);

for (int iEntry = 0; dst1->LoadTree(iEntry) >= 0; ++iEntry)
    {
        dst1->GetEntry(iEntry);

        if (s == 0) {
            h11st->Fill(m);
            h21st->Fill(c);
        } else if (s == 1) {
            h12st->Fill(m);
            h22st->Fill(c);
        }
    }

    h11st->SetFillColor(38);
    h11st->SetFillStyle(1001);

    h21st->SetFillColor(38);
    h21st->SetFillStyle(1001);

    h12st->SetFillColor(46);
    h12st->SetFillStyle(1001);

    h22st->SetFillColor(46);
    h22st->SetFillStyle(1001);

    hs1->Add(h11st,"");
    hs1->Add(h12st,"");

    hs2->Add(h21st,"");
    hs2->Add(h22st,"");
//--------------------------------------------------------------------------------------------------//
    // This lower section sets up the plot and gives it its labels.
    c1->cd(1);
    hs1->Draw();
    gPad->BuildLegend(0.78,0.695,0.980,0.935,"","f");
    
    c1->cd(2);
    hs2->Draw();
    gPad->BuildLegend(0.78,0.695,0.980,0.935,"","f");
 
    hs1->GetXaxis()->SetTitle("mass [GeV]");
    hs1->GetYaxis()->SetTitle("Entry");

    hs2->GetXaxis()->SetTitle("Candidates");
    hs2->GetYaxis()->SetTitle("Entry");
 
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