void Multi_Plot()
{
    TTree* dst1 = _file0->Get<TTree>("KMCtree");
    TTree* dst2 = _file0->Get<TTree>("PiMCtree");
    TTree* dst3 = _file0->Get<TTree>("EMCtree");

    auto c1 = new TCanvas("c1", "Transverse Momentum", 600, 600);
    c1->Divide(2, 2);

    double m1;
    double m2;
    double m3;
    //double m4;

    dst1->SetBranchAddress("pt", &m1);
    dst2->SetBranchAddress("pt", &m2);
    dst3->SetBranchAddress("pt", &m3);
    //dst1->SetBranchAddress("K_theta", &m4);

    TH1F* draw1 = new TH1F("h1", "Kaon Transverse Momentum (Generator Level: D_s+ -> [D0 -> ^K- pi+] e+ nu_e)", 100, 0, 5);
    TH1F* draw2 = new TH1F("h2", "Pion Transverse Momentum (Generator Level: D_s+ -> [D0 -> K- ^pi+] e+ nu_e)", 100, 0, 5);
    TH1F* draw3 = new TH1F("h3", "Electron Transverse Momentum (Generator Level: D_s+ -> [D0 -> K- pi+] ^e+ nu_e)", 100, 0, 0.40);
    //TH1F* draw4 = new TH1F("h4", "Kaon Reconsturction Polar Angle", 100, 0, 4);

    for (int i = 0; dst1->LoadTree(i) >= 0; i++)
    {
        dst1->GetEntry(i);
        draw1->Fill(m1);
    }

    for (int i = 0; dst2->LoadTree(i) >= 0; i++)
    {
        dst2->GetEntry(i);
        draw2->Fill(m2);
    }

    for (int i = 0; dst3->LoadTree(i) >= 0; i++)
    {
        dst3->GetEntry(i);
        draw3->Fill(m3);
    }

    //for (int i = 0; dst2->LoadTree(i) >= 0; i++)
    //{
    //    dst2->GetEntry(i);
    //    draw2->Fill(m2);

    //    dst2->GetEntry(i);
    //    draw4->Fill(m4);
    //}

    c1->cd(1);
    draw1->GetXaxis()->SetTitle("momentum [GeV/c]");
    draw1->GetYaxis()->SetTitle("Events");
    draw1->Draw();

    c1->cd(2);
    draw2->GetXaxis()->SetTitle("momentum [GeV/c]");
    draw2->GetYaxis()->SetTitle("Events");
    draw2->Draw();

    c1->cd(3);
    draw3->GetXaxis()->SetTitle("momentum [GeV/c]");
    draw3->GetYaxis()->SetTitle("Events");
    draw3->Draw();

    //c1->cd(4);
    //draw4->GetXaxis()->SetTitle("Radians");
    //draw4->GetYaxis()->SetTitle("Events");
    //draw4->Draw();

}

//void Multi_Plot()
//{
//    TTree* dst1 = _file0->Get<TTree>("KMCtree");
//    TTree* dst2 = _file0->Get<TTree>("PiMCtree");
//    TTree* dst3 = _file0->Get<TTree>("EMCtree");
//
//    auto c1 = new TCanvas("c1", "Mass", 600, 600);
//    c1->Divide(2, 2);
//
//    double m1;
//    double m2;
//    double m3;
//    //double m4;
//
//    dst1->SetBranchAddress("useCMSFrame__bop__bc", &m1);
//    dst2->SetBranchAddress("useCMSFrame__bop__bc", &m2);
//    dst3->SetBranchAddress("useCMSFrame__bop__bc", &m3);
//    //dst1->SetBranchAddress("K_theta", &m4);
//
//    TH1F* draw1 = new TH1F("h1", "Kaon CMS Momentum (Generator Level: D+ -> [D0 -> ^K- pi+] e+ nu_e)", 100, 0, 5);
//    TH1F* draw2 = new TH1F("h2", "Pion CMS Momentum (Generator Level: D+ -> [D0 -> K- ^pi+] e+ nu_e)", 100, 0, 5);
//    TH1F* draw3 = new TH1F("h3", "Electron CMS Momentum (Generator Level: D+ -> [D0 -> K- pi+] ^e+ nu_e)", 100, 0, 0.025);
//    //TH1F* draw4 = new TH1F("h4", "Kaon Reconsturction Polar Angle", 100, 0, 4);
//
//    for (int i = 0; dst1->LoadTree(i) >= 0; i++)
//    {
//        dst1->GetEntry(i);
//        draw1->Fill(m1);
//    }
//
//    for (int i = 0; dst2->LoadTree(i) >= 0; i++)
//    {
//        dst2->GetEntry(i);
//        draw2->Fill(m2);
//    }
//
//    for (int i = 0; dst3->LoadTree(i) >= 0; i++)
//    {
//        dst3->GetEntry(i);
//        draw3->Fill(m3);
//    }
//
//    //for (int i = 0; dst2->LoadTree(i) >= 0; i++)
//    //{
//    //    dst2->GetEntry(i);
//    //    draw2->Fill(m2);
//
//    //    dst2->GetEntry(i);
//    //    draw4->Fill(m4);
//    //}
//
//    c1->cd(1);
//    draw1->GetXaxis()->SetTitle("momentum [GeV/c]");
//    draw1->GetYaxis()->SetTitle("Events");
//    draw1->Draw();
//
//    c1->cd(2);
//    draw2->GetXaxis()->SetTitle("momentum [GeV/c]");
//    draw2->GetYaxis()->SetTitle("Events");
//    draw2->Draw();
//
//    c1->cd(3);
//    draw3->GetXaxis()->SetTitle("momentum [GeV/c]");
//    draw3->GetYaxis()->SetTitle("Events");
//    draw3->Draw();
//
//    //c1->cd(4);
//    //draw4->GetXaxis()->SetTitle("Radians");
//    //draw4->GetYaxis()->SetTitle("Events");
//    //draw4->Draw();
//
//}