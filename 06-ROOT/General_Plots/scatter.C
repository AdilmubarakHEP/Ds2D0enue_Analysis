void scatter()
{
  auto c1 = new TCanvas("c1","Profile histogram example",200,10,700,500);

  TProfile* m_Profile[2] = {nullptr};
  m_Profile[1] = new TProfile("hprof","Profile of pz versus px",100,-4,4,0,20);
  Float_t px, py, pz;
  for ( Int_t i=0; i<25000; i++) {
    gRandom->Rannor(px,py);
    pz = px*px + py*py;
    m_Profile[1]->Fill(px,pz,1);
  }

  m_Profile[1]->GetXaxis()->SetRangeUser(0,2);

  m_Profile[1]->Draw();

  cout << m_Profile[1]->GetEntries() << endl;

  c1->SaveAs("test.png");
}