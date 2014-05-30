void dEdx()
{
  const Int_t nev = 10000;
  const Double_t density = 8.960;//g cm3
  const Double_t mass = 0.1057;
  const Double_t length = 0.2;
  app->InitMC("geometry/cubox");
  app->SetPrimaryPDG(-13);

  TH1F* hloss = new TH1F("hloss","; -dE [MeV]",100,1,5);
  std::vector<double> Eloss(nev);

  TGraph* gdEdx =  new TGraph(100);
  Int_t i = 0;
  TFolder* histfolder=(TFolder*)gROOT->FindObjectAny("/Geant4/Histograms");
  TProfile* hprim=(TProfile*)gROOT->FindObjectAny("/Geant4/Histograms/hPrimaryEnergy"); 
  TCanvas* c = new TCanvas("c");
  double momentum = 1;
  for(int i = 0 ; i < nev ; ++i) {
    Double_t beta = momentum/sqrt(momentum*momentum + mass*mass);
    Double_t gamma = 1/sqrt(1-beta*beta);
    app->SetPrimaryMomentum(momentum);
    hprim->Reset();
    app->RunMC(1,!i);
    double Eafter = hprim->GetBinContent(51);
    std::cout <<  sqrt(mass*mass+momentum*momentum) << ", " << Eafter << '\n';
    double loss = sqrt(mass*mass+momentum*momentum)-Eafter;
    loss *= 1000; //MeV
    Eloss[i] = loss;
    hloss->Fill(loss);
    loss = loss / density /length; //MeV g-1 cm2
    //Eloss[i] = loss;
    std::cout << "betagamma:" << momentum/mass << "     -dE/dex:" << loss << '\n';
    //gdEdx->SetPoint(i,momentum/mass,-loss);
  }  
  TCanvas* c1 = new TCanvas("c1");
  hloss->Draw();
  TCanvas* c2 = new TCanvas("c2");
  gdEdx->Set(1);
  c2->Clear();
  c2->cd();
  c2->SetLogx();
  gdEdx->Draw("A*");
  gdEdx->GetXaxis()->SetTitle("#beta#gamma");
  gdEdx->GetYaxis()->SetTitle("-dE/dx [MeV g^{-1} cm^{2}]");
  c2->Modified();
  c2->Update();
}
