void dEdx_solve()
{
  const Int_t nev = 400;
  const Double_t density = 8.960;//g cm3
  const Double_t mass = 0.1057;
  const Double_t length = 0.2;
  app->InitMC("geometry/cubox");
  app->SetPrimaryPDG(13);

  TH1F* hloss;
  std::vector<double> Eloss(nev); 
  TGraphErrors* gdEdx =  new TGraphErrors(1000);
  TCanvas* c2 = new TCanvas("c2");
  //gdEdx->Set(1);
  c2->Clear();
  c2->cd();
  TH1F *h1= new TH1F("h1",";#beta#gamma;-dE/dx [MeV g^{-1} cm^{2}]",100,pow(10,-3),pow(10,8));
  h1->Draw();
  h1->SetMaximum(100);
  h1->SetMinimum(0.1);
  c2->SetLogx();
  c2->SetLogy();
  gdEdx->Draw("P");
  Int_t gdEdx_i = 0;
  TFolder* histfolder=(TFolder*)gROOT->FindObjectAny("/Geant4/Histograms");
  TProfile* hprim=(TProfile*)gROOT->FindObjectAny("/Geant4/Histograms/hPrimaryEnergy"); 
  for(double x = -3; x < 7 ; x+=0.2) {
    double momentum = pow(10,x);
    Eloss.clear();
    double minDE=100000000;
    double maxDE=0;
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
      Eloss.push_back(loss);
      if(loss < minDE) minDE=loss;
      if(loss > maxDE) maxDE=loss;
      if(momentum == 1) hloss->Fill(loss);
      loss = loss / density /length; //MeV g-1 cm2
      //Eloss[i] = loss;
    }
    delete hloss;
    if(minDE < 0.1) minDE = 0.1;
    hloss = new TH1F("hloss","; -dE [MeV]",10000,minDE,maxDE);
    for(int j = 0 ; j < Eloss.size() ; ++j) {
	hloss->Fill(Eloss[j]);
    }
    double dEdx = hloss->GetMean() / density /length; //MeV g-1 cm2	  
    std::cout << "betagamma:" << momentum/mass << "     -dE/dex:" << dEdx << '\n';
    gdEdx->SetPoint(gdEdx_i,momentum/mass,dEdx);
    gdEdx->SetPointError(gdEdx_i,0,hloss->GetMeanError() / density /length);
    c2->Modified();
    c2->Update();
    ++gdEdx_i;
  }
  TCanvas* c1 = new TCanvas("c1");
  hloss->Draw();  
 
}
