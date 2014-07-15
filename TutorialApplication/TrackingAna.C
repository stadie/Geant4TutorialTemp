void TrackingAna() {
  gROOT->ProcessLine(".L trackingSolv.C+O");
  
  TutorialApplication* app = (TutorialApplication*)TutorialApplication::Instance();

  // initialize tracker geometry: position of silicon layers, pitch, magnetic field, additional material
  double layerpos[] = {-45,-30,45};//[cm]position of 1st, 2nd and 3rd layer
  double pitch = 0.0150;//[cm] pitch of Si strip
  double materialLength= 0.6;//[cm] lenght of iron as additional material per layer
  double bfield = 2;//[T] strength of magentic field
  
  TString geom("geometry/tracker2(");
  geom+=layerpos[0]; geom.Append(",");
  geom+=layerpos[1]; geom.Append(",");
  geom+=layerpos[2]; geom.Append(",");
  geom+=pitch; geom.Append(",");
  geom+=materialLength; geom.Append(",");
  geom+=bfield; geom.Append(")");
  app->InitMC(geom);

  //book histograms for hit residuals

  //book histograms for momentum measurement analysis
  TH1F *hpt[6];
  TH1F *hptpull = new TH1F("hptpull","; (p_{T}^{meas} - p_{T}^{true})/#sigma",100,-10,10);
 


  // define particle and control parameters of loop   
  double pts[] = {1,2,4,7,10,20};
  // generate  some events
  TObjArray* clusters = new TObjArray();
  clusters->SetOwner(true);
  
  
  for(int ipt = 0 ; ipt < 6 ;  ++ipt) {
    unsigned int nevt = 200;
    TVector3 p;
    double pt = pts[ipt];//[GeV]
    hpt[ipt] = new TH1F("hpt","; p_{T} [GeV]",25,pt-.15*pt,pt+0.15*pt);
    //app->SetPrimaryPDG(13);    // +/-11: PDG code of e+/- 
    /* other PDG codes     22: Photon    +-13: muon   
       +/-211: pion   +/-2212: proton     */
  
    for(unsigned int i=0;i<nevt;++i) {
      bool draw = !(i%((nevt/20) ? (nevt/20) : 1));
      
      double theta = TMath::Pi()/2+gRandom->Uniform(-.15,0.15);//phi
      double z = gRandom->Uniform(-0.5,0.5);
      p.SetMagThetaPhi(pt,theta,0);
      app->SetPrimaryPDG(gRandom->Uniform(-1,1) < 0 ? -13 : 13);
      app->SetPrimaryMomentum(p);
      app->SetPrimaryVertex(-50,0,z);
      
      removeAllHelices();
      app->RunMC(1, draw);
      updateClusters(clusters);
      reconstructHitsBinary(clusters);
      plotResdiuals(clusters);
      if(clusters->GetEntriesFast() >=3) {
	Track *t = fitTrack(clusters);
	if(t->chi2() < 1) { 
	  if(draw) {t->helix()->Draw(); gPad->Modified(); gPad->Update();}
	  hpt[ipt]->Fill(t->pt());
	  std::cout << " Chi2:" << t->chi2() << " Pt:" << t->pt() << " +- " << t->ptErr() << std::endl;
	  hptpull->Fill((t->pt()-pt)/t->ptErr());
	}
      }
    }
  }
  TCanvas* c = new TCanvas("c");
  c->cd();
  c->Divide(2,3);
  for(i = 0 ; i < 6 ; ++i) {
    c->cd(i+1);
    hpt[i]->Fit("gaus");
  }
    
  //hlayer3->Draw();
  //hresid2->Draw();
  //hpt->Draw();
  //hpt->Fit("gaus");
}

