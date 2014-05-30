void scint()
{  
  // identifier of tracking media
  Int_t indVac = gGeoManager->GetMedium("Vacuum")->GetId();
  Int_t indPb  = gGeoManager->GetMedium("Pb")->GetId();
  Int_t indFe  = gGeoManager->GetMedium("Fe")->GetId();
  Int_t indAr  = gGeoManager->GetMedium("ArgonGas")->GetId();
  Int_t indScint = gGeoManager->GetMedium("Scintillator")->GetId();

  //half lengths!!! Defines a 10 cm x 4 cm x 4 cm  box
  Double_t dx = 50/2;
  Double_t dr = 10/2;

  //number of divisions
  Int_t nx = 1;
  Int_t nz = 1;

  Double_t offsetX     = 0;
  Double_t exphXYZ[3]  = {2 * dx + offsetX,dr,dr};
  Double_t caloXYZ[3]  = {dx,dr,dr};
  Double_t layerXYZ[3] = {dx/nx,dr,dr};
  
  Double_t* ubuf(0);
  
  // experimental hall, "world volume"
  TGeoVolume* top = gGeoManager->Volume("EXPH","BOX",indVac,exphXYZ,3);
  gGeoManager->SetTopVolume(top);
  
  // calorimeter block 
  gGeoManager->Volume("CALB","BOX",indScint,caloXYZ,3);
  gGeoManager->Node("CALB",1,"EXPH",0.0,0.0,0.0,0,kTRUE,ubuf);
  
  // calorimeter cells
  gGeoManager->Volume("LAYB","BOX",indScint,layerXYZ,3);
  for (Int_t ix=0;ix<nx;ix++)
    gGeoManager->Node("LAYB",ix,"CALB",
		      (ix-nx/2)*2*dx/nx,0.,0.,0,kTRUE,ubuf);
  
  gGeoManager->CloseGeometry();
  gMC->SetRootGeometry();
}
