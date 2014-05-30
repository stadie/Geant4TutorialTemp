void g2()
{  
  // identifier of tracking media
  Int_t indVac = gGeoManager->GetMedium("Vacuum")->GetId();
  Int_t indPb  = gGeoManager->GetMedium("Pb")->GetId();
  Int_t indAr  = gGeoManager->GetMedium("ArgonGas")->GetId();

  //half lengths!!! Defines a 50 cm x 50 cm x 50 cm  box
  Double_t dx = 50/2;
  Double_t dr = 50/2;
  Double_t dz = 50/2;

  //number of divisions
  Int_t nx = 25;
  Int_t nz = 25;

  Double_t offsetX     = 0;
  Double_t exphXYZ[3]  = {2 * dx + offsetX,dr,dr};
  Double_t caloXYZ[3]  = {dx,dr,dr};
  Double_t layerXYZ[3] = {dx/nx,dr,dr};
  
  Double_t* ubuf(0);
  
  // experimental hall, "world volume"
  TGeoVolume* top = gGeoManager->Volume("EXPH","BOX",indVac,exphXYZ,3);
  gGeoManager->SetTopVolume(top);
  
  // calorimeter block 
  gGeoManager->Volume("CALB","BOX",indPb,caloXYZ,3);
  gGeoManager->Node("CALB",1,"EXPH",0.0,0.0,0.0,0,kTRUE,ubuf);
  
  // calorimeter cells
  gGeoManager->Volume("LAYB","BOX",indPb,layerXYZ,3);
  for (Int_t ix=0;ix<nx;ix++)
    for (Int_t iz=0;iz<nz;iz++)
      gGeoManager->Node("LAYB",iz*nx+ix,"CALB",
		      (ix-(nx-1)/2)*2*dx/nx,
		      0.,
		      (iz-(nz-1)/2)*2*dz/nz,
		      0,kTRUE,ubuf);
  
  gGeoManager->CloseGeometry();
  gMC->SetRootGeometry();

  /*
    TGeomWrapper wrapper;
    
    //material, other choices:wrapper.ArGas(),Pb(),Fe(),U(),Scint(),Vacuum()
    //                                 2       3    4   5     6       1
    Int_t mat = wrapper.Pb();
    
    //half lengths!!! Defines a 50cm x 50cm x 50cm  box
    Double_t dx = 50/2;
    Double_t dr = 50/2;
    
    //number of divisions
    Int_t nx = 25;
    Int_t nz = 25;
    
    Double_t offsetX = 0;
    Double_t block[3] ={2 * dx + offsetX/2,dr,dr};
    Double_t calo[3] = {dx,dr,dr};
    wrapper.Gsvolu("EXPH","BOX",wrapper.Vacuum(),block,3);
    wrapper.Gsvolu("TMP1","BOX",mat,calo,3);
    wrapper.Gsdvn("TMP2","TMP1",nx,1);
    wrapper.Gsdvn("CALO","TMP2",nz,3);
    wrapper.Gspos("TMP1",1,"EXPH",dx + offsetX,0.,0.,0,"ONLY");
  */
}
