//
// Script to define a simple geometry, a lead box, for TGeantApplication
//



void buildLayer(Double_t x) 
{
  static int nstrips = 0;
  static int nstructs = 0;
  const TGeoMedium *si  = gGeoManager->GetMedium("Si");  
  const double thickness = 0.0400;
  const double pitch     = 0.0150;
  const double length    = 50.0;

  Double_t* ubuf(0);

  TGeoVolume* top = gGeoManager->GetVolume("EXPH");
  TGeoVolume* struc = gGeoManager->GetVolume("Structure");
  
  top->AddNode(struc,nstructs++,new TGeoTranslation(x-0.16-pitch/2,0,0));
  double minz = -(x+50.0) * 0.5;
  double maxz = -minz;
  for(double z = minz - pitch/2 ; z < maxz ; z += pitch) {
    stringstream ss1; 
    ss1<<"Strip" << nstructs << "_" << nstrips++;
    TGeoVolume* v = gGeoManager->MakeBox(ss1.str().c_str(),si,thickness/2,length,pitch/2);
    top->AddNode(v,1,new TGeoTranslation(x,0,z));
  }
  top->AddNode(struc,nstructs++,new TGeoTranslation(x+pitch/2+0.16,0,0));
}

void tracker()
{
  TVirtualMagField *B = new TGeoUniformMagField(0.0,20,0);
  gMC->SetMagField(B);
 
  //
  // TRACKING MEDIA
  //
  Int_t    ifield =     2;  // magnetic field
  Double_t fieldm =    20;  //
  Double_t epsil  =  .001;  // Tracking precision,
  Double_t stemax = -0.01;  // Maximum displacement for multiple scat
  Double_t tmaxfd =  -90.;  // Maximum angle due to field deflection
  Double_t deemax =   -.3;  // Maximum fractional energy loss, DLS
  Double_t stmin  =   -.0001;
  
  gGeoManager->Medium("Vacuum2", 10, 1,
		      0, ifield, fieldm, tmaxfd, stemax, deemax, epsil, stmin); 
  epsil  =  .0001;
  stemax = -0.1;
  Double_t a = 28.0855;
  Double_t z = 14;
  Double_t density = 2.33;
  Double_t radl = 9.36;
  Double_t absl = 0;
  gGeoManager->Material("Si", a, z, density, 11, radl, absl);
  
  gGeoManager->Medium("Si",11, 11,
		      0, ifield, fieldm, tmaxfd, stemax, deemax, epsil, stmin); 
  gGeoManager->Medium("Fe2",12, 4,
		      0, ifield, fieldm, tmaxfd, stemax, deemax, epsil, stmin); 

  // identifier of tracking media
  Int_t indVac = gGeoManager->GetMedium("Vacuum2")->GetId();
  Int_t indSi  = gGeoManager->GetMedium("Si")->GetId();
  Int_t indFe  = gGeoManager->GetMedium("Fe2")->GetId();

  //half lengths!!! Defines a 22 cm x 20 cm x 20 cm box
  Double_t exphXYZ[3] = { 51.,50.,50. };
  Double_t structureXYZ[3] = { 0.3/2,50.,50. };
  

  Double_t* ubuf(0);
  
  // experimental hall, "world volume"
  TGeoVolume* top = gGeoManager->Volume("EXPH","BOX",indVac,exphXYZ,3);
  gGeoManager->SetTopVolume(top);
  
  // box made of Pb 
  gGeoManager->Volume("Structure","BOX",indFe,structureXYZ,3); 

  buildLayer(-45.0);
  buildLayer(-30.0);
  buildLayer(45.0);
}
    
