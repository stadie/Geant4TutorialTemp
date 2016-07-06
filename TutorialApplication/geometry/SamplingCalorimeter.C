/* Definition of a sampling calorimetr 
    consisting of an absorber (Fe or Pb) and and an active material
    scaled size in terms of average hadronic interaction lengh
*/
void SamplingCalorimeter
(Double_t width1 = 2., Double_t width2 = 1.0, Double_t f=0.7, Int_t IAbs=1)
// width1: width of absorber in cm
// width2: width of active material in cm
// f: size scaling factor (depth = 2*f*lambda_i, width = f*lambda_i)
// IAbs: Absorber type 1: Pb  2: Fe
//      default values: a typical EM calorimeter
{
 // define properties of known media
  // identifier of tracking media
  Int_t indVac = gGeoManager->GetMedium("Vacuum")->GetId();
  Int_t indPb  = gGeoManager->GetMedium("Pb")->GetId();
  Int_t indFe  = gGeoManager->GetMedium("Fe")->GetId();
  Int_t indSci = gGeoManager->GetMedium("Scintillator")->GetId();
  // some hadronic interaction lengths 
  Double_t lambda_i_Pb = 17.59; // Lead
  Double_t lambda_i_Fe = 16.77; // Iron
  Double_t lambda_i_Scint = 78.80; // PVT

  //define material to be used here 
  Double_t lambda_i_Abs=lambda_i_Pb;
  Int_t indAbs=indPb;
  If(IAbs==2){
    lambda_iAbs=lambda_i_Fe;
    indAbs=indFe; 
    }

  Double_t x = width1/lambda_i_Abs + (width2/lambda_i_Scint);
  Double_t comb_lambda_i = (width1+width2)/x; //combined interaction length

  // calculate size of calorimeter in terms of combined interaction length
  Double_t offsetX      = 0.0;
  Double_t exphXYZ[3]   = {comb_lambda_i*2*f,comb_lambda_i*f,comb_lambda_i*f}; 

  Double_t layer1XYZ[3] = {width1 / 2, exphXYZ[1], exphXYZ[2]};
  Double_t layer2XYZ[3] = {width2 / 2, exphXYZ[1], exphXYZ[2]};
  
  Double_t*ubuf(0);
  
  // experimental hall, "world volume"  
  TGeoVolume* top = gGeoManager->Volume("EXPH","BOX",indVac,exphXYZ,3);
  gGeoManager->SetTopVolume(top);

  // calorimeter layers
  Int_t ilayer(1);
  Double_t x=2.*offsetX;
  for (;x<exphXYZ[0]-2*layer1XYZ[0]-2*layer2XYZ[0];++ilayer) {
    cout<<"ilayer="<<ilayer<<endl;
    x+=layer1XYZ[0];
    stringstream ss1; ss1<<"ABS"<<ilayer;
    gGeoManager->Volume(ss1.str().c_str(),"BOX",indPb,layer1XYZ,3);
    gGeoManager->Node(ss1.str().c_str(),1,"EXPH",x,0.,0.,0,kTRUE,ubuf);
    x+=layer1XYZ[0];
    x+=layer2XYZ[0];
    stringstream ss2; ss2<<"SCI"<<ilayer;
    gGeoManager->Volume(ss2.str().c_str(),"BOX",indSci,layer2XYZ,3);
    gGeoManager->Node(ss2.str().c_str(),1,"EXPH",x,0.,0.,0,kTRUE,ubuf);
    x+=layer2XYZ[0];
  }  
}
