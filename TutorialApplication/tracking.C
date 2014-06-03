#include "include/TutorialApplication.hh"
#include "TH1.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TObjArray.h"
#include "TGeoTrack.h"
#include "TROOT.h"
#include "TGeoMatrix.h"
#include "TMatrixTSym.h"
#include "THelix.h"
#include "Minuit2/FCNBase.h"
#include "TFitterMinuit.h"
#include "TList.h"
#include "TPad.h"


typedef TMatrixTSym<double>	TMatrixDSym;

TH1F *hlayer1 = new TH1F("hlayer1","layer1;z [cm]; counts",100/0.0150,-50,50);
TH1F *hlayer2 = new TH1F("hlayer2","layer2;z [cm]; counts",100/0.0150,-50,50);
TH1F *hlayer3 = new TH1F("hlayer3","layer3;z [cm]; counts",100/0.0150,-50,50);

TH1F *hresid2 = new TH1F("hresid2","resid2; z_{hit}-z_{orig} [cm]; events",100,-0.1,0.1);


class Cluster : public TVector3 {
public:
  Cluster(double x = 0, double y = 0, double startz = 0, double pitch = 0,unsigned char* strips = 0, 
	  int nstrips = 0) : 
    TVector3(x,y, startz), fStartz(startz), fPitch(pitch), fNstrips(nstrips), 
    fErrX(0), fErrY(0), fErrZ(0) {
    for(int i = 0 ; i < fNstrips ; ++i) {
      fStrips[i] = strips[i];
    }
  } 
  
  int nStrips() const { return fNstrips;}
  unsigned int signal(int istrip) const {return fStrips[istrip];}
  double  ZofFirstStrip() const { return fStartz;}
  double pitch() const { return fPitch;}

  double errX() const {return fErrX;}
  void setErrX(double ex) { fErrX = ex;} 
  double errY() const {return fErrY;}
  void setErrY(double ey) { fErrY = ey;} 
  double errZ() const {return fErrZ;}
  void setErrZ(double ez) { fErrZ = ez;} 
private:
  double fStartz;
  double fPitch;
  unsigned char fStrips[10];
  int fNstrips;
  double fErrX,fErrY,fErrZ;
  ClassDef(Cluster,0)
};

ClassImp(Cluster)

unsigned char getSignal(int volid) 
{ 
  TutorialApplication* app = (TutorialApplication*)TutorialApplication::Instance();
  
  int c = app->depEinVol(volid) * 600000;
  //add noise
  //c += gRandom->Gaus(0,3);
  //noise cut
  int noisecut = 0;
  if( c < noisecut ) return 0;
  if(c > 255) return 255;
  return c;
}



int updateClusters(TObjArray* clusters) 
{ 
  hlayer1->Reset();
  hlayer2->Reset();
  hlayer3->Reset();
  clusters->Clear();
  TObjArray* vollist = gGeoManager->GetListOfUVolumes();
  TGeoVolume* top = gGeoManager->GetMasterVolume();
  //top->GetNodes()->Print();
  unsigned char strips[10];
  TH1F* hcurrent = 0;
  for(int j = 3 ; j < vollist->GetEntries() ; ++j) {
    strips[0] = getSignal(j);
    if(strips[0]) {
      //std::cout << "volid = " << j << " Edep:" << app->depEinVol(j) << '\n';
      TGeoVolume* vol = (TGeoVolume*)vollist->At(j);
      TString name(vol->GetName());
      if(name.BeginsWith("Strip")) {
	int nstrips = 0;
	name += "_1";
	if(name.BeginsWith("Strip1")) hcurrent = hlayer1; 
	if(name.BeginsWith("Strip3")) hcurrent = hlayer2;
	if(name.BeginsWith("Strip5")) hcurrent = hlayer3;
 
	TGeoNode *node = top->FindNode(name);
	double zmin,zmax;
	vol->GetShape()->GetAxisRange(3,zmin,zmax);
	double pitch = zmax-zmin;
	const double* pos = node->GetMatrix()->GetTranslation();
	//std::cout << pos[0] << "," << pos[1] << ", " << pos[2]  << ":" << (int)strips[0] << '\n';
	hcurrent->Fill(pos[2],strips[0]);
	for(nstrips = 1 ; nstrips < 10 ; ++nstrips) {
	  strips[nstrips] = getSignal(j+nstrips);
	  if(! strips[nstrips]) break;
	  hcurrent->Fill(pos[2]+nstrips*pitch,strips[nstrips]);
	}
	std::cout << "cluster:" << pos[0] << "," << pos[1] << ", " << pos[2]  << ":" << (int)strips[0] << ", " << (int)strips[1] << " nstrips = " << nstrips << '\n';
	clusters->Add(new Cluster(pos[0],pos[1],pos[2],pitch,strips,nstrips));
	j+= nstrips;
      }
    }
  } 
  return clusters->GetEntriesFast();
}

int reconstructHits(TObjArray* clusters) 
{
  for(int i = 0 ; i < clusters->GetEntriesFast() ; ++i) {
    Cluster* c = (Cluster*)clusters->At(i);
    //use first strip for position
    c->SetZ(c->ZofFirstStrip());
    c->setErrZ(c->nStrips()*c->pitch());
  }
  return clusters->GetEntriesFast();
}


double getTrueZ(double detx) {
  //get primary track
  TObjArray* tracks = gGeoManager->GetListOfTracks();
  TGeoTrack* track = (TGeoTrack*)tracks->At(0);
  //track->Print();
  Double_t x,y,z,t;
  Int_t j=0;
  do {
    track->GetPoint(j,x,y,z,t);
    j++;
  } while(detx > x);
  double x1 = x;
  double z1 = z;
  j -= 2;
  track->GetPoint(j,x,y,z,t);
  double x2 = x;
  double z2 = z;
  //std::cout << x1 << "," << z1 << "; " << x2 << "," << z2 << '\n';
  //std::cout << (z2-z1)/(x2-x1)*(detx-x1) + z1 << '\n';
  return (z2-z1)/(x2-x1)*(detx-x1) + z1;
}


void plotResdiuals(TObjArray* clusters) {
  for(int i = 0 ; i < clusters->GetEntriesFast() ; ++i) {
    Cluster* c = (Cluster*)clusters->At(i);
    double x = c->X();
    //only use second layer
    if(x != -30) continue;
    double zorig = getTrueZ(x);    
    hresid2->Fill(zorig-c->Z());
  }
}




void tracking()
{
  TutorialApplication* app = (TutorialApplication*)TutorialApplication::Instance();

  // initialize calorimeter volumes and material   
  //app->InitMC("geometry/tracker"); 


  // define particle and control parameters of loop   
  unsigned int nevt = 1;
  double p = 1.0;
  app->SetPrimaryPDG(-13);    // +/-11: PDG code of e+/- 
  /* other PDG codes     22: Photon    +-13: muon   
                     +/-211: pion   +/-2212: proton     */
  app->SetPrimaryMomentum(p);
  // generate  some events
  TObjArray* clusters = new TObjArray();
  clusters->SetOwner(true);
  hresid2->Reset();
  for(unsigned int i=0;i<nevt;++i) {
    bool draw = !i;
    // p = gRandom->Gaus(2.0,0.1);
    //app->SetPrimaryMomentum(p);
    app->RunMC(1, draw);
    updateClusters(clusters);
    reconstructHits(clusters);
    plotResdiuals(clusters);
  }
  
  TCanvas* c = new TCanvas("c");
  c->cd();
  hlayer3->Draw();
  //hresid2->Draw();
}
