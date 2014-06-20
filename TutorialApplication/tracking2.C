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
TH1F *hpt = new TH1F("hpt","; p_{T} [GeV]",100,0,10);
TH1F *hptpull = new TH1F("hptpull","; (p_{T}^{meas} - p_{T}^{true})/#sigma",100,-10,10);

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

class Track {
public:
  Track(double R, double x0, double z0) :
    fR(R), fX0(x0), fZ0(z0),fCov(3) {}
  
  THelix* helix() const;
  
  int charge() const { return fR > 0 ? 1 : -1;}
  double r() const { return std::abs(fR);}
  double x0() const { return fX0;}
  double z0() const { return fZ0;}
  double phi0() const { return 0;}
  


  double pt() const { return 0;}//needs changes

  double rErr() const { return sqrt(fCov(0,0));}
  double ptErr() const { return 0;}//needs changes


  double cov(int i, int j) const { return fCov(i,j);}
  
  void setParameters(double a, double b, double c) {
    fR  = a;
    fX0 = b;
    fZ0 = c;
  }
  
  void setCov(int i, int j, double c) { fCov(i,j) = c;}
  
  double x(double lambda) const { return 0;}//needs changes
  double z(double lambda) const { return 0;}//needs changes
  double y(double) const { return 0; }
  
  double lambdaFromX(double posx) const { //needs changes
    return 0;
  }

  
private:
  double fR, fX0, fZ0;
  TMatrixDSym fCov;
};

THelix* Track::helix() const {
  double xyz[3],v[3],range[2],axis[3];
  xyz[0] = x(0);
  xyz[1] = y(0);
  xyz[2] = z(0);
  v[0] = charge() * r();
  v[1] = 0;
  v[2] = 0;
  range[0] = 50;//std::min(50.0,x0()+r());
  range[1] = -50;//std::max(-50.0,x0()-r());
  axis[0] = 0;
  axis[1] = 1;
  axis[2] = 0;
  THelix* h = new THelix(xyz,v,1,range,kHelixY,axis);
  //h->Print();
  h->SetLineColor(7);
  return h;
}


ClassImp(Cluster)

unsigned char getSignal(int volid) 
{ 
  TutorialApplication* app = (TutorialApplication*)TutorialApplication::Instance();
  
  int c = app->depEinVol(volid) * 600000;
  //add noise
  //c += gRandom->Gaus(0,3);
  //noise cut
  int noisecut = 12;
  //if( c < noisecut ) return 0;
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

int reconstructHitsSimple(TObjArray* clusters) 
{
  for(int i = 0 ; i < clusters->GetEntriesFast() ; ++i) {
    Cluster* c = (Cluster*)clusters->At(i);
    //use first strip for position
    c->SetZ(c->ZofFirstStrip());
    c->setErrZ(c->nStrips()*c->pitch());
  }
  return clusters->GetEntriesFast();
}

int reconstructHitsBinary(TObjArray* clusters) 
{
  for(int i = 0 ; i < clusters->GetEntriesFast() ; ++i) {
    Cluster* c = (Cluster*)clusters->At(i);
    unsigned int maxc = c->signal(0);
    int maxstrip = 0;
    for(int j = 0 ; j < c->nStrips() ; ++j) {
      if(maxc < c->signal(j)) {
	maxstrip = j;
	maxc = c->signal(maxstrip);
      }
    }
    //std::cout << "maxstrip " << maxstrip << ", " << maxc << '\n';
    //use strip with max charge for position
    c->SetZ(c->ZofFirstStrip()+maxstrip*c->pitch());
    c->setErrZ(c->pitch());
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


// taken from http://root.cern.ch/viewvc/trunk/math/minuit2/test/testMinimize.cxx?view=markup
class TrackFCN : public ROOT::Minuit2::FCNBase { 

public: 

  TrackFCN(Track *t, TObjArray *clusters) : fTrack(t), fClusters(clusters) {}

  double operator() (const std::vector<double> & x) const;
  
  double Up() const { return 1.; }

private: 
  Track* fTrack;
  TObjArray *fClusters;
};

double TrackFCN::operator()(const std::vector<double> & x) const {
  //set new track parameters:
  fTrack->setParameters(x[0],x[1],x[2]);
  //std::cout << "track: r,x0,z0 = " << fTrack->r() << ", " << fTrack->x0() << ", " << fTrack->z0() << '\n'; 
  double chi2 = 0;
  for(int i = 0 ; i < fClusters->GetEntriesFast() ; ++i) {
    Cluster* c = (Cluster*)fClusters->At(i);
    double x = c->X();
    double z = c->Z();
    double lambda = fTrack->lambdaFromX(x);
    
    //std::cout << "hit:" << x << ", " << z << "   track:" << fTrack->x(lambda) << ", " << fTrack->z(lambda) 
    //	      << "  lambda = " << lambda << '\n';
    
    double dZ = z - fTrack->z(lambda);
    double dX = x - fTrack->x(lambda);
    chi2 += (dZ*dZ)/c->errZ()/c->errZ() + dX*dX*1e6;
  }
  //std::cout << "Chi2:  " << chi2 << '\n';
  return chi2;
}



Track* fitTrack(TObjArray* clusters) {
  Track *t = new Track(100,-50,100);

  TrackFCN fcn(t,clusters);
  
  TFitterMinuit * minuit = new TFitterMinuit(3);
  
  minuit->SetMinuitFCN(&fcn);
  //minuit->SetParameter(0,"Curvature",t->curvature(),0.00001,0,0);
  //minuit->SetParameter(1,"Phi0",t->phi0(),0.001,0,0);
  //minuit->SetParameter(2,"D0",t->d0(),0.1,0,0); 
  minuit->SetParameter(0,"R",(t->charge() * t->r()),1,0,0);
  minuit->SetParameter(1,"X0",t->x0(),1,0,0);
  minuit->SetParameter(2,"Z0",t->z0(),1,0,0);
  minuit->SetPrintLevel(1);
  // create Minimizer (default is Migrad)
  minuit->CreateMinimizer();
  int iret = minuit->Minimize();
  minuit->PrintResults(1,0);
  if(iret != 0) {
    std::cout << "track fit failed.\n";
    t->setParameters(1.0e10,0,0);
  } else {
    for(int i = 0 ; i < 3 ; ++i) {
      for(int j = i ; j < 3; ++j) {
	t->setCov(i,j,minuit->GetCovarianceMatrixElement(i,j));
      }
    }
  }
  return t;
}

void removeAllHelices() {
  TObjLink *lnk = gPad->GetListOfPrimitives()->FirstLink();
  while (lnk) {
    TObject* to = lnk->GetObject();
    if(to->InheritsFrom(THelix::Class())) to->Delete(); 
    lnk = lnk->Next();    
  }
}


void tracking2()
{
  TutorialApplication* app = (TutorialApplication*)TutorialApplication::Instance();
  // initialize calorimeter volumes and material   
  //app->InitMC("geometry/tracker"); 


  // define particle and control parameters of loop   
  unsigned int nevt = 1;
  double p = 5.0;
  app->SetPrimaryPDG(-13);    // +/-11: PDG code of e+/- 
  /* other PDG codes     22: Photon    +-13: muon   
                     +/-211: pion   +/-2212: proton     */
  app->SetPrimaryMomentum(p);
  // generate  some events
  TObjArray* clusters = new TObjArray();
  clusters->SetOwner(true);
  hresid2->Reset();
  hpt->Reset();
  for(unsigned int i=0;i<nevt;++i) {
    bool draw = !i;
    // p = gRandom->Gaus(2.0,0.1);
    //app->SetPrimaryMomentum(p);
    removeAllHelices();
    app->RunMC(1, draw);
    updateClusters(clusters);
    reconstructHitsBinary(clusters);
    plotResdiuals(clusters);
    if(clusters->GetEntriesFast() >= 3) {
      Track *t = fitTrack(clusters);
      if(draw) t->helix()->Draw();
      hpt->Fill(t->pt());
      hptpull->Fill((t->pt()-p)/t->ptErr());
    } else {
      std::cout << "Warning: Not enough hits for track fit.\n";
    }
  }
  /*
    TCanvas* c = new TCanvas("c");
    c->cd();
    //hlayer3->Draw();
    hresid2->Draw();
  */
}
