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
  Track(double c, double phi, double d) :
    fCurv(c), fPhi0(phi), fD0(d),fCov(3) {}
  
  THelix* helix() const;
  
  int charge() const { return fCurv > 0 ? 1 : -1;} 
  double r() const { return std::abs(1/fCurv);}
  double curvature() const { return fCurv ? std::abs(fCurv) : 1e-10;}
  double d0() const { return fD0;}
  double phi0() const { return fPhi0;}
  double cov(int i, int j) const { return fCov(i,j);}

  void setCurvature(double c) { fCurv = c;};
  void setD0(double d0) { fD0 = d0;};
  void setPhi0(double phi0) { fPhi0 = phi0;};
  void setParameters(double a, double b, double c) {
    fCurv  = a;
    fPhi0  = b;
    fD0    = c;
  }
  double pt() const { return 1.49898e-04 * 2 * 20*r();}

  double rErr() const { return sqrt(fCov(0,0))/fCurv/fCurv;}
  double ptErr() const { return 1.49898e-04 * 2 *  20* rErr();}

  
  void setCov(int i, int j, double c) { fCov(i,j) = c;}
  
  double x(double lambda) const { return x0() + r() * charge() * sin(charge()*lambda + phi0());}
  double z(double lambda) const { return z0() - r() * charge() * cos(charge()*lambda + phi0());}
  double y(double lambda) const { return 0; }
  
  double lambdaFromX(double posx) const {return (charge()*asin( (posx-x0()) /charge()/r() ) - fPhi0);}
  
  double x0() const {return -sin(fPhi0) * (d0()+charge()*r()) -50;}
  double z0() const {return  cos(fPhi0) * (d0()+charge()*r());}
  
private:
  double fCurv, fPhi0, fD0;
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
  range[0] = std::min(50.0,x0()+r());
  range[1] = std::max(-50.0,x0()-r());
  axis[0] = 0;
  axis[1] = 1;
  axis[2] = 0;
  THelix* h = new THelix(xyz,v,1,range,kHelixY,axis);
  //h->Print();
  h->SetLineColor(7);
  return h;
}


ClassImp(Cluster)

unsigned char getSignal(int n) 
{ 
  if(! n) return 0;
  TutorialApplication* app = (TutorialApplication*)TutorialApplication::Instance();
  
  int c = app->depEinNode(n) * 600000;
  //add noise
  c += gRandom->Gaus(0,3);
  //noise cut
  int noisecut = 20;
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
  gGeoManager->CdTop();
  //gGeoManager->CdDown(0);
  //TGeoVolume* top = gGeoManager->GetMasterVolume();  
  
  unsigned char strips[10];
  TH1F* hcurrent = 0;
  //loop over layers
  for(int i = 0, ml = gGeoManager->GetCurrentNode()->GetNdaughters() ; i < ml ; ++i) {
    hcurrent = 0;
    if(i == 0) hcurrent = hlayer1;
    if(i == 1) hcurrent = hlayer2;
    if(i == 2) hcurrent = hlayer3;
    //gGeoManager->GetCurrentNode()->Print();
    gGeoManager->CdDown(i);
    //gGeoManager->GetCurrentNode()->Print();
    for(int j = 0, mj =  gGeoManager->GetCurrentNode()->GetNdaughters() ; j < mj ; ++j) {
      gGeoManager->CdDown(j);
      //gGeoManager->GetCurrentNode()->Print();
      TString name(gGeoManager->GetCurrentNode()->GetName());
      if(! name.BeginsWith("SiLayer")) { gGeoManager->CdUp(); continue;}
      //loop over strips
      for(int k = 0, mk = gGeoManager->GetCurrentNode()->GetNdaughters(); k < mk ; ++k) {
	gGeoManager->CdDown(k);
	//gGeoManager->GetCurrentNode()->Print();
	strips[0] = getSignal(gGeoManager->GetCurrentNodeId());
	if(strips[0]) {
	  int nstrips = 0;
	  double zmin,zmax;
	  gGeoManager->GetCurrentVolume()->GetShape()->GetAxisRange(3,zmin,zmax);
	  double pitch = zmax-zmin;
	  double local[3]={0,0,0};
	  double pos[3]={0,0,0};
	  gGeoManager->LocalToMaster(local,pos);
	  //std::cout << "Node:" << gGeoManager->GetCurrentNodeId() << ", " << gGeoManager->GetCurrentNode() << '\n';
	  //std::cout << pos[0] << "," << pos[1] << ", " << pos[2]  << ":" << (int)strips[0] << '\n';
	  if(hcurrent) hcurrent->Fill(pos[2],strips[0]);
	  for(nstrips = 1 ; nstrips < 10 ; ++nstrips) {
	    gGeoManager->CdUp();
	    gGeoManager->CdDown(k+nstrips);
	    strips[nstrips] = getSignal(gGeoManager->GetCurrentNodeId());
	    if(! strips[nstrips]) break;
	    if(hcurrent) hcurrent->Fill(pos[2]+nstrips*pitch,strips[nstrips]);
	  }
	  std::cout << "cluster:" << pos[0] << "," << pos[1] << ", " << pos[2]  << ":" << (int)strips[0] << ", " << (int)strips[1] << " nstrips = " << nstrips << '\n';
	  clusters->Add(new Cluster(pos[0],pos[1],pos[2],pitch,strips,nstrips));
	  k+= nstrips;
	}//strip with charge
	gGeoManager->CdUp();
      }//loop within silayer
      gGeoManager->CdUp();
    }//loop within layer
    gGeoManager->CdUp();
  }//loop over layers
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
    //use first strip for position
    c->SetZ(c->ZofFirstStrip()+maxstrip*c->pitch());
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
    
    // std::cout << "hit:" << x << ", " << z << "   track:" << fTrack->x(lambda) << ", " << fTrack->z(lambda) 
    //	      << "  lambda = " << lambda << '\n';
    
    double dZ = z - fTrack->z(lambda);
    double dX = x - fTrack->x(lambda);
    chi2 += (dZ*dZ)/c->errZ()/c->errZ() + dX*dX*1e6;
  }
  //std::cout << "Chi2:  " << chi2 << '\n';
  return chi2;
}



Track* fitTrack(TObjArray* clusters) {
  Track *t = new Track(0,0,0);

  TrackFCN fcn(t,clusters);
  
  TFitterMinuit * minuit = new TFitterMinuit(3);
  
  minuit->SetMinuitFCN(&fcn);
  minuit->SetParameter(0,"Curvature",t->curvature(),0.001,0,0);
  minuit->SetParameter(1,"Phi0",t->phi0(),0.1,0,0);
  minuit->SetParameter(2,"D0",t->d0(),0.1,0,0); 
  //minuit->SetParameter(0,"R",(t->charge() * t->r()),1,0,0);
  //minuit->SetParameter(1,"X0",t->x0(),1,0,0);
  //minuit->SetParameter(2,"Z0",t->z0(),1,0,0);
  minuit->SetPrintLevel(1);
  // create Minimizer (default is Migrad)
  minuit->CreateMinimizer();
  int iret = minuit->Minimize();
  minuit->PrintResults(1,0);
  if(iret != 0) {
    std::cout << "track fit failed.\n";
    t->setParameters(1e-10,0,0);
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


void trackingSolv()
{
  TutorialApplication* app = (TutorialApplication*)TutorialApplication::Instance();

  //std::cout << "app:" << app << '\n';
  //if(! app) return;
  //gMC->SetCut("CUTELE",0.0000000005); 
  //gMC->SetCut("CUTGAM",0.0000000005);
 
  // initialize calorimeter volumes and material   
  //app->InitMC("geometry/tracker"); 


  // define particle and control parameters of loop   
  unsigned int nevt = 1000;
  double p = 3.0;
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
    if(clusters->GetEntriesFast() >=3) {
      Track *t = fitTrack(clusters);
      if(t->r() < 1e06) { 
	if(draw) t->helix()->Draw();
	hpt->Fill(t->pt());
	std::cout << "Pt:" << t->pt() << " +- " << t->ptErr() << std::endl;
	hptpull->Fill((t->pt()-p)/t->ptErr());
      }
    }
  }
  /*
    TCanvas* c = new TCanvas("c");
    c->cd();
    //hlayer3->Draw();
    hresid2->Draw();
  */
}
