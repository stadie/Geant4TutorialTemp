////////////////////////////////////////////////////////////////////////////////
//
// tutorial application using geant4_vmc via TVirutalMCApplication
// ---------------------------------------------------------------
//
// This is an adaption of the TGeantApplication class written by
// Hartmut Stadie, which is an adaption of the Ex01MCApplication class
// written by Ivana Hrivnacova.
//
//            11/03/2008 Philipp Schieferdecker <philipp.schieferdecker@cern.ch>
////////////////////////////////////////////////////////////////////////////////

#ifndef TUTORIALAPPLICATION_H
#define TUTORIALAPPLICATION_H

#include <TVirtualMCApplication.h>
#include <TVirtualMCStack.h>
#include <TGeoManager.h>
#include <TFolder.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TVirtualMCGeometry.h>
#include <TGeoMCGeometry.h>
#include <TVirtualPad.h>
#include <TDatabasePDG.h>

#include <iostream>
#include <map>


class TVirtualGeoTrack;
class TGeoNode;

class TutorialApplication : public TVirtualMCApplication
{
public:
  //
  // construction / destruction
  //
  TutorialApplication(const char* name="TutorialApplication",
		      const char* title="MC using Geant4");
  virtual ~TutorialApplication();

public:
  //
  // member functions
  //
  void InitMC(const char *setup="geometry/g1");
  void RunMC(Int_t nofEvents = 1, bool draw=true);
  void FinishRun();
  
  // static access method
  static TutorialApplication* Instance();
  
  // construct materials & volumes
  void ConstructMaterials();
  void ConstructVolumes();
  
  // functions from TVirtualMCApplication
  virtual void ConstructGeometry();
  virtual void InitGeometry();
  virtual void GeneratePrimaries();
  virtual void BeginEvent();
  virtual void BeginPrimary(){;}
  virtual void PreTrack();
  virtual void Stepping();
  virtual void PostTrack();
  virtual void FinishPrimary(){;}
  virtual void FinishEvent();
  virtual void Field(const Double_t* x, Double_t* b) const;
  
  Double_t TrackingRmax() const { return DBL_MAX; }
  Double_t TrackingZmax() const { return DBL_MAX; }
  
  
  void SetPrimaryMomentum(Double_t p) { fPrimaryMomentum = p; }
  void SetPrimaryPDG(Int_t pdg) { fPrimaryPDG = pdg; }
  void SetPrimaryGeant(Int_t geant);
  void SetDrawPad(TVirtualPad* pad) { fPad = pad; }
  
  void DrawEvent();
  void Help();
  
  
  double depEinVol(int voluid) const;
  double depE(TGeoNode* n) const { std::map<TGeoNode*,double>::const_iterator i =  fDepEinVol.find(n); return i != fDepEinVol.end() ? i->second : 0;}
  const std::map<TGeoNode*,double>& depEMap() const { return fDepEinVol;}

  //
  // member data
  //
private:
  TVirtualMCStack* fStack;
  TFolder*         fTopFolder;
  TFolder*         fHistFolder;
  
  Double_t         fPrimaryMomentum;
  Int_t            fPrimaryPDG;
  
  TVirtualGeoTrack*fCurrentTrack; // !
  TString          fGeoName; // !
  Bool_t           fInit; // !
  TVirtualPad*     fPad; //!
  
  //histograms
  TH1F*     hEdepLong;
  TH1F*     hEdepTrans;
  TProfile* hPrimaryEnergy;

  std::map<TGeoNode*,double> fDepEinVol;
  
  
  ClassDef(TutorialApplication,1);
};


//
// implementation of inline functions
//

//______________________________________________________________________________
inline
TutorialApplication* TutorialApplication::Instance()
{ 
  return (TutorialApplication*)(TVirtualMCApplication::Instance());
}


//______________________________________________________________________________
inline
void TutorialApplication::BeginEvent()
{
  fStack->Clear();
  fDepEinVol.clear();
}


//______________________________________________________________________________
inline
void TutorialApplication::Field(const Double_t* /* x */, Double_t* b) const
{
  b[0] = 0.;
  b[1] = 0.;
  b[2] = 0.;
}


//______________________________________________________________________________
inline
void TutorialApplication::PostTrack()
{
  fCurrentTrack = NULL;
}


//______________________________________________________________________________
inline
void TutorialApplication::SetPrimaryGeant(Int_t geant)
{
  fPrimaryPDG = TDatabasePDG::Instance()->ConvertGeant3ToPdg(geant);
}


#endif
