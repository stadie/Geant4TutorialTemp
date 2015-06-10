#include "TGeoManager.h"
#include "TGeoTrack.h"
#include "TDatabasePDG.h"

#include <iostream>

Int_t CountChargedinScint()
{
  Int_t ncharged = 0;

  TObjArray* tracks = gGeoManager->GetListOfTracks();

  for(Int_t i=0,l=tracks->GetEntriesFast();i<l;++i) {
    TGeoTrack* track = (TGeoTrack*)tracks->At(i);
    //track->Print();
    
    Int_t pdg = track->GetPDG();
    Double_t charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
    if( charge == 0 ) continue;
    
    Double_t x,y,z,t;
    TGeoNode* lastnode = NULL;
    for(Int_t j=0,lp=track->GetNpoints();j<lp;++j) {
      track->GetPoint(j,x,y,z,t);
      TGeoNode* node = gGeoManager->FindNode(x,y,z);
      if(! node ) continue;
      //node->Print();
      
      if( lastnode == node ) continue;
      lastnode = node;
      //is scintillator ?
      //std::cout << node->GetMedium()->GetMaterial()->GetZ() << std::endl;
      
      if(node->GetMedium()->GetMaterial()->GetZ() == 1)
	++ncharged;
    }
    //std::cout << "charge:" << ncharged << std::endl;
  }
  //std::cout << ncharged << std::endl;
  return ncharged;
}
