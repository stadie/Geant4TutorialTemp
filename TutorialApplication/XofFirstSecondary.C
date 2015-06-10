#include "TGeoManager.h"
#include "TGeoTrack.h"
#include "TDatabasePDG.h"

#include <iostream>

Float_t XofFirstSecondary()
{
  TObjArray* tracks = gGeoManager->GetListOfTracks();
  
  if(tracks->GetEntriesFast() < 2) return 999;
  
  //get first secondary
  TGeoTrack* track = (TGeoTrack*)tracks->At(1);
  
  //get the first point of this track
  Double_t x,y,z,t;
  track->GetPoint(0,x,y,z,t);
  
  return x;
}
