#ifndef TUTORIALMAINFRAME_HH
#define TUTORIALMAINFRAME_HH
#include <RQ_OBJECT.h>
#include <TQObject.h>

class TutorialApplication;
class TGMainFrame;
class TGTextEntry;
class TGHSlider;

class TutorialMainFrame {
  RQ_OBJECT("TutorialMainFrame")

 public:
  // construction / destruction
  TutorialMainFrame(TutorialApplication* app);
  virtual ~TutorialMainFrame();

  // member functions
 public:
  void HandleMomentumSlider(int pos);
  void HandleMomentumText();
  void HandleType(int kf);
  void HandleNextEvent();
  void HandleExit();

  // member data
 private:
  TGMainFrame* fMain;
  TutorialApplication* fApp;
  TGTextEntry* fMomentumText;
  TGHSlider* fMomentumSlider;
};

#endif
