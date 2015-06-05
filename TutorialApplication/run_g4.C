{
  // Load basic libraries
  gROOT->LoadMacro("$G4PATH/geant4_vmc/examples/macro/basiclibs.C");
  basiclibs();
  
  // Load Geant4 libraries
  gROOT->LoadMacro("$G4PATH/geant4_vmc/examples/macro/g4libs.C");
  g4libs();
  
  // Load the tutorial application library
  //gSystem->Load("libTutorialApplication");
  gSystem->SetIncludePath("-I/usr/include -I$G4PATH/geant4_vmc/include/geant4vmc -I$G4PATH/geant4.9.5.p02-install/include/Geant4 -Iinclude -I/usr/include/geant4");
  gROOT->LoadMacro("src/TutorialStack.cxx+g");
  gROOT->LoadMacro("src/TGeomWrapper.cc+g");
  gROOT->LoadMacro("src/TutorialApplication.cxx+g");
  gROOT->LoadMacro("src/TutorialMainFrame.cxx+g");


  // MC application
  TutorialApplication* app 
    = new TutorialApplication("TutorialApplication",
			      "Tutorial Application for HEP Lecture @KIT/UHH");
  
  // configure Geant4
  gROOT->LoadMacro("g4Config.C");
  Config();

  // instantiate graphical user interface for tutorial application
  new TutorialMainFrame(app);
}
