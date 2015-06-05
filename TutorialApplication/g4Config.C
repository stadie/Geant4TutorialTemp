//Modified Geant4 config file for inclusion of hadronic 
//interactions using the physics list "QGSP_BERT"
//Benjamin Treiber, Karlsruhe 2011-10-27

void Config()
{
  // RunConfiguration for Geant4 
  TG4RunConfiguration* runConfiguration = new
    //TG4RunConfiguration("geomRootToGeant4"); // only em interactions
    //TG4RunConfiguration("geomRootToGeant4","QGSP_BERT","stepLimiter+specialCuts"); 
                       //incl. hdronic interactions also, for Geant3 energy cuts-offs 
    TG4RunConfiguration("geomRootToGeant4","QGSP_BERT"); 
                       // incl. hdronic interactions also, Geant 4 range cuts
    //TG4RunConfiguration("geomRootToGeant4","FTFP_BERT");
  // TGeant4
  TGeant4* geant4
    = new TGeant4("TGeant4", "The Geant4 Monte Carlo", runConfiguration);

  cout << "Geant4 has been created." << endl;
 
  // Customise Geant4 setting
  // (verbose level, global range cut, ..)
  // geant4->ProcessGeantMacro("g4config.in");
  geant4->ProcessGeantCommand("/mcPhysics/rangeCuts 0.1 mm");
  // geant4->ProcessGeantCommand("/mcPhysics/rangeCutForGamma 1 mm");
 
}
