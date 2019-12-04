#include "globals.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

//#include "RootGeometryManager.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
 : fConstruction(Det)
{
  fPienuDir = new G4UIdirectory("/pienu/");
  fPienuDir -> SetGuidance(" Pienu Simulation ");

  fUpdateCmd = new G4UIcmdWithoutParameter("/pienu/Update",this);
  fUpdateCmd->SetGuidance("Update geometry.");
  fUpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  fUpdateCmd->SetGuidance("if you changed geometrical value(s)");
  fUpdateCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

//  fExportCmd = new G4UIcmdWithAString("/pienu/export",this);
//  fExportCmd->SetGuidance("Export the geometry to a file.");
//  fExportCmd->SetParameterName("RootFile", false);

  verboseCmd = new G4UIcmdWithAnInteger("/pienu/Verbose",this);
  verboseCmd->SetGuidance("Set verbose for ");
  verboseCmd->SetParameterName("verb",false);
  verboseCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

DetectorMessenger::~DetectorMessenger()
{
  delete fPienuDir;
  delete fUpdateCmd;
//  delete fExportCmd;
  delete verboseCmd;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if( command == fUpdateCmd )
    fConstruction->UpdateGeometry();
  else if( command == verboseCmd )
    fConstruction->SetVerbose(verboseCmd->GetNewIntValue(newValue));
//  else if( cmd == fExportCmd )
//    RootGeometryManager::Get()->Export(newValue);
}
