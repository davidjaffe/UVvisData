#include "RunActionMessenger.hh"

#include "RunAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

#include "G4ios.hh"
#include "globals.hh"

#include "Randomize.hh"

RunActionMessenger::RunActionMessenger(RunAction* RA)
  : runAction (RA)
{
  RndmDir = new G4UIdirectory("/rndm/");
  RndmDir->SetGuidance("Rndm status control.");
  
  SetAutoSeedCmd = new G4UIcmdWithABool("/rndm/autoSeed",this);
  SetAutoSeedCmd->SetGuidance("Switch on/off time-based random seeds");
  SetAutoSeedCmd->SetGuidance(" true: run seeds determined by system time");
  SetAutoSeedCmd->SetGuidance("false: use command 'random/resetEngineFrom'");
  SetAutoSeedCmd->SetGuidance("Default = false");
  SetAutoSeedCmd->SetParameterName("autoSeed", false);
  SetAutoSeedCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
}

RunActionMessenger::~RunActionMessenger()
{
  delete RndmDir; delete SetAutoSeedCmd;
}

void RunActionMessenger::SetNewValue(G4UIcommand* command,G4String newValues)
{ 
  if(command == SetAutoSeedCmd)
      runAction->SetAutoSeed(SetAutoSeedCmd->GetNewBoolValue(newValues)); 
}
