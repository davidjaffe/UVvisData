#include "globals.hh"

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

#include "EventAction.hh"
#include "EventActionMessenger.hh"

EventActionMessenger::EventActionMessenger(EventAction* EvAct)
     :eventAction(EvAct)
{ 
  drawCmd = new G4UIcmdWithAString("/event/drawTracks",this);
  drawCmd->SetGuidance("Draw the tracks in the event");
  drawCmd->SetGuidance("  Choice : none, neutral, charged, charged+n, all");
  drawCmd->SetParameterName("choice",true);
  drawCmd->SetDefaultValue("all");
  drawCmd->SetCandidates("none neutral charged charged+n all");
  drawCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  printCmd = new G4UIcmdWithAnInteger("/event/printModulo",this);
  printCmd->SetGuidance("Print events modulo n");
  printCmd->SetParameterName("EventNb",false);
  printCmd->SetRange("EventNb>0");
  printCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  debugCmd = new G4UIcmdWithAnInteger("/event/debugEvent",this);
  debugCmd->SetGuidance("D event to debug");
  debugCmd->SetParameterName("fNb",false);
  debugCmd->SetRange("fNb>=0");
  debugCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

EventActionMessenger::~EventActionMessenger()
{
  delete drawCmd;
  delete printCmd;
  delete debugCmd;
}

void EventActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if(command == drawCmd)
    {eventAction->SetDrawFlag(newValue);}

  if(command == printCmd)
    {eventAction->SetPrintModulo(printCmd->GetNewIntValue(newValue));}

  if(command == debugCmd)
    {eventAction->AddEventToDebug(debugCmd->GetNewIntValue(newValue));}
}
