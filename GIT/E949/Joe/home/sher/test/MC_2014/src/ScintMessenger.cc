#include "globals.hh"

#include "ScintMessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

ScintMessenger::ScintMessenger(ScintConstructor* c)
                : GeoConstructorMessenger(c,"Control the scint. geometry."),
                  fConstructor(c)
{

  fTypeCMD = new G4UIcmdWithAString(CommandName("type"),this);
  fTypeCMD->SetGuidance("Set the type of the scintillator: Box or Tubs.");
  fTypeCMD->SetParameterName("type",true);
  fTypeCMD->SetDefaultValue("Box");
  fTypeCMD->SetCandidates("Box" "Tubs");

  fMaterialCMD = new G4UIcmdWithAString(CommandName("material"),this);
  fMaterialCMD->SetGuidance("Set the material name of the scintillator");
  fMaterialCMD->SetParameterName("material",true);
  fMaterialCMD->SetDefaultValue("Scintillator");

  fWidthCMD = new G4UIcmdWithADoubleAndUnit(CommandName("width"),this);
  fWidthCMD->SetGuidance("Set the width of the scintillator.");
  fWidthCMD->SetParameterName("Width",false);
  fWidthCMD->SetUnitCategory("Length");

  fHeightCMD = new G4UIcmdWithADoubleAndUnit(CommandName("height"),this);
  fHeightCMD->SetGuidance("Set the height of the scintillator.");
  fHeightCMD->SetParameterName("Height",false);
  fHeightCMD->SetUnitCategory("Length");

  fLengthCMD = new G4UIcmdWithADoubleAndUnit(CommandName("length"),this);
  fLengthCMD->SetGuidance("Set the length of the scintillator.");
  fLengthCMD->SetParameterName("Length",false);
  fLengthCMD->SetUnitCategory("Length");

  fPositionCMD = new G4UIcmdWith3VectorAndUnit(CommandName("position"),this);
  fPositionCMD->SetGuidance("Set the position of the scint. counter.");
  fPositionCMD->SetParameterName("X","Y","Z",true,true);
  fPositionCMD->SetDefaultUnit("cm");
  fPositionCMD->SetUnitCandidates("microm mm cm m km");
  fPositionCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
}

ScintMessenger::~ScintMessenger()
{
  delete fTypeCMD;
  delete fWidthCMD;
  delete fHeightCMD;
  delete fLengthCMD;
  delete fPositionCMD;
  delete fMaterialCMD;
}

void ScintMessenger::SetNewValue(G4UIcommand *cmd, G4String val)
{
  if (cmd==fTypeCMD) {
     fConstructor->SetType(val);
  }
  else if (cmd==fMaterialCMD) {
     fConstructor->SetMaterialName(val);
  }
  else if (cmd==fWidthCMD) {
     fConstructor->SetWidth(fWidthCMD->GetNewDoubleValue(val));
  }
  else if (cmd==fHeightCMD) {
     fConstructor->SetHeight(fHeightCMD->GetNewDoubleValue(val));
  }
  else if (cmd==fLengthCMD) {
     fConstructor->SetLength(fLengthCMD->GetNewDoubleValue(val));
  }
  else if (cmd==fPositionCMD) {
     fConstructor->SetPosition(fPositionCMD->GetNew3VectorValue(val));
  }
  else {
     GeoConstructorMessenger::SetNewValue(cmd,val);
  }
}
