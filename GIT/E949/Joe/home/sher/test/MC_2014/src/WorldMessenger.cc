#include "globals.hh"

#include "WorldMessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

WorldMessenger::WorldMessenger(WorldConstructor* c)
                : GeoConstructorMessenger(c,"Control the world geometry."),
                  fConstructor(c)
{
  fWidthCMD = new G4UIcmdWithADoubleAndUnit(CommandName("width"),this);
  fWidthCMD->SetGuidance("Set the width of the world.");
  fWidthCMD->SetParameterName("Width",false);
  fWidthCMD->SetUnitCategory("Length");

  fHeightCMD = new G4UIcmdWithADoubleAndUnit(CommandName("height"),this);
  fHeightCMD->SetGuidance("Set the height of the world.");
  fHeightCMD->SetParameterName("Height",false);
  fHeightCMD->SetUnitCategory("Length");

  fLengthCMD = new G4UIcmdWithADoubleAndUnit(CommandName("length"),this);
  fLengthCMD->SetGuidance("Set the length of the world.");
  fLengthCMD->SetParameterName("Length",false);
  fLengthCMD->SetUnitCategory("Length");

  fUseStandAloneModuleCMD =
          new G4UIcmdWithABool(CommandName("UseStandAloneModule"),this);
  fUseStandAloneModuleCMD->SetGuidance(
                                      "Set true to have only a single module");
  fUseStandAloneModuleCMD->SetParameterName("UseStandAloneModule",false);
  fUseStandAloneModuleCMD->SetDefaultValue(false);

  fStandAloneModuleNameCMD =
          new G4UIcmdWithAString(CommandName("StandAloneModuleName"),this);
  fStandAloneModuleNameCMD->SetGuidance("Name of Stand Alone Module");
  fStandAloneModuleNameCMD->SetParameterName("StandAloneModuleName",false);
  fStandAloneModuleNameCMD->SetDefaultValue("NoValueSet");
}

WorldMessenger::~WorldMessenger()
{
  delete fWidthCMD;
  delete fHeightCMD;
  delete fLengthCMD;

  delete fUseStandAloneModuleCMD;
  delete fStandAloneModuleNameCMD;
}

void WorldMessenger::SetNewValue(G4UIcommand *cmd, G4String val)
{
  if (cmd==fWidthCMD) {
     fConstructor->SetWidth(fWidthCMD->GetNewDoubleValue(val));
  }
  else if (cmd==fHeightCMD) {
     fConstructor->SetHeight(fHeightCMD->GetNewDoubleValue(val));
  }
  else if (cmd==fLengthCMD) {
     fConstructor->SetLength(fLengthCMD->GetNewDoubleValue(val));
  }
  else if (cmd==fUseStandAloneModuleCMD) {
     fConstructor->SetUseStandAloneModule
                              (fUseStandAloneModuleCMD->GetNewBoolValue(val));
  }
  else if (cmd==fStandAloneModuleNameCMD) {
                              fConstructor->SetStandAloneModuleName(val);
  }
  else {
     GeoConstructorMessenger::SetNewValue(cmd,val);
  }
}
