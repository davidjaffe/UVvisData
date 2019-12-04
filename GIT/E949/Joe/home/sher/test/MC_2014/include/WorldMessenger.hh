#ifndef WorldMessenger_h
#define WorldMessenger_h 1

#include "WorldConstructor.hh"
#include "GeoConstructorMessenger.hh"

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

class WorldMessenger : public GeoConstructorMessenger
{
  private:

    WorldConstructor* fConstructor;

    G4UIcmdWithADoubleAndUnit* fWidthCMD;
    G4UIcmdWithADoubleAndUnit* fHeightCMD;
    G4UIcmdWithADoubleAndUnit* fLengthCMD;

    G4UIcmdWithABool*   fUseStandAloneModuleCMD;
    G4UIcmdWithAString* fStandAloneModuleNameCMD;

  public:

    WorldMessenger(WorldConstructor* c);
    virtual ~WorldMessenger();

    void SetNewValue(G4UIcommand *cmd, G4String val);

};

#endif
