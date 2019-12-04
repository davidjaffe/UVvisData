#ifndef ScintMessenger_h
#define ScintMessenger_h 1

#include "ScintConstructor.hh"
#include "GeoConstructorMessenger.hh"

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

class ScintMessenger : public GeoConstructorMessenger
{
  private:

    ScintConstructor* fConstructor;

    G4UIcmdWithAString*        fTypeCMD;
    G4UIcmdWithADoubleAndUnit* fWidthCMD;
    G4UIcmdWithADoubleAndUnit* fHeightCMD;
    G4UIcmdWithADoubleAndUnit* fLengthCMD;
    G4UIcmdWith3VectorAndUnit* fPositionCMD;
    G4UIcmdWithAString*        fMaterialCMD;

  public:

    ScintMessenger(ScintConstructor* c);
    virtual ~ScintMessenger();

    void SetNewValue(G4UIcommand *cmd, G4String val);

};
#endif
