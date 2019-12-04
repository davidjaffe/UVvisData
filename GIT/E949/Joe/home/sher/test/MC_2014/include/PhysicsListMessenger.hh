#ifndef PhysicsListMessenger_h
#define PhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "G4DecayTable.hh"
#include "G4VDecayChannel.hh"

class PhysicsList;

class G4UIdirectory;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;

/// Provide control of the physics list and cut parameters

class PhysicsListMessenger : public G4UImessenger
{
public:

    PhysicsListMessenger(PhysicsList* );
    virtual ~PhysicsListMessenger();

    void SetNewValue(G4UIcommand*, G4String);

private:

  PhysicsList* fPhysicsList;

  G4UIdirectory* fDirectory;
  G4UIdirectory* fDecayDirectory;

  G4UIcmdWithADoubleAndUnit* fGammaCutCMD;
  G4UIcmdWithADoubleAndUnit* fElectCutCMD;
  G4UIcmdWithADoubleAndUnit* fPosCutCMD;
  G4UIcmdWithADoubleAndUnit* fAllCutCMD;
//  G4UIcmdWithAString*        fPhysicsListCMD;
//  G4UIcmdWithoutParameter*   fListCMD;

  G4UIcmdWithADouble* fPienuCMD;
  G4UIcmdWithADouble* fPimunuCMD;

  G4ParticleTable* particleTable;
  G4ParticleDefinition* particleDef;

  G4DecayTable* table;
  G4VDecayChannel* mode;

};
#endif
