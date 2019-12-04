//Author: L.Doria 2011

#ifndef ParticleGunMessenger_h
#define ParticleGunMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class PrimaryGeneratorAction;

class G4UIdirectory;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;

class ParticleGunMessenger : public G4UImessenger
{

public:
  ParticleGunMessenger(PrimaryGeneratorAction*);
  virtual ~ParticleGunMessenger();

  void SetNewValue(G4UIcommand*, G4String);

private:
  PrimaryGeneratorAction* fPrimaryGeneratorAction;

  G4UIdirectory* fPileUpDirectory;
  G4UIdirectory* fBeamDirectory;

  G4UIcmdWithABool* fActive;
  G4UIcmdWithADouble* fRate;
  G4UIcmdWithADouble* fTmin;
  G4UIcmdWithADouble* fTmax;

  G4UIcmdWithABool* fDatabeam;

};
#endif
