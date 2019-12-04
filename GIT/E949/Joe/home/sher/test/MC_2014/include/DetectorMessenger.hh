#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;

class DetectorConstruction;

class DetectorMessenger : public G4UImessenger
{
public:

  DetectorMessenger(DetectorConstruction* );
  virtual ~DetectorMessenger();

  void SetNewValue(G4UIcommand*, G4String);

private:

  DetectorConstruction* fConstruction;

  G4UIdirectory*             fPienuDir;

  G4UIcmdWithoutParameter*   fUpdateCmd;
  G4UIcmdWithAnInteger*      verboseCmd;

};

#endif
