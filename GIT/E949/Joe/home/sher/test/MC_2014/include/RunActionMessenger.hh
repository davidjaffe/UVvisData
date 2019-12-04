#ifndef RunActionMessenger_h
#define RunActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4ios.hh"

class RunAction;
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

class RunActionMessenger : public G4UImessenger
{
  public:

   RunActionMessenger(RunAction* );
  ~RunActionMessenger();

   void SetNewValue(G4UIcommand* ,G4String );

  private:

   RunAction*              runAction;
   
   G4UIdirectory*             RndmDir;
   G4UIcmdWithAnInteger*      RndmSaveCmd;
   G4UIcmdWithABool*          SetAutoSeedCmd;
 
};

#endif
