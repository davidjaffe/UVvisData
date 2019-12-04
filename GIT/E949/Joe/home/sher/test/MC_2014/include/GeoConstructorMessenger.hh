#ifndef GeoConstructorMessenger_h
#define GeoConstructorMessenger_h 1

class GeoConstructor;

#include "G4UImessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"

class GeoConstructorMessenger: public G4UImessenger
{
  private:

    GeoConstructor* fConstructor;

    /// The UI directory for this messenger.
    G4UIdirectory* fDirectory;

    /// The directory name for this messenger
    G4String fDirName;

    G4UIcmdWithABool* fVisibleCMD;
    G4UIcmdWithABool* fVisibleDaughtersCMD;

    G4UIcommand*      fSensitiveCMD;

  public:

    GeoConstructorMessenger(GeoConstructor* c, const char* guide=NULL);
    virtual ~GeoConstructorMessenger();

    /// Return the name of the directory.
    G4String GetDirectory(void) {return fDirName;};

    /// Build a command name with the directory prefix.
    G4String CommandName(G4String cmd) {
        G4String name = GetDirectory() + cmd;
        return name;
    };

    void SetNewValue(G4UIcommand *cmd, G4String val);

};

#endif
