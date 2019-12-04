#include "GeoConstructor.hh"
#include "GeoConstructorMessenger.hh"

GeoConstructorMessenger::GeoConstructorMessenger(GeoConstructor* c,
                                                 const char* guide)
{
    fConstructor = c;
    fDirName = c->GetName() + "/";
    fDirectory = new G4UIdirectory(fDirName);
    if (guide) fDirectory->SetGuidance(guide);
    else fDirectory->SetGuidance("Commands for a PIENU sub-detector.");

    fVisibleCMD = new G4UIcmdWithABool(CommandName("visible"),this);
    fVisibleCMD->SetGuidance("The object is drawn if this is true.");
    fVisibleCMD->SetParameterName("Visibility",false);

    fVisibleDaughtersCMD = new G4UIcmdWithABool(
                                         CommandName("visibleDaughters"),this);
    fVisibleDaughtersCMD->SetGuidance(
                    "The daughters of this object are drawn if this is true.");
    fVisibleDaughtersCMD->SetParameterName("Visibility",false);

    fSensitiveCMD = new G4UIcommand(CommandName("sensitive"),this);
    fSensitiveCMD->SetGuidance("Set the name of the sensitive detector");
    G4UIparameter* nameParam = new G4UIparameter('s');
    nameParam->SetParameterName("Name");
    fSensitiveCMD->SetParameter(nameParam);
    G4UIparameter* typeParam = new G4UIparameter('s');
    typeParam->SetParameterName("Type");
    fSensitiveCMD->SetParameter(typeParam);
}

void GeoConstructorMessenger::SetNewValue(G4UIcommand *cmd, G4String val)
{
    if (cmd == fVisibleCMD) {
        fConstructor->SetVisible(fVisibleCMD->GetNewBoolValue(val));
    }
    else if (cmd == fVisibleDaughtersCMD) {
        fConstructor->SetVisibleDaughters(
            fVisibleDaughtersCMD->GetNewBoolValue(val));
    }
    else if (cmd==fSensitiveCMD) {
        std::istringstream buf(val.c_str());
        std::string name, type;
        buf >> name;
        buf >> type;
        fConstructor->SetSensitiveDetector(name,type);
    }
}

GeoConstructorMessenger::~GeoConstructorMessenger()
{
    delete fDirectory;
    delete fVisibleCMD;
    delete fVisibleDaughtersCMD;

    delete fSensitiveCMD;
}
