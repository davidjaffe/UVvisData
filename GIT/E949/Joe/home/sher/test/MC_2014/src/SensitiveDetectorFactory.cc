#include "SensitiveDetectorFactory.hh"

#include "G4SDManager.hh"

#include "ScintSD.hh"
#include "CrystalSD.hh"
#include "WireChamberSD.hh"
#include "ZirconiumSD.hh"

SensitiveDetectorFactory::SensitiveDetectorFactory(G4String type)
  : fType(type) {}

SensitiveDetectorFactory::~SensitiveDetectorFactory() {}

G4VSensitiveDetector* SensitiveDetectorFactory::GetSD(G4String name)
{
  return G4SDManager::GetSDMpointer()->FindSensitiveDetector(name,false);
}

G4VSensitiveDetector* SensitiveDetectorFactory::MakeSD(G4String name)
{
  G4VSensitiveDetector* sd = GetSD(name);

  if (!sd)
  {
    if (fType == "segment") {
      if (name == "scintSD") {
	sd = new ScintSD(name);
      } else if (name == "crystalSD") {
	sd = new CrystalSD(name);
      } else if (name == "wireChamberSD") {
	sd = new WireChamberSD(name);

      } else if (name == "zirconiumSD") {
	sd = new ZirconiumSD(name);          // DvB's addition

      } else {
         G4cout << "No such sensitive detector " << name << G4endl;
         G4Exception("sensitive detector type not implemented","",FatalErrorInArgument,"");
      }
      G4SDManager::GetSDMpointer()->AddNewDetector(sd);
    }
  }
  return sd;
}
