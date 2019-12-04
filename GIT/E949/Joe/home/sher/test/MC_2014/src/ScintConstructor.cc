#include "globals.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"

#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4TwoVector.hh"

#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"

#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "ScintSD.hh"

#include "G4UserLimits.hh"

#include "SensitiveDetectorFactory.hh"

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "Materials.hh"

#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4GeometryManager.hh"

#include "G4VisAttributes.hh"

#include "GeoConstructor.hh"
#include "ScintConstructor.hh"

#include "ScintMessenger.hh"

ScintConstructor::~ScintConstructor() {}

void ScintConstructor::Init(void) {

  SetMessenger(new GeoConstructorMessenger(this));

  SetVisible(false);
  SetSensitiveDetector("scintSD","segment");

  fType   = Box;
  fWidth  = 100.00*mm;
  fHeight = 100.00*mm;
  fLength =   6.60*mm;
  fRmin   =   0.0 *mm;
  fRmax   = 171.45*mm;
  fPosition.setX(  15.*cm);
  fPosition.setY( -50.*cm);
  fPosition.setZ(-176.*cm);
  fMaterialName = "Scintillator";
}

G4LogicalVolume *ScintConstructor::GetPiece(void) {

  G4LogicalVolume* logVolume;

  switch (fType) {
    case Box:

      logVolume
	= new G4LogicalVolume(new G4Box(GetName(),
                                        GetWidth()/2.0,
                                        GetHeight()/2.0,
                                        GetLength()/2.0),
                              FindMaterial(GetMaterialName()),
                              GetName());

      break;
    case Tubs:

      logVolume
        = new G4LogicalVolume(new G4Tubs(GetName(),
                                         GetRmin(),
                                         GetRmax(),
                                         GetLength()/2.0,
                                         0.*deg,360.*deg),
                              FindMaterial(GetMaterialName()),
                              GetName());
      break;
    default:

      logVolume
        = new G4LogicalVolume(new G4Box(GetName(),
                                        GetWidth()/2.0,
                                        GetHeight()/2.0,
                                        GetLength()/2.0),
                              FindMaterial(GetMaterialName()),
                              GetName());
  }

  if (GetSensitiveDetector()) {
     logVolume->SetSensitiveDetector(GetSensitiveDetector());
  }

  if (GetVisible()) {
      G4VisAttributes *visual = new G4VisAttributes();
      visual->SetColor(0.2,0.5,0.2,1); // green
      logVolume->SetVisAttributes(visual);
  } else {
      logVolume->SetVisAttributes(G4VisAttributes::Invisible);
  }

  return logVolume;

}
