#include "globals.hh"

#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4Material.hh"
#include "G4RunManager.hh"
#include "G4GeometryManager.hh"

#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"

#include "Materials.hh"
#include "DetectorConstruction.hh"

#include "WorldConstructor.hh"
//#include "RootGeometryManager.hh"

DetectorConstruction::DetectorConstruction()
{
  fDetectorMessenger = new DetectorMessenger(this);
  fWorldConstructor = new WorldConstructor("/pienu",this);
}

DetectorConstruction::~DetectorConstruction()
{
  delete fDetectorMessenger;
  delete fWorldConstructor;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  DefineMaterials();
  G4VPhysicalVolume* world = ConstructDetector();

//  RootGeometryManager::Get()->Update(world);

//    G4VPersistencyManager *pMan 
//        = G4VPersistencyManager::GetPersistencyManager();
//    if (pMan) pMan->Store(world);

  return world;
}

void DetectorConstruction::DefineMaterials()
{
  Materials* materials = Materials::GetInstance();

  //define materials for the experiment

  if(!(materials->GetMaterial("G4_Galactic")))G4Exception("Vacuum not defined","",FatalErrorInArgument,"");
  if(!(materials->GetMaterial("G4_AIR")))G4Exception("Air not defined","",FatalErrorInArgument,"");
  if(!(materials->GetMaterial("Scintillator")))G4Exception("Sci not defined","",FatalErrorInArgument,"");
  if(!(materials->GetMaterial("G4_SODIUM_IODIDE")))G4Exception("NaI not defined","",FatalErrorInArgument,"");
  if(!(materials->GetMaterial("G4_CESIUM_IODIDE")))G4Exception("CsI not defined","",FatalErrorInArgument,"");
  if(!(materials->GetMaterial("G10")))G4Exception("G10 not defined","",FatalErrorInArgument,"");
  if(!(materials->GetMaterial("Rohacell")))G4Exception("Rohacell not defined","",FatalErrorInArgument,"");
  if(!(materials->GetMaterial("G4_Al")))G4Exception("Aluminum not defined","",FatalErrorInArgument,"");
  if(!(materials->GetMaterial("Iron")))G4Exception("Iron not defined","",FatalErrorInArgument,"");
  if(!(materials->GetMaterial("G4_MYLAR")))G4Exception("Mylar not defined","",FatalErrorInArgument,"");
  if(!(materials->GetMaterial("Isobutane")))G4Exception("Isobutane not defined","",FatalErrorInArgument,"");
  if(!(materials->GetMaterial("CF4")))G4Exception("CF4 not defined","",FatalErrorInArgument,"");
  if(!(materials->GetMaterial("WCGas")))G4Exception("WC Gas not defined","",FatalErrorInArgument,"");

  G4cout << G4endl << "The materials defined are: " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

}

G4VPhysicalVolume* DetectorConstruction::ConstructDetector()
{

  G4String name = "None";
  G4LogicalVolume* vol = NULL;

  name = fWorldConstructor->GetName();
  vol = fWorldConstructor->GetPiece();

  //------------------------------
  // World
  //------------------------------

  //  Must place the World Physical volume unrotated at (0,0,0).
  G4VPhysicalVolume* physWorld
        = new G4PVPlacement(0,               // no rotation
                            G4ThreeVector(), // position (0,0,0)
                            name,            // name
                            vol,             // logical volume
                            0,               // mother  volume
                            false,
                            0);

  return physWorld;

}

void DetectorConstruction::UpdateGeometry()
{
  // clean-up previous geometry

  G4GeometryManager::GetInstance()->OpenGeometry();

  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  //define new one

  G4VPhysicalVolume* world = Construct();

  G4RunManager::GetRunManager()->DefineWorldVolume(world);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  //RootGeometryManager::Get()->Update(world);

//  G4VPersistencyManager *pMan
//      = G4VPersistencyManager::GetPersistencyManager();
//  if (pMan) pMan->Store(world);
}

void DetectorConstruction::SetVerbose(G4int val)
{
  verbose = val;
}
