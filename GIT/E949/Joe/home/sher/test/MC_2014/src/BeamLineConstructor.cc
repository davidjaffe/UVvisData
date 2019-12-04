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
#include "WireChamberSD.hh"

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
#include "BeamLineConstructor.hh"

#include "GeoConstructorMessenger.hh"

#include "G4SubtractionSolid.hh"

BeamLineConstructor::~BeamLineConstructor() {}

void BeamLineConstructor::Init(void) {

  SetVisible(false);
  //  SetSensitiveDetector("wireChamberSD","segment");

  SetMessenger(new GeoConstructorMessenger(this));
}

G4LogicalVolume *BeamLineConstructor::GetPiece(void) {


  G4VisAttributes *visual = new G4VisAttributes();
  visual->SetColor(0.5,0.1,0.1,1); // dark red
  visual->SetForceWireframe(true);


  G4double beamline8inchPipeLength          = 47.35*cm;
  G4double beamline8inchRadiusIn            = 98*mm;
  G4double beamline8inchRadiusOut           = 101.6*mm;

  G4double beamline4inchPipeLength          = 391.262*mm;
  G4double beamline4inchRadiusIn            = 48.514*mm;
  G4double beamline4inchRadiusOut           = 51.054*mm;

  G4double beamline8to4inchDiskLength       = 38.1*mm;
  
  G4double beamlineVacuumWindowLength       = 0.0762*mm;

  G4double beamlineEnvelopeLength   = beamline8inchPipeLength + beamline4inchPipeLength + beamlineVacuumWindowLength;
  
 

  G4Tubs* beamline_8inch_pipe = new G4Tubs("beamline_8inch_pipe",
					   beamline8inchRadiusIn,
					   beamline8inchRadiusOut,
					   beamline8inchPipeLength/2.,
					   0.*deg,360.*deg);
  
  G4LogicalVolume* beamline_8inch_pipe_log = new G4LogicalVolume(beamline_8inch_pipe,
								FindMaterial("Iron"),
								"beamline_8inch_pipe",0,0,0);

  
  G4Tubs* beamline_8inch_vacuum = new G4Tubs("beamline_8inch_vacuum",
					   0.,
					   beamline8inchRadiusIn,
					   beamline8inchPipeLength/2.,
					   0.*deg,360.*deg);
  
  G4LogicalVolume* beamline_8inch_vacuum_log = new G4LogicalVolume(beamline_8inch_vacuum,
								FindMaterial("Vac"),
								"beamline_8inch_vacuum",0,0,0);
  
  G4Tubs* beamline_4inch_pipe = new G4Tubs("beamline_4inch_pipe",
					   beamline4inchRadiusIn,
					   beamline4inchRadiusOut,
					   beamline4inchPipeLength/2.,
					   0.*deg,360.*deg);
  
  G4LogicalVolume* beamline_4inch_pipe_log = new G4LogicalVolume(beamline_4inch_pipe,
								FindMaterial("Iron"),
								"beamline_4inch_pipe",0,0,0);
  
  G4Tubs* beamline_4inch_vacuum = new G4Tubs("beamline_4inch_vacuum",
					   0.,
					   beamline4inchRadiusIn,
					   beamline4inchPipeLength/2.,
					   0.*deg,360.*deg);
  
  G4LogicalVolume* beamline_4inch_vacuum_log = new G4LogicalVolume(beamline_4inch_vacuum,
								FindMaterial("Vac"),
								"beamline_4inch_vacuum",0,0,0);


   G4Tubs* beamline_8to4inch_disk = new G4Tubs("beamline_8to4inch_disk",
					   beamline4inchRadiusOut,
					   beamline8inchRadiusOut,
					   beamline8to4inchDiskLength/2.,
					   0.*deg,360.*deg);
  
  G4LogicalVolume* beamline_8to4inch_disk_log = new G4LogicalVolume(beamline_8to4inch_disk,
								FindMaterial("Iron"),
								"beamline_8to4inch_disk",0,0,0);


   G4Tubs* beamline_vacuum_window = new G4Tubs("beamline_vacuum_window",
					   0,
					   beamline4inchRadiusOut,
					   beamlineVacuumWindowLength/2.,
					   0.*deg,360.*deg);
  
  G4LogicalVolume* beamline_vacuum_window_log = new G4LogicalVolume(beamline_vacuum_window,
								FindMaterial("Iron"),
								"beamline_vacuum_window",0,0,0);


  G4Tubs* beamline_envelope =  new G4Tubs("beamline_envelope",
					   0,
					   beamline8inchRadiusOut,
					   beamlineEnvelopeLength/2.,
					   0.*deg,360.*deg);
  
  G4LogicalVolume* beamline_envelope_log = new G4LogicalVolume(beamline_envelope,
								FindMaterial("G4_AIR"),
								"beamline_envelope",0,0,0);


  
  G4double z_8inch = (beamline8inchPipeLength-beamlineEnvelopeLength)*(0.5);
  
  new G4PVPlacement(0,G4ThreeVector(0.,0.,z_8inch),
		    beamline_8inch_pipe_log,"beamline8inchpipe",beamline_envelope_log,false,0);
  new G4PVPlacement(0,G4ThreeVector(0.,0.,z_8inch),
		    beamline_8inch_vacuum_log,"beamline8inchvacuum",beamline_envelope_log,false,0);
  
  G4double z_4inch = beamlineEnvelopeLength/2.-beamlineVacuumWindowLength-beamline4inchPipeLength/2.;
  
  new G4PVPlacement(0,G4ThreeVector(0.,0.,z_4inch),
                    beamline_4inch_pipe_log,"beamline4inchpipe",beamline_envelope_log,false,0);
  
  new G4PVPlacement(0,G4ThreeVector(0.,0.,z_4inch),
                    beamline_4inch_vacuum_log,"beamline4inchvacuum",beamline_envelope_log,false,0);
  
  G4double z_8to4inch = z_8inch + beamline8inchPipeLength/2.+beamline8to4inchDiskLength/2.;
  
  
  new G4PVPlacement(0,G4ThreeVector(0.,0.,z_8to4inch),
		    beamline_8to4inch_disk_log,"beamline8to4inchdisk",beamline_envelope_log,false,0);
  
  G4double z = beamlineEnvelopeLength/2.-beamlineVacuumWindowLength/2.;
  
  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
		    beamline_vacuum_window_log,"beamlinevacuumwindow",beamline_envelope_log,false,0);
  beamline_8inch_pipe_log->SetVisAttributes(visual);
  beamline_4inch_pipe_log->SetVisAttributes(visual);
  beamline_8to4inch_disk_log->SetVisAttributes(visual);
  beamline_vacuum_window_log->SetVisAttributes(visual);
  
  
  return beamline_envelope_log;
  
}
