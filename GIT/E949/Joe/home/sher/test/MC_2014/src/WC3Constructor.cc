#include "globals.hh"

#include "G4Box.hh"
#include "G4Trd.hh"
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
#include "WC3Constructor.hh"

#include "GeoConstructorMessenger.hh"

WC3Constructor::~WC3Constructor() {}

void WC3Constructor::Init(void) {

  SetVisible(false);
  SetSensitiveDetector("wireChamberSD","segment");

  SetMessenger(new GeoConstructorMessenger(this));
}

G4LogicalVolume *WC3Constructor::GetPiece(void) {

  G4double wc3AlCoverRmax     = (314.0/2.)*mm;
  G4double wc3AlCoverRmin     = (247.2/2.)*mm;
  G4double wc3AlCoverLength   =   4.0*mm;

  G4double wc3G10CoverRmax    = (314.0/2.)*mm;
  G4double wc3G10CoverRmin    = wc3G10CoverRmax - 9.25*mm;
  G4double wc3G10CoverLength  = 4.0*mm;

  G4double wc3AlSpacerRmax    = (314.0/2. - 9.25)*mm;
  G4double wc3AlSpacerRmin    = wc3AlSpacerRmax - 12.0*mm;
  G4double wc3AlSpacerLength  = 1.75*mm;

  G4double wc3G10SpacerARmax   = wc3G10CoverRmin;
  G4double wc3G10SpacerARmin   = (247.2/2.)*mm;
  G4double wc3G10SpacerLength = 2.0*mm;

  G4double wc3G10SpacerBRmax   = wc3G10CoverRmin;
  G4double wc3G10SpacerBRmin   = wc3AlSpacerRmin;

  G4double wc3Radius          = wc3AlCoverRmax;
  G4double wc3Length          = 2.*wc3AlCoverLength+4.*wc3G10CoverLength;

  G4double wc3GasRadius       = wc3AlCoverRmax;
  G4double wc3GasLength       = 4.*wc3G10CoverLength;

  G4double wc3FoilRadius      = wc3AlSpacerRmin;
  G4double wc3FoilLength      = 0.0508*mm;

  G4double wc3PlaneRadius      = wc3AlCoverRmin;
  G4double wc3PlaneLength      = wc3G10CoverLength - wc3FoilLength;

  G4Tubs* wc3_envelope
    = new G4Tubs("wc3_envelope",0.0,wc3Radius,wc3Length/2.,0.*deg,360.*deg);

  G4LogicalVolume* wc3_envelope_log
    = new G4LogicalVolume(wc3_envelope,FindMaterial("G4_AIR"),"wc3_envelope",0,0,0);

  G4VisAttributes *visual = new G4VisAttributes();
  visual->SetColor(0.5,0.1,0.1,1); // dark red
  visual->SetForceWireframe(true);

  wc3_envelope_log->SetVisAttributes(visual);

  G4Tubs* wc3_gas
    = new G4Tubs("wc3_gas",0.0,wc3GasRadius,wc3GasLength/2.,0.*deg,360.*deg);

  G4LogicalVolume* wc3_gas_log
    = new G4LogicalVolume(wc3_gas,FindMaterial("WCGas"),"wc3_gas",0,0,0);

  wc3_gas_log->SetVisAttributes(visual);

  new G4PVPlacement(0,G4ThreeVector(),
                    wc3_gas_log,"wc3_gas",wc3_envelope_log,false,0);

  G4Tubs* wc3_AlCover = new G4Tubs("wc3_AlCover",
                                   wc3AlCoverRmin,
                                   wc3AlCoverRmax,
                                   wc3AlCoverLength/2.,
                                   0.*deg,360.*deg);

  G4LogicalVolume* wc3_AlCover_log =
   new G4LogicalVolume(wc3_AlCover,FindMaterial("G4_Al"),"wc3_AlCover",0,0,0);

  wc3_AlCover_log->SetVisAttributes(visual);

  G4double z = -wc3Length/2. + wc3AlCoverLength/2.;

  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    wc3_AlCover_log,"wc3_AlCover",wc3_envelope_log,false,0);

  z =  wc3Length/2. - wc3AlCoverLength/2.;

  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    wc3_AlCover_log,"wc3_AlCover",wc3_envelope_log,false,1);

  G4Tubs* wc3_AlSpacer = new G4Tubs("wc3_AlSpacer",
                                    wc3AlSpacerRmin,
                                    wc3AlSpacerRmax,
                                    wc3AlSpacerLength/2.,
                                    0.*deg,360.*deg);

  G4LogicalVolume* wc3_AlSpacer_log =
   new G4LogicalVolume(wc3_AlSpacer,FindMaterial("G4_Al"),"wc3_AlSpacer",0,0,0);

  wc3_AlSpacer_log->SetVisAttributes(visual);

  z = -wc3GasLength/2. + 
          (wc3G10CoverLength/2. - wc3AlSpacerLength)/2. + wc3AlSpacerLength/2.;

  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    wc3_AlSpacer_log,"wc3_AlSpacer",wc3_gas_log,false,0);

  z = z +  wc3G10CoverLength;

  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    wc3_AlSpacer_log,"wc3_AlSpacer",wc3_gas_log,false,1);

  z = z +  wc3G10CoverLength;

  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    wc3_AlSpacer_log,"wc3_AlSpacer",wc3_gas_log,false,2);

  z = z +  wc3G10CoverLength;

  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    wc3_AlSpacer_log,"wc3_AlSpacer",wc3_gas_log,false,3);

  G4Tubs* wc3_G10Cover = new G4Tubs("wc3_G10Cover",
                                    wc3G10CoverRmin,
                                    wc3G10CoverRmax,
                                    2.*wc3G10CoverLength,
                                    0.*deg,360.*deg);

  G4LogicalVolume* wc3_G10Cover_log =
   new G4LogicalVolume(wc3_G10Cover,FindMaterial("G10"),"wc3_G10Cover",0,0,0);

  wc3_G10Cover_log->SetVisAttributes(visual);

  new G4PVPlacement(0,G4ThreeVector(),
                    wc3_G10Cover_log,"wc3_G10Cover",wc3_gas_log,false,0);

  G4Tubs* wc3_G10SpacerA = new G4Tubs("wc3_G10SpacerA",
                                      wc3G10SpacerARmin,
                                      wc3G10SpacerARmax,
                                      wc3G10SpacerLength/2.,
                                      0.*deg,360.*deg);

  G4LogicalVolume* wc3_G10SpacerA_log =
  new G4LogicalVolume(wc3_G10SpacerA,FindMaterial("G10"),"wc3_G10SpacerA",0,0,0);

  G4Tubs* wc3_G10SpacerB = new G4Tubs("wc3_G10SpacerB",
                                     wc3G10SpacerBRmin,
                                     wc3G10SpacerBRmax,
                                     wc3G10SpacerLength/2.,
                                     0.*deg,360.*deg);

  G4LogicalVolume* wc3_G10SpacerB_log =
  new G4LogicalVolume(wc3_G10SpacerB,FindMaterial("G10"),"wc3_G10SpacerB",0,0,0);

  wc3_G10SpacerA_log->SetVisAttributes(visual);
  wc3_G10SpacerB_log->SetVisAttributes(visual);

  z = -wc3GasLength/2. + wc3G10CoverLength - wc3G10SpacerLength/2.;

  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    wc3_G10SpacerA_log,"wc3_G10SpacerA",wc3_gas_log,false,0);

  z = z +  wc3G10CoverLength;

  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    wc3_G10SpacerA_log,"wc3_G10SpacerA",wc3_gas_log,false,1);

  z = z +  wc3G10CoverLength;

  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    wc3_G10SpacerA_log,"wc3_G10SpacerA",wc3_gas_log,false,2);

  z = z +  wc3G10CoverLength;

  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    wc3_G10SpacerB_log,"wc3_G10SpacerB",wc3_gas_log,false,3);

  G4Tubs* wc3_foil
    = new G4Tubs("wc3_foil",0.0,wc3FoilRadius,wc3FoilLength/2.,0.*deg,360.*deg);

  G4LogicalVolume* wc3_foil_log
    = new G4LogicalVolume(wc3_foil,FindMaterial("G4_MYLAR"),"wc3_foil",0,0,0);

  wc3_foil_log->SetVisAttributes(visual);

  G4Tubs* wc3_plane
    = new G4Tubs("wc3_plane",0.0,wc3PlaneRadius,wc3PlaneLength/2.,0.*deg,360.*deg);

  G4LogicalVolume* wc3_plane_log
    = new G4LogicalVolume(wc3_plane,FindMaterial("WCGas"),"wc3_plane",0,0,0);

  wc3_plane_log->SetVisAttributes(visual);

  z = -wc3GasLength/2. + wc3FoilLength/2.;

  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    wc3_foil_log,"wc3_foil",wc3_gas_log,false,0);

  z = -wc3GasLength/2. + wc3G10CoverLength/2. - wc3FoilLength/2.;

  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    wc3_foil_log,"wc3_foil",wc3_gas_log,false,1);

  z = z + wc3G10CoverLength/2.;

  G4RotationMatrix* g4rot = new G4RotationMatrix();
  g4rot->rotateZ(180.*degree);

  new G4PVPlacement(g4rot,G4ThreeVector(0.,0.,z),
                    wc3_plane_log,"wc3_plane",wc3_gas_log,false,0);

  z = z + wc3G10CoverLength/2.;

  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    wc3_foil_log,"wc3_foil",wc3_gas_log,false,2);

  z = z + wc3G10CoverLength/2.;

  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(60*degree);

  new G4PVPlacement(g4rot,G4ThreeVector(0.,0.,z),
                    wc3_plane_log,"wc3_plane",wc3_gas_log,false,1);

  z = z + wc3G10CoverLength/2.;

  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    wc3_foil_log,"wc3_foil",wc3_gas_log,false,3);

  z = z + wc3G10CoverLength/2.;

  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(-60*degree);

  new G4PVPlacement(g4rot,G4ThreeVector(0.,0.,z),
                    wc3_plane_log,"wc3_plane",wc3_gas_log,false,2);

  z = z + wc3G10CoverLength/2.;

  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    wc3_foil_log,"wc3_foil",wc3_gas_log,false,4);

  z =  wc3GasLength/2. - wc3FoilLength/2.;

  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    wc3_foil_log,"wc3_foil",wc3_gas_log,false,5);

//
// --> Fill these planes with trapezoidal drift cells
//
  G4int nwires = 48;
  G4double wire_spacing = 4.8*mm;
  G4double dy = wc3PlaneLength/2.;
  G4double dz = wire_spacing/2.;

  g4rot = new G4RotationMatrix();
  *g4rot = stringToRotationMatrix("Z90,Y90");
  *g4rot = g4rot->inverse();
  if (*g4rot == G4RotationMatrix()) g4rot = NULL;

  for (G4int i = 1 ; i <= nwires ; i++) {

    G4double ds1 = static_cast<G4double>(i-nwires/2-1)*wire_spacing;
    G4double ds2 = static_cast<G4double>(i-nwires/2  )*wire_spacing;

    G4double dx1 = std::sqrt(wc3PlaneRadius*wc3PlaneRadius - ds1*ds1);
    G4double dx2 = std::sqrt(wc3PlaneRadius*wc3PlaneRadius - ds2*ds2);

    G4double x = static_cast<G4double>(2*i-nwires-1)*wire_spacing/2.;

    G4Trd* wc3_trd = new G4Trd("Cell",dx1,dx2,dy,dy,dz);
    G4LogicalVolume* wc3_trd_log =
              new G4LogicalVolume(wc3_trd,FindMaterial("WCGas"),"Cell",0,0,0);

    wc3_trd_log->SetSensitiveDetector(GetSensitiveDetector());

    new G4PVPlacement(g4rot,G4ThreeVector(x,0.,0.),
                      wc3_trd_log,"Cell",wc3_plane_log,false,i-1);
  }

  return wc3_envelope_log;

}
