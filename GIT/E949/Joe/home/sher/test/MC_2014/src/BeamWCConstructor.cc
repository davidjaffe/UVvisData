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
#include "BeamWCConstructor.hh"

#include "GeoConstructorMessenger.hh"

BeamWCConstructor::~BeamWCConstructor() {}

void BeamWCConstructor::Init(void) {

  SetVisible(false);
  SetSensitiveDetector("wireChamberSD","segment");

  SetMessenger(new GeoConstructorMessenger(this));
}

G4LogicalVolume *BeamWCConstructor::GetPiece(void) {

  G4double wcAlWindowR = 101.6*mm/2. + 15.875*mm;
  G4double wcAlWindowL = 6.35*mm;

  G4double wcMylarWindowR = 101.6*mm/2.;
  G4double wcMylarWindowL = 0.025*mm;

  G4double wcAlWindowAirR = wcMylarWindowR;
  G4double wcAlWindowAirL = wcAlWindowL - wcMylarWindowL;

  G4double wcAlSpacerRmin   = wcAlWindowAirR;
  G4double wcAlSpacerRmax   = wcAlWindowR;
  G4double wcAlSpacerLength = 11.125*mm;

  G4double wcAlFlangeRmin   = wcAlWindowAirR;
  G4double wcAlFlangeRmax   = wcAlFlangeRmin + 7.512*mm;
  G4double wcAlFlangeLength = 3.175*mm;

  G4double v1VetoScintRmin   = 80.0*mm/2.;
  G4double v1VetoScintRmax   = v1VetoScintRmin + 12.0*mm;
  G4double v1VetoScintLength = 3.175*mm;

  G4double wcG10RingRmin   = wcMylarWindowR;
  G4double wcG10RingRmax   = wcAlWindowR;
  G4double wcG10RingLength = 1.588*mm;

  G4double wcMylarPlaneR = wcMylarWindowR;
  G4double wcMylarPlaneL = 0.006*mm;

  G4double wcPlaneR = wcMylarWindowR;
  G4double wcPlaneL = 2.*wcG10RingLength - 2.*wcMylarPlaneL;

  G4double beamWC_envelopeR = wcAlWindowR;
  G4double beamWC_envelopeL = v1VetoScintLength + 1.829*mm +
                              wcAlFlangeLength + 2.*wcAlWindowL +
                              7.*wcAlSpacerLength + 6.*wcG10RingLength;

  G4Tubs* beamWC_envelope
    = new G4Tubs("beamWC_envelope",
                 0.0,beamWC_envelopeR,beamWC_envelopeL/2.,0.*deg,360.*deg);

  G4LogicalVolume* beamWC_envelope_log
    = new G4LogicalVolume(beamWC_envelope,FindMaterial("G4_AIR"),"beamWC_envelope",0,0,0);

  G4VisAttributes *visual = new G4VisAttributes();
  visual->SetColor(0.5,0.1,0.1,1); // dark red
  visual->SetForceWireframe(true);

  beamWC_envelope_log->SetVisAttributes(visual);

  G4double beamWC_gasR = wcAlWindowR;
  G4double beamWC_gasL = 7.*wcAlSpacerLength + 6.*wcG10RingLength;

  G4Tubs* beamWC_gas
    = new G4Tubs("beamWC_gas",
                 0.0,beamWC_gasR,beamWC_gasL/2.,0.*deg,360.*deg);

  G4LogicalVolume* beamWC_gas_log
    = new G4LogicalVolume(beamWC_gas,FindMaterial("WCGas"),"beamWC_gas",0,0,0);

  beamWC_gas_log->SetVisAttributes(visual);

  G4double z = beamWC_envelopeL/2. - wcAlWindowL - beamWC_gasL/2.;

  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    beamWC_gas_log,"beamWC_gas",beamWC_envelope_log,false,0);

  G4Tubs* v1VetoScint = new G4Tubs("v1VetoScint",
                                   v1VetoScintRmin,
                                   v1VetoScintRmax,
                                   v1VetoScintLength/2.,
                                   0.*deg,360.*deg);

  G4LogicalVolume* v1VetoScint_log =
   new G4LogicalVolume(v1VetoScint,FindMaterial("Scintillator"),"v1VetoScint",0,0,0);

  v1VetoScint_log->SetVisAttributes(visual);

  z = -beamWC_envelopeL/2. + v1VetoScintLength/2.;

  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    v1VetoScint_log,"v1VetoScint",beamWC_envelope_log,false,0);

  G4Tubs* beamWC_AlFlange = new G4Tubs("beamWC_AlFlange",
                                       wcAlFlangeRmin,
                                       wcAlFlangeRmax,
                                       wcAlFlangeLength/2.,
                                       0.*deg,360.*deg);

  G4LogicalVolume* beamWC_AlFlange_log =
   new G4LogicalVolume(beamWC_AlFlange,FindMaterial("G4_Al"),"beamWC_AlFlange",0,0,0);

  beamWC_AlFlange_log->SetVisAttributes(visual);

  z = -beamWC_envelopeL/2. + v1VetoScintLength + 1.829*mm + wcAlFlangeLength/2.;

  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    beamWC_AlFlange_log,"beamWC_AlFlange",beamWC_envelope_log,false,0);

  G4Tubs* beamWC_AlWindow = new G4Tubs("beamWC_AlWindow",
                                       0.0,
                                       wcAlWindowR,
                                       wcAlWindowL/2.,
                                       0.*deg,360.*deg);

  G4LogicalVolume* beamWC_AlWindow_log =
   new G4LogicalVolume(beamWC_AlWindow,FindMaterial("G4_Al"),"beamWC_AlWindow",0,0,0);

  beamWC_AlWindow_log->SetVisAttributes(visual);

  G4Tubs* beamWC_AlWindowAir = new G4Tubs("beamWC_AlWindowAir",
                                          0.0,
                                          wcAlWindowAirR,
                                          wcAlWindowAirL/2.,
                                          0.*deg,360.*deg);

  G4LogicalVolume* beamWC_AlWindowAir_log =
   new G4LogicalVolume(beamWC_AlWindowAir,FindMaterial("G4_AIR"),"beamWC_AlWindowAir",0,0,0);

  beamWC_AlWindowAir_log->SetVisAttributes(visual);

  z = wcAlWindowL/2. - wcAlWindowAirL/2.;

  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    beamWC_AlWindowAir_log,"beamWC_AlWindowAir",beamWC_AlWindow_log,false,0);

  G4Tubs* beamWC_MylarWindow = new G4Tubs("beamWC_MylarWindow",
                                       0.0,
                                       wcMylarWindowR,
                                       wcMylarWindowL/2.,
                                       0.*deg,360.*deg);

  G4LogicalVolume* beamWC_MylarWindow_log =
   new G4LogicalVolume(beamWC_MylarWindow,FindMaterial("G4_MYLAR"),"beamWC_MylarWindow",0,0,0);

  beamWC_MylarWindow_log->SetVisAttributes(visual);

  z = -wcAlWindowL/2. + wcMylarWindowL/2.;

  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    beamWC_MylarWindow_log,"beamWC_MylarWindow",beamWC_AlWindow_log,false,0);

  G4RotationMatrix* g4rotate = new G4RotationMatrix();
  g4rotate->rotateY(180.*degree);

  z = -beamWC_envelopeL/2. + v1VetoScintLength + 1.829*mm
      + wcAlFlangeLength + wcAlWindowL/2.;

  new G4PVPlacement(g4rotate,G4ThreeVector(0.,0.,z),
                    beamWC_AlWindow_log,"beamWC_AlWindow",beamWC_envelope_log,false,0);

  z = beamWC_envelopeL/2. - wcAlWindowL/2.;

  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    beamWC_AlWindow_log,"beamWC_AlWindow",beamWC_envelope_log,false,1);

  G4Tubs* beamWC_AlSpacer = new G4Tubs("beamWC_AlSpacer",
                                      wcAlSpacerRmin,
                                      wcAlSpacerRmax,
                                      wcAlSpacerLength/2.,
                                      0.*deg,360.*deg);

  G4LogicalVolume* beamWC_AlSpacer_log =
   new G4LogicalVolume(beamWC_AlSpacer,FindMaterial("G4_Al"),"beamWC_AlSpacer",0,0,0);

  beamWC_AlSpacer_log->SetVisAttributes(visual);

  z = -beamWC_gasL/2. + wcAlSpacerLength/2.;

  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    beamWC_AlSpacer_log,"beamWC_AlSpacer",beamWC_gas_log,false,0);

  G4Tubs* beamWC_G10Ring = new G4Tubs("beamWC_G10Ring",
                                      wcG10RingRmin,
                                      wcG10RingRmax,
                                      wcG10RingLength/2.,
                                      0.*deg,360.*deg);

  G4LogicalVolume* beamWC_G10Ring_log =
   new G4LogicalVolume(beamWC_G10Ring,FindMaterial("G10"),"beamWC_G10Ring",0,0,0);

  beamWC_G10Ring_log->SetVisAttributes(visual);

  G4Tubs* beamWC_plane
    = new G4Tubs("beamWC_plane",
                 0.0,wcPlaneR,wcPlaneL/2.,0.*deg,360.*deg);

  G4LogicalVolume* beamWC_plane_log
    = new G4LogicalVolume(beamWC_plane,FindMaterial("WCGas"),"beamWC_plane",0,0,0);

  beamWC_plane_log->SetVisAttributes(visual);

  G4Tubs* beamWC_mylar
    = new G4Tubs("beamWC_mylar",
                 0.0,wcMylarPlaneR,wcMylarPlaneL/2.,0.*deg,360.*deg);

  G4LogicalVolume* beamWC_mylar_log
    = new G4LogicalVolume(beamWC_mylar,FindMaterial("G4_MYLAR"),"beamWC_mylar",0,0,0);

  beamWC_mylar_log->SetVisAttributes(visual);

  G4RotationMatrix* g4rot[6];
  g4rot[0] = new G4RotationMatrix();
  g4rot[0]->rotateZ(180.*degree);
  g4rot[1] = new G4RotationMatrix();
  g4rot[1]->rotateZ(-60*degree);
  g4rot[2] = new G4RotationMatrix();
  g4rot[2]->rotateZ(60*degree);
  g4rot[3] = g4rot[0];
  g4rot[4] = g4rot[1];
  g4rot[5] = g4rot[2];

  for (G4int j = 0 ; j < 6 ; j++) {

      z = z + wcAlSpacerLength/2. + wcG10RingLength/2.;

      new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                        beamWC_G10Ring_log,"beamWC_G10Ring",beamWC_gas_log,false,j);

      G4double zz = z - wcG10RingLength/2. + wcMylarPlaneL/2.;

      new G4PVPlacement(0,G4ThreeVector(0.,0.,zz),
                        beamWC_mylar_log,"beamWC_mylar",beamWC_gas_log,false,2*j);

      zz = z + 3.*wcG10RingLength/2. - wcMylarPlaneL/2.;

      new G4PVPlacement(0,G4ThreeVector(0.,0.,zz),
                        beamWC_mylar_log,"beamWC_mylar",beamWC_gas_log,false,2*j+1);

      zz = z + wcG10RingLength/2.;

      new G4PVPlacement(g4rot[j],G4ThreeVector(0.,0.,zz),
                        beamWC_plane_log,"beamWC_plane",beamWC_gas_log,false,j);

      z = z + wcG10RingLength/2. + wcAlSpacerLength/2.;

      new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                        beamWC_AlSpacer_log,"beamWC_AlSpacer",beamWC_gas_log,false,1);
  }

//
// --> Fill these planes with trapezoidal drift cells
//
  G4int nwires = 40;
  G4double wire_spacing = 2.4*mm;
  G4double dy = wcPlaneL/2.;
  G4double dz = wire_spacing/2.;

  g4rotate = new G4RotationMatrix();
  *g4rotate = stringToRotationMatrix("Z90,Y90");
  *g4rotate = g4rotate->inverse();
  if (*g4rotate == G4RotationMatrix()) g4rotate = NULL;

  for (G4int i = 1 ; i <= nwires ; i++) {

    G4double ds1 = static_cast<G4double>(i-nwires/2-1)*wire_spacing;
    G4double ds2 = static_cast<G4double>(i-nwires/2  )*wire_spacing;

    G4double dx1 = std::sqrt(wcPlaneR*wcPlaneR - ds1*ds1);
    G4double dx2 = std::sqrt(wcPlaneR*wcPlaneR - ds2*ds2);

    G4double x = static_cast<G4double>(2*i-nwires-1)*wire_spacing/2.;

    G4Trd* wc_trd = new G4Trd("Cellb",dx1,dx2,dy,dy,dz);
    G4LogicalVolume* wc_trd_log =
              new G4LogicalVolume(wc_trd,FindMaterial("WCGas"),"Cellb",0,0,0);

    wc_trd_log->SetSensitiveDetector(GetSensitiveDetector());

    new G4PVPlacement(g4rotate,G4ThreeVector(x,0.,0.),
                      wc_trd_log,"Cellb",beamWC_plane_log,false,i-1);
  }

  return beamWC_envelope_log;

}
