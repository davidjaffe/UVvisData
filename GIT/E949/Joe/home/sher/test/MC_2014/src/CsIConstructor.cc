#include "globals.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"

#include "G4ExtrudedSolid.hh"
#include "G4DisplacedSolid.hh"

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
#include "CsIConstructor.hh"

#include "GeoConstructorMessenger.hh"

CsIConstructor::~CsIConstructor() {}

void CsIConstructor::Init(void) {

  SetVisible(false);
  SetSensitiveDetector("crystalSD","segment");

  SetMessenger(new GeoConstructorMessenger(this));
}

G4LogicalVolume* CsIConstructor::GetPiece(void) {

  G4LogicalVolume* logVolume = ConstructCsI_Frame();
  ConstructCsI_Upstream();
  ConstructCsI_Downstream();
  
  return logVolume;

}
G4LogicalVolume* CsIConstructor::ConstructCsI_Frame(void)
{
  inch = 2.54*cm;

  G4double csi_flange_thickness = 0.25*inch;

  G4double csi_ring1_ir        = 22.113*inch/2.;
  G4double csi_ring1_thickness =  0.074*inch;
  G4double csi_ring1_length    = 22.185*inch;

  G4double csi_ring2_ir        = 29.138*inch/2.;
  G4double csi_ring2_thickness =  0.074*inch;
  G4double csi_ring2_length    = 21.175*inch;

  G4double csi_ring3_ir        = 36.207*inch/2.;
  G4double csi_ring3_thickness =  0.120*inch;
  G4double csi_ring3_length    = csi_ring1_length;

  G4double rmin = csi_ring1_ir - csi_flange_thickness;
  G4double rmax = csi_ring3_ir + csi_ring3_thickness + csi_flange_thickness;
  G4double dz   = csi_ring3_length;

  G4Tubs* csi_frame
    = new G4Tubs("csi_frame",rmin,rmax,dz/2.,0.*deg,360.*deg);

  G4LogicalVolume* csi_frame_log
    = new G4LogicalVolume(csi_frame,FindMaterial("Iron"),"csi_frame",0,0,0);

  G4VisAttributes* visual = new G4VisAttributes();
  visual->SetColor(0.5,0.1,0.1,1); // dark red
  visual->SetForceWireframe(true);
  csi_frame_log->SetVisAttributes(visual);

  rmin = csi_ring1_ir + csi_ring1_thickness;
  rmax = csi_ring3_ir;
  dz   = csi_ring3_length;

  G4Tubs* csi_air
    = new G4Tubs("csi_air",rmin,rmax,dz/2.,0.*deg,360.*deg);

  G4LogicalVolume* csi_air_log
    = new G4LogicalVolume(csi_air,FindMaterial("G4_AIR"),"csi_air",0,0,0);

  csi_air_log->SetVisAttributes(visual);

  new G4PVPlacement(0,G4ThreeVector(),
                    csi_air_log,"csi_air",csi_frame_log,false,0);

  rmin = csi_ring1_ir + csi_ring1_thickness;
  rmax = csi_ring2_ir;
  dz   = csi_ring3_length/2.;

  G4Tubs* csi_space1
    = new G4Tubs("csi_space1",rmin,rmax,dz/2.,0.*deg,360.*deg);

  rmin = csi_ring2_ir + csi_ring2_thickness;
  rmax = csi_ring3_ir;
  dz   = csi_ring3_length/2.;

  G4Tubs* csi_space2
    = new G4Tubs("csi_space2",rmin,rmax,dz/2.,0.*deg,360.*deg);

  csi_space1_log
   = new G4LogicalVolume(csi_space1,FindMaterial("G4_AIR"),"csi_space1",0,0,0);
  csi_space2_log
   = new G4LogicalVolume(csi_space2,FindMaterial("G4_AIR"),"csi_space2",0,0,0);
  csi_space3_log
   = new G4LogicalVolume(csi_space1,FindMaterial("G4_AIR"),"csi_space3",0,0,0);
  csi_space4_log
   = new G4LogicalVolume(csi_space2,FindMaterial("G4_AIR"),"csi_space4",0,0,0);

  csi_space1_log->SetVisAttributes(visual);
  csi_space2_log->SetVisAttributes(visual);
  csi_space3_log->SetVisAttributes(visual);
  csi_space4_log->SetVisAttributes(visual);

  new G4PVPlacement(0,G4ThreeVector(0.,0.,-dz/2.),
                    csi_space1_log,"csi_space1",csi_air_log,false,0);
  new G4PVPlacement(0,G4ThreeVector(0.,0.,-dz/2.),
                    csi_space2_log,"csi_space2",csi_air_log,false,1);

  G4RotationMatrix* g4rot = new G4RotationMatrix();
  *g4rot = stringToRotationMatrix("Y180");
  *g4rot = g4rot->inverse();
  if (*g4rot == G4RotationMatrix()) g4rot = NULL;

  new G4PVPlacement(g4rot,G4ThreeVector(0.,0., dz/2.),
                    csi_space3_log,"csi_space3",csi_air_log,false,2);
  new G4PVPlacement(g4rot,G4ThreeVector(0.,0., dz/2.),
                    csi_space4_log,"csi_space4",csi_air_log,false,3);

  rmin = csi_ring2_ir;
  rmax = csi_ring2_ir + csi_ring2_thickness;
  dz = csi_ring2_length;

  G4Tubs* csi_shell
    = new G4Tubs("csi_shell",rmin,rmax,dz/2.,0.*deg,360.*deg);

  G4LogicalVolume* csi_shell_log
    = new G4LogicalVolume(csi_shell,FindMaterial("Iron"),"csi_shell",0,0,0);

  csi_shell_log->SetVisAttributes(visual);

  new G4PVPlacement(0,G4ThreeVector(),
                    csi_shell_log,"csi_shell",csi_air_log,false,0);

  // Flange
  G4double csi_flange_length = 1.*inch;
  
  rmin = csi_ring3_ir + csi_ring3_thickness;
  rmax = rmin + csi_flange_thickness;
  dz = csi_ring3_length - 2*csi_flange_length;

  G4Tubs* csi_flange1
    = new G4Tubs("csi_flange1",rmin,rmax,dz/2.,0.*deg,360.*deg);

  rmax = csi_ring1_ir;
  rmin = rmax - csi_flange_thickness;
  dz = csi_ring1_length - 2*csi_flange_length;

  G4Tubs* csi_flange2
    = new G4Tubs("csi_flange2",rmin,rmax,dz/2.,0.*deg,360.*deg);

  G4LogicalVolume* csi_flange1_log
   = new G4LogicalVolume(csi_flange1,
                         FindMaterial("G4_AIR"),
                         "csi_flange1",0,0,0);
  G4LogicalVolume* csi_flange2_log
   = new G4LogicalVolume(csi_flange2,
                         FindMaterial("G4_AIR"),
                         "csi_flange2",0,0,0);

  csi_flange1_log->SetVisAttributes(visual);
  csi_flange2_log->SetVisAttributes(visual);

  new G4PVPlacement(0,G4ThreeVector(),
                    csi_flange1_log,"csi_flange1",csi_frame_log,false,0);
  new G4PVPlacement(0,G4ThreeVector(),
                    csi_flange2_log,"csi_flange2",csi_frame_log,false,1);

  // Inner Rib
  G4double csi_irib_x = 3.438*inch;
  G4double csi_irib_y = 0.074*inch;
  G4double csi_irib_z = 9.750*inch;

  // Outer Rib
  G4double csi_orib_x = 3.460*inch;
  G4double csi_orib_y = 0.074*inch;
  G4double csi_orib_z = 9.750*inch;

  G4Box* csi_irib
    = new G4Box("csi_irib",csi_irib_x/2.,csi_irib_y/2.,csi_irib_z/2.);
  G4Box* csi_orib
    = new G4Box("csi_orib",csi_orib_x/2.,csi_orib_y/2.,csi_orib_z/2.);

  G4ThreeVector transVector1(((csi_ring1_ir+csi_ring1_thickness)+csi_ring2_ir)/2.,0,0);
  G4ThreeVector transVector2(((csi_ring2_ir+csi_ring2_thickness)+csi_ring3_ir)/2.,0,0);

  G4DisplacedSolid* csi_irib_displaced = new G4DisplacedSolid("csi_irib",
                                                               csi_irib,
                                                               0,
                                                               transVector1);
  G4DisplacedSolid* csi_orib_displaced = new G4DisplacedSolid("csi_orib",
                                                               csi_orib,
                                                               0,
                                                               transVector2);

  G4LogicalVolume* csi_irib_log =
  new G4LogicalVolume(csi_irib_displaced,FindMaterial("Iron"),"csi_irib",0,0,0);
  G4LogicalVolume* csi_orib_log =
  new G4LogicalVolume(csi_orib_displaced,FindMaterial("Iron"),"csi_orib",0,0,0);

  visual = new G4VisAttributes();
  visual->SetColor(0.5,0.1,0.1,1); // dark red
  visual->SetForceWireframe(true);
  csi_irib_log->SetVisAttributes(visual);
  csi_orib_log->SetVisAttributes(visual);

  G4double angle = -90.0*degree;
  G4double z_irib = csi_ring3_length/4. - csi_irib_z/2.;
  G4double z_orib = csi_ring3_length/4. - csi_orib_z/2.;

  for (G4int i = 0; i < 7; ++i) {

      G4RotationMatrix* g4rot = new G4RotationMatrix();
      g4rot->rotateZ(angle);

      new G4PVPlacement(g4rot,G4ThreeVector(0.0,0.0,z_irib),
                        csi_irib_log,"csi_irib",csi_space1_log,false,i);
      new G4PVPlacement(g4rot,G4ThreeVector(0.0,0.0,z_orib),
                        csi_orib_log,"csi_orib",csi_space2_log,false,i);

      angle += (360./7.)*degree;
  }

  G4int ii = 0;
  G4double angle0 = -90.0*degree - 9.58*degree;

  for (G4int i = 0; i < 3; ++i) {

      G4RotationMatrix* g4rot = new G4RotationMatrix();
      g4rot->rotateZ(angle0);

      new G4PVPlacement(g4rot,G4ThreeVector(0.0,0.0,z_irib),
                        csi_irib_log,"csi_irib",csi_space3_log,false,ii);

      ii++;
      angle = angle0 + 51.73*degree;

      g4rot = new G4RotationMatrix();
      g4rot->rotateZ(angle);

      new G4PVPlacement(g4rot,G4ThreeVector(0.0,0.0,z_irib),
                        csi_irib_log,"csi_irib",csi_space3_log,false,ii);

      ii++;
      angle0 += 120.00*degree;
  }

  ii = 0;
  angle0 = -90.0*degree - 9.58*degree;

  for (G4int i = 0; i < 3; ++i) {

      G4RotationMatrix* g4rot = new G4RotationMatrix();
      g4rot->rotateZ(angle0);

      new G4PVPlacement(g4rot,G4ThreeVector(0.0,0.0,z_orib),
                        csi_orib_log,"csi_orib",csi_space4_log,false,ii);

      ii++;
      angle = angle0 + 53.52*degree;

      g4rot = new G4RotationMatrix();
      g4rot->rotateZ(angle);

      new G4PVPlacement(g4rot,G4ThreeVector(0.0,0.0,z_orib),
                        csi_orib_log,"csi_orib",csi_space4_log,false,ii);

      ii++;
      angle0 += 120.00*degree;
  }

  return csi_frame_log;

}

void CsIConstructor::ConstructCsI_Upstream(void)
{
  inch = 2.54*cm;

  G4double crystal_length = 9.843*inch;

  // Crystal14

  G4double x = 3.199 * std::cos((180.-102.86)*deg);
  G4double y = 3.199 * std::sin((180.-102.86)*deg);

  std::vector<G4TwoVector> polygon;
  polygon.push_back(G4TwoVector((-2.978/2.)*inch, (-3.369/2.)*inch));
  polygon.push_back(G4TwoVector((-2.978/2.-x)*inch, (-3.369/2.+y)*inch));
  polygon.push_back(G4TwoVector((-0.3/2.)*inch, (3.369/2.)*inch));
  polygon.push_back(G4TwoVector(( 0.3/2.)*inch, (3.369/2.)*inch));
  polygon.push_back(G4TwoVector(( 2.978/2.+x)*inch, (-3.369/2.+y)*inch));
  polygon.push_back(G4TwoVector(( 2.978/2.)*inch, (-3.369/2.)*inch));

  G4ExtrudedSolid* crystal14_solid = new G4ExtrudedSolid("Crystal14", 
                                                         polygon,
                                                         crystal_length/2.,
                                                         G4TwoVector(),
                                                         1.0,
                                                         G4TwoVector(),
                                                         1.0);

  // Crystal21

  x = 3.259 * std::cos((180.-98.57)*deg);
  y = 3.259 * std::sin((180.-98.57)*deg);

  polygon.clear();
  polygon.push_back(G4TwoVector((-3.010/2.)*inch, (-3.370/2.)*inch));
  polygon.push_back(G4TwoVector((-3.010/2.-x)*inch, (-3.370/2.+y)*inch));
  polygon.push_back(G4TwoVector((-0.3/2.)*inch, (3.370/2.)*inch));
  polygon.push_back(G4TwoVector(( 0.3/2.)*inch, (3.370/2.)*inch));
  polygon.push_back(G4TwoVector(( 3.010/2.+x)*inch, (-3.370/2.+y)*inch));
  polygon.push_back(G4TwoVector(( 3.010/2.)*inch, (-3.370/2.)*inch));

  G4ExtrudedSolid* crystal21_solid = new G4ExtrudedSolid("Crystal21",
                                                         polygon,
                                                         crystal_length/2.,
                                                         G4TwoVector(),
                                                         1.0,
                                                         G4TwoVector(),
                                                         1.0);

  // Crystal27

  x = 3.277 * std::cos((180.-96.67)*deg);
  y = 3.277 * std::sin((180.-96.67)*deg);

  polygon.clear();
  polygon.push_back(G4TwoVector((-3.146/2.)*inch, (-3.370/2.)*inch));
  polygon.push_back(G4TwoVector((-3.146/2.-x)*inch, (-3.370/2.+y)*inch));
  polygon.push_back(G4TwoVector((-0.3/2.)*inch, (3.370/2.)*inch));
  polygon.push_back(G4TwoVector(( 0.3/2.)*inch, (3.370/2.)*inch));
  polygon.push_back(G4TwoVector(( 3.146/2.+x)*inch, (-3.370/2.+y)*inch));
  polygon.push_back(G4TwoVector(( 3.146/2.)*inch, (-3.370/2.)*inch));

  G4ExtrudedSolid* crystal27_solid = new G4ExtrudedSolid("Crystal27U",
                                                         polygon,
                                                         crystal_length/2.,
                                                         G4TwoVector(),
                                                         1.0,
                                                         G4TwoVector(),
                                                         1.0);


  G4RotationMatrix* g4rot = new G4RotationMatrix();
  *g4rot = stringToRotationMatrix("Z-90");
  *g4rot = g4rot->inverse();
  if (*g4rot == G4RotationMatrix()) g4rot = NULL;

  G4double csi_ir = 28.1*cm;
  G4double csi_shell_thickness = 0.16*cm;
  G4double csi_crystal_dr = 8.6*cm;

  G4double CsI1_rmin = csi_ir+csi_shell_thickness;
  G4double CsI1_rmax = CsI1_rmin+csi_crystal_dr;

  G4double CsI2_rmin = CsI1_rmax+csi_shell_thickness;
  G4double CsI2_rmax = CsI2_rmin+csi_crystal_dr;

  G4ThreeVector transVector1((CsI1_rmin+CsI1_rmax)/2.,0,0);
  G4ThreeVector transVector2((CsI2_rmin+CsI2_rmax)/2.,0,0);

  G4DisplacedSolid* crystal14_displaced = new G4DisplacedSolid("Crystal14",
                                                               crystal14_solid,
                                                               g4rot,
                                                               transVector1);

  G4DisplacedSolid* crystal21i_displaced = new G4DisplacedSolid("Crystal21",
                                                               crystal21_solid,
                                                               g4rot,
                                                               transVector1);

  G4DisplacedSolid* crystal21o_displaced = new G4DisplacedSolid("Crystal21",
                                                               crystal21_solid,
                                                               g4rot,
                                                               transVector2);

  G4DisplacedSolid* crystal27_displaced = new G4DisplacedSolid("Crystal27U",
                                                               crystal27_solid,
                                                               g4rot,
                                                               transVector2);

  G4LogicalVolume* crystal14_log
    = new G4LogicalVolume(crystal14_displaced,FindMaterial("G4_CESIUM_IODIDE"),
                                                            "Crystal14",0,0,0);
  G4LogicalVolume* crystal21i_log
    = new G4LogicalVolume(crystal21i_displaced,FindMaterial("G4_CESIUM_IODIDE"),
                                                            "Crystal21",0,0,0);
  G4LogicalVolume* crystal21o_log
    = new G4LogicalVolume(crystal21o_displaced,FindMaterial("G4_CESIUM_IODIDE"),
                                                            "Crystal21",0,0,0);
  G4LogicalVolume* crystal27_log
    = new G4LogicalVolume(crystal27_displaced,FindMaterial("G4_CESIUM_IODIDE"),
                                                           "Crystal27U",0,0,0);

  G4VisAttributes* visual = new G4VisAttributes();
  visual->SetColor(0.0,0.5,1.0,1); // blue
  visual->SetForceWireframe(true);
  crystal14_log->SetVisAttributes(visual);
  crystal21i_log->SetVisAttributes(visual);
  crystal21o_log->SetVisAttributes(visual);
  crystal27_log->SetVisAttributes(visual);

//------------------------------------------------
// Sensitive detectors
//------------------------------------------------

  crystal14_log->SetSensitiveDetector(GetSensitiveDetector());
  crystal21i_log->SetSensitiveDetector(GetSensitiveDetector());
  crystal21o_log->SetSensitiveDetector(GetSensitiveDetector());
  crystal27_log->SetSensitiveDetector(GetSensitiveDetector());

  // *** Build CsI Upstream Outer Ring ***

  G4double csi_ring3_length    = 22.185*inch;
  G4double z_position = csi_ring3_length/4. - crystal_length/2.;
  G4ThreeVector position(0,0,z_position);

  G4double angle = -90.0*degree - 6.43*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27U",
                           csi_space2_log,false,0);

  angle += 13.037*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal21o_log,"Crystal21",
                           csi_space2_log,false,1);

  angle += 12.957*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27U",
                           csi_space2_log,false,2);

  angle += 12.718*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27U",
                           csi_space2_log,false,3);

  angle += 12.717*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27U",
                           csi_space2_log,false,4);

  angle += 13.037*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal21o_log,"Crystal21",
                           csi_space2_log,false,5);

  angle += 12.957*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27U",
                           csi_space2_log,false,6);

  angle += 12.718*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27U",
                           csi_space2_log,false,7);

  angle += 12.717*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27U",
                           csi_space2_log,false,8);

  angle += 13.037*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal21o_log,"Crystal21",
                           csi_space2_log,false,9);

  angle += 12.957*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27U",
                           csi_space2_log,false,10);

  angle += 12.718*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27U",
                           csi_space2_log,false,11);

  angle += 12.717*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27U",
                           csi_space2_log,false,12);

  angle += 13.037*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal21o_log,"Crystal21",
                           csi_space2_log,false,13);

  angle += 12.957*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27U",
                           csi_space2_log,false,14);

  angle += 12.718*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27U",
                           csi_space2_log,false,15);

  angle += 12.717*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27U",
                           csi_space2_log,false,16);

  angle += 13.037*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal21o_log,"Crystal21",
                           csi_space2_log,false,17);

  angle += 12.957*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27U",
                           csi_space2_log,false,18);

  angle += 12.718*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27U",
                           csi_space2_log,false,19);

  angle += 12.717*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27U",
                           csi_space2_log,false,20);

  angle += 13.037*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal21o_log,"Crystal21",
                           csi_space2_log,false,21);

  angle += 12.957*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27U",
                           csi_space2_log,false,22);

  angle += 12.718*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27U",
                           csi_space2_log,false,23);

  angle += 12.717*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27U",
                           csi_space2_log,false,24);

  angle += 13.037*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal21o_log,"Crystal21",
                           csi_space2_log,false,25);

  angle += 12.957*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27U",
                           csi_space2_log,false,26);

  angle += 12.718*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27U",
                           csi_space2_log,false,27);

  // *** Build CsI Upstream Inner Ring ***

  angle = -90.0*degree - 8.57*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal21i_log,"Crystal21",
                           csi_space1_log,false,0);

  angle += 17.803*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal14_log,"Crystal14",
                           csi_space1_log,false,1);

  angle += 17.199*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal21i_log,"Crystal21",
                           csi_space1_log,false,2);

  angle += 16.426*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal21i_log,"Crystal21",
                           csi_space1_log,false,3);

  angle += 17.803*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal14_log,"Crystal14",
                           csi_space1_log,false,4);

  angle += 17.199*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal21i_log,"Crystal21",
                           csi_space1_log,false,5);

  angle += 16.426*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal21i_log,"Crystal21",
                           csi_space1_log,false,6);

  angle += 17.803*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal14_log,"Crystal14",
                           csi_space1_log,false,7);

  angle += 17.199*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal21i_log,"Crystal21",
                           csi_space1_log,false,8);

  angle += 16.426*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal21i_log,"Crystal21",
                           csi_space1_log,false,9);

  angle += 17.803*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal14_log,"Crystal14",
                           csi_space1_log,false,10);

  angle += 17.199*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal21i_log,"Crystal21",
                           csi_space1_log,false,11);

  angle += 16.426*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal21i_log,"Crystal21",
                           csi_space1_log,false,12);

  angle += 17.803*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal14_log,"Crystal14",
                           csi_space1_log,false,13);

  angle += 17.199*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal21i_log,"Crystal21",
                           csi_space1_log,false,14);

  angle += 16.426*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal21i_log,"Crystal21",
                           csi_space1_log,false,15);

  angle += 17.803*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal14_log,"Crystal14",
                           csi_space1_log,false,16);

  angle += 17.199*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal21i_log,"Crystal21",
                           csi_space1_log,false,17);

  angle += 16.426*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal21i_log,"Crystal21",
                           csi_space1_log,false,18);

  angle += 17.803*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal14_log,"Crystal14",
                           csi_space1_log,false,19);

  angle += 17.199*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal21i_log,"Crystal21",
                           csi_space1_log,false,20);

}

void CsIConstructor::ConstructCsI_Downstream(void)
{
  inch = 2.54*cm;

  G4double crystal_length = 9.843*inch;

  // Crystal19

  G4double x = 3.245 * std::cos((180.-99.47)*deg);
  G4double y = 3.245 * std::sin((180.-99.47)*deg);

  std::vector<G4TwoVector> polygon;
  polygon.push_back(G4TwoVector((-3.082/2.)*inch, (-3.370/2.)*inch));
  polygon.push_back(G4TwoVector((-3.082/2.-x)*inch, (-3.370/2.+y)*inch));
  polygon.push_back(G4TwoVector((-0.3/2.)*inch, (3.370/2.)*inch));
  polygon.push_back(G4TwoVector(( 0.3/2.)*inch, (3.370/2.)*inch));
  polygon.push_back(G4TwoVector(( 3.082/2.+x)*inch, (-3.370/2.+y)*inch));
  polygon.push_back(G4TwoVector(( 3.082/2.)*inch, (-3.370/2.)*inch));

  G4ExtrudedSolid* crystal19_solid = new G4ExtrudedSolid("Crystal19",
                                                         polygon,
                                                         crystal_length/2.,
                                                         G4TwoVector(),
                                                         1.0,
                                                         G4TwoVector(),
                                                         1.0);

  // Crystal25

  x = 3.271 * std::cos((180.-97.20)*deg);
  y = 3.271 * std::sin((180.-97.20)*deg);

  polygon.clear();
  polygon.push_back(G4TwoVector((-3.212/2.)*inch, (-3.370/2.)*inch));
  polygon.push_back(G4TwoVector((-3.212/2.-x)*inch, (-3.370/2.+y)*inch));
  polygon.push_back(G4TwoVector((-0.3/2.)*inch, (3.370/2.)*inch));
  polygon.push_back(G4TwoVector(( 0.3/2.)*inch, (3.370/2.)*inch));
  polygon.push_back(G4TwoVector(( 3.212/2.+x)*inch, (-3.370/2.+y)*inch));
  polygon.push_back(G4TwoVector(( 3.212/2.)*inch, (-3.370/2.)*inch));

  G4ExtrudedSolid* crystal25_solid = new G4ExtrudedSolid("Crystal25",
                                                         polygon,
                                                         crystal_length/2.,
                                                         G4TwoVector(),
                                                         1.0,
                                                         G4TwoVector(),
                                                         1.0);

  // Crystal27

  x = 3.277 * std::cos((180.-96.67)*deg);
  y = 3.277 * std::sin((180.-96.67)*deg);

  polygon.clear();
  polygon.push_back(G4TwoVector((-3.146/2.)*inch, (-3.370/2.)*inch));
  polygon.push_back(G4TwoVector((-3.146/2.-x)*inch, (-3.370/2.+y)*inch));
  polygon.push_back(G4TwoVector((-0.3/2.)*inch, (3.370/2.)*inch));
  polygon.push_back(G4TwoVector(( 0.3/2.)*inch, (3.370/2.)*inch));
  polygon.push_back(G4TwoVector(( 3.146/2.+x)*inch, (-3.370/2.+y)*inch));
  polygon.push_back(G4TwoVector(( 3.146/2.)*inch, (-3.370/2.)*inch));

  G4ExtrudedSolid* crystal27_solid = new G4ExtrudedSolid("Crystal27D",
                                                         polygon,
                                                         crystal_length/2.,
                                                         G4TwoVector(),
                                                         1.0,
                                                         G4TwoVector(),
                                                         1.0);

  G4RotationMatrix* g4rot = new G4RotationMatrix();
  *g4rot = stringToRotationMatrix("Z-90");
  *g4rot = g4rot->inverse();
  if (*g4rot == G4RotationMatrix()) g4rot = NULL;

  G4double csi_ir = 28.1*cm;
  G4double csi_shell_thickness = 0.16*cm;
  G4double csi_crystal_dr = 8.6*cm;

  G4double CsI1_rmin = csi_ir+csi_shell_thickness;
  G4double CsI1_rmax = CsI1_rmin+csi_crystal_dr;

  G4double CsI2_rmin = CsI1_rmax+csi_shell_thickness;
  G4double CsI2_rmax = CsI2_rmin+csi_crystal_dr;

  G4ThreeVector transVector1((CsI1_rmin+CsI1_rmax)/2.,0,0);
  G4ThreeVector transVector2((CsI2_rmin+CsI2_rmax)/2.,0,0);

  G4DisplacedSolid* crystal19i_displaced = new G4DisplacedSolid("Crystal19",
                                                               crystal19_solid,
                                                               g4rot,
                                                               transVector1);

  G4DisplacedSolid* crystal19o_displaced = new G4DisplacedSolid("Crystal19",
                                                               crystal19_solid,
                                                               g4rot,
                                                               transVector2);

  G4DisplacedSolid* crystal25_displaced = new G4DisplacedSolid("Crystal25",
                                                               crystal25_solid,
                                                               g4rot,
                                                               transVector2);

  G4DisplacedSolid* crystal27_displaced = new G4DisplacedSolid("Crystal27D",
                                                               crystal27_solid,
                                                               g4rot,
                                                               transVector1);

  G4LogicalVolume* crystal19i_log
    = new G4LogicalVolume(crystal19i_displaced,FindMaterial("G4_CESIUM_IODIDE"),
                                                            "Crystal19",0,0,0);
  G4LogicalVolume* crystal19o_log
    = new G4LogicalVolume(crystal19o_displaced,FindMaterial("G4_CESIUM_IODIDE"),
                                                            "Crystal19",0,0,0);
  G4LogicalVolume* crystal25_log
    = new G4LogicalVolume(crystal25_displaced,FindMaterial("G4_CESIUM_IODIDE"),
                                                           "Crystal25",0,0,0);
  G4LogicalVolume* crystal27_log
    = new G4LogicalVolume(crystal27_displaced,FindMaterial("G4_CESIUM_IODIDE"),
                                                           "Crystal27D",0,0,0);

  G4VisAttributes* visual = new G4VisAttributes();
  visual->SetColor(0.0,0.5,1.0,1); // blue
  visual->SetForceWireframe(true);
  crystal19i_log->SetVisAttributes(visual);
  crystal19o_log->SetVisAttributes(visual);
  crystal25_log->SetVisAttributes(visual);
  crystal27_log->SetVisAttributes(visual);

//------------------------------------------------
// Sensitive detectors
//------------------------------------------------

  crystal19i_log->SetSensitiveDetector(GetSensitiveDetector());
  crystal19o_log->SetSensitiveDetector(GetSensitiveDetector());
  crystal25_log->SetSensitiveDetector(GetSensitiveDetector());
  crystal27_log->SetSensitiveDetector(GetSensitiveDetector());

  // *** Build CsI Downstream Outer Ring ***

  G4double csi_ring3_length    = 22.185*inch;
  G4double z_position = csi_ring3_length/4. - crystal_length/2.;
  G4ThreeVector position(0,0,z_position);

  G4double angle = -90.0*degree - 9.58*degree + 20.3*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal25_log,"Crystal25",
                           csi_space4_log,false,0);

  angle -= 13.450*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal19o_log,"Crystal19",
                           csi_space4_log,false,1);

  angle -= 13.450*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal25_log,"Crystal25",
                           csi_space4_log,false,2);

  angle -= 13.300*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal25_log,"Crystal25",
                           csi_space4_log,false,3);

  angle -= 13.300*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal25_log,"Crystal25",
                           csi_space4_log,false,4);

  angle -= 13.300*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal25_log,"Crystal25",
                           csi_space4_log,false,5);

  angle -= 13.300*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal25_log,"Crystal25",
                           csi_space4_log,false,6);

  angle -= 13.300*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal25_log,"Crystal25",
                           csi_space4_log,false,7);

  angle -= 13.300*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal25_log,"Crystal25",
                           csi_space4_log,false,8);

  angle -= 13.300*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal25_log,"Crystal25",
                           csi_space4_log,false,9);

  angle -= 13.450*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal19o_log,"Crystal19",
                           csi_space4_log,false,10);

  angle -= 13.450*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal25_log,"Crystal25",
                           csi_space4_log,false,11);

  angle -= 13.300*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal25_log,"Crystal25",
                           csi_space4_log,false,12);

  angle -= 13.300*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal25_log,"Crystal25",
                           csi_space4_log,false,13);

  angle -= 13.300*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal25_log,"Crystal25",
                           csi_space4_log,false,14);

  angle -= 13.300*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal25_log,"Crystal25",
                           csi_space4_log,false,15);

  angle -= 13.300*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal25_log,"Crystal25",
                           csi_space4_log,false,16);

  angle -= 13.300*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal25_log,"Crystal25",
                           csi_space4_log,false,17);

  angle -= 13.300*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal25_log,"Crystal25",
                           csi_space4_log,false,18);

  angle -= 13.450*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal19o_log,"Crystal19",
                           csi_space4_log,false,19);

  angle -= 13.450*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal25_log,"Crystal25",
                           csi_space4_log,false,20);

  angle -= 13.300*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal25_log,"Crystal25",
                           csi_space4_log,false,21);

  angle -= 13.300*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal25_log,"Crystal25",
                           csi_space4_log,false,22);

  angle -= 13.300*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal25_log,"Crystal25",
                           csi_space4_log,false,23);

  angle -= 13.300*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal25_log,"Crystal25",
                           csi_space4_log,false,24);

  angle -= 13.300*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal25_log,"Crystal25",
                           csi_space4_log,false,25);

  angle -= 13.300*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal25_log,"Crystal25",
                           csi_space4_log,false,26);

  // *** Build CsI Downstream Inner Ring ***

  angle = -90.0*degree - 9.58*degree + 25.8*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal19i_log,"Crystal19",
                           csi_space3_log,false,0);

  angle -= 16.915*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal19i_log,"Crystal19",
                           csi_space3_log,false,1);

  angle -= 17.616*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27D",
                           csi_space3_log,false,2);

  angle -= 17.219*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal19i_log,"Crystal19",
                           csi_space3_log,false,3);

  angle -= 16.794*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal19i_log,"Crystal19",
                           csi_space3_log,false,4);

  angle -= 16.794*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal19i_log,"Crystal19",
                           csi_space3_log,false,5);

  angle -= 17.500*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27D",
                           csi_space3_log,false,6);

  angle -= 17.162*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal19i_log,"Crystal19",
                           csi_space3_log,false,7);

  angle -= 16.915*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal19i_log,"Crystal19",
                           csi_space3_log,false,8);

  angle -= 17.616*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27D",
                           csi_space3_log,false,9);

  angle -= 17.219*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal19i_log,"Crystal19",
                           csi_space3_log,false,10);

  angle -= 16.794*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal19i_log,"Crystal19",
                           csi_space3_log,false,11);

  angle -= 16.794*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal19i_log,"Crystal19",
                           csi_space3_log,false,12);

  angle -= 17.500*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27D",
                           csi_space3_log,false,13);

  angle -= 17.162*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal19i_log,"Crystal19",
                           csi_space3_log,false,14);

  angle -= 16.915*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal19i_log,"Crystal19",
                           csi_space3_log,false,15);

  angle -= 17.616*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27D",
                           csi_space3_log,false,16);

  angle -= 17.219*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal19i_log,"Crystal19",
                           csi_space3_log,false,17);

  angle -= 16.794*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal19i_log,"Crystal19",
                           csi_space3_log,false,18);

  angle -= 16.794*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal19i_log,"Crystal19",
                           csi_space3_log,false,19);

  angle -= 17.500*degree;
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(angle);

  new G4PVPlacement(g4rot,position,crystal27_log,"Crystal27D",
                           csi_space3_log,false,20);

}
