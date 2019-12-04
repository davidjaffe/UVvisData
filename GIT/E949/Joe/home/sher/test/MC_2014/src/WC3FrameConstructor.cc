#include "globals.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"

#include "G4DisplacedSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4TwoVector.hh"

#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"

#include "G4UserLimits.hh"

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "Materials.hh"

#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4GeometryManager.hh"

#include "G4VisAttributes.hh"

#include "GeoConstructor.hh"
#include "WC3FrameConstructor.hh"

#include "GeoConstructorMessenger.hh"

#include "G4SubtractionSolid.hh"

WC3FrameConstructor::~WC3FrameConstructor() {}

void WC3FrameConstructor::Init(void) {

  SetVisible(false);

  SetMessenger(new GeoConstructorMessenger(this));
}

G4LogicalVolume *WC3FrameConstructor::GetPiece(void) {

  G4double wc3G10BoardWidth  = 73.66*mm;
  G4double wc3G10BoardHeight = 71.8*mm;
  G4double wc3G10BoardLength = 1.75*mm;

  G4double wc3G10PieceOffset = 6.41*mm;

  G4double wc3G10PieceWidth   =  1.6*mm;
  G4double wc3G10PieceHeight  = 71.8*mm;
  G4double wc3G10Piece1Length = 26.0*mm;
  G4double wc3G10Piece2Length = 22.0*mm;
  G4double wc3G10Piece3Length = 18.0*mm;

  G4double wc3G10RingRmin   = 314.0*mm/2.;
  G4double wc3G10RingRmax   = wc3G10RingRmin+10.0*mm;
  G4double wc3G10RingLength = 2.0*mm;

  G4double wc3FrameRmin   = wc3G10RingRmin;
  G4double wc3FrameRmax   = wc3G10RingRmax - wc3G10PieceOffset + wc3G10BoardWidth;
  G4double wc3FrameLength = wc3G10Piece1Length + wc3G10RingLength;

  G4double z_1st = wc3FrameLength/2-wc3G10RingLength/2.;
  G4double z_2nd = z_1st - 4*mm;
  G4double z_3rd = z_2nd - 4*mm;

  G4double wc3G10PieceRadius = wc3G10RingRmax-wc3G10PieceOffset-wc3G10PieceWidth/2.;
  G4double wc3G10BoardRadius = wc3G10RingRmax-wc3G10PieceOffset+wc3G10BoardWidth/2.;

  G4Tubs* wc3_frame_intermediate = new G4Tubs("wc3_frame_intermediate",
                                 wc3FrameRmin,
                                 1.05*wc3FrameRmax,
                                 wc3FrameLength/2.,
                                 0.*deg,360.*deg);

  G4Tubs* wc3_frame_delete = new G4Tubs("wc3_frame_delete",
                                 wc3FrameRmin*1.1,
                                 1.05*wc3FrameRmax,
                                 wc3FrameLength/2.-2*mm,
                                 0.*deg,360.*deg);

  G4ThreeVector zTrans(0, 0, 4*mm);
  
  G4SubtractionSolid* wc3_frame =
    new G4SubtractionSolid("wc3_frame", wc3_frame_intermediate,wc3_frame_delete,0,zTrans);
  
  G4LogicalVolume* wc3_frame_log
    = new G4LogicalVolume(wc3_frame,FindMaterial("G4_AIR"),"wc3_frame",0,0,0);

  G4VisAttributes *visual0 = new G4VisAttributes();
  visual0->SetColor(0.5,0.1,0.1,1); // dark red
  visual0->SetForceWireframe(true);

  wc3_frame_log->SetVisAttributes(visual0);

  G4Tubs* wc3_ring = new G4Tubs("wc3_ring",
                                 wc3G10RingRmin,
                                 wc3G10RingRmax,
                                 wc3G10RingLength/2.,
                                  35.*degree,
                                 110.*degree);

  G4LogicalVolume* wc3_ring_log
    = new G4LogicalVolume(wc3_ring,FindMaterial("G10"),"wc3_ring",0,0,0);

  G4VisAttributes *visual = new G4VisAttributes();
  visual->SetColor(1.,0.,0.);
  visual->SetForceSolid(true);

  wc3_ring_log->SetVisAttributes(visual);

  new G4PVPlacement(0    ,G4ThreeVector(0.,0., z_1st),  //4*mm
                    wc3_ring_log,"wc3_ring",wc3_frame_log,false,0);

  G4RotationMatrix* g4rot = new G4RotationMatrix();
  g4rot->rotateZ(-120.*degree);

  new G4PVPlacement(g4rot,G4ThreeVector(0.,0., z_2nd),
                    wc3_ring_log,"wc3_ring",wc3_frame_log,false,1);

  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(120.*degree);
  
  new G4PVPlacement(g4rot,G4ThreeVector(0.,0.,z_3rd),
                    wc3_ring_log,"wc3_ring",wc3_frame_log,false,2);

  G4Box* wc3_piece1 = new G4Box("wc3_piece1",
                                wc3G10PieceWidth/2.,
                                wc3G10PieceHeight/2.,
                                wc3G10Piece1Length/2.);

  G4Box* wc3_piece2 = new G4Box("wc3_piece2",
                                wc3G10PieceWidth/2.,
                                wc3G10PieceHeight/2.,
                                wc3G10Piece2Length/2.);

  G4Box* wc3_piece3 = new G4Box("wc3_piece3",
                                wc3G10PieceWidth/2.,
                                wc3G10PieceHeight/2.,
                                wc3G10Piece3Length/2.);

  G4Box* wc3_board = new G4Box("wc3_board",
                                wc3G10BoardWidth/2.,
                                wc3G10BoardHeight/2.,
                                wc3G10BoardLength/2.);

  G4ThreeVector pieceVector(wc3G10PieceRadius,0,0);

  G4DisplacedSolid* wc3_piece1_displaced =
                            new G4DisplacedSolid("wc3_piece1",
                                                 wc3_piece1,
                                                 0,
                                                 pieceVector);

  G4DisplacedSolid* wc3_piece2_displaced =
                            new G4DisplacedSolid("wc3_piece2",
                                                 wc3_piece2,
                                                 0,
                                                 pieceVector);

  G4DisplacedSolid* wc3_piece3_displaced =
                            new G4DisplacedSolid("wc3_piece3",
                                                 wc3_piece3,
                                                 0,
                                                 pieceVector);

  G4ThreeVector boardVector(wc3G10BoardRadius,0,0);

  G4DisplacedSolid* wc3_board_displaced =
                            new G4DisplacedSolid("wc3_board",
                                                 wc3_board,
                                                 0,
                                                 boardVector);

  G4LogicalVolume* wc3_piece1_log =
    new G4LogicalVolume(wc3_piece1_displaced,FindMaterial("G10"),
                                             "wc3_piece1",0,0,0);

  G4LogicalVolume* wc3_piece2_log =
    new G4LogicalVolume(wc3_piece2_displaced,FindMaterial("G10"),
                                             "wc3_piece2",0,0,0);

  G4LogicalVolume* wc3_piece3_log =
    new G4LogicalVolume(wc3_piece3_displaced,FindMaterial("G10"),
                                             "wc3_piece3",0,0,0);

  G4LogicalVolume* wc3_board_log =
    new G4LogicalVolume(wc3_board_displaced,FindMaterial("G10"),
                                            "wc3_board",0,0,0);

  G4double z  = z_1st - wc3G10RingLength/2. - wc3G10Piece1Length/2.;
  G4double zz = -wc3FrameLength/2. + wc3G10BoardLength/2.;

  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(-50.0*degree);

  new G4PVPlacement(g4rot,G4ThreeVector(0,0,z),wc3_piece1_log,"wc3_piece1",
                    wc3_frame_log,false,0);
  new G4PVPlacement(g4rot,G4ThreeVector(0,0,zz),wc3_board_log,"wc3_board",
                    wc3_frame_log,false,0);

  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(-90.0*degree);

  new G4PVPlacement(g4rot,G4ThreeVector(0,0,z),wc3_piece1_log,"wc3_piece1",
                    wc3_frame_log,false,1);
  new G4PVPlacement(g4rot,G4ThreeVector(0,0,zz),wc3_board_log,"wc3_board",
                    wc3_frame_log,false,1);

  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(-130.0*degree);

  new G4PVPlacement(g4rot,G4ThreeVector(0,0,z),wc3_piece1_log,"wc3_piece1",
                    wc3_frame_log,false,2);
  new G4PVPlacement(g4rot,G4ThreeVector(0,0,zz),wc3_board_log,"wc3_board",
                    wc3_frame_log,false,2);

  z = z_2nd - wc3G10RingLength/2. - wc3G10Piece2Length/2.;

  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(-170.0*degree);

  new G4PVPlacement(g4rot,G4ThreeVector(0,0,z),wc3_piece2_log,"wc3_piece2",
                    wc3_frame_log,false,0);
  new G4PVPlacement(g4rot,G4ThreeVector(0,0,zz),wc3_board_log,"wc3_board",
                    wc3_frame_log,false,3);

  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(-210.0*degree);

  new G4PVPlacement(g4rot,G4ThreeVector(0,0,z),wc3_piece2_log,"wc3_piece2",
                    wc3_frame_log,false,1);
  new G4PVPlacement(g4rot,G4ThreeVector(0,0,zz),wc3_board_log,"wc3_board",
                    wc3_frame_log,false,4);

  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(-250.0*degree);

  new G4PVPlacement(g4rot,G4ThreeVector(0,0,z),wc3_piece2_log,"wc3_piece2",
                    wc3_frame_log,false,2);
  new G4PVPlacement(g4rot,G4ThreeVector(0,0,zz),wc3_board_log,"wc3_board",
                    wc3_frame_log,false,5);

  z = z_3rd - wc3G10RingLength/2. - wc3G10Piece3Length/2.;

  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(-290.0*degree);

  new G4PVPlacement(g4rot,G4ThreeVector(0,0,z),wc3_piece3_log,"wc3_piece3",
                    wc3_frame_log,false,0);
  new G4PVPlacement(g4rot,G4ThreeVector(0,0,zz),wc3_board_log,"wc3_board",
                    wc3_frame_log,false,6);

  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(-330.0*degree);

  new G4PVPlacement(g4rot,G4ThreeVector(0,0,z),wc3_piece3_log,"wc3_piece3",
                    wc3_frame_log,false,1);
  new G4PVPlacement(g4rot,G4ThreeVector(0,0,zz),wc3_board_log,"wc3_board",
                    wc3_frame_log,false,7);

  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(-10.0*degree);

  new G4PVPlacement(g4rot,G4ThreeVector(0,0,z),wc3_piece3_log,"wc3_piece3",
                    wc3_frame_log,false,2);
  new G4PVPlacement(g4rot,G4ThreeVector(0,0,zz),wc3_board_log,"wc3_board",
                    wc3_frame_log,false,8);

  return wc3_frame_log;

}
