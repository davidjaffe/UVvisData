#include "globals.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"

#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4TwoVector.hh"

#include "G4PVReplica.hh"

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
#include "SiStripConstructor.hh"

#include "GeoConstructorMessenger.hh"

#include "G4SubtractionSolid.hh"

SiStripConstructor::~SiStripConstructor() {}

void SiStripConstructor::Init(void) {

  SetVisible(false);
  SetSensitiveDetector("wireChamberSD","segment");

  SetMessenger(new GeoConstructorMessenger(this));
}

G4LogicalVolume *SiStripConstructor::GetPiece(void) {


  G4VisAttributes *visual = new G4VisAttributes();
  visual->SetColor(0.5,0.1,0.1,1); // dark red
  visual->SetForceWireframe(true);


  G4double sistripG10PCBWidth                = 150.0*mm;
  G4double sistripG10PCBHeight               = 150.0*mm;
  G4double sistripG10PCBLength               = 1.6*mm;

  G4double sistripG10PCBFrontWindowWidth     = 65.0*mm;
  G4double sistripG10PCBFrontWindowHeight    = 65.0*mm;
  G4double sistripG10PCBFrontWindowLength    = 0.8*mm;

  G4double sistripG10PCBRearWindowWidth      = 85.0*mm;
  G4double sistripG10PCBRearWindowHeight     = 85.0*mm;
  G4double sistripG10PCBRearWindowLength     = 0.8*mm;

  G4double sistripSiliconWidth               = 61.44*mm;//64*mm;
  G4double sistripSiliconHeight              = 61.44*mm;//64*mm;
  G4double sistripSiliconLength              = 0.285*mm;
  G4double Nstrips = 48;

  G4double sistripSubstrateWidth             = 84*mm;
  G4double sistripSubstrateHeight            = 84*mm;
  G4double sistripSubstrateLength            = 1*mm;
  
  G4double sistripSubstrateWindowHeight      = 62*mm;
  G4double sistripSubstrateWindowWidth       = 62*mm;
  G4double sistripSubstrateWindowLength      = 1*mm;
  
  G4double sistripAluPlateWidth              = 270.5*mm;
  G4double sistripAluPlateHeight             = 130*mm;
  G4double sistripAluPlateLength             = 1.29*mm;
  
  G4double sistripAluPlateWindowWidth        = 66*mm;
  G4double sistripAluPlateWindowHeight       = 66*mm;
  G4double sistripAluPlateWindowLength       = 1.29*mm;


  G4Box* sistrip_alu_plate = new G4Box("sistrip_alu_plate_intermediate",
					       sistripAluPlateWidth/2.,
					       sistripAluPlateHeight/2,
					       sistripAluPlateLength/2.);
  G4Box* sistrip_alu_plate_window = new G4Box("sistrip_alu_plate_window",
					 sistripAluPlateWindowWidth/2.,
					 sistripAluPlateWindowHeight/2,
					 sistripAluPlateWindowLength/2.);
  // Al plate window is Aluminum, and then a window of Air is placed 
  // inside it to create an aluminum plate with a window
  G4LogicalVolume* sistrip_alu_plate_window_log
    = new G4LogicalVolume(sistrip_alu_plate_window,FindMaterial("G4_AIR"),"sistrip_alu_plate_window",0,0,0);
  G4LogicalVolume* sistrip_alu_plate_log
    = new G4LogicalVolume(sistrip_alu_plate,FindMaterial("G4_Al"),"sistrip_alu_plate",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),
                    sistrip_alu_plate_window_log,"sistrip_alu_plate_window",sistrip_alu_plate_log,false,0);
  
  
  G4double sistripPreampWidth                = 60*mm;
  G4double sistripPreampHeight               = 130*mm;
  G4double sistripPreampLength               = 1.6*mm;
  
  G4Box* sistrip_preamp = new G4Box("sistrip_preamp",
			       sistripPreampWidth/2.,
			       sistripPreampHeight/2,
			       sistripPreampLength/2.);
  G4LogicalVolume* sistrip_preamp_log
    = new G4LogicalVolume(sistrip_preamp,FindMaterial("G10"),"sistrip_preamp",0,0,0);


  G4double sistripPreampCoverWidth           = 60*mm;
  G4double sistripPreampCoverHeight          = 130*mm;
  G4double sistripPreampCoverLength          = 1.29*mm;

  G4Box* sistrip_preamp_cover = new G4Box("sistrip_preamp_cover",
				     sistripPreampCoverWidth/2.,
				     sistripPreampCoverHeight/2,
				     sistripPreampCoverLength/2.);
  G4LogicalVolume* sistrip_preamp_cover_log
    = new G4LogicalVolume(sistrip_preamp_cover,FindMaterial("G4_Al"),"sistrip_preamp_cover",0,0,0);


  G4double zBeetweenSiPCB                    = 1*mm;
  G4double zBeetweenSiPCBCovers              = 2*mm;

  G4double sistripEnvelopeWidth              = sistripAluPlateWidth;
  G4double sistripEnvelopeHeight             = sistripAluPlateWidth;
  G4double sistripEnvelopeLength             = sistripG10PCBLength*2+zBeetweenSiPCB+zBeetweenSiPCBCovers*2+sistripPreampCoverLength*2;

  G4double sistripSDEnvelopeWidth            = sistripG10PCBWidth;
  G4double sistripSDEnvelopeHeight           = sistripG10PCBHeight;
  G4double sistripSDEnvelopeLength           = sistripG10PCBLength*2+zBeetweenSiPCB;

  

    
  
  
  
  G4Box* sistrip_envelope = new G4Box("sistrip_envelope",
                                 sistripEnvelopeWidth/2.,
                                 sistripEnvelopeHeight/2.,
				 sistripEnvelopeLength/2.);
  G4LogicalVolume* sistrip_envelope_log
    = new G4LogicalVolume(sistrip_envelope,FindMaterial("G4_AIR"),"sistrip_envelope",0,0,0);
  
    G4Box* sistrip_sd_envelope = new G4Box("sistrip_sd_envelope",
                                 sistripSDEnvelopeWidth/2.,
                                 sistripSDEnvelopeHeight/2.,
				 sistripSDEnvelopeLength/2.);
  G4LogicalVolume* sistrip_sd_envelope_log
    = new G4LogicalVolume(sistrip_sd_envelope,FindMaterial("G4_AIR"),"sistrip_envelope",0,0,0);

  G4Box* sistrip_silicon = new G4Box("sistrip_silicon",
				sistripSiliconWidth/2.,
				sistripSiliconHeight/2.,
				sistripSiliconLength/2.);
  G4LogicalVolume* sistrip_silicon_log
    = new G4LogicalVolume(sistrip_silicon,FindMaterial("G4_Si"),"sistrip_silicon",0,0,0);
  
  
  G4Box* sistrip_sistrip = new G4Box("sistrip_sistrip",
				sistripSiliconWidth/2./Nstrips,
				sistripSiliconHeight/2.,
				sistripSiliconLength/2.);
  G4LogicalVolume* sistrip_sistrip_log
    = new G4LogicalVolume(sistrip_sistrip,FindMaterial("G4_Si"),"sistrip_sistrip",0,0,0);
  

  G4Box* sistrip_substrate_intermediate = new G4Box("sistrip_substrate_intermediate",
					       sistripSubstrateWidth/2.,
					       sistripSubstrateHeight/2.,
					       sistripSubstrateLength/2.);
  G4Box* sistrip_substrate_window = new G4Box("sistrip_substrate_window",
					 sistripSubstrateWindowWidth/2.,
					 sistripSubstrateWindowHeight/2.,
					 sistripSubstrateWindowLength/2.);
  G4SubtractionSolid* sistrip_substrate =
    new G4SubtractionSolid("sistrip_substrate", sistrip_substrate_intermediate, sistrip_substrate_window);
  
  G4LogicalVolume* sistrip_substrate_log
    = new G4LogicalVolume(sistrip_substrate,FindMaterial("Ceramic"),"sistrip_substrate",0,0,0);
  
  G4Box* sistrip_g10_pcb = new G4Box("sistrip_g10_pcb",
                                sistripG10PCBWidth/2.,
                                sistripG10PCBHeight/2.,
                                sistripG10PCBLength/2.);  
  G4Box* sistrip_g10_pcb_front_window = new G4Box("sistrip_g10_pcb_front_window",
					     sistripG10PCBFrontWindowWidth/2.,
					     sistripG10PCBFrontWindowHeight/2.,
					     sistripG10PCBFrontWindowLength/2.);
  G4Box* sistrip_g10_pcb_rear_window = new G4Box("sistrip_g10_pcb_rear_window",
					    sistripG10PCBRearWindowWidth/2.,
					    sistripG10PCBRearWindowHeight/2.,
					    sistripG10PCBRearWindowLength/2.); 
  
  
  // ---------------------------------------------------------------------------------------------------
  // --   Making a G10 PCB with the help of boolean solids  --
  // ---------------------------------------------------------
   // ----------   Aleksey Sher, November, 2010 ---------------
  
  G4RotationMatrix* g4rot = new G4RotationMatrix();
  
  G4ThreeVector zTrans(0, 0, sistripG10PCBFrontWindowLength/(-2.));
  
  G4SubtractionSolid* sistrip_g10_pcb_intermediate =
    new G4SubtractionSolid("sistrip_g10_pcb_final", sistrip_g10_pcb, sistrip_g10_pcb_front_window,0,zTrans);
  
  G4ThreeVector zTrans1(0, 0, sistripG10PCBRearWindowLength/(2.));
  
  G4SubtractionSolid* sistrip_g10_pcb_final =
    new G4SubtractionSolid("sistrip_g10_pcb_final", sistrip_g10_pcb_intermediate, sistrip_g10_pcb_rear_window,0,zTrans1);
  
  
  G4LogicalVolume* sistrip_g10_pcb_log
    = new G4LogicalVolume(sistrip_g10_pcb_final,FindMaterial("G10"),"sistrip_g10_pcb",0,0,0);
  
  // -------------         G10 PCB construction is complete ---------------------------------------------
  
  
  
  
  
  // -------------- P L A C E M E N T S -----------------------------------------------------------------

  


  
  G4double z_pcb = zBeetweenSiPCB/2.+sistripG10PCBLength/2.;
  
  G4double z = - z_pcb;
  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    sistrip_g10_pcb_log,"sistripG10PCB",sistrip_sd_envelope_log,false,0);
  z = z_pcb;
  g4rot->rotateX(180.*degree);  
  new G4PVPlacement(g4rot,G4ThreeVector(0.,0.,z),
                    sistrip_g10_pcb_log,"sistripG10PCB",sistrip_sd_envelope_log,false,1);
  

  z = -z_pcb - sistripSiliconLength/2.;
  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    sistrip_silicon_log,"sistripSilicon",sistrip_sd_envelope_log,false,0);

  new G4PVReplica("sistripSistrip",sistrip_sistrip_log,sistrip_silicon_log,kXAxis,48,1.28*mm,0);
  
  
  z = -z;
  g4rot->rotateZ(90.*degree); 
  new G4PVPlacement(g4rot,G4ThreeVector(0.,0.,z),
                    sistrip_silicon_log,"sistripSilicon",sistrip_sd_envelope_log,false,1);
  


  z = z_pcb - sistripSubstrateLength/2.;
  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    sistrip_substrate_log,"sistripSubstrate",sistrip_sd_envelope_log,false,0);
  z = -z;
  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    sistrip_substrate_log,"sistripSubstrate",sistrip_sd_envelope_log,false,1);


  new G4PVPlacement(0,G4ThreeVector(0.,0.,0),
                    sistrip_sd_envelope_log,"sistripSDEnvelope",sistrip_envelope_log,false,0);
  

  z= z_pcb + sistripG10PCBLength/2.+ zBeetweenSiPCBCovers+sistripAluPlateLength/2;
  new G4PVPlacement(0,G4ThreeVector(0.,0.,z),
                    sistrip_alu_plate_log,"sistripAluPlate",sistrip_envelope_log,false,0);
  
  g4rot = new G4RotationMatrix();
  g4rot->rotateZ(90.*degree); 
  z = -z;
  new G4PVPlacement(g4rot,G4ThreeVector(0.,0.,z),
                    sistrip_alu_plate_log,"sistripAluPlate",sistrip_envelope_log,false,1);

  

  z = z_pcb;
  G4double x = sistripAluPlateWidth/2.-sistripPreampWidth/2.;
  new G4PVPlacement(0,G4ThreeVector(x,0.,z),
                    sistrip_preamp_log,"sistripPreamp",sistrip_envelope_log,false,0); 
  x = -x;
  new G4PVPlacement(0,G4ThreeVector(x,0.,z),
                    sistrip_preamp_log,"sistripPreamp",sistrip_envelope_log,false,1);  
  
  G4double y = x;
  z = -z;
  new G4PVPlacement(g4rot,G4ThreeVector(0.,y,z),
                    sistrip_preamp_log,"sistripPreamp",sistrip_envelope_log,false,2); 
  y = -y;
  new G4PVPlacement(g4rot,G4ThreeVector(0.,y,z),
                    sistrip_preamp_log,"sistripPreamp",sistrip_envelope_log,false,3); 
  
  

  z = z_pcb - sistripG10PCBLength/2.- zBeetweenSiPCBCovers - sistripAluPlateLength/2;
  new G4PVPlacement(0,G4ThreeVector(x,0.,z),
                    sistrip_preamp_cover_log,"sistripPreampCover",sistrip_envelope_log,false,0); 
  x = -x;
  new G4PVPlacement(0,G4ThreeVector(x,0.,z),
                    sistrip_preamp_cover_log,"sistripPreampCover",sistrip_envelope_log,false,1); 
  z = -z;
  
  new G4PVPlacement(g4rot,G4ThreeVector(0.,y,z),
                    sistrip_preamp_cover_log,"sistripPreampCover",sistrip_envelope_log,false,2); 
  y = -y;
  new G4PVPlacement(g4rot,G4ThreeVector(0.,y,z),
		    sistrip_preamp_cover_log,"sistripPreampCover",sistrip_envelope_log,false,3); 
  
  
  sistrip_sistrip_log->SetSensitiveDetector(GetSensitiveDetector());
  
  sistrip_silicon_log->SetVisAttributes(visual);
  sistrip_g10_pcb_log->SetVisAttributes(visual);
  sistrip_substrate_log->SetVisAttributes(visual);
  sistrip_alu_plate_log->SetVisAttributes(visual);
  sistrip_preamp_log->SetVisAttributes(visual);      
  sistrip_preamp_cover_log->SetVisAttributes(visual);      
  
    return sistrip_envelope_log;
    
}
