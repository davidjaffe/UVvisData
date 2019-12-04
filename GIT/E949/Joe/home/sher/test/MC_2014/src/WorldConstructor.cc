#include "globals.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"

#include "G4RotationMatrix.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VPhysicalVolume.hh"

#include "G4VisAttributes.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"

#include "SensitiveDetectorFactory.hh"

#include "GeoConstructor.hh"
#include "WorldConstructor.hh"
#include "WorldMessenger.hh"

// To debug a particular volume define
// DETECTOR_CONSTRUCTOR, DETECTOR_NAME
// and include the required file.
#undef DETECTOR_CONSTRUCTOR
#ifndef DETECTOR_CONSTRUCTOR
#include "CsIConstructor.hh"
#define DETECTOR_CONSTRUCTOR CsIConstructor
#define DETECTOR_NAME        "CsI"
#include "WC3Constructor.hh"
//#define DETECTOR_CONSTRUCTOR WC3Constructor
//#define DETECTOR_NAME        "WC3"

#include "WC3FrameConstructor.hh"
//#define DETECTOR_CONSTRUCTOR WC3FrameConstructor
//#define DETECTOR_NAME        "WC3Frame"
//#endif

#include "BeamWCConstructor.hh"
//#define DETECTOR_CONSTRUCTOR BeamWCConstructor
//#define DETECTOR_NAME        "BeamWC"

#include "SiStripConstructor.hh"
//#define DETECTOR_CONSTRUCTOR SiStripConstructor
//#define DETECTOR_NAME        "SiStrip"
#endif

#include "ScintConstructor.hh"


#include "BeamLineConstructor.hh"

void WorldConstructor::Init(void)
{
  SetMessenger(new WorldMessenger(this));
  AddConstructor(new CsIConstructor("CsI",this));
  AddConstructor(new BeamLineConstructor("BeamLine",this));
  AddConstructor(new SiStripConstructor("SiStrip",this));
  AddConstructor(new BeamWCConstructor("BeamWC",this));
  AddConstructor(new WC3Constructor("WC3",this));
  AddConstructor(new WC3FrameConstructor("WC3Frame",this));
  AddConstructor(new ScintConstructor("B1",this));
  AddConstructor(new ScintConstructor("B2",this));
  AddConstructor(new ScintConstructor("SW",this));
  AddConstructor(new ScintConstructor("Target",this));
  AddConstructor(new ScintConstructor("T1",this));
  AddConstructor(new ScintConstructor("T2",this));
  AddConstructor(new ScintConstructor("V2",this));
  AddConstructor(new ScintConstructor("V3",this));
  
  fWidth =  1.0*m;
  fHeight = 1.0*m;
  fLength = 3.0*m;
  
  SetUseStandAloneModule(false);
  SetStandAloneModuleName("NoNameSpecified");
}

WorldConstructor::~WorldConstructor() {}

G4LogicalVolume *WorldConstructor::GetPiece(void)
{
  // Make sure the region-store is empty before
  // we start constructing the new detector.
  G4RegionStore* regionStore = G4RegionStore::GetInstance();
  
  // Create a region outside of the detector to define cuts.
  G4Region* Region = regionStore->GetRegion("Region",false);
  if (Region) regionStore->DeRegister(Region);
  Region = new G4Region("Region");
//   G4ProductionCuts* Cuts = new G4ProductionCuts();
//   Cuts->SetProductionCut(10*mm);
//   Region->SetProductionCuts(Cuts);
  
   G4cout << "##############################################" << G4endl;
  G4cout << "#     CONSTRUCT THE PIENU WORLD GEOMETRY     #" << G4endl;
  G4cout << "##############################################" << G4endl;


  G4double  BeamLine_lngth = 86.484*cm;
  G4double  BeamWC_lngth = 108.282*mm;
  G4double  B1_lngth =  6.604*mm;
  //G4double  B1_lngth =  67*mm;     // experimental to find out where to put Zr foil
  G4double  B2_lngth = 3.07*mm;
  G4double  SiStrip_lngth = 10.78*mm;
  G4double  Target_lngth = 8.05*mm;
  G4double  T1_lngth = 3.04*mm;
  //G4double T1_lngth = 18*mm;       // experimental to find out where to put Zr foil
  G4double  WC3_lngth = 24*mm;
  G4double  T2_lngth = 6.6*mm;
  G4double  NaI_lngth =  48.6*cm; //48*cm;
  G4double  CsI_crystal_lngth = 25*cm;
  G4double  V2_lngth = 6.35*mm;
  G4double  V3_lngth = 6.35*mm;
  //G4double  Degrader_lngth = 3mm;

  G4double  Target_z = 0*cm;


  G4double  BeamLine_BeamWC_gap = 2*mm-1.2076*mm; //2.591*mm
  G4double  BeamWC_B1_gap = 1.887*mm;  // 2*mm
  G4double  B1_B2_gap = 4.175*mm; //4.175*mm
//  G4double  B2_SiStrip1_gap = 1.5*mm-0.19075*mm;   // 1*mm
  G4double  B2_SiStrip1_gap = 1.*mm;   // 1*mm
  G4double  SiStrip1_SiStrip2_gap = 1.*mm; // 1*mm
//  G4double  SiStrip2_Target_gap = 1.55*mm+0.19075*mm;  // 0.959  
  G4double  SiStrip2_Target_gap = 0.9*mm;  // 0.959  



//  G4double  Target_SiStrip3_gap =  1.6*mm+0.1675*mm+1*mm;  // 2.52*mm
  G4double  Target_SiStrip3_gap =  2.525*mm+0.1425*mm;  // 2.52*mm -0.1*mm 
  G4double  SiStrip3_T1_gap = 2*mm-0.1675*mm-1*mm+0.2375-0.1425*mm; //1*mm
  G4double  T1_WC3_gap = 24*mm-0.2696*mm-2.25*mm+1.96*mm-1*mm;// +2.34*mm; 
// 
  G4double  WC3_T2_gap = 1*mm;
  G4double  T2_NaI_gap = 1*mm;

  G4double  CsIedge_BinaFace_gap = 12.4968*cm;


  G4double  SiStrip3_z = Target_z+Target_lngth/2.+Target_SiStrip3_gap+SiStrip_lngth/2.;
  
  G4double  T1_z  = SiStrip3_z+SiStrip_lngth/2.+SiStrip3_T1_gap+T1_lngth/2.;

  G4double  WC3_z  = T1_z+T1_lngth/2.+T1_WC3_gap+WC3_lngth/2.;

  G4double  SiStrip2_z = Target_z-Target_lngth/2.-SiStrip2_Target_gap-SiStrip_lngth/2.;
  
  G4double  SiStrip1_z = SiStrip2_z - SiStrip_lngth/2.-SiStrip1_SiStrip2_gap-SiStrip_lngth/2.;

  G4double  B2_z = SiStrip1_z - SiStrip_lngth/2.-B2_SiStrip1_gap - B2_lngth/2.;
 
  G4double  B1_z = B2_z - B2_lngth/2.-B1_B2_gap - B1_lngth/2.;

  G4double  BeamWC_z = B1_z - B1_lngth/2. - BeamWC_B1_gap - BeamWC_lngth/2.;
  
  G4double  T2_z = WC3_z + WC3_lngth/2.+WC3_T2_gap+T2_lngth/2.;

  G4double V2_z = WC3_z - WC3_lngth/2.-V2_lngth/2.-2*mm;
  G4double V3_z = WC3_z -4*V3_lngth/2.-2*mm;

  G4double  NaI_z = T2_z +T2_lngth/2.+T2_NaI_gap + NaI_lngth/2;
  
  G4double  CsI_z = NaI_z-NaI_lngth/2.+(CsI_crystal_lngth-CsIedge_BinaFace_gap);
  

  G4double  BeamLine_z = BeamWC_z - BeamWC_lngth/2. - BeamLine_BeamWC_gap - BeamLine_lngth/2.;
  
  G4cout << " B1 " << B1_z << " B2 " << B2_z << " S1 " << SiStrip1_z << " S2 " << SiStrip2_z << " Tg " << Target_z << " S3 " << SiStrip3_z << " T1 " <<  T1_z  << " WC3 " << WC3_z << " T2 " << T2_z << " Bina " << NaI_z << G4endl;

  G4ThreeVector detectorCenter;
  //    detectorCenter.setZ(15*cm);
  detectorCenter.setZ(130*cm);
  
  G4LogicalVolume *logWorld
    = new G4LogicalVolume(new G4Box(GetName(),
				    fWidth/2,   // horizontal
				    fHeight/2,  // vertical
				    fLength/2), // beam direction
			  FindMaterial("G4_AIR"),
			  GetName());
  
  logWorld->SetVisAttributes(G4VisAttributes::Invisible);
  
  //------------------------------
  // Stand alone sections
  //------------------------------
  
    if (fUseStandAloneModule) {
        GeoConstructor* theContructor = 
                               Find<GeoConstructor>(fStandAloneModuleName);
        if (!theContructor) {
           G4cout << "WARNING Name of standalone module not known: " 
                     << fStandAloneModuleName << G4endl;
        }

        assert(theContructor); //should throw exception here,really

        G4LogicalVolume* theVolume = theContructor->GetPiece();
        new G4PVPlacement(0,        // rotation
                          G4ThreeVector(0,0,0),
                          theVolume, // logical volume
                          theContructor->GetName(), // name
                          logWorld, // mother  volume
                          false,    // no boolean operations
                          0);       // not a copy.
        return logWorld;
    }

    //------------------------------ 
    // Detector --- Construct the toplevel detector volume.
    //------------------------------ 
    
     GeoConstructor& CsI = Get<GeoConstructor>("CsI");
    
     G4LogicalVolume* logCsI = CsI.GetPiece(); // logical volume
     detectorCenter.setZ(CsI_z);
         new G4PVPlacement(0,                   // rotation
                           detectorCenter,      // position
                     logCsI,         // logical volume
                   CsI.GetName(),  // name
                   logWorld,            // mother  volume
                     false,               // no boolean operations
                     0);
    
         Region->AddRootLogicalVolume(logCsI);
    

    GeoConstructor& SiStrip = Get<GeoConstructor>("SiStrip");

    G4LogicalVolume* logSiStrip = SiStrip.GetPiece(); // logical volume

    //    detectorCenter.setZ(-20*cm);
    detectorCenter.setZ(SiStrip1_z);
    
    new G4PVPlacement(0,                   // rotation
		      detectorCenter,      // position
                      logSiStrip,         // logical volume
                      "SS1",  // name
                      logWorld,            // mother  volume
                      false,               // no boolean operations
                      1);
       detectorCenter.setZ(SiStrip2_z); 
       new G4PVPlacement(0,                   // rotation
		      detectorCenter,      // position
                      logSiStrip,         // logical volume
                      "SS2",  // name
                      logWorld,            // mother  volume
                      false,               // no boolean operations
                      2);
       detectorCenter.setZ(SiStrip3_z); 
       new G4PVPlacement(0,                   // rotation
		      detectorCenter,      // position
                      logSiStrip,         // logical volume
                      "SS3",  // name
                      logWorld,            // mother  volume
                      false,               // no boolean operations
                      3);       
    Region->AddRootLogicalVolume(logSiStrip);


    
    GeoConstructor& BeamLine = Get<GeoConstructor>("BeamLine");

    G4LogicalVolume* logBeamLine = BeamLine.GetPiece(); // logical volume

    detectorCenter.setZ(BeamLine_z);
    new G4PVPlacement(0,                   // rotation
		      detectorCenter,      // position
                      logBeamLine,         // logical volume
                      BeamLine.GetName(),  // name
                      logWorld,            // mother  volume
                      false,               // no boolean operations
                      0);
    
    Region->AddRootLogicalVolume(logBeamLine);


    
     GeoConstructor& BeamWC = Get<GeoConstructor>("BeamWC");

     G4LogicalVolume* logBeamWC = BeamWC.GetPiece(); // logical volume

     detectorCenter.setZ(BeamWC_z);
     new G4PVPlacement(0,                   // rotation
 		      detectorCenter,      // position
                       logBeamWC,         // logical volume
                       BeamWC.GetName(),  // name
                       logWorld,            // mother  volume
                       false,               // no boolean operations
                       0);
    
     Region->AddRootLogicalVolume(logBeamWC);

    
    GeoConstructor& WC3 = Get<GeoConstructor>("WC3");

    G4LogicalVolume* logWC3 = WC3.GetPiece(); // logical volume

    detectorCenter.setZ(WC3_z);
    new G4PVPlacement(0,                   // rotation
		      detectorCenter,      // position
                      logWC3,         // logical volume
                      WC3.GetName(),  // name
                      logWorld,            // mother  volume
                      false,               // no boolean operations
                      0);
    
    Region->AddRootLogicalVolume(logWC3);

    GeoConstructor& WC3Frame = Get<GeoConstructor>("WC3Frame");

    G4LogicalVolume* logWC3Frame = WC3Frame.GetPiece(); // logical volume

    detectorCenter.setZ(WC3_z-10*mm);
    new G4PVPlacement(0,                   // rotation
		      detectorCenter,      // position
                      logWC3Frame,         // logical volume
                      WC3Frame.GetName(),  // name
                      logWorld,            // mother  volume
                      false,               // no boolean operations
                      0);
    
    Region->AddRootLogicalVolume(logWC3Frame);

    

    ScintConstructor& B1 = Get<ScintConstructor>("B1");
    B1.SetWidth(100.0*mm);
    B1.SetHeight(100.0*mm);
    //B1.SetLength(6.604*mm);
    B1.SetLength(B1_lngth);
    G4ThreeVector position(0.,0., B1_z);
    B1.SetPosition(position);

    G4LogicalVolume* logB1 = B1.GetPiece();

    G4VisAttributes* visual = new G4VisAttributes();
    visual->SetColor(0.2,0.5,0.2,1); // green
    visual->SetForceWireframe(true);
    logB1->SetVisAttributes(visual);

    G4ThreeVector B1Position = B1.GetPosition();

    new G4PVPlacement(0,               // rotation
                      B1Position,      // position
                      logB1,           // logical volume
                      B1.GetName(),    // name
                      logWorld,        // mother  volume
                      false,           // no boolean operations
                      1);

    Region->AddRootLogicalVolume(logB1);
    
    
    
    ScintConstructor& B2 = Get<ScintConstructor>("B2");

    B2.SetWidth(45.0*mm);
    B2.SetHeight(45.0*mm);
    B2.SetLength(3.07*mm);
    position.setZ(B2_z);
    B2.SetPosition(position);

    G4LogicalVolume* logB2 = B2.GetPiece();

    logB2->SetVisAttributes(visual);

    G4ThreeVector B2Position = B2.GetPosition();

    new G4PVPlacement(0,               // rotation
                      B2Position,      // position
                      logB2,           // logical volume
                      B2.GetName(),    // name
                      logWorld,        // mother  volume
                      false,           // no boolean operations
                      2);

    Region->AddRootLogicalVolume(logB2);



    ScintConstructor& TA = Get<ScintConstructor>("Target");

    TA.SetWidth(70.0*mm);
    TA.SetHeight(70.0*mm);
    TA.SetLength(8.05*mm);
    position.setZ(Target_z);
    TA.SetPosition(position);

    G4LogicalVolume* logTarget = TA.GetPiece();

    visual->SetColor(0.5,0.5,0.5,1); // white
    logTarget->SetVisAttributes(visual);

    //logTarget->SetUserLimits(new G4UserLimits(0.1*mm));

    G4ThreeVector TargetPosition = TA.GetPosition();

    G4RotationMatrix* g4rot = new G4RotationMatrix();
    g4rot->rotateZ(45.0*degree);

    new G4PVPlacement(g4rot,           // rotation
                      TargetPosition,  // position
                      logTarget,       // logical volume
                      TA.GetName(),    // name
                      logWorld,        // mother  volume
                      false,           // no boolean operations
                      3);

    Region->AddRootLogicalVolume(logTarget);

    ScintConstructor& T1 = Get<ScintConstructor>("T1");

    T1.SetWidth(80.0*mm);
    T1.SetHeight(80.0*mm);
    //T1.SetLength(3.04*mm);
    T1.SetLength(T1_lngth);
    position.setZ(T1_z);
    T1.SetPosition(position);

    G4LogicalVolume* logT1 = T1.GetPiece();

    visual->SetColor(0.5,0.5,1.0,0); // blue-green
    logT1->SetVisAttributes(visual);

    G4ThreeVector T1Position = T1.GetPosition();

    new G4PVPlacement(g4rot,           // rotation
                      T1Position,      // position
                      logT1,           // logical volume
                      T1.GetName(),    // name
                      logWorld,        // mother  volume
                      false,           // no boolean operations
                      4);

    Region->AddRootLogicalVolume(logT1);

    ScintConstructor& T2 = Get<ScintConstructor>("T2");

    T2.SetType("Tubs");
    T2.SetRmin(  0.0 *mm);
    T2.SetRmax(171.45*mm);
    T2.SetLength(6.6*mm);
    position.setZ(T2_z);
    T2.SetPosition(position);

    G4LogicalVolume* logT2 = T2.GetPiece();

    visual->SetColor(0.5,0.5,1.0,0); // blue-green
    logT2->SetVisAttributes(visual);

    G4ThreeVector T2Position = T2.GetPosition();

    new G4PVPlacement(0,               // rotation
                      T2Position,      // position
                      logT2,           // logical volume
                      T2.GetName(),    // name
                      logWorld,        // mother  volume
                      false,           // no boolean operations
                      5);

    Region->AddRootLogicalVolume(logT2);





//---------------- adding V2

    ScintConstructor& V2 = Get<ScintConstructor>("V2");

    V2.SetType("Tubs");
    V2.SetRmin(107.95*mm);
    V2.SetRmax(150.65*mm);
    V2.SetLength(6.35*mm);
    position.setZ(V2_z);
    V2.SetPosition(position);

    G4LogicalVolume* logV2 = V2.GetPiece();

    visual->SetColor(0.5,0.5,1.0,0); // blue-green
    logV2->SetVisAttributes(visual);

    G4ThreeVector V2Position = V2.GetPosition();

    new G4PVPlacement(0,               // rotation
                      V2Position,      // position
                      logV2,           // logical volume
                      V2.GetName(),    // name
                      logWorld,        // mother  volume
                      false,           // no boolean operations
                      6);

    Region->AddRootLogicalVolume(logV2);
// ------------adding V3
    ScintConstructor V3 = Get<ScintConstructor>("V3");

    V3.SetType("Tubs");
    V3.SetRmin(177.8 *mm);
    V3.SetRmax(241.3*mm);
    V3.SetLength(6.35*mm);
    position.setZ(V3_z);
    V3.SetPosition(position);

    G4LogicalVolume* logV3 = V3.GetPiece();

    visual->SetColor(0.5,0.5,1.0,0); // blue-green
    logV3->SetVisAttributes(visual);

    G4ThreeVector V3Position = V3.GetPosition();

    new G4PVPlacement(0,               // rotation
                      V3Position,      // position
                      logV3,           // logical volume
                      V3.GetName(),    // name
                      logWorld,        // mother  volume
                      false,           // no boolean operations
                      7);

    Region->AddRootLogicalVolume(logV3);

    //------------------------------ NaI

    G4double NaI_rmin = 0.0*cm;
    G4double NaI_rmax = 241.3*mm;
    G4double NaI_dz   = 481*mm;

    G4Tubs* NaI_tub
            = new G4Tubs("NaI",NaI_rmin,NaI_rmax,NaI_dz/2.,0.*deg,360.*deg);

    G4LogicalVolume* NaI_log = 
        new G4LogicalVolume(NaI_tub,
                            FindMaterial("G4_SODIUM_IODIDE"),
                            "NaI",0,0,0);

    visual->SetColor(0.5,0.5,0.5,1); // white
    NaI_log->SetVisAttributes(visual);

    SensitiveDetectorFactory factory("segment");
    NaI_log->SetSensitiveDetector(factory.MakeSD("crystalSD"));

    new G4PVPlacement(0,G4ThreeVector(0.0,0.0,NaI_z),
                      NaI_log,"NaI",logWorld,false,0);

    //------------------------------ NaI enclosure

    G4double bina_shell_rmin = NaI_rmax+0.1*cm;
    G4double bina_shell_rmax = NaI_rmax+0.4*cm;
    G4double bina_shell_dz   = 481*mm;///52.3*cm;

    G4Tubs* bina_shell_tub
      = new G4Tubs("shell",bina_shell_rmin,bina_shell_rmax,bina_shell_dz/2.,
                                                           0.*deg,360.*deg);

    G4LogicalVolume* bina_shell_log
      = new G4LogicalVolume(bina_shell_tub,
                            FindMaterial("G4_Al"),
                            "shell",0,0,0);

    bina_shell_log->SetVisAttributes(visual);

    new G4PVPlacement(0,G4ThreeVector(0.0,0.0,NaI_z),
                      bina_shell_log,"shell",logWorld,false,0);

    G4double bina_front_rmin = 0.0*cm;
    G4double bina_front_rmax = 22.86*cm;
    G4double bina_front_dz   = 0.5*mm;

    G4Tubs* bina_front_tub
      = new G4Tubs("front",bina_front_rmin,bina_front_rmax,bina_front_dz/2.,
                                                           0.*deg,360.*deg);

    G4LogicalVolume* bina_front_log
      = new G4LogicalVolume(bina_front_tub,
                            FindMaterial("G4_Al"),
                            "front",0,0,0);
    bina_front_log->SetSensitiveDetector(factory.MakeSD("crystalSD")); //Chloe's addition
    bina_front_log->SetVisAttributes(visual);

    new G4PVPlacement(0,G4ThreeVector(0.0,0.0,NaI_z-NaI_dz/2-0.25*mm),
                      bina_front_log,"front",logWorld,false,0);

    G4double bina_flange1_rmin = NaI_rmax+.4*cm;
    G4double bina_flange1_rmax = bina_flange1_rmin+2.38*cm;
    G4double bina_flange1_dz   = 5.55*cm+13.462*mm;

    G4Tubs* bina_flange1_tub
      = new G4Tubs("flange1",bina_flange1_rmin,bina_flange1_rmax,
                                       bina_flange1_dz/2.,0.*deg,360.*deg);

    G4LogicalVolume* bina_flange1_log
      = new G4LogicalVolume(bina_flange1_tub,
                            FindMaterial("G4_Al"),
                            "flange1",0,0,0);

    bina_flange1_log->SetVisAttributes(visual);

    new G4PVPlacement(0,
               G4ThreeVector(0.0,0.0,NaI_z-NaI_dz/2-0.5*mm+bina_flange1_dz/2.-3.1962*cm),
               bina_flange1_log,"flange1_u",logWorld,false,0);

    new G4PVPlacement(0,
		      G4ThreeVector(0.0,0.0,NaI_z+NaI_dz/2-bina_flange1_dz/2.+3.1962*cm),
               bina_flange1_log,"flange1_d",logWorld,false,0);

    G4double bina_flange2_rmin = 17.8*cm;
    G4double bina_flange2_rmax = NaI_rmax+0.4*cm;
    G4double bina_flange2_dz   = 3.1962*cm;
    
    G4Tubs* bina_flange2_tub
      = new G4Tubs("flange2",bina_flange2_rmin,bina_flange2_rmax,
		    bina_flange2_dz/2.,0.*deg,360.*deg);
    
    G4LogicalVolume* bina_flange2_log
      = new G4LogicalVolume(bina_flange2_tub,
                            FindMaterial("G4_Al"),
                            "flange2",0,0,0);
    
    bina_flange2_log->SetVisAttributes(visual);
    
    new G4PVPlacement(0,G4ThreeVector(0.0,0.0,NaI_z-NaI_dz/2-0.5*mm-15.981*mm),
		      bina_flange2_log,"flange2_u",logWorld,false,0);

    new G4PVPlacement(0,G4ThreeVector(0.0,0.0,NaI_z+NaI_dz/2+15.981*mm),
		      bina_flange2_log,"flange2_d",logWorld,false,0);

    

//     //   ------------------------------
    
   //  //Dorothea'as addition of Zirconium foil
    
    G4double Zr_rmin = 0*cm;
    G4double Zr_rmax = 7*cm;
    G4double Zr_dz   = 0.25*mm;  // default: 0.25mm for Zr; 0.71mm for mylar
    
    G4Tubs* Zr_tub
      = new G4Tubs("Zr",Zr_rmin,Zr_rmax,
		   Zr_dz/2.,0.*deg,360.*deg);
    
    G4LogicalVolume* Zr_log
      = new G4LogicalVolume(Zr_tub,
			    FindMaterial("G4_Zr"),   // G4_Zr for Zirconium foil, G4_MYLAR for mylar foil
			    "Zr",0,0,0);
    
    Zr_log->SetVisAttributes(visual);
    
    new G4PVPlacement(0,G4ThreeVector(0.0,0.0,T1_z+T1_lngth/2+2*mm),
		      Zr_log,"Zr",logWorld,false,0);

    G4cout << "T1 Z " << T1_z << " Zr Z " << T1_z+T1_lngth/2+2*mm  <<  G4endl;
    
    // make sensitive detector for Zr
    Zr_log->SetSensitiveDetector(factory.MakeSD("zirconiumSD"));
    
    SetSensitiveDetector("zirconiumSD","segment");
    Zr_log->SetSensitiveDetector(GetSensitiveDetector());
    
    //================================================================
    // Make five zirconium foils to track the momentum in each slice
    //================================================================
        
   //  G4double z_coord = T1_z+T1_lngth/2+2*mm;
//     G4double Zr_rmin = 0*cm;
//     G4double Zr_rmax = 7*cm;
//     G4double Zr_dz   = 0.05*mm;  // for five slices of 0.05 mm each, making up 0.25 mm in total
//     G4Tubs* Zr_tub
//       = new G4Tubs("Zr",Zr_rmin,Zr_rmax,
// 		   Zr_dz/2.,0.*deg,360.*deg);
    
//     // First slice 
//     G4double Zr1_coord = z_coord - 3*Zr_dz;
    
//     G4LogicalVolume* Zr1_log
//       = new G4LogicalVolume(Zr_tub,
// 			    FindMaterial("G4_Zr"),   // G4_Zr for Zirconium foil, G4_MYLAR for mylar foil
// 			    "Zr1",0,0,0);
//     Zr1_log->SetVisAttributes(visual);
//     new G4PVPlacement(0,G4ThreeVector(0.0,0.0,Zr1_coord),
// 		      Zr1_log,"Zr1",logWorld,false,0);
//     // make sensitive detector for Zr
//     Zr1_log->SetSensitiveDetector(factory.MakeSD("zirconiumSD"));
//     SetSensitiveDetector("zirconiumSD","segment");
//     Zr1_log->SetSensitiveDetector(GetSensitiveDetector());

//     // Second Slice
//     G4double Zr2_coord = z_coord - 2*Zr_dz;
    
//     G4LogicalVolume* Zr2_log
//       = new G4LogicalVolume(Zr_tub,
// 			    FindMaterial("G4_Zr"),   // G4_Zr for Zirconium foil, G4_MYLAR for mylar foil
// 			    "Zr2",0,0,0);
//     Zr2_log->SetVisAttributes(visual);
//     new G4PVPlacement(0,G4ThreeVector(0.0,0.0,Zr2_coord),
// 		      Zr2_log,"Zr2",logWorld,false,0);
//     Zr2_log->SetSensitiveDetector(factory.MakeSD("zirconiumSD"));
//     SetSensitiveDetector("zirconiumSD","segment");
//     Zr2_log->SetSensitiveDetector(GetSensitiveDetector());

//     // Third Slice
//     G4double Zr3_coord = z_coord - 1*Zr_dz;
    
//     G4LogicalVolume* Zr3_log
//       = new G4LogicalVolume(Zr_tub,
// 			    FindMaterial("G4_Zr"),   // G4_Zr for Zirconium foil, G4_MYLAR for mylar foil
// 			    "Zr3",0,0,0);
//     Zr3_log->SetVisAttributes(visual);
//     new G4PVPlacement(0,G4ThreeVector(0.0,0.0,Zr3_coord),
// 		      Zr3_log,"Zr3",logWorld,false,0);
//     Zr3_log->SetSensitiveDetector(factory.MakeSD("zirconiumSD"));
//     SetSensitiveDetector("zirconiumSD","segment");
//     Zr3_log->SetSensitiveDetector(GetSensitiveDetector());

//     // Fourth Slice
//     G4double Zr4_coord = z_coord;
    
//     G4LogicalVolume* Zr4_log
//       = new G4LogicalVolume(Zr_tub,
// 			    FindMaterial("G4_Zr"),   // G4_Zr for Zirconium foil, G4_MYLAR for mylar foil
// 			    "Zr4",0,0,0);
//     Zr4_log->SetVisAttributes(visual);
//     new G4PVPlacement(0,G4ThreeVector(0.0,0.0,Zr4_coord),
// 		      Zr4_log,"Zr4",logWorld,false,0);
//     Zr4_log->SetSensitiveDetector(factory.MakeSD("zirconiumSD"));
//     SetSensitiveDetector("zirconiumSD","segment");
//     Zr4_log->SetSensitiveDetector(GetSensitiveDetector());

//     // Fifth Slice
//     G4double Zr5_coord = z_coord + Zr_dz;
    
//     G4LogicalVolume* Zr5_log
//       = new G4LogicalVolume(Zr_tub,
// 			    FindMaterial("G4_Zr"),   // G4_Zr for Zirconium foil, G4_MYLAR for mylar foil
// 			    "Zr5",0,0,0);
//     Zr5_log->SetVisAttributes(visual);
//     new G4PVPlacement(0,G4ThreeVector(0.0,0.0,Zr5_coord),
// 		      Zr5_log,"Zr5",logWorld,false,0);
//     Zr5_log->SetSensitiveDetector(factory.MakeSD("zirconiumSD"));
//     SetSensitiveDetector("zirconiumSD","segment");
//     Zr5_log->SetSensitiveDetector(GetSensitiveDetector());


//     G4cout << "Zr1 " << Zr1_coord << " Zr2 " << Zr2_coord  <<  " Zr3 " << Zr3_coord << " Zr4 " << Zr4_coord << " Zr5 " << Zr5_coord << G4endl;

     //    ------------------------------
    return logWorld;
}
