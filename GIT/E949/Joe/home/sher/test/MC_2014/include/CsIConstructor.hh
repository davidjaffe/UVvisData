#ifndef CsIConstructor_h
#define CsIConstructor_h 1

#include "GeoConstructor.hh"
#include "G4ThreeVector.hh"

class G4LogicalVolume;

class CsIConstructor : public GeoConstructor {

  public:

    CsIConstructor(G4String n, DetectorConstruction* c)
       : GeoConstructor(n,c) {Init();};

    CsIConstructor(G4String n, GeoConstructor* p)
       : GeoConstructor(n,p) {Init();};

    virtual ~CsIConstructor();

    /// Construct and return a G4 volume for the object.  This is a pure
    /// virtual function, which means it must be implemented by the inheriting
    /// classes.  This returns an unplaced logical volume which faces along
    /// the Z axis.
    virtual G4LogicalVolume *GetPiece(void);

  private:

    void Init(void);

    void ConstructCsI_Upstream(void);
    void ConstructCsI_Downstream(void);

    G4LogicalVolume* ConstructCsI_Frame(void);

    G4LogicalVolume* csi_space1_log;
    G4LogicalVolume* csi_space2_log;
    G4LogicalVolume* csi_space3_log;
    G4LogicalVolume* csi_space4_log;

    G4double inch;

};

#endif
