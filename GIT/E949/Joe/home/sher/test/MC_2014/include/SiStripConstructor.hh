#ifndef SiStripConstructor_h
#define SiStripConstructor_h 1

#include "GeoConstructor.hh"
#include "G4ThreeVector.hh"

class G4LogicalVolume;

class SiStripConstructor : public GeoConstructor {

  public:

    SiStripConstructor(G4String n, DetectorConstruction* c)
        : GeoConstructor(n,c) {Init();};

    SiStripConstructor(G4String n, GeoConstructor* p)
        : GeoConstructor(n,p) {Init();};

    virtual ~SiStripConstructor();

    /// Construct and return a G4 volume for the object.  This is a pure
    /// virtual function, which means it must be implemented by the inheriting
    /// classes.  This returns an unplaced logical volume which faces along
    /// the Z axis.
    virtual G4LogicalVolume *GetPiece(void);

  private:

    void Init(void);

};
#endif
