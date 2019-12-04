#ifndef WorldConstructor_h
#define WorldConstructor_h 1

#include "GeoConstructor.hh"

class G4LogicalVolume;

/// Construct the PIENU world geometry.  The origin is located at the target
/// position.  The beam runs in the Y-Z plane and +Y is pointing up with X 
/// completing the right-handed coordinate system.
class WorldConstructor : public GeoConstructor
{
  public:

    WorldConstructor(G4String n, DetectorConstruction* c)
        : GeoConstructor(n,c) {Init();};
    virtual ~WorldConstructor();

    /// Construct and return a G4 volume for the object.  This is a pure
    /// virtual function, which means it must be implemented by the inheriting
    /// classes.  This returns an unplaced logical volume which faces along
    /// the Z axis.
    virtual G4LogicalVolume *GetPiece(void);

    /// Set the width of the world
    void SetWidth(G4double v) {fWidth = v;}
    
    /// Set the height of the world
    void SetHeight(G4double v) {fHeight = v;}
    
    /// Set the length of the world
    void SetLength(G4double v) {fLength = v;}

    /// Get the width of the world
    G4double GetWidth() {return fWidth;}

    /// Get the height of the world
    G4double GetHeight() {return fHeight;}

    /// Get the length of the world
    G4double GetLength() {return fLength;}

    /// Choose to have nothing but one stand alone module
    virtual void SetUseStandAloneModule(bool UseSAM)
        { fUseStandAloneModule=UseSAM; }
    virtual void SetStandAloneModuleName(G4String name)
        { fStandAloneModuleName=name; }

  private:

    void Init(void);

    /// The distance along the Z axis to simulate.
    G4double fLength;

    /// The side to side dimension to be simulated
    G4double fWidth;

    /// The height of the world coordinate system. 
    G4double fHeight;

    G4bool fUseStandAloneModule;
    G4String fStandAloneModuleName;

};

#endif
