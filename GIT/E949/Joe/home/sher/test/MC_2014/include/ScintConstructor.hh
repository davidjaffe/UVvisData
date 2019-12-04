#ifndef ScintConstructor_h
#define ScintConstructor_h 1

#include "GeoConstructor.hh"
#include "G4ThreeVector.hh"

class G4LogicalVolume;

class ScintConstructor : public GeoConstructor {

  public:

    ScintConstructor(G4String n, DetectorConstruction* c)
        : GeoConstructor(n,c) {Init();};

    ScintConstructor(G4String n, GeoConstructor* p)
        : GeoConstructor(n,p) {Init();};

    virtual ~ScintConstructor();

    enum Shape {Box, Tubs};

    /// Construct and return a G4 volume for the object.  This is a pure
    /// virtual function, which means it must be implemented by the inheriting
    /// classes.  This returns an unplaced logical volume which faces along
    /// the Z axis.
    virtual G4LogicalVolume *GetPiece(void);

    /// Set the type of the scintillator.
    void SetType(G4String t) {
      if(t == "Box")fType = Box;
      if(t == "Tubs")fType = Tubs;
    }

    /// Get the type of the scintillator.
    Shape GetType(void) {return fType;}

    /// Set the width of the scintillator.
    void SetWidth(G4double w) {fWidth = w;}

    /// Get the width of the scintillator.
    G4double GetWidth(void) {return fWidth;}

    /// Set the height of the scintillator.
    void SetHeight(G4double w) {fHeight = w;}

    /// Get the height of the scintillator.
    G4double GetHeight(void) {return fHeight;}

    /// Set the length of the scintillator along the beam axis.
    void SetLength(G4double w) {fLength = w;}

    /// Get the length of the scintillator along the beam axis.
    G4double GetLength(void) {return fLength;}

    /// Set the inner radius of a circular scintillator.
    void SetRmin(G4double r) {fRmin = r;}

    /// Get the inner radius of a scintillator.
    G4double GetRmin(void) {return fRmin;}

    /// Set the outer radius of a circular scintillator.
    void SetRmax(G4double r) {fRmax = r;}

    /// Get the outer radius of a scintillator.
    G4double GetRmax(void) {return fRmax;}

    /// Set the position of the scintillator.
    void SetPosition(G4ThreeVector v) {fPosition = v;}

    /// Get the position of the scintillator.
    G4ThreeVector GetPosition() {return fPosition;}

    /// Set the material name of the scintillator.
    void SetMaterialName(G4String n) {fMaterialName = n;}

    /// Get the position of the scintillator.
    G4String GetMaterialName() {return fMaterialName;}

  private:

    void Init(void);

    /// The type of scintillator - Box; Tubs
    Shape fType;

    /// The width of the scintillator 
    /// (horizontal size perpendicular to the beam axis)
    G4double fWidth;

    /// The height of the scintillator 
    /// (vertical size perpendicular to the beam axis).
    G4double fHeight;

    /// The length of the scintillator 
    /// (parallel to the beam axis).
    G4double fLength;

    /// The inner radius of a circular scintillator
    G4double fRmin;

    /// The outer radius of a circular scintillator
    G4double fRmax;

    /// The position of the scintillator
    G4ThreeVector fPosition;

    /// The material of the scintillator
    G4String fMaterialName;
};
#endif
