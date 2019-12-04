#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "DetectorMessenger.hh"
#include "G4VUserDetectorConstruction.hh"

class GeoConstructor;
class G4VPhysicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DetectorConstruction();
    virtual ~DetectorConstruction();

    /// Construct the detector and define the world volume.
    virtual G4VPhysicalVolume* Construct();

    /// This really constructs the detector, but doesn't define
    /// materials before it's constructed.  This is called by Construct()
    G4VPhysicalVolume* ConstructDetector();

    /// Return the detector construction messenger
    DetectorMessenger* GetMessenger(void) {return fDetectorMessenger;}

    void SetVerbose(G4int);

    /// Update the geometry
    void UpdateGeometry();

  protected:

    /// Define the materials used in the detector.
    void DefineMaterials(void);

    /// A messenger for the DetectorConstruction object.
    DetectorMessenger* fDetectorMessenger;

    /// A constructor to create the world.
    GeoConstructor* fWorldConstructor;

  private:

    G4int verbose;

};

#endif
