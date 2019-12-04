#ifndef SensitiveDetectorFactory_h
#define SensitiveDetectorFactory_h 1

#include "G4String.hh"

class G4VSensitiveDetector;

/// Create an object to build a sensitive detector (SD) based on a name and
/// type.  The name of the SD will be used by the G4 SD manager to see if the
/// detector already exists.  If it doesn't then a new SD of "type" will be
/// built.

class SensitiveDetectorFactory {

public:

    /// Build a factory to build sensitive detectors specified by "type".
    SensitiveDetectorFactory(G4String type);
    virtual ~SensitiveDetectorFactory();

    /// Get pointer to a sensitive detector built by this factory, but return
    /// null if the detector doesn't exist.
    G4VSensitiveDetector* GetSD(G4String name);

    /// Get pointer to a sensitive detector built by this factory, and create
    /// a new sensitive detector if required.
    G4VSensitiveDetector* MakeSD(G4String name);

private:

    /// The type of sensitive detector that this will build.
    G4String fType;

};
#endif
