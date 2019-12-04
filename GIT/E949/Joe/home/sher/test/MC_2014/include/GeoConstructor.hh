#ifndef GeoConstructor_h
#define GeoConstructor_h 1

#include "globals.hh"
#include "G4Material.hh"
#include "G4RotationMatrix.hh"
#include "SensitiveDetectorFactory.hh"

#include "GeoConstructorMessenger.hh"

class G4LogicalVolume;
class G4VSensitiveDetector;

class DetectorConstruction;

/// Base for classes that construct elements of the detector.  The method 
/// get GetPiece() should construct a new unplaced element that is facing
/// up along the Z axis.  The calling code should then rotate the object into
/// position and and place the object into a large mother volume.

class GeoConstructor
{
  public:

    GeoConstructor(G4String n, DetectorConstruction* c);
    GeoConstructor(G4String n, GeoConstructor* parent);
    virtual ~GeoConstructor();

  private:

    /// The short local name of the constructor.
    G4String fLocalName;

    /// The name of the constructor.
    G4String fName;

    /// The G4VUserDetectorConstruction class that this is cooperating with.
    DetectorConstruction* fConstruction;

    /// The parent of this constructor
    GeoConstructor* fParent;

    /// The user interface messenger that will provide values for this class. 
    GeoConstructorMessenger* fMessenger;

    /// If this is true, then draw constructed object.
    bool fVisible;

    /// The sub constructors that might be used in this class.
    std::map<G4String,GeoConstructor*> fSubConstructors;

    /// The sensitive detector for this component.
    G4VSensitiveDetector* fSensitiveDetector;

  protected:

    G4Material* FindMaterial(G4String m);

  public:

    /// Construct and return a G4 volume for the object.  This is a pure
    /// virtual function, which means it must be implemented by the inheriting
    /// classes.  This returns an unplaced logical volume which faces along
    /// the Z axis.
    virtual G4LogicalVolume *GetPiece(void) = 0;

    /// Return the base name of the object that this builds.
    G4String GetName(void);

    /// Return the base name of the object that this builds.
    G4String GetLocalName(void) {return fLocalName;};

    /// Set the base name of the logical volume that this builds.
    void SetLocalName(const G4String& name);

    /// Set the visibility of the constructed object.
    void SetVisible(bool v);
    
    /// Set the visibility of the object daughters
    void SetVisibleDaughters(bool v);
    
    /// Get the visibility of the constructed object.
    bool GetVisible(void) const {return fVisible;}

    /// Return the detector construction that is building this piece.
    DetectorConstruction* GetConstruction(void) {
        return fConstruction;
    }
    
    /// Return the messenger for this constructor
    GeoConstructorMessenger* GetMessenger(void) {return fMessenger;}

    /// Set the messenger for this constructor.
    void SetMessenger(GeoConstructorMessenger* m) {
        fMessenger = m;
    }

    /// Add a new sub constructor to the current geometry constructor.  This
    /// should be done in the derived class constructor.  All geometry 
    /// constructors that will be used by the derived class should be added
    /// using this method.
    void AddConstructor(GeoConstructor* c) {
        if (fSubConstructors.find(c->GetLocalName()) 
            != fSubConstructors.end()) {
            G4cout << "Multiply defined constructor name " 
                   << c->GetName()
                   << G4endl;
            G4Exception("GeoConstructor::AddConstructor():", " Multiple defines",FatalErrorInArgument,"");
        }
        fSubConstructors[c->GetLocalName()] = c;
    }

    /// Get a sub constructor by name and do the dynamic cast.  This will
    /// abort with an error message if you request an undefined name.
    template <class T> T& Get(G4String n) {
        std::map<G4String,GeoConstructor*>::iterator p 
            = fSubConstructors.find(n);
        if (p==fSubConstructors.end()) {
            G4cout << "Error in " << GetName() << G4endl;
            G4cout << "   Undefined constructor requested "
                   << n << G4endl;
            G4cout << "   The choices are: " << G4endl;
            for (p = fSubConstructors.begin();
                 p != fSubConstructors.end();
                 ++p) {
                G4cout << "        \"" << (*p).first << "\""
                       << " for constructor: " << (*p).second->GetName()
                       << G4endl;
            }
            G4Exception("GeoConstructor::Get<>:", " Undefined constructor",FatalErrorInArgument,"");
        }
        T* c = dynamic_cast<T*>((*p).second);
        if (!c) {
            G4cout << "Error in " << GetName() << G4endl;
            G4cout << "  Cannot type cast " << n 
                      << " to requested class" << G4endl;
            G4Exception("GeoConstructor::Get<>:", " Bad typecast",FatalErrorInArgument,"");
        }
        return *c;
    }

    /// Find a sub constructor by name and do the dynamic cast.  This returns
    /// a pointer that will be NULL if you request an undefined name.
    template <class T> T* Find(G4String n) {
        std::map<G4String,GeoConstructor*>::iterator p 
            = fSubConstructors.find(n);
        if (p==fSubConstructors.end()) return NULL;
        T* c = dynamic_cast<T*>((*p).second);
        return c;
    }

    /// Set the sensitive detector for this component.
    virtual void SetSensitiveDetector(G4VSensitiveDetector* s) {
        fSensitiveDetector = s;
    }

    /// Set the sensitive detector for this component by name.
    void SetSensitiveDetector(G4String name, G4String type) {
        SensitiveDetectorFactory factory(type);
        SetSensitiveDetector(factory.MakeSD(name));
    }

    /// Get the sensitive detector for this component.
    G4VSensitiveDetector* GetSensitiveDetector(void) {
        return fSensitiveDetector;
    }

    /// stringToRotationMatrix() converts a string "X90,Y45" into a
    /// G4RotationMatrix.
    /// This is an active rotation, in that the object is first rotated
    /// around the parent's X axis by 90 degrees, then the object is
    /// further rotated around the parent's Y axis by 45 degrees.
    /// The return value points to a G4RotationMatrix on the heap, so
    /// it is persistent. Angles are in degrees, can have decimals,
    /// and can be negative. Axes are X, Y, Z.
    static G4RotationMatrix stringToRotationMatrix(G4String rotation);
};

#endif
