#include "GeoConstructor.hh"
#include "DetectorConstruction.hh"

GeoConstructor::GeoConstructor(G4String n, DetectorConstruction* c) 
    : fLocalName(n), fName(n), fConstruction(c), fParent(NULL),
      fMessenger(NULL), fVisible(false), fSensitiveDetector(NULL) { }

GeoConstructor::GeoConstructor(G4String n, GeoConstructor* p) 
    : fLocalName(n), fName(n), fConstruction(NULL), fParent(p), 
      fMessenger(NULL), fVisible(false), fSensitiveDetector(NULL)
{
    fName = fParent->GetName() + "/" + fLocalName;
    fConstruction = fParent->GetConstruction();
    fMessenger = fParent->GetMessenger();
}

GeoConstructor::~GeoConstructor() {}

G4LogicalVolume *GeoConstructor::GetPiece(void)
{
  G4Exception("GeoConstructor::GetPiece(): Not implemented","",FatalErrorInArgument,"");
    return NULL;
}

G4Material* GeoConstructor::FindMaterial(G4String name)
{
    G4Material* material = G4Material::GetMaterial(name,true);
    return material;
}

G4String GeoConstructor::GetName(void) {
    return fName;
}

/// Set the local name of the object to be built.  The local name of the
/// logical volume being set, and also sets the full name of the volume.
void GeoConstructor::SetLocalName(const G4String& name)
{
    fLocalName = name;
    fName = name;
    if (fParent) fName = fParent->GetName() + "/" + fName;
}

/// Set the visibility of the constructed object.
void GeoConstructor::SetVisible(bool v)
{
    if (fVisible != v) {
        G4cout << "Set Visibility for " << GetName()
                     << " from " << fVisible << " to " << v << G4endl;
    }
    fVisible = v;
}

/// Set the visibility of the constructed object.
void GeoConstructor::SetVisibleDaughters(bool v)
{
    if (fVisible != v) {
        G4cout << "Set daughter visibility for " << GetName()
                     << " to " << v << G4endl;
    }

    SetVisible(false);

    for (std::map<G4String,GeoConstructor*>::iterator p 
             = fSubConstructors.begin();
         p!=fSubConstructors.end();
         ++p) {
        (*p).second->SetVisible(v);
    }
}

G4RotationMatrix GeoConstructor::stringToRotationMatrix(G4String rotation)
{
  // We apply successive rotations OF THE OBJECT around the FIXED
  // axes of the parent's local coordinates; rotations are applied
  // left-to-right (rotation="r1,r2,r3" => r1 then r2 then r3).

  G4RotationMatrix rot;

  unsigned int place = 0;

  while (place < rotation.size()) {

        G4double angle;
        char* p(0);
        G4String current=rotation.substr(place+1);
        angle = strtod(current.c_str(),&p) * deg;

        if (!p || (*p != ',' && *p != '\0')) {
          G4cout << "Invalid rotation specification: "
                     << rotation.c_str() << G4endl;

           return rot;
        }

        G4RotationMatrix thisRotation;

        switch(rotation.substr(place,1).c_str()[0]) {
              case 'X': case 'x':
                thisRotation = G4RotationMatrix(CLHEP::HepRotationX(angle));
                break;
              case 'Y': case 'y':
                thisRotation = G4RotationMatrix(CLHEP::HepRotationY(angle));
                break;
              case 'Z': case 'z':
                thisRotation = G4RotationMatrix(CLHEP::HepRotationZ(angle));
                break;
              default:
                G4cout << " Invalid rotation specification: "
                           << rotation << G4endl;
                return rot;
        }

       rot = thisRotation * rot;
       place = rotation.find(',',place);
       if (place > rotation.size()) break;
       ++place;
  }

  return rot;

}
