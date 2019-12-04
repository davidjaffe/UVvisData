#ifndef VolumeID_h
#define VolumeID_h 1

#include <vector>

#include "G4VPhysicalVolume.hh"
#include "G4TouchableHandle.hh"

/// A class to provide a unique identifier for all physical volumes.
/// One of the features of G4 is that physical volumes are identified by their
/// position in the volume tree.  This has many advantages, but it means that
/// it can be quite tricky to figure out if two physical volume pointers refer
/// to the same physical volume, or different replicas of the same physical
/// volume.  (Consider the case where you have an unreplicated sub-volume of
/// a replicated parent).  This class provides an unique comparable object that
/// allows an equality test between volumes which takes into account the
/// full position in the hierarchy.  It is similar in function to
/// G4TouchableHandle, but provides better comparison operators.

class VolumeID {

  public:

    /// Construct a new volume ID.
    VolumeID(const G4TouchableHandle& handle);
    /// Construct an empty volume ID.
    VolumeID();
    ~VolumeID();

    /// Explicitly add a new volume to the volume ID.
    void AddVolume(G4VPhysicalVolume* fPhysVol,G4int fReplica);

    G4bool operator !() {
        return (fVolumes.size()==0);
    };

    VolumeID& operator = (const G4TouchableHandle& handle);
    VolumeID& operator = (const VolumeID& id);

    friend G4bool operator ==(const VolumeID& x, 
                            const G4TouchableHandle& y);
    friend G4bool operator ==(const G4TouchableHandle& x,
                            const VolumeID& y);
    friend G4bool operator ==(const VolumeID& x,
                            const VolumeID& y);

    friend G4bool operator !=(const VolumeID& x, 
                            const G4TouchableHandle& y);
    friend G4bool operator !=(const G4TouchableHandle& x,
                            const VolumeID& y);
    friend G4bool operator !=(const VolumeID& x,
                            const VolumeID& y);

    friend G4bool operator <(const VolumeID& x,
                           const VolumeID& y);

    G4VPhysicalVolume* GetVolume(G4int depth=0) const {
      return fVolumes[depth].fPhysVol;
    }
    G4int GetReplicaNumber(G4int depth=0) const {
      return fVolumes[depth].fReplica;
    }

  private:

    class Volume {
      public:
        G4VPhysicalVolume* fPhysVol;
        G4int fReplica;
    };
    typedef std::vector<Volume> Volumes;

    /// This is a vector of information required to uniquely ID the volume.
    /// It is arranges so that the most local information is at the front of
    /// the vector and the world volume will be the last element in the
    /// vector. 

    Volumes fVolumes;

};
#endif
