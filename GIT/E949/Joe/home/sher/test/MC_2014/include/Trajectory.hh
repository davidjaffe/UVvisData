#ifndef Trajectory_h_seen
#define Trajectory_h_seen 1

#include "G4ios.hh"

#include "globals.hh"
#include "G4VTrajectory.hh"
#include "G4Allocator.hh"

#include "G4ParticleDefinition.hh"
#include "G4TrajectoryPoint.hh"
#include "G4Track.hh"
#include "G4Step.hh"

typedef std::vector<G4VTrajectoryPoint*> TrajectoryPointContainer;
class Trajectory : public G4VTrajectory
{

//--------
   public: // without description
//--------

// Constructor/Destructor

     Trajectory();
     Trajectory(const G4Track* aTrack);
     Trajectory(Trajectory &);
     virtual ~Trajectory();

// Operators

     inline void* operator new(size_t);
     inline void  operator delete(void*);
     inline int operator == (const Trajectory& right) const
     { return (this==&right); }

// Get/Set functions

     inline G4int GetTrackID() const { return fTrackID; }
     inline G4int GetParentID() const { return fParentID; }
     inline G4String GetParticleName() const { return ParticleName; }
     inline G4double GetCharge() const { return PDGCharge; }
     inline G4int GetPDGEncoding() const { return PDGEncoding; }
     inline G4ThreeVector GetInitialMomentum() const { return initialMomentum; }

// Other member functions

     virtual void ShowTrajectory(std::ostream& os=G4cout) const;
     virtual void DrawTrajectory(G4int i_mode=0) const;
     virtual void AppendStep(const G4Step* aStep);
     virtual void MergeTrajectory(G4VTrajectory* secondTrajectory);

     G4ParticleDefinition* GetParticleDefinition();

     virtual int GetPointEntries() const
     { return fpPointsContainer->size(); }
     virtual G4VTrajectoryPoint* GetPoint(G4int i) const
     { return (*fpPointsContainer)[i]; }

    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
    virtual std::vector<G4AttValue>* CreateAttValues() const;

//---------
   private:
//---------

// Member data

     TrajectoryPointContainer* fpPointsContainer;

     G4int fTrackID;
     G4int fParentID;
     G4double PDGCharge;
     G4int    PDGEncoding;
     G4String ParticleName;
     G4ThreeVector initialMomentum;

};

#if defined G4TRACKING_ALLOC_EXPORT
  extern G4DLLEXPORT G4Allocator<Trajectory> aTrajectoryAllocator;
#else
  extern G4DLLIMPORT G4Allocator<Trajectory> aTrajectoryAllocator;
#endif

inline void* Trajectory::operator new(size_t) {
    void* aTrajectory = (void*) aTrajectoryAllocator.MallocSingle();
    return aTrajectory;
}

inline void Trajectory::operator delete(void* aTrajectory) {
    aTrajectoryAllocator.FreeSingle((Trajectory*)aTrajectory);
}

#endif
