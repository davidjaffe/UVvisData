#ifndef TrajectoryPoint_h
#define TrajectoryPoint_h 1

#include "globals.hh"

#include "G4TrajectoryPoint.hh"
#include "G4ThreeVector.hh"
#include "G4Allocator.hh"

#include "G4StepStatus.hh"

class G4Track;
class G4Step;
class G4VProcess;

class TrajectoryPoint : public G4TrajectoryPoint {

//--------
  public: // without description
//--------

// Constructor/Destructor

    TrajectoryPoint();
    TrajectoryPoint(const G4Track* aTrack);
    TrajectoryPoint(const G4Step* aStep);
    TrajectoryPoint(const TrajectoryPoint &right);
    virtual ~TrajectoryPoint();

// Operators

    inline void *operator new(size_t);
    inline void operator delete(void *aTrajectoryPoint);
    inline int operator==(const TrajectoryPoint& right) const
    { return (this==&right); };

// Get/Set functions

    inline G4double GetTime() const { return fTime; };
    inline const G4ThreeVector GetMomentum() const { return fMomentum; };
    inline G4StepStatus GetStepStatus() const { return fStepStatus; };
    inline G4String GetVolumeName() const { return fVolumeName; };

// Get method for HEPRep style attributes

   virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
   virtual std::vector<G4AttValue>* CreateAttValues() const;

//---------
  private:
//---------

// Member data

    G4double fTime;
    G4ThreeVector fMomentum;
    G4StepStatus fStepStatus;
    G4String fVolumeName;

};

#if defined G4TRACKING_ALLOC_EXPORT
  extern G4DLLEXPORT G4Allocator<TrajectoryPoint> aTrajPointAllocator;
#else
  extern G4DLLIMPORT G4Allocator<TrajectoryPoint> aTrajPointAllocator;
#endif

inline void* TrajectoryPoint::operator new(size_t)
{
    void *aTrajectoryPoint = (void *) aTrajPointAllocator.MallocSingle();
    return aTrajectoryPoint;
}

inline void TrajectoryPoint::operator delete(void *aTrajectoryPoint)
{
    aTrajPointAllocator.FreeSingle(
        (TrajectoryPoint *) aTrajectoryPoint);
}

#endif
