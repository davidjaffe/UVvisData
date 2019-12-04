#ifndef SegmentSD_h
#define SegmentSD_h 1

#include "G4VSensitiveDetector.hh"

#include "HitSegment.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

/// A sensitive detector to create HitSegment based hits.
class SegmentSD : public G4VSensitiveDetector
{
    
  public:

    SegmentSD(G4String name);
    virtual ~SegmentSD();
    
    void Initialize(G4HCofThisEvent*);
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    void EndOfEvent(G4HCofThisEvent*);

    virtual HitSegment* GetNewHit() { return new HitSegment(); }

  protected:

    /// The collection of hits that is being filled in the current event.  It
    /// is constructed in Initialize, filled in ProcessHits, and added the the
    /// event in EndOfEvent.

    HitSegmentCollection* fHits;

    /// The hit collection id of fHits
    int fHCID;

    /// The last hit that was found.
    int fLastHit;
};

#endif
