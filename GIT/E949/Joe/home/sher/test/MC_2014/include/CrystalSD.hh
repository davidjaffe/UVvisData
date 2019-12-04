#ifndef CrystalSD_h
#define CrystalSD_h 1

#include "SegmentSD.hh"
#include "CrystalHitSegment.hh"

class G4Step;
class G4TouchableHistory;

class CrystalSD : public SegmentSD
{

public:

  CrystalSD(G4String name) : SegmentSD(name) {
    collectionName.insert("crystalCollection");
  }

  virtual ~CrystalSD() {}

  virtual CrystalHitSegment* GetNewHit() {
    return new CrystalHitSegment(); }

private:

};

#endif
