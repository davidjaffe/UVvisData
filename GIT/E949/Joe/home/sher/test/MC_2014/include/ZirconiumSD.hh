#ifndef ZirconiumSD_h
#define ZirconiumSD_h 1

#include "SegmentSD.hh"
#include "ZirconiumHitSegment.hh"

class G4Step;
class G4TouchableHistory;

class ZirconiumSD : public SegmentSD
{

public:

  ZirconiumSD(G4String name) : SegmentSD(name) {
    collectionName.insert("zirconiumCollection");
  }

  virtual ~ZirconiumSD() {}

  virtual ZirconiumHitSegment* GetNewHit() {
    return new ZirconiumHitSegment(); }

private:

};

#endif
