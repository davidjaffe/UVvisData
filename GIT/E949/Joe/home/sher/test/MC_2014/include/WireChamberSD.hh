#ifndef WireChamberSD_h
#define WireChamberSD_h 1

#include "SegmentSD.hh"
#include "WireChamberHitSegment.hh"

class WireChamberSD : public SegmentSD
{

public:

  WireChamberSD(G4String name) : SegmentSD(name) {
    collectionName.insert("wireChamberCollection");
  }

  virtual ~WireChamberSD() {}

  virtual WireChamberHitSegment* GetNewHit() {
    return new WireChamberHitSegment();
  }

private:

};

#endif
