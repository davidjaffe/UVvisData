#ifndef ScintSD_h
#define ScintSD_h 1

#include "SegmentSD.hh"
#include "ScintHitSegment.hh"

class ScintSD : public SegmentSD
{

public:

  ScintSD(G4String name) : SegmentSD(name) {
    collectionName.insert("scintCollection");
  }
  virtual ~ScintSD() {}

  virtual ScintHitSegment* GetNewHit() {return new ScintHitSegment();}

private:

};

#endif
