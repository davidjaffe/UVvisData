#ifndef WireChamberHitSegment_h
#define WireChamberHitSegment_h 1

#include "HitSegment.hh"


class WireChamberHitSegment : public HitSegment
{

public:

  WireChamberHitSegment() {}
  virtual ~WireChamberHitSegment() {}

  virtual G4bool CheckMergeCondition(HitSegment* ) { return false; }

};

#endif
