#ifndef ZirconiumHitSegment_h
#define ZirconiumHitSegment_h 1

#include "HitSegment.hh"

class ZirconiumHitSegment : public HitSegment
{

public:

  ZirconiumHitSegment() {}
  virtual ~ZirconiumHitSegment() {}

  virtual G4bool CheckMergeCondition(HitSegment* ) { return true; }

};
#endif
