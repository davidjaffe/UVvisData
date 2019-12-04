#ifndef CrystalHitSegment_h
#define CrystalHitSegment_h 1

#include "HitSegment.hh"

class CrystalHitSegment : public HitSegment
{

public:

  CrystalHitSegment() {}
  virtual ~CrystalHitSegment() {}

  virtual G4bool CheckMergeCondition(HitSegment* ) { return true; }

};
#endif
