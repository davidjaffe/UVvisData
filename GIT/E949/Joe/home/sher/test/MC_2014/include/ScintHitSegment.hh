#ifndef ScintHitSegment_h
#define ScintHitSegment_h 1

#include "HitSegment.hh"

class ScintHitSegment : public HitSegment
{

public:

  ScintHitSegment() {}
  virtual ~ScintHitSegment() {}

  virtual G4bool CheckMergeCondition(HitSegment* rhs) {
    if ( fParentID == rhs->GetTrackID() ) {
         //Only merge delta rays into their parent
         //particle if they overlap in time
         if (fStartT < rhs->GetStopT() ||
            (fStopT  > rhs->GetStartT() && fStartT < rhs->GetStartT())) {
             return true;
         } else {
             return false;
         }
    } else {
         return false;
    }
  }

};

#endif
