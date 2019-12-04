
#ifndef TrackingAction_h
#define TrackingAction_h 1

#include "globals.hh"
#include "G4UserTrackingAction.hh"

class G4Track;
class RunAction;

class TrackingAction : public G4UserTrackingAction {

    public:
        TrackingAction(RunAction* action);
        ~TrackingAction();

        void PreUserTrackingAction(const G4Track*);
        void PostUserTrackingAction(const G4Track*);

    private:

	RunAction* runAction;

};
#endif
