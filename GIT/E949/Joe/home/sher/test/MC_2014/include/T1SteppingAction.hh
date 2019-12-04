
#ifndef T1SteppingAction_h
#define T1SteppingAction_h 1

#include "globals.hh"

#include "G4UserSteppingAction.hh"

class RunAction;

class T1SteppingAction : public G4UserSteppingAction {

    public:

        T1SteppingAction(RunAction* run);
        ~T1SteppingAction();

        void UserSteppingAction(const G4Step*);

    private:

        // G4int DecayInFlightCounter;

//         RunAction* runAction;
//         G4bool savedT1Once;
//         G4bool savedT3Once;

};

#endif
