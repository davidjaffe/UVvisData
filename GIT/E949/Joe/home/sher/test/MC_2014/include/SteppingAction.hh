
#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "globals.hh"

#include "G4UserSteppingAction.hh"

#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

class RunAction;

class SteppingAction : public G4UserSteppingAction {

    public:

        SteppingAction(RunAction* run);
        ~SteppingAction();

        void UserSteppingAction(const G4Step*);
        void InvokedProcesses(const G4Step*);

    private:

        G4int DecayInFlightCounter;

        RunAction* runAction;
        G4bool savedT1Once;
        G4bool savedT3Once;

};

#endif
