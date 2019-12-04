#ifndef EventAction_h
#define EventAction_h 1

#include "globals.hh"
#include "G4UserEventAction.hh"

#include <vector>

class G4Event;

class RunAction;
class EventActionMessenger;

class G4UImanager;

class HitSegment;

// This class provides user hooks at the beginning and end of an event.
//
class EventAction : public G4UserEventAction {

    public:

        // Constructor
        //
        EventAction(RunAction*);

        // Destructor
        //
        virtual ~EventAction();

        // Called at beginning of each event.
        //
        void  BeginOfEventAction(const G4Event*);

        // Called at end of each event.
        //
        void  EndOfEventAction(const G4Event*);

        // Used to sort the hit collection
        //
        static G4bool SortPredicate(const HitSegment* a,
                                    const HitSegment* b);

        void SetDrawFlag(G4String val) { drawFlag = val; };
        void SetPrintModulo(G4int val) { printModulo = val; };

        void AddEventToDebug(G4int val) {
             selectedEvents.push_back(val);
             nSelected++; };

    private:

        RunAction* runAction;
        EventActionMessenger* eventMessenger;

        G4int printModulo;

        // drawFlags = all, charged, neutral, charged+n
        G4String drawFlag;

        G4int nSelected;
        G4int iSelected;
        G4bool debugStarted;

        G4UImanager*          UI;
        std::vector<G4int>    selectedEvents;

};
#endif
