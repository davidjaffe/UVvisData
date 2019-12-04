
#include "TrackingAction.hh"

#include "RunAction.hh"

#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4RunManager.hh"

// Default Constructor
//
TrackingAction::TrackingAction(RunAction* action) {

    runAction = action;

}

TrackingAction::~TrackingAction() {

}

// This function is called when a new track is taken off of the stack
//
void TrackingAction::PreUserTrackingAction(const G4Track* theTrack) {

    G4ThreeVector momentum =
    theTrack->GetDynamicParticle()->GetMomentum();

    G4ThreeVector position = theTrack->GetPosition();

    const G4String theParticleName =
    theTrack->GetDefinition()->GetParticleName();

    G4String theCreatorProcess;
    if ( theTrack->GetCreatorProcess() ) {
       theCreatorProcess = theTrack->GetCreatorProcess()->GetProcessName();
    }

    // If the particle is the primary pion...
    if (theParticleName == "pi+" && theTrack->GetParentID()== 0) {
        // Save information about this pion for diagnostics
        runAction->SavePionStart(position, momentum);
    }

    // If the particle is a muon created by the decay of the primary pion...
    if (theParticleName == "mu+" && 
        theTrack->GetParentID() == 1 &&
        theCreatorProcess == "Decay") {
          // Save information about this muon for diagnostics
          runAction->SaveMuonStart(position, momentum);
    }

    // If the particle is a positron that came from the primary pion
    // or the secondary muon...
    if (theParticleName == "e+" &&
       ((theCreatorProcess == "Decay" && theTrack->GetParentID()==1) ||
         theCreatorProcess == "DecayWithSpin")) {
          // Save information about this positron for diagnostics
          runAction->SavePositronStart(position, momentum);
    }
}

// This function is called when a particle decays or goes below
// threshold.
//
void TrackingAction::PostUserTrackingAction(const G4Track* theTrack) {

    G4String theVolume = theTrack->GetVolume()->GetName();

    const G4String theParticleName =
    theTrack->GetDefinition()->GetParticleName();

    G4ThreeVector position = theTrack->GetPosition();

    G4double kineticEnergy = theTrack->GetDynamicParticle()->GetKineticEnergy();

    G4ThreeVector momentum =
    theTrack->GetDynamicParticle()->GetMomentum();

    G4String theCreatorProcess;
    if ( theTrack->GetCreatorProcess() ) {
       theCreatorProcess = theTrack->GetCreatorProcess()->GetProcessName();
    }

    // If the particle is the primary pi+...
    if (theParticleName == "pi+" && theTrack->GetParentID()== 0) {
       // Did it stop in the target?
       if (theVolume == "Target" && kineticEnergy == 0) {
          // Save information about this pion for diagnostics
          runAction->SavePionStop(position);
        }
    }

    // If the particle is a mu+ that was created from the primary
    // pi+...
    if (theParticleName == "mu+" && theTrack->GetParentID() == 1  &&
        theCreatorProcess == "Decay" ) {
        // Did it stop in the target?
        if (theVolume == "Target" && kineticEnergy == 0) {
        }
    }

}
