
#include "StackingAction.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4ios.hh"

StackingAction::StackingAction()
{}

StackingAction::~StackingAction()
{}

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  //keep primary particle
  if (aTrack->GetParentID() == 0) return fUrgent;

  // kill all secondaries
  // if(aTrack->GetParentID() != 0) return fKill;

  G4ParticleDefinition * particleType = aTrack->GetDefinition();

  if(particleType==G4NeutrinoE::NeutrinoEDefinition())return fKill;
  if(particleType==G4NeutrinoMu::NeutrinoMuDefinition())return fKill;
  if(particleType==G4NeutrinoTau::NeutrinoTauDefinition())return fKill;
  if(particleType==G4AntiNeutrinoE::AntiNeutrinoEDefinition())return fKill;
  if(particleType==G4AntiNeutrinoMu::AntiNeutrinoMuDefinition())return fKill;
  if(particleType==G4AntiNeutrinoTau::AntiNeutrinoTauDefinition())return fKill;

  return fUrgent;
}

void StackingAction::NewStage() {}

void StackingAction::PrepareNewEvent() {}

