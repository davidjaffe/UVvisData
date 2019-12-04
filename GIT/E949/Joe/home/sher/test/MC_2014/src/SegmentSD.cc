#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VProcess.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

#include "SegmentSD.hh"
#include "HitSegment.hh"

SegmentSD::SegmentSD(G4String name)
  : G4VSensitiveDetector(name), fHits(NULL), fHCID(-1), fLastHit(-1) {}

SegmentSD::~SegmentSD() {}

void SegmentSD::Initialize(G4HCofThisEvent* HCE)
{
  fHits =
     new HitSegmentCollection(GetName(), GetCollectionName(0));

  if (fHCID<0) {
      G4String hcName = GetName() + "/" + GetCollectionName(0);
      fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(hcName);
  }

  HCE->AddHitsCollection(fHCID, fHits);
}

G4bool SegmentSD::ProcessHits(G4Step* theStep, G4TouchableHistory*)
{
  G4double energyDeposit = theStep->GetTotalEnergyDeposit();

  if (energyDeposit <= 0.) return true;

  HitSegment* currentHit = NULL;

  G4StepPoint* preStepPoint = theStep->GetPreStepPoint();
  const G4VProcess* preProcess = preStepPoint->GetProcessDefinedStep();

  G4ProcessType preProcessType = fNotDefined;
  if (preProcess) preProcessType = preProcess->GetProcessType();

  // If the presStepPoint is of type fTransporation, then we need a new
  // hit since we are crossing a geometry boundary. On the other hand,
  // we look at the most recent hit stored in the collection to see if
  // it is the 'same hit' (same track ID and volume tree).

  if (preProcessType != fTransportation) {
     if (fHits) {
        fLastHit = fHits->entries()-1;
        if (fLastHit>-1) {
           HitSegment* tmpHit =
              dynamic_cast<HitSegment*>((*fHits)[fLastHit]);
           if (tmpHit->SameHit(theStep)){
              currentHit = tmpHit;
           }
        }
     }
  }

  if (!currentHit) {
     currentHit = GetNewHit();
     fHits->insert(currentHit);
     fLastHit = fHits->entries()-1;
  }

  currentHit->AddStep(theStep);

  return true;
}

void SegmentSD::EndOfEvent(G4HCofThisEvent* ) { }
