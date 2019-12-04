#include "G4Step.hh"
#include "G4VProcess.hh"

#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Polyline.hh"
#include "G4Color.hh"

#include "G4TouchableHandle.hh"
#include "G4VPhysicalVolume.hh"

#include "G4EmSaturation.hh"
#include "G4LossTableManager.hh"

#ifdef G4VIS_USE
#include "G4TouchableHistory.hh"
#include "G4AffineTransform.hh"
#endif

#include "VolumeID.hh"
#include "HitSegment.hh"

G4Allocator<HitSegment> HitSegmentAllocator;

HitSegment::HitSegment()
  : fTrackID(-1), fParentID(-1),
    fEnergyDeposit(0), fPosition(0,0,0), fTrackLength(0),
    fStartX(0), fStartY(0), fStartZ(0), fStartT(0),
    fStopX(0),  fStopY(0),  fStopZ(0),  fStopT(0),
    fPDG(-1), fParticleName(""),fCreatorFlag(0),
    fCreatorProcessName(""), fEbirk(0),
    fMerged(false) {}

HitSegment::HitSegment(const HitSegment &right)
    : G4VHit(),
    fHitVolume(right.fHitVolume),
    fTrackID(right.fTrackID),
    fParentID(right.fParentID),
    fEnergyDeposit(right.fEnergyDeposit),
    fPosition(right.fPosition),
    fTrackLength(right.fTrackLength),
    fStartX(right.fStartX), fStartY(right.fStartY),
    fStartZ(right.fStartZ), fStartT(right.fStartT),
    fStopX(right.fStopX), fStopY(right.fStopY),
    fStopZ(right.fStopZ), fStopT(right.fStopT),
    fPDG(right.fPDG), fParticleName(right.fParticleName),
    fCreatorFlag(right.fCreatorFlag),
    fCreatorProcessName(right.fCreatorProcessName),
    fEbirk(right.fEbirk),
    fMerged(right.fMerged) {
#ifdef G4VIS_USE
    fRotation = right.fRotation;
    fTranslation = right.fTranslation;
#endif
    }

HitSegment::~HitSegment() {}

HitSegment HitSegment::operator=(const HitSegment& op2)
{
    fHitVolume = op2.fHitVolume;
    fTrackID = op2.fTrackID;
    fParentID = op2.fParentID;
    fEnergyDeposit = op2.fEnergyDeposit;
    fPosition = op2.fPosition;
    fTrackLength = op2.fTrackLength;
    fStartX = op2.fStartX;
    fStartY = op2.fStartY;
    fStartZ = op2.fStartZ;
    fStartT = op2.fStartT;
    fStopX  = op2.fStopX;
    fStopY  = op2.fStopY;
    fStopZ  = op2.fStopZ;
    fStopT  = op2.fStopT;
    fPDG    = op2.fPDG;
    fParticleName   = op2.fParticleName;
    fCreatorFlag   = op2.fCreatorFlag;
    fCreatorProcessName = op2.fCreatorProcessName;
    fEbirk  = op2.fEbirk;
    fMerged = op2.fMerged;
#ifdef G4VIS_USE
    fRotation    = op2.fRotation;
    fTranslation = op2.fTranslation;
#endif
    return *this;
}

void HitSegment::Clear()
{
  fTrackID = -1;
  fParentID = -1;
  fEnergyDeposit = 0;
  fStartT = 0;
  fStopT  = 0;
  fPDG    = -1;
  fEbirk  = 0;
}

HitSegment HitSegment::operator+=(const HitSegment& op2)
{
   fPosition *= fEnergyDeposit;
   fPosition += op2.GetEnergyDeposit()*op2.GetEnergyWeightedPosition();

   fEnergyDeposit += op2.GetEnergyDeposit();

   fPosition /= fEnergyDeposit;

   fTrackLength += op2.GetTrackLength();
   fEbirk       += op2.GetEbirk();

   if(op2.GetStopT() > fStopT )fStopT = op2.GetStopT();

   return *this;
}

G4bool HitSegment::SameHit(G4Step* theStep)
{
  // Check that the track IDs are the same

  G4int trackId = theStep->GetTrack()->GetTrackID();
  if (fTrackID != trackId) return false;

  // Check that the hit and new step are in the same volume

  const G4TouchableHandle& touchable
      = theStep->GetPreStepPoint()->GetTouchableHandle();

  if (fHitVolume != touchable) return false;

  return true;
}

void HitSegment::AddStep(G4Step* theStep)
{
  G4TouchableHandle touchable
      = theStep->GetPreStepPoint()->GetTouchableHandle();  

  G4ThreeVector prePos  = theStep->GetPreStepPoint()->GetPosition();
  G4ThreeVector postPos = theStep->GetPostStepPoint()->GetPosition();

  G4StepStatus stepStatus = theStep->GetPostStepPoint()->GetStepStatus();
  G4ParticleDefinition* particle =  theStep->GetTrack()->GetDefinition();

  G4double energyDeposit = theStep->GetTotalEnergyDeposit();
  G4double stepLength = theStep->GetStepLength();

  // Occasionally, a neutral particle will produce a particle below
  // threshold, and it will be recorded as generating the hit.  All of the
  // energy should be deposited at the stopping point of the track.

  if (stepStatus == fPostStepDoItProc
        && std::abs(particle->GetPDGCharge()) == 0.0) {
     prePos = postPos;
     stepLength = 0.0*cm;
  }

  /// First make sure the hit has been initialized.

  if (!fHitVolume) {

     fHitVolume = touchable;

     fTrackID = theStep->GetTrack()->GetTrackID();
     fParentID = theStep->GetTrack()->GetParentID();

     fStartX = prePos.x();
     fStartY = prePos.y();
     fStartZ = prePos.z();
     fStartT = theStep->GetPreStepPoint()->GetGlobalTime();
     fStopX = postPos.x();
     fStopY = postPos.y();
     fStopZ = postPos.z();
     fStopT = theStep->GetPostStepPoint()->GetGlobalTime();

     fPDG = theStep->GetTrack()->GetDefinition()->GetPDGEncoding();
     fParticleName = theStep->GetTrack()->GetDefinition()->GetParticleName();

     const G4VProcess* CreatorProcess =
                                   theStep->GetTrack()->GetCreatorProcess();
 
     if (CreatorProcess) {
       fCreatorProcessName = CreatorProcess->GetProcessName();
       G4ProcessType CreatorProcessType = CreatorProcess->GetProcessType();
       fCreatorFlag = CreatorProcessType + 1;
     } else{
       fCreatorFlag = 0;
     }

#ifdef G4VIS_USE
     const G4TouchableHistory* aHist = 
           (G4TouchableHistory*)(theStep->GetPreStepPoint()->GetTouchable());
     G4AffineTransform aTrans = aHist->GetHistory()->GetTopTransform();
     fRotation = aTrans.NetRotation();
     aTrans.Invert();
     fTranslation = aTrans.NetTranslation();
#endif

  } else {

     // Record the track that contributes to this hit.
     G4int trackId = theStep->GetTrack()->GetTrackID();
     G4int parentId = theStep->GetTrack()->GetParentID();

     if (trackId == fTrackID) {
        fStopX = postPos.x();
        fStopY = postPos.y();
        fStopZ = postPos.z();
        fStopT = theStep->GetPostStepPoint()->GetGlobalTime();
     } else {
        fTrackID = trackId;
        fParentID = parentId;
     }
  }

  fPosition *= fEnergyDeposit;

  fEnergyDeposit += energyDeposit;

  G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
  fEbirk         += emSaturation->VisibleEnergyDeposition(theStep);

  fTrackLength   += stepLength;

  fPosition += energyDeposit*((postPos+prePos)/2.);
  fPosition /= fEnergyDeposit;

}

#ifdef G4VIS_USE
void HitSegment::Draw(void)
{
  G4VPhysicalVolume* physVolume = fHitVolume.GetVolume();

  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if (pVVisManager) {
     G4VisAttributes attribs(G4Color(1.,0.,0.));
     attribs.SetForceSolid(true);
     G4Transform3D trans(G4Transform3D (fRotation, fTranslation));
     pVVisManager->Draw(*physVolume,attribs,trans);
  }
}
#endif

void HitSegment::Print(void)
{
    G4cout << "Hit Energy Deposit: "
           << G4BestUnit(fEnergyDeposit,"Energy")
           << "Hit Visible Energy: "
           << G4BestUnit(fEbirk,"Energy")
           << G4endl;
}
