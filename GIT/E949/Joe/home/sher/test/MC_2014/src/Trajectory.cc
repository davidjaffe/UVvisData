#include "Trajectory.hh"
#include "TrajectoryPoint.hh"
#include "G4ParticleTable.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"

//#define G4ATTDEBUG
#ifdef G4ATTDEBUG
#include "G4AttCheck.hh"
#endif

G4Allocator<Trajectory> aTrajectoryAllocator;

Trajectory::Trajectory()
    : fpPointsContainer(0), fTrackID(0), fParentID(0),
      PDGCharge(0.0), PDGEncoding(0), ParticleName(""),
      initialMomentum(G4ThreeVector()) {;}

Trajectory::Trajectory(const G4Track* aTrack)
{
    G4ParticleDefinition * fpParticleDefinition = aTrack->GetDefinition();
    ParticleName = fpParticleDefinition->GetParticleName();
    PDGCharge = fpParticleDefinition->GetPDGCharge();
    PDGEncoding = fpParticleDefinition->GetPDGEncoding();
    fTrackID = aTrack->GetTrackID();
    fParentID = aTrack->GetParentID();
    initialMomentum = aTrack->GetMomentum();
    fpPointsContainer = new TrajectoryPointContainer();
    // Following is for the first trajectory point
    fpPointsContainer->push_back(new TrajectoryPoint(aTrack));
}

Trajectory::Trajectory(Trajectory & right) : G4VTrajectory() 
{
    ParticleName = right.ParticleName;
    PDGCharge = right.PDGCharge;
    PDGEncoding = right.PDGEncoding;
    fTrackID = right.fTrackID;
    fParentID = right.fParentID;
    initialMomentum = right.initialMomentum;
    fpPointsContainer = new TrajectoryPointContainer();

    for(size_t i=0;i<right.fpPointsContainer->size();++i) {
        TrajectoryPoint* rightPoint
            = (TrajectoryPoint*)((*(right.fpPointsContainer))[i]);
        fpPointsContainer->push_back(new TrajectoryPoint(*rightPoint));
    }
}

Trajectory::~Trajectory()
{
    for(size_t i=0;i<fpPointsContainer->size();++i){
        delete  (*fpPointsContainer)[i];
    }
    fpPointsContainer->clear();

    delete fpPointsContainer;
}

void Trajectory::ShowTrajectory(std::ostream& os) const
{
    // Invoke the default implementation in G4VTrajectory...
    G4VTrajectory::ShowTrajectory(os);
    // ... or override with your own code here.
}

void Trajectory::DrawTrajectory(G4int i_mode) const
{
    // Invoke the default implementation in G4VTrajectory...
    G4VTrajectory::DrawTrajectory(i_mode);
    // ... or override with your own code here.
}

void Trajectory::AppendStep(const G4Step* aStep)
{
    fpPointsContainer->push_back(new TrajectoryPoint(aStep));
}

G4ParticleDefinition* Trajectory::GetParticleDefinition()
{
    return (G4ParticleTable::GetParticleTable()->FindParticle(ParticleName));
}

void Trajectory::MergeTrajectory(G4VTrajectory* secondTrajectory)
{
    if(!secondTrajectory) return;

    Trajectory* second = (Trajectory*)secondTrajectory;
    G4int ent = second->GetPointEntries();
    // initial point of the second trajectory should not be merged
    for(G4int i=1; i<ent; ++i) {
        fpPointsContainer->push_back((*(second->fpPointsContainer))[i]);
    }
    delete (*second->fpPointsContainer)[0];
    second->fpPointsContainer->clear();
}

const std::map<G4String,G4AttDef>* Trajectory::GetAttDefs() const
{
    G4bool isNew;
    std::map<G4String,G4AttDef>* store
        = G4AttDefStore::GetInstance("Trajectory",isNew);

    if (isNew) {

      G4String ID("ID");
      (*store)[ID] = G4AttDef(ID,"Track ID","Bookkeeping","","G4int");

      G4String PID("PID");
      (*store)[PID] = G4AttDef(PID,"Parent ID","Bookkeeping","","G4int");

      G4String PN("PN");
      (*store)[PN] = G4AttDef(PN,"Particle Name","Physics","","G4String");

      G4String Ch("Ch");
      (*store)[Ch] = G4AttDef(Ch,"Charge","Physics","e+","G4double");

      G4String PDG("PDG");
      (*store)[PDG] = G4AttDef(PDG,"PDG Encoding","Physics","","G4int");

      G4String IMom("IMom");
      (*store)[IMom] = G4AttDef(IMom,
                       "Momentum of track at start of trajectory",
                       "Physics","G4BestUnit","G4ThreeVector");

      G4String IMag("IMag");
      (*store)[IMag] = G4AttDef(IMag,
                       "Magnitude of momentum of track at start of trajectory",
                       "Physics","G4BestUnit","G4double");

        G4String NTP("NTP");
        (*store)[NTP] = G4AttDef(NTP,"No. of points","Bookkeeping","","G4int");

    }
    return store;
}

std::vector<G4AttValue>* Trajectory::CreateAttValues() const
{
  std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

  values->push_back
    (G4AttValue("ID",G4UIcommand::ConvertToString(fTrackID),""));

  values->push_back
    (G4AttValue("PID",G4UIcommand::ConvertToString(fParentID),""));

  values->push_back(G4AttValue("PN",ParticleName,""));

  values->push_back
    (G4AttValue("Ch",G4UIcommand::ConvertToString(PDGCharge),""));

  values->push_back
    (G4AttValue("PDG",G4UIcommand::ConvertToString(PDGEncoding),""));

  values->push_back
    (G4AttValue("IMom",G4BestUnit(initialMomentum,"Energy"),""));

  values->push_back
    (G4AttValue("IMag",G4BestUnit(initialMomentum.mag(),"Energy"),""));

  values->push_back
    (G4AttValue("NTP",G4UIcommand::ConvertToString(GetPointEntries()),""));

#ifdef G4ATTDEBUG
  G4cout << G4AttCheck(values,GetAttDefs());
#endif
    return values;
}
