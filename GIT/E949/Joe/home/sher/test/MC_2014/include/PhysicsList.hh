#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "globals.hh"
#include "G4VModularPhysicsList.hh"

class G4VPhysicsConstructor;
class PhysicsListMessenger;

class StepMax;

class PhysicsList : public G4VModularPhysicsList
{
public:

    PhysicsList(G4String);
    virtual ~PhysicsList();

    void SetCuts();
    void SetCutForGamma(G4double);
    void SetCutForElectron(G4double);
    void SetCutForPositron(G4double);

    void SetStepMax(G4double);
    StepMax* GetStepMaxProcess();
    void AddStepMax();

    void ConstructParticle();
    void ConstructProcess();

    /// Set the radiative muon decay branching ratio.
    void SetRadMuonDecayBR(G4double);

private:

    G4double fCutForGamma;
    G4double fCutForElectron;
    G4double fCutForPositron;

    G4double MaxChargedStep;
    StepMax* stepMaxProcess;
 
    PhysicsListMessenger* fMessenger;

    G4double radMuonDecayBR;

};

#endif
