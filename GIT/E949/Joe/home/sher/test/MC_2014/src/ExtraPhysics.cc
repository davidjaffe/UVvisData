#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4StepLimiter.hh"
#include "MaxTimeCuts.hh"
#include "MinEkineCuts.hh"

#include "ExtraPhysics.hh"

ExtraPhysics::ExtraPhysics() 
    : G4VPhysicsConstructor("Extra") { }

ExtraPhysics::~ExtraPhysics() { }

void ExtraPhysics::ConstructParticle() { }

void ExtraPhysics::ConstructProcess()
{
    G4cout << "ExtraPhysics:: Add Extra Physics Processes" 
              << G4endl;

    theParticleIterator->reset();

    while ((*theParticleIterator)()) {
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();
        G4double charge = particle->GetPDGCharge();

        if (!pmanager) {
            std::ostringstream o;
            o << "Particle " << particleName << "without a Process Manager";
            G4Exception("ExtraPhysics::ConstructProcess()","",
                         FatalException,o.str().c_str());
        }

        // All charged particles should have a step limiter to make sure that
        // the steps do not get too long.
        if (charge != 0.0) {
            pmanager->AddDiscreteProcess(new G4StepLimiter("Step Limit"));
        }
        // time cuts for ONLY neutrons:
        if(particleName == "neutron") {
            pmanager->AddDiscreteProcess(new MaxTimeCuts());
        }
        // Energy cuts to kill charged (embedded in method) particles:
        // pmanager->AddDiscreteProcess(new MinEkineCuts());
    }
}
