#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

#include "G4LossTableManager.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4PhysListFactory.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "StepMax.hh"
#include "ExtraPhysics.hh"
//#include "G4OpticalPhysics.hh"

#include "G4ProcessTable.hh"
#include "G4ParticleTable.hh"

#include "G4PionDecayMakeSpin.hh"
#include "G4DecayWithSpinDepol.hh"

#include "G4DecayTable.hh"

#include "G4MuonDecayChannelWithSpin.hh"
#include "G4MuonDecayChannel.hh"
#include "G4MuonDecayChannelM.hh"

#include "G4RadiativeMuonCaptureProcess.hh"

#include "G4MuonRadiativeDecayChannelWithSpin.hh"
#include "RadiativeMuCapCrossSection.hh"


PhysicsList::PhysicsList(G4String physicsList) : G4VModularPhysicsList()
{
    G4LossTableManager::Instance();

    defaultCutValue  = 1.*mm;
    fCutForGamma     = 1*mm;
    fCutForElectron  = 0.1*mm;
    fCutForPositron  = 0.1*mm;

    radMuonDecayBR   = 0.014;

    fMessenger = new PhysicsListMessenger(this);

    SetVerboseLevel(0);

    G4PhysListFactory factory;
    G4VModularPhysicsList* phys = 0;

    if(factory.IsReferencePhysList(physicsList))
      phys =factory.GetReferencePhysList(physicsList);

    // Physics List is defined via environment variable PHYSLIST
    if(!phys) phys = factory.ReferencePhysList();

    for (G4int i = 0; ; ++i) {
       G4VPhysicsConstructor* elem =
                  const_cast<G4VPhysicsConstructor*> (phys->GetPhysics(i));
       if (elem == NULL) break;
       G4cout << "RegisterPhysics: " << elem->GetPhysicsName() << G4endl;
       RegisterPhysics(elem);
    }

    RegisterPhysics(new ExtraPhysics());
    //RegisterPhysics(new G4OpticalPhysics());

    stepMaxProcess = new StepMax();
}

PhysicsList::~PhysicsList()
{
    if (fMessenger) delete fMessenger;

    if (stepMaxProcess) delete stepMaxProcess;
}

void PhysicsList::ConstructParticle()
{
    G4VModularPhysicsList::ConstructParticle();

    G4DecayTable* MuonPlusDecayTable = new G4DecayTable();
    G4DecayTable* MuonMinusDecayTable = new G4DecayTable();

    G4double br = 1. - radMuonDecayBR;
    if (br > 0.) {

      
      MuonPlusDecayTable -> Insert(new
                           G4MuonDecayChannelWithSpin("mu+",br));
      MuonMinusDecayTable -> Insert(new
				    G4MuonDecayChannelWithSpin("mu-",br));
      
      G4cout << "****************** NOT RADIATIVE" << G4endl;
            
      /*
      //Massive Nu Decay
      MuonPlusDecayTable -> Insert(new
				    G4MuonDecayChannelM("mu+",br));
      MuonMinusDecayTable -> Insert(new
				    G4MuonDecayChannelM("mu-",br));
      */

      // // Radiative muon capture signal
//       MuonMinusDecayTable -> Insert(new 
// 				    G4RadiativeMuonCaptureProcess("mu-",br));

    }
    
    br = radMuonDecayBR;
    if (br > 0.) {
      
      MuonPlusDecayTable -> Insert(new
				   G4MuonRadiativeDecayChannelWithSpin("mu+",br));
      MuonMinusDecayTable -> Insert(new
				    G4MuonRadiativeDecayChannelWithSpin("mu-",br));      

      G4cout << "****************** RADIATIVE" << G4endl;      

      
     //  //Massive Nu Decay
//       MuonPlusDecayTable -> Insert(new
// 				   G4MuonDecayChannelM("mu+",br));
//       MuonMinusDecayTable -> Insert(new
// 				    G4MuonDecayChannelM("mu-",br));

      // // Radiative muon capture signal
//       MuonMinusDecayTable->Insert(new G4RadiativeMuonCaptureChannel("mu-",br));
      
    }
    
    G4MuonPlus::MuonPlusDefinition()   -> SetDecayTable(MuonPlusDecayTable);
    G4MuonMinus::MuonMinusDefinition() -> SetDecayTable(MuonMinusDecayTable);
}

void PhysicsList::ConstructProcess()
{
    G4VModularPhysicsList::ConstructProcess();

    G4DecayWithSpinDepol* decayWithSpin = new G4DecayWithSpinDepol();

    G4ProcessTable* processTable = G4ProcessTable::GetProcessTable();

    G4VProcess* decay;
    decay = processTable->FindProcess("Decay",G4MuonPlus::MuonPlus());

    G4ProcessManager* fManager;
    fManager = G4MuonPlus::MuonPlus()->GetProcessManager();

    if (fManager) {
      if (decay) fManager->RemoveProcess(decay);
      fManager->AddProcess(decayWithSpin);
      // set ordering for PostStepDoIt and AtRestDoIt
      fManager ->SetProcessOrdering(decayWithSpin, idxPostStep);
      fManager ->SetProcessOrdering(decayWithSpin, idxAtRest);
    }

    decay = processTable->FindProcess("Decay",G4MuonMinus::MuonMinus());

    fManager = G4MuonMinus::MuonMinus()->GetProcessManager();

   //  G4HadronicProcess *RadiativeCapture = new G4RadiativeMuonCaptureProcess();
    
//     RadiativeCapture->AddDataSet(RadiativeCrossSection);

    if (fManager) {
      if (decay) fManager->RemoveProcess(decay);
      fManager->AddProcess(decayWithSpin);
      // set ordering for PostStepDoIt and AtRestDoIt
      fManager ->SetProcessOrdering(decayWithSpin, idxPostStep);
      fManager ->SetProcessOrdering(decayWithSpin, idxAtRest);
      //if (G4RadiativeMuonCaptureProcess->IsApplicable(const G4ParticleDefinition&) {
	   //  fManager->AddProcess(RadiativeCapture);
// 	    fManager->SetProcessOrdering(RadiativeCapture,idxPostStep);
					 //  }
    }

    G4PionDecayMakeSpin* poldecay = new G4PionDecayMakeSpin();

    decay = processTable->FindProcess("Decay",G4PionPlus::PionPlus());

    fManager = G4PionPlus::PionPlus()->GetProcessManager();

    if (fManager) {
      if (decay) fManager->RemoveProcess(decay);
      fManager->AddProcess(poldecay);
      // set ordering for PostStepDoIt and AtRestDoIt
      // fManager ->SetProcessOrdering(poldecay, idxPostStep);
//       fManager ->SetProcessOrdering(poldecay, idxAtRest);
    }

    decay = processTable->FindProcess("Decay",G4PionMinus::PionMinus());

    fManager = G4PionMinus::PionMinus()->GetProcessManager();

    if (fManager) {
      if (decay) fManager->RemoveProcess(decay);
      fManager->AddProcess(poldecay);
      // set ordering for PostStepDoIt and AtRestDoIt
      fManager ->SetProcessOrdering(poldecay, idxPostStep);
      fManager ->SetProcessOrdering(poldecay, idxAtRest);
    }

    AddStepMax();
}

void PhysicsList::SetCuts()
{
    if (verboseLevel >0) {
        G4cout << "PhysicsList::SetCuts:";
        G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length")
               << G4endl;
    }

    // set cut values for gamma at first and for e- second and next for e+,
    // because some processes for e+/e- need cut values for gamma
    SetCutValue(fCutForGamma, "gamma");
    SetCutValue(fCutForElectron, "e-");
    SetCutValue(fCutForPositron, "e+");

    if (verboseLevel>0) DumpCutValuesTable();
}

void PhysicsList::SetCutForGamma(G4double cut)
{
    fCutForGamma = cut;
    SetParticleCuts(fCutForGamma, G4Gamma::Gamma());
}

void PhysicsList::SetCutForElectron(G4double cut)
{
    fCutForElectron = cut;
    SetParticleCuts(fCutForElectron, G4Electron::Electron());
}

void PhysicsList::SetCutForPositron(G4double cut)
{
    fCutForPositron = cut;
    SetParticleCuts(fCutForPositron, G4Positron::Positron());
}

void PhysicsList::SetStepMax(G4double step)
{
  MaxChargedStep = step ;
  stepMaxProcess->SetStepMax(MaxChargedStep);
}

StepMax* PhysicsList::GetStepMaxProcess()
{
  return stepMaxProcess;
}

void PhysicsList::AddStepMax()
{
  // Step limitation seen as a process

  theParticleIterator->reset();
  while ((*theParticleIterator)()){
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();

      if (stepMaxProcess->IsApplicable(*particle) && !particle->IsShortLived())
      {
         if (pmanager) pmanager ->AddDiscreteProcess(stepMaxProcess);
      }
  }
}

void PhysicsList::SetRadMuonDecayBR(G4double br)
{
  radMuonDecayBR = br;
}
