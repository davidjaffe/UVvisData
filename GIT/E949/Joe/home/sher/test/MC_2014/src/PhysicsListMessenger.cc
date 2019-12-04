
#include "PhysicsListMessenger.hh"
#include "PhysicsList.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "G4PionRadiativeDecayChannelMine.hh"

#include "G4ParticleTable.hh"

PhysicsListMessenger::PhysicsListMessenger(PhysicsList* pPhys)
    : fPhysicsList(pPhys)
{
    fDirectory = new G4UIdirectory("/exp/phys/");
    fDirectory->SetGuidance("Control the physics lists");

    fGammaCutCMD = new G4UIcmdWithADoubleAndUnit("/exp/phys/gammaCut",this);
    fGammaCutCMD->SetGuidance("Set gamma cut.");
    fGammaCutCMD->SetParameterName("Gcut",false);
    fGammaCutCMD->SetUnitCategory("Length");
    fGammaCutCMD->SetRange("Gcut>0.0");
    fGammaCutCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    fElectCutCMD = new G4UIcmdWithADoubleAndUnit("/exp/phys/electronCut",
                                                 this);
    fElectCutCMD->SetGuidance("Set electron cut.");
    fElectCutCMD->SetParameterName("Ecut",false);
    fElectCutCMD->SetUnitCategory("Length");
    fElectCutCMD->SetRange("Ecut>0.0");
    fElectCutCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    fPosCutCMD = new G4UIcmdWithADoubleAndUnit("/exp/phys/positronCut",
                                               this);
    fPosCutCMD->SetGuidance("Set positron cut.");
    fPosCutCMD->SetParameterName("Pcut",false);
    fPosCutCMD->SetUnitCategory("Length");
    fPosCutCMD->SetRange("Pcut>0.0");
    fPosCutCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    fAllCutCMD = new G4UIcmdWithADoubleAndUnit("/exp/phys/allCuts",this);
    fAllCutCMD->SetGuidance("Set cut for all.");
    fAllCutCMD->SetParameterName("cut",false);
    fAllCutCMD->SetUnitCategory("Length");
    fAllCutCMD->SetRange("cut>0.0");
    fAllCutCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    fDecayDirectory = new G4UIdirectory("/decay/");
    fDecayDirectory->SetGuidance("Decay chain control commands.");

    fPienuCMD = new G4UIcmdWithADouble("/decay/Pienu", this);
    //    fPienuCMD->SetGuidance("Sets the pi to decay into enu and enug ");
    fPienuCMD->SetGuidance("Set the branching ratio for the non-radiative channel. Default == 1.");
    fPienuCMD->SetParameterName("br",false);
    fPienuCMD->SetRange("br>0.&&br<=1.");
    fPienuCMD->SetDefaultValue((G4double)1.);
    fPienuCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    fPimunuCMD = new G4UIcmdWithADouble("/decay/Pimunu", this);
    //    fPimunuCMD->SetGuidance("Sets the pi to decay into munu and munug");
    fPimunuCMD->SetGuidance("Set the branching ratio for the non-radiative channel. Default == 1.");
    fPimunuCMD->SetParameterName("br",false);
    fPimunuCMD->SetRange("br>0.&&br<=1.");
    fPimunuCMD->SetDefaultValue((G4double)1.);
    fPimunuCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

}

PhysicsListMessenger::~PhysicsListMessenger()
{
    delete fGammaCutCMD;
    delete fElectCutCMD;
    delete fPosCutCMD;
    delete fAllCutCMD;
//    delete fPhysicsListCMD;
//    delete fListCMD;
}

void PhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{

    if (command == fPienuCMD) {
       // br = 0.999983     
       G4double br = fPienuCMD->GetNewDoubleValue(newValue);
       particleTable = G4ParticleTable::GetParticleTable();
       particleDef = particleTable->FindParticle("pi+");
       mode = new G4PhaseSpaceDecayChannel("pi+",br,2,"e+","nu_e");
       table=new G4DecayTable();
       table->Insert(mode);
       // br = 0.000017
       br = 1. - br;
       if (br >0.) {
          mode = new G4PionRadiativeDecayChannelMine("pi+",br);
          table->Insert(mode);
       }
       particleDef->SetDecayTable(table);
    }

    if (command == fPimunuCMD) {
       particleTable = G4ParticleTable::GetParticleTable();
       particleDef = particleTable->FindParticle("pi+");
       mode = new G4PhaseSpaceDecayChannel("pi+",1.0,2,"mu+","nu_mu");
       table=new G4DecayTable();
       table->Insert(mode);
       particleDef->SetDecayTable(table);
       G4double br = fPimunuCMD->GetNewDoubleValue(newValue);
       br = 1. - br;
       fPhysicsList->SetRadMuonDecayBR(br);
       fPhysicsList->ConstructParticle();
    }

    if (command == fGammaCutCMD) {
        fPhysicsList->SetCutForGamma(fGammaCutCMD
                                     ->GetNewDoubleValue(newValue));
    }
    else if (command == fElectCutCMD) {
        fPhysicsList->SetCutForElectron(fElectCutCMD
                                     ->GetNewDoubleValue(newValue));
    }
    else if (command == fPosCutCMD) {
        fPhysicsList->SetCutForPositron(fPosCutCMD
                                     ->GetNewDoubleValue(newValue));
    }
    else if (command == fAllCutCMD) {
        G4double cut = fAllCutCMD->GetNewDoubleValue(newValue);
        fPhysicsList->SetCutForGamma(cut);
        fPhysicsList->SetCutForElectron(cut);
        fPhysicsList->SetCutForPositron(cut);
    }
}
