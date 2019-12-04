#include "RunAction.hh"
#include "SteppingAction.hh"
//#include "T1SteppingAction.hh"

// #include "G4SteppingManager.hh"
// #include "G4Track.hh"
// #include "G4Step.hh"
// #include "G4ios.hh"
// #include "G4UnitsTable.hh"

#include "G4VProcess.hh"
#include "G4StepStatus.hh"
#include "G4Trajectory.hh"


SteppingAction::SteppingAction(RunAction* run) {

    runAction = run;
    savedT1Once = false;
    savedT3Once = false;
    
}



SteppingAction::~SteppingAction() { }

// This function is called each time Geant takes a step
//
void SteppingAction::UserSteppingAction(const G4Step* theStep) {
      
    G4String thePreVolume          = " ";
    G4String thePostVolume         = " ";
    G4String theProcessName        = " ";
    G4String theParticleName       = " ";
    G4String theCreatorProcessName = " ";

    G4ThreeVector preMomentum(0,0,0);
    G4ThreeVector prePosition(0,0,0);

    G4double preTime = 0;
    G4double preEnergy = 0;

    G4ThreeVector postMomentum(0,0,0);
    G4ThreeVector postPosition(0,0,0);

    G4double postTime = 0;
    G4double postEnergy = 0;
    
    G4double EnergyDeposit = 0;
    //G4double AngleT2 = 0;
    //G4double AngleB1 = 0;
    G4ThreeVector AngleB1;
    G4ThreeVector AngleT2;

    G4int anhlT1=0;
    G4int anhlT2=0;
    G4int anhlTg=0;
    G4int inelB1=0;
    G4int inelB1gamma=0;
    G4int inelB2=0;
    G4int inelB2gamma=0;
    G4int inelTg=0;
    G4int inelTggamma=0;
    G4int mucap=0;
    G4int pid=0;
    
    G4double ZrPreP;
    G4double ZrPostP;

    G4StepPoint* pPrePoint = theStep->GetPreStepPoint();
    if (pPrePoint) {
        G4VPhysicalVolume* prePhysVol = pPrePoint->GetPhysicalVolume();
        if (prePhysVol) thePreVolume = prePhysVol->GetName();
        prePosition = pPrePoint->GetPosition();
        preMomentum = pPrePoint->GetMomentum();
	preTime     = pPrePoint->GetGlobalTime();
	preEnergy   = pPrePoint->GetKineticEnergy();
    } else { return; }

    G4StepPoint* pPostPoint = theStep->GetPostStepPoint();
    if (pPostPoint) {
        G4VPhysicalVolume* postPhysVol = pPostPoint->GetPhysicalVolume();
        if (postPhysVol) thePostVolume = postPhysVol->GetName();
        postPosition = pPostPoint->GetPosition();
        postMomentum = pPostPoint->GetMomentum();
	postTime     = pPostPoint->GetGlobalTime();
	postEnergy   = pPostPoint->GetKineticEnergy();
        const G4VProcess* pProcess = pPostPoint->GetProcessDefinedStep();
        if (pProcess) theProcessName = pProcess->GetProcessName();
    } else { return; }

    G4Track* theTrack = theStep->GetTrack();
    if (theTrack) {
      theParticleName = theTrack->GetDefinition()->GetParticleName();
      const G4VProcess* theCreatorProcess = theTrack->GetCreatorProcess();
      if (theCreatorProcess) theCreatorProcessName = theCreatorProcess->GetProcessName();

    } else { return; }
    
    //Record the position of a photon-capture with neutron emission
    if (theParticleName == "neutron" && 
	theCreatorProcessName == "PhotonInelastic" &&
	theTrack->GetCurrentStepNumber() == 1 &&
	thePostVolume == "NaI" && thePreVolume == "NaI"
	){
      runAction->SPhotonuclear(postPosition);
    } 

    if (theParticleName == "pi+") {
      if (theTrack->GetCurrentStepNumber() == 1)
	runAction->SPionStart(prePosition,preTime,preMomentum,preEnergy);
      if (theProcessName == "Decay" | theProcessName == "DecayWithSpin")
	runAction->SPionDecay(postPosition,postTime,postMomentum,postEnergy);
    }
    if (theParticleName == "mu+" | theParticleName == "mu-") {
      if (theTrack->GetCurrentStepNumber() == 1) {
	runAction->SMuonStart(prePosition,preTime,preMomentum,preEnergy);
      }
      if (theProcessName == "Decay" | theProcessName == "DecayWithSpin"){
	runAction->SMuonDecay(postPosition,postTime,postMomentum,postEnergy);
    	runAction->SMuPolarization(theTrack->GetPolarization()); //Save the polarization of the muon
      }
    }
    
   
    if (theParticleName == "mu-" && theTrack->GetParentID() == 0 ){
      if (theTrack->GetTrackStatus() == fStopAndKill) {
       runAction->SMuonStop(postPosition,postTime,postMomentum,postEnergy);
       G4cout << "Muon Stop Process " << theProcessName << G4endl;
     }
   }

    // Muon energy
    // in TG
    if(theParticleName == "mu+" && thePreVolume ==  "/pienu/Target") {
      runAction->SMuETG(postEnergy);
     //  G4cout << "**********************" << G4endl;
//       G4cout << postEnergy << G4endl;
    }

    //Chloe's addition
    if (theParticleName == "e+" 
	&& thePostVolume == "/pienu/Target" 
	&& theProcessName == "eIoni") {
      EnergyDeposit = preEnergy - postEnergy;
      runAction->SBhabha(EnergyDeposit);
    }

    if (theParticleName == "e-"
	&& thePreVolume == "/pienu/Target"
	&& theCreatorProcessName == "eIoni"
	&& theTrack->GetCurrentStepNumber() == 1) {
      EnergyDeposit = preEnergy - postEnergy;
      runAction->SBhabhaEl(EnergyDeposit);
    }
  
    // Annihilations in flight in T1
    if (theParticleName == "e+" 
	&& thePostVolume == "/pienu/T1"
	&& theProcessName == "annihil") {
      anhlT1 = 1;
      runAction->SAnnihilT1(anhlT1);
    }
 
    // Annihilations in flight in T2
    if (theParticleName == "e+" 
	&& thePostVolume == "/pienu/T2"
	&& theProcessName == "annihil") {
      anhlT2 = 1;
      runAction->SAnnihilT2(anhlT2);
    }

    // Annihilations in flight in the target
    if (theParticleName == "e+" 
	&& thePostVolume == "/pienu/Target"
	&& theProcessName == "annihil"
	&& theTrack->GetCurrentStepNumber() == 1) {
      anhlTg = 1;
      runAction->SAnnihilTg(anhlTg);
    }
    
    // Inelastic scattering in B1
    if (theParticleName == "pi+"
	&& thePostVolume == "/pienu/B1"
	&& theProcessName == "PionPlusInelastic") {
      inelB1 = 1;
      runAction->SInelasticB1(inelB1);
    }

   //  // Inelastic scattering in B1
//     if (theParticleName == "neutron"
// 	&& thePostVolume == "/pienu/NaI"
// 	&& theProcessName == "NeutronInelastic") {
// 	//&& postEnergy>0) {
//       inelB1 = 1;
//       runAction->SInelasticB1(inelB1);
//     }

    if(theParticleName == "gamma"
       && thePreVolume == "/pienu/B1"
       && theCreatorProcessName == "PionPlusInelastic" 
       && theTrack->GetCurrentStepNumber() == 1) {
      inelB1gamma = 1;
      //InvokedProcesses(theStep);
      runAction->SInelasticB1gamma(postEnergy);
    }
       
    // Inelastic scattering in B2
    if (theParticleName == "pi+"
	&& thePostVolume == "/pienu/B2"
	&& theProcessName == "PionPlusInelastic") {
      inelB2 = 1;
      runAction->SInelasticB2(inelB2);
    }
    
    if(theParticleName == "gamma"
       && thePreVolume == "/pienu/B2"
       && theCreatorProcessName == "PionPlusInelastic"
       && theTrack->GetCurrentStepNumber() == 1) {
      inelB2gamma = 1;
      runAction->SInelasticB2gamma(postEnergy);
    }

    // Inelastic scattering in the target
    if (theParticleName == "pi+"
	&& thePostVolume == "/pienu/Target"
	&& theProcessName == "PionPlusInelastic") {
      inelTg = 1;
      runAction->SInelasticTg(inelTg);
    }

    if(theParticleName == "gamma"
       && thePostVolume == "/pienu/Target"
       && theCreatorProcessName == "PionPlusInelastic"
       && theTrack->GetCurrentStepNumber() == 1) {
      inelTggamma = 1;
      //InvokedProcesses(theStep);
      runAction->SInelasticTggamma(postEnergy);
      //G4cout << theCreatorProcessName << G4endl;
    }
    
    // Multiple Scattering in T2
    if (theParticleName == "e+"
	&& thePreVolume == "/pienu/T2"
	&& theProcessName == "msc") { 
      AngleT2[0] = preMomentum.theta() - postMomentum.theta();
      AngleT2[1] = atan2(postMomentum.x(),postMomentum.z());
      AngleT2[2] = acos(postMomentum.y()/sqrt(pow(postMomentum.x(),2)+pow(postMomentum.y(),2)+pow(postMomentum.z(),2)));
      runAction->SMscT2Angle(AngleT2);
    }

    //Multiple Scattering in B1
    if (theParticleName == "e+"
	&& thePostVolume == "/pienu/B1"
	&& theProcessName == "msc") {
      // InvokedProcesses(theStep);
      // AngleB1[0]: using deltaPhi() method of CHLEP ThreeVector class for azimuthal angle difference before and after scattering
      AngleB1[0] = preMomentum.theta() - postMomentum.theta();
      //AngleB1[0] = preMomentum.deltaPhi(postMomentum);
      //AngleB1[1]: azimuthal angle after scattering
      AngleB1[1] = atan2(postMomentum.x(),postMomentum.z());
      //AngleB1[2]: polar angle after scattering in [-pi/2,pi/2]
      AngleB1[2] = acos(postMomentum.y()/sqrt(pow(postMomentum.x(),2)+pow(postMomentum.y(),2)+pow(postMomentum.z(),2)));// - pi/2;
      //AngleB1[2] = atan(sqrt(pow(postMomentum.x(),2)+pow(postMomentum.y(),2)));
      runAction->SMscB1Angle(AngleB1);
      //G4cout << "x " << postMomentum.x() << " y " << postMomentum.y() << " z " << postMomentum.z() << G4endl;
    }
    
    //========================
    // Muon capture start
    //========================

    // Muon energy in Bina from beam muons (do they make it?)
    if (theParticleName == "mu-" && thePostVolume == "NaI" && theTrack->GetParentID() == 0) {
      runAction->SBeamMuonBina(prePosition,postTime,preMomentum,preEnergy);
    }
    
    // How often does it occur? 
    if (theParticleName == "mu-") {
      if(theProcessName == "muMinusCaptureAtRest") {
	mucap = 1;
      }
      if(theProcessName == "DecayWithSpin") {
	mucap = 2;}
      runAction->SMuCap(mucap);
    }
         
    //Which particles are produced by muon capture?
    //if (theCreatorProcessName=="muMinusCaptureAtRest" && theTrack->GetCurrentStepNumber() == 1) {
    if (theProcessName=="muMinusCaptureAtRest"){
      InvokedProcesses(theStep); 
      pid = theTrack->GetDefinition()->GetPDGEncoding();
      //G4cout << "PID " << pid << G4endl;
      runAction->SMuCapPID(pid);
    }
    
    // Position, time, momentum, energy right before capture...
    if (theParticleName == "mu-") {
      if(theProcessName == "muMinusCaptureAtRest") {
	runAction->SMuonCapture(prePosition,postTime,preMomentum,preEnergy);
      }
    }
    
    // Position, time, momentum, energy right before capture in zirconium...
    if (theParticleName == "mu-") {
      if(theProcessName == "muMinusCaptureAtRest") {
	if(thePreVolume == "Zr" && thePostVolume == "Zr") {
	  runAction->SMuonCaptureZr(prePosition,postTime,preMomentum,preEnergy);
	}
      }
    }

    // What does the energy spectrum of neutrons look like in Bina?
    if (theCreatorProcessName=="muMinusCaptureAtRest") {
      if (theParticleName == "neutron"){
	if (thePreVolume == "NaI" && thePostVolume == "NaI") {
	  runAction->SCaptureNeutron(prePosition,postTime,preMomentum,preEnergy);
	}
      }
    }

    // What is the multiplicity of neutrons after muon capture at rest?
    if (theCreatorProcessName=="muMinusCaptureAtRest") {
      if (theParticleName == "neutron"){
	
	
      }
    }

    // What does the energy spectrum of protons look like in T1?
    if (theCreatorProcessName=="muMinusCaptureAtRest") {
      if (theParticleName == "proton"){
	if (thePreVolume == "T1" && thePostVolume == "T1") {
	  runAction->SCaptureProton(prePosition,postTime,preMomentum,preEnergy);
	}
      }
    }

    //What does the energy spectrum of gammas look like in Bina?
    if (theCreatorProcessName=="muMinusCaptureAtRest") {
      if (theParticleName == "gamma"){
	if (thePreVolume == "NaI" && thePostVolume == "NaI") {
	  runAction->SCaptureGamma(prePosition,postTime,preMomentum,preEnergy);
	}
      }
    }
    
    // check momentum at the beginning and end of the zirconium foil
    if (theParticleName == "mu-" && theTrack->GetParentID() == 0 ){
      if (thePostVolume == "Zr" && thePreVolume != "Zr") {
	ZrPreP = preMomentum.z();
	runAction->SZrPreP(ZrPreP);}
      if (thePreVolume == "Zr" && thePostVolume != "Zr") {
	ZrPostP = postMomentum.z();
	runAction->SZrPostP(ZrPostP); }
    }

//     // check momentum and energy in five slices of the zirconium foil
//     if (theParticleName == "mu-" && theTrack->GetParentID() == 0 ){
//       G4double z = postPosition.z();
//       //G4cout << "Z coordinate " << z << G4endl;
//       // First slice
//       if (thePostVolume != "Zr1" && thePreVolume == "Zr1") {
// 	runAction->SZrFirstSlice(prePosition,preMomentum,preEnergy);
//       }
//       // Second slice
//       if (thePostVolume != "Zr2" && thePreVolume == "Zr2") {
// 	runAction->SZrSecondSlice(prePosition,preMomentum,preEnergy);
//       }
//       // Third slice
//       if (thePostVolume != "Zr3" && thePreVolume == "Zr3") {
// 	runAction->SZrThirdSlice(prePosition,preMomentum,preEnergy);
//       }
//       // Fourth slice
//       if (thePostVolume != "Zr4" && thePreVolume == "Zr4") {
// 	runAction->SZrFourthSlice(prePosition,preMomentum,preEnergy);
//       }
//       // Fifth slice
//       if (thePostVolume != "Zr5" && thePreVolume == "Zr5") {
// 	runAction->SZrFifthSlice(prePosition,preMomentum,preEnergy);
//       }
//     }

  //   // Number of muons stopping in Zr
//     if (theParticleName == "mu-" && theTrack->GetParentID() == 0 ){
//       if (theTrack->GetTrackStatus() == fStopAndKill) {
// 	if (thePostVolume == "Zr1" || thePostVolume == "Zr2" || thePostVolume == "Zr3" || thePostVolume == "Zr4" || thePostVolume == "Zr5") {
// 	  runAction->SMuStopZr(1);
// 	}
//       }
//     }

    
    // Number of muons stopping in Zr
    if (theParticleName == "mu-" && theTrack->GetParentID() == 0 ){
      if (theTrack->GetTrackStatus() == fStopAndKill) {
      	if (thePostVolume == "Zr") {
	  runAction->SMuStopZr(1);
	}
      }
    }

    // Number of muons stopping in T1
     if (theParticleName == "mu-" && theTrack->GetParentID() == 0 ){
       if (theTrack->GetTrackStatus() == fStopAndKill) {
       //if ( postEnergy == 0) {
	if (thePostVolume == "T1") {
	  runAction->SMuStopT1(1);
	}
      }
    }
    
     // Number of muons stopping in S3
    if (theParticleName == "mu-" && theTrack->GetParentID() == 0 ){
      if (theTrack->GetTrackStatus() == fStopAndKill) {
	if (thePostVolume == "SS3") {
	  runAction->SMuStopS3(1);
	}
      }
    }
    
    
    // Number of muons stopping in or after Zr
    if (theParticleName == "mu-" && theTrack->GetParentID() == 0 ){
      if (theTrack->GetTrackStatus() == fStopAndKill) {
	if (thePostVolume == "Zr" || thePostVolume == "T2" || thePostVolume == "NaI") {
	  runAction->SMuStopAfterZr(1);
	}
      }
    }

    
     // Energy of Muon Bremsstrahlung in Bina
    if (theCreatorProcessName=="muBrems") {
      runAction->SMuBrems(prePosition,postTime,preMomentum,preEnergy);
    }

    // Study signal process
    if (theParticleName == "gamma" && theCreatorProcessName=="RadiativeMuonCapture" && thePreVolume == "Zr" && thePostVolume == "Zr" ) {
      runAction->SGammaSignal(prePosition,postTime,preMomentum,preEnergy);
    }
      


    //***************************
    // Simulate Gamma ray signal
    //***************************

    // Where are the gammas generated? in the zirconium? -> check
    if (theParticleName == "gamma" && theTrack->GetCurrentStepNumber() == 1 && theTrack->GetParentID() == 0) {
      runAction->SGammaStart(prePosition,postTime,preMomentum,preEnergy);
    }

    // Gamma properties for those generated in the zirconium
    if (theParticleName == "gamma" && theTrack->GetCurrentStepNumber() == 1 && thePostVolume == "Zr") {
      runAction->SGammaStartZr(prePosition,postTime,preMomentum,preEnergy);
    }

    // What does the energy spectrum look like in Bina?
    if (theParticleName == "gamma" && theTrack->GetParentID() == 0 
	&& thePreVolume != "NaI" && thePostVolume == "NaI") {
      runAction->SGammaBina(prePosition,postTime,preMomentum,preEnergy);
      //G4cout << "CREATOR " << theCreatorProcessName << " kinE " << preEnergy <<  G4endl;
    }

    // What does it look like for electrons originating from the photons?
    if (theParticleName == "e-" || theParticleName == "e+"  
	&& thePreVolume == "NaI" && thePostVolume == "NaI") {
      runAction->SElectronBina(prePosition,postTime,preMomentum,preEnergy);
    }
    
    // using cross section & kinetic energy of muons
    if (theProcessName == "RadiativeMuonCapture") {
      runAction->SRadCap(prePosition,postTime,postMomentum,postEnergy);
    }
    
    // Muon Captre Background
    if (theParticleName == "gamma" && theCreatorProcessName == "RadiativeMuonCapture" && theTrack->GetCurrentStepNumber() == 1) {
      runAction->SGammaRadCap(postPosition,postTime,postMomentum,postEnergy);
    }


    //========================
    // Muon capture end
    //========================
    

    if (theParticleName == "pi+"
	&& thePostVolume == "/pienu/B1" ){
      InvokedProcesses(theStep);      
    }


    if (theParticleName == "e-" ) {
      if (theTrack->GetCurrentStepNumber() == 1 && theTrack->GetParentID() == 0){
	runAction->SElectronStart(prePosition,postTime,preMomentum,preEnergy);
      }
      if (theTrack->GetTrackStatus() == fStopAndKill) 
	runAction->SElectronStop(postPosition,postTime,postMomentum,postEnergy);
    }

    
    if (theParticleName == "e+" ) {
      //if (theCreatorProcessName=="DecayWithSpin" |
      //theCreatorProcessName=="Decay" ) {
      if (theTrack->GetCurrentStepNumber() == 1){
	runAction-> SPositronStart(prePosition,postTime,preMomentum,preEnergy);
	runAction->SEmomentum(preMomentum);//Save the decay positron polarization
      }
      if (theTrack->GetTrackStatus() == fStopAndKill) 
	runAction->SPositronStop(postPosition,postTime,postMomentum,postEnergy);
      //}
    }
       

    // If the particle is a pion...
    if (theParticleName == "pi+") {
        // If the particle is the primary pion...
        if (theTrack->GetParentID() == 0) {
	  // If the pi+ is stopped...
          if (theTrack->GetTrackStatus() == fStopAndKill) {
             savedT1Once = false;
             savedT3Once = false;
          }
        }
    }
    // If the particle is a positron...
    if (theParticleName == "e+") {
        // If the positron is in T3 
        if (thePreVolume == "T2") {
            if (savedT3Once == false) {
               runAction->SaveT3(prePosition, preMomentum);
               savedT3Once = true;
            }
        }
        // If the positron is in T1...
        if (thePreVolume == "T1") {
            if (savedT1Once == false) {
               runAction->SaveT1(prePosition, preMomentum);
               savedT1Once = true;
            }
        }
    }
    
}

void SteppingAction::InvokedProcesses(const G4Step* theStep) {
  G4cout << "**************** new Step ****************************" << G4endl;
  
  G4cout << "Step is limited by "
   << theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
   << G4endl;
  G4cout << "Processes involved to the step" << G4endl;
  G4StepStatus stepStatus = fpSteppingManager->GetfStepStatus();

  if(stepStatus==fAtRestDoItProc)
  {
    G4ProcessVector* procAtRest = fpSteppingManager->GetfAtRestDoItVector();
    G4SelectedAtRestDoItVector* selProcAtRest
     = fpSteppingManager->GetfSelectedAtRestDoItVector();
    size_t MAXofAtRestLoops = fpSteppingManager->GetMAXofAtRestLoops();
    for(size_t i1=0;i1<MAXofAtRestLoops;i1++)
    {
      if((*selProcAtRest)[MAXofAtRestLoops-i1-1]==2)
      { G4cout << "  At rest : " << (*procAtRest)[i1]->GetProcessName() << " (forced)" << G4endl; }
      else if((*selProcAtRest)[MAXofAtRestLoops-i1-1]==1)
      { G4cout << "  At rest : " << (*procAtRest)[i1]->GetProcessName() << G4endl; }
    }
  }

  if(stepStatus!=fExclusivelyForcedProc && stepStatus!=fAtRestDoItProc)
  {
    G4ProcessVector* procAlong = fpSteppingManager->GetfAlongStepDoItVector();
    size_t MAXofAlongStepLoops = fpSteppingManager->GetMAXofAlongStepLoops();
    for(size_t i2=0;i2<MAXofAlongStepLoops;i2++)
    {
      if((*procAlong)[i2]!=0)
      G4cout << "  Along step : " << (*procAlong)[i2]->GetProcessName() << G4endl;
    }
  }

  if(stepStatus!=fAtRestDoItProc)
  {
    G4ProcessVector* procPost = fpSteppingManager->GetfPostStepDoItVector();
    G4SelectedPostStepDoItVector* selProcPost
     = fpSteppingManager->GetfSelectedPostStepDoItVector();
    size_t MAXofPostStepLoops = fpSteppingManager->GetMAXofPostStepLoops();
    for(size_t i3=0;i3<MAXofPostStepLoops;i3++)
    {
      if((*selProcPost)[MAXofPostStepLoops-i3-1]==2)
      { G4cout << "  Post step : " << (*procPost)[i3]->GetProcessName() << " (forced)" << G4endl; }
      else if((*selProcPost)[MAXofPostStepLoops-i3-1]==1)
      { G4cout << "  Post step : " << (*procPost)[i3]->GetProcessName() << G4endl; }
    }
  }

  G4int nSecAtRest = fpSteppingManager->GetfN2ndariesAtRestDoIt();
  G4int nSecAlong  = fpSteppingManager->GetfN2ndariesAlongStepDoIt();
  G4int nSecPost   = fpSteppingManager->GetfN2ndariesPostStepDoIt();
  G4int nSecTotal  = nSecAtRest+nSecAlong+nSecPost;
  G4TrackVector* secVec = fpSteppingManager->GetfSecondary();

  if(nSecTotal>0)
  {
    G4cout << "  :----- List of 2ndaries - " << std::setw(3) << nSecTotal
           << " (Rest=" << std::setw(2) << nSecAtRest
           << ",Along=" << std::setw(2) << nSecAlong
           << ",Post="  << std::setw(2) << nSecPost << ")" << G4endl;

    for(size_t lp1=(*secVec).size()-nSecTotal; lp1<(*secVec).size(); lp1++)
    {
      G4cout << "    : "
             << G4BestUnit((*secVec)[lp1]->GetPosition(), "Length") << " "
             << std::setw( 9) << G4BestUnit((*secVec)[lp1]->GetKineticEnergy() , "Energy") << " "
             << std::setw(18) << (*secVec)[lp1]->GetDefinition()->GetParticleName()
             << " generated by " << (*secVec)[lp1]->GetCreatorProcess()->GetProcessName() << G4endl;
    }
  }
}



