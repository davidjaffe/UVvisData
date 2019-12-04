
#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "G4ios.hh"

#include "TFile.h"
#include "TNtuple.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TTree.h"

#include "RunActionMessenger.hh"

#include "Randomize.hh"

#include <iomanip>
#include <string>
#include <iostream>
#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
// for lock file
#include <fcntl.h> // for open()
#include <cerrno> // for errno 
#include <cstdio> // for perror()

// Default Constructor
//
RunAction::RunAction()
  : autoSeed(false)
{
    runActionMessenger = new RunActionMessenger(this);

    t1energyDeposit = 0;
    t1pathLength = 0;
    t3energyDeposit = 0;
    t3pathLength = 0;

    pionStopPosition.setX(0);      // DvB: changed to match version 2.1.3.1 of CLHEP
    pionStopPosition.setY(0);
    pionStopPosition.setZ(0);
}

// Default Destructor
//
RunAction::~RunAction() {

    delete runActionMessenger;

}

// Open a ROOT file and create an Ntuple
//

void RunAction::OpenRoot() {


  int fd;

  do{
  fd=open("password.lck", O_WRONLY | O_CREAT | O_EXCL);
  
  //  G4cout << "We are opening the lock file"<<G4endl;
  //  G4cout << "fd = " << fd<<G4endl;
  }
  while (fd<0);

  
  // Intelligently select filenames to use for output
  std::ifstream fin;
  bool isOpen = false;
    int i = 1;
    std::string filename = "";
    std::string extension=".root";
    //char env;
    
    while (isOpen == false) {
      
        filename.clear();
        filename.insert(0,itoa(i,10));
	char *env = getenv("MCOUTPUT");
	filename.insert(0,env);
        filename.append(extension);
        fin.open(filename.c_str(),std::ios::in);

        if (fin.fail()) {

            std::cout << "Saving to: " << filename << std::endl;
            isOpen = true;

        } else {

            i++;

        }

        fin.close();
    }

    // Create a new Root file
    file = new TFile(filename.c_str(),"RECREATE");

    remove("password.lck");
    // Create TNtuple to store the hits
    hitTuple = new TNtuple("hits","Hits","eventID:volumeID0:volumeID1:volumeID2:energyDeposit:startX:startY:startZ:startT:stopX:stopY:stopZ:stopT:Ebirk:PID:CreatorProcess",40000);

    aTree =new TTree("tree","pienu");
    aTree = ((TTree*)gROOT->FindObject("tree"));
    aTree->Branch("eventID",&eventID,"eventID/I");
    aTree->Branch("PiStartX",PStartX,"PiStartX[4]/D");
    aTree->Branch("PiStartP",PStartP,"PiStartP[4]/D");
    aTree->Branch("PiDecayX",PDecayX,"PiDecayX[4]/D");
    aTree->Branch("PiDecayP",PDecayP,"PiDecayP[4]/D");
    aTree->Branch("MuStartX",MStartX,"MuStartX[4]/D");
    aTree->Branch("MuStartP",MStartP,"MuStartP[4]/D");
    aTree->Branch("MuDecayX",MDecayX,"MuDecayX[4]/D");
    aTree->Branch("MuDecayP",MDecayP,"MuDecayP[4]/D");
    aTree->Branch("MuStopX",MStopX,"MuStopX[4]/D");
    aTree->Branch("MuStopP",MStopP,"MuStopP[4]/D");
    aTree->Branch("PoStartX",EStartX,"PoStartX[4]/D");
    aTree->Branch("PoStartP",EStartP,"PoStartP[4]/D");
    aTree->Branch("PoStopX",EStopX,"PoStopX[4]/D");
    aTree->Branch("PoStopP",EStopP,"PoStopP[4]/D");

    aTree->Branch("PhotonuclearX",EPhNX,"PhotonuclearX[3]/D");
    aTree->Branch("BhabhaE",&EBh,"BhabhaE/D");
    aTree->Branch("MuETG",&EMuTG,"MuETG/D");
    aTree->Branch("BhabhaEEl",&EBhEl,"BhaBhaEEl/D");
    aTree->Branch("AnnihilT1",&AnhlT1,"AnnihilT1/I");
    aTree->Branch("AnnihilT2",&AnhlT2,"AnnihilT2/I");
    aTree->Branch("AnnihilTg",&AnhlTg,"AnnihilTg/I");
    aTree->Branch("InelasticB1",&InelB1,"InelasticB1/I");
    aTree->Branch("InelasticB1gammaE",&InelB1gammaE,"InelasticB1gammaE/D");
    aTree->Branch("InelasticB2",&InelB2,"InelasticB2/I");
    aTree->Branch("InelasticB2gammaE",&InelB2gammaE,"InelasticB2gammaE/D");
    aTree->Branch("InelasticTg",&InelTg,"InelasticTg/I");
    aTree->Branch("InelasticTggammaE",&InelTggammaE,"InelastigTggammaE/D");
    aTree->Branch("MscT2Angle",&MscT2Angle,"MscT2Angle[3]/D");
    aTree->Branch("MscB1Angle",&MscB1Angle,"MscB1Angle[3]/D");
    aTree->Branch("MuPolarization",MuPol,"MuPolarization[3]/D");
    aTree->Branch("EMomentum",EMom,"EMomentum[3]/D");
    aTree->Branch("MuCap",&MuCap,"MuCap/I");
    aTree->Branch("MuCapPID",&MuCapPID,"MuCapPID/I");
    aTree->Branch("MuStopZr",&MuStopZr,"MuStopZr/I");
    aTree->Branch("MuStopT1",&MuStopT1,"MuStopT1/I");
    aTree->Branch("MuStopS3",&MuStopS3,"MuStopS3/I");
    aTree->Branch("MuStopAfterZr",&MuStopAfterZr,"MuStopAfterZr/I");
    aTree->Branch("ZrPreP",&MuZrPreP,"ZrPreP/D");
    aTree->Branch("ZrPostP",&MuZrPostP,"ZrPostP/D");
    aTree->Branch("ZrFirstSliceP",ZrFirstSliceP,"ZrFirstSliceP[4]/D");
    aTree->Branch("ZrFirstSliceX",ZrFirstSliceX,"ZrFirstSliceX[4]/D");
    aTree->Branch("ZrSecondSliceP",ZrSecondSliceP,"ZrSecondSliceP[4]/D");
    aTree->Branch("ZrSecondSliceX",ZrSecondSliceX,"ZrSecondSliceX[4]/D");
    aTree->Branch("ZrThirdSliceP",ZrThirdSliceP,"ZrThirdSliceP[4]/D");
    aTree->Branch("ZrThirdSliceX",ZrThirdSliceX,"ZrThirdSliceX[4]/D");
    aTree->Branch("ZrFourthSliceP",ZrFourthSliceP,"ZrFourthSliceP[4]/D");
    aTree->Branch("ZrFourthSliceX",ZrFourthSliceX,"ZrFourthSliceX[4]/D");
    aTree->Branch("ZrFifthSliceP",ZrFifthSliceP,"ZrFifthSliceP[4]/D");
    aTree->Branch("ZrFifthSliceX",ZrFifthSliceX,"ZrFifthSliceX[4]/D");
    aTree->Branch("MuCapX",MuCapX,"MuCapX[4]/D");
    aTree->Branch("MuCapP",MuCapP,"MuCapP[4]/D");
    aTree->Branch("BeamMuonBinaX",BeamMuonBinaX,"BeamMuonBinaX[4]/D");
    aTree->Branch("BeamMuonBinaP",BeamMuonBinaP,"BeamMuonBinaP[4]/D");
    aTree->Branch("MuCapZrX",MuCapZrX,"MuCapZrX[4]/D");
    aTree->Branch("MuCapZrP",MuCapZrP,"MuCapZrP[4]/D");
    aTree->Branch("CapNeutronX",CapNeutronX,"CapNeutronX[4]/D");
    aTree->Branch("CapNeutronP",CapNeutronP,"CapNeutronP[4]/D");
    aTree->Branch("CapProtonX",CapProtonX,"CapProtonX[4]/D");
    aTree->Branch("CapProtonP",CapProtonP,"CapProtonP[4]/D");
    aTree->Branch("CapGammaX",CapGammaX,"CapGammaX[4]/D");
    aTree->Branch("CapGammaP",CapGammaP,"CapGammaP[4]/D");
    aTree->Branch("ElectronStartX",ElectronStartX,"ElectronStartX[4]/D");
    aTree->Branch("ElectronStartP",ElectronStartP,"ElectronStartP[4]/D");
    aTree->Branch("ElectronStopX",ElectronStopX,"ElectronStopX[4]/D");
    aTree->Branch("ElectronStopP",ElectronStopP,"ElectronStopP[4]/D");
    aTree->Branch("GammaSignalX",GammaSignalX,"GammaSignalX[4]/D");
    aTree->Branch("GammaSignalP",GammaSignalP,"GammaSignalP[4]/D");
    aTree->Branch("MuBremsX",MuBremsX,"MuBremsX[4]/D");
    aTree->Branch("MuBremsP",MuBremsP,"MuBremsP[4]/D");
    aTree->Branch("GammaStartX",GammaStartX,"GammaStartX[4]/D");
    aTree->Branch("GammaStartP",GammaStartP,"GammaStartP[4]/D");
    aTree->Branch("GammaStartZrX",GammaStartZrX,"GammaStartZrX[4]/D");
    aTree->Branch("GammaStartZrP",GammaStartZrP,"GammaStartZrP[4]/D");
    aTree->Branch("GammaBinaX",GammaBinaX,"GammaBinaX[4]/D");
    aTree->Branch("GammaBinaP",GammaBinaP,"GammaBinaP[4]/D");
    aTree->Branch("ElectronBinaX",ElectronBinaX,"ElectronBinaX[4]/D");
    aTree->Branch("ElectronBinaP",ElectronBinaP,"ElectronBinaP[4]/D");
    aTree->Branch("RadCapX",RadCapX,"RadCapX[4]/D");
    aTree->Branch("RadCapP",RadCapP,"RadCapP[4]/D");
    aTree->Branch("GammaRadCapX",GammaRadCapX,"GammaRadCapX[4]/D");
    aTree->Branch("GammaRadCapP",GammaRadCapP,"GammaRadCapP[4]/D");


    //Diagnostic Histograms

    // Particle Start Position/Momentum and Stop Position
    pionStartR = new TH1D("pionStartR", "Pion Start Position R Coordinate", 40, 30, 30);
    pionStartZ = new TH1D("pionStartZ", "Pion Start Position Z Coordinate", 20, -175, -125);
    pionStartPR = new TH1D("pionStartPR", "Pion Start Momentum R Component", 300, -75, 75);
    pionStartPZ = new TH1D("pionStartPZ", "Pion Start Momentum Z Component", 200, 50, 100);
    pionStartPRP = new TH1D("pionStartPRP", "Pion Start Momentum Pr/P", 100, -1, 1);
    pionStartPZP = new TH1D("pionStartPZP", "Pion Start Momentum Pz/P", 100, -1, 1);
    pionStopR = new TH1D("pionStopR", "Pion Stop Position R Coordinate", 50, 0, 50);
    pionStopZ = new TH1D("pionStopZ", "Pion Stop Position Z Coordinate", 20, -0.4, 0.6);
    muonStartR = new TH1D("muonStartR", "Muon Start Position R Coordinate", 50, 0, 50);
    muonStartZ = new TH1D("muonStartZ", "Muon Start Position Z Coordinate", 20, -0.4, 0.6);
    muonStartPR = new TH1D("muonStartPR", "Muon Start Momentum R Component", 300, -75, 75);
    muonStartPZ = new TH1D("muonStartPZ", "Muon Start Momentum Z Component", 300, -75, 75);
    muonStartPRP = new TH1D("muonStartPRP", "Muon Start Momentum Pr/P", 100, -1, 1);
    muonStartPZP = new TH1D("muonStartPZP", "Muon Start Momentum Pz/P", 100, -1, 1);
    positronStartR = new TH1D("positronStartR", "Positron Start Position R Coordinate", 50, 0, 50);
    positronStartZ = new TH1D("positronStartZ", "Positron Start Position Z Coordinate", 20, -0.4, 0.6);
    positronStartPR = new TH1D("positronStartPR", "Positron Start Momentum R Component", 300, -75, 75);
    positronStartPZ = new TH1D("positronStartPZ", "Positron Start Momentum Z Component", 300, -75, 75);
    positronStartPRP = new TH1D("positronStartPRP", "Positron Start Momentum Pr/P", 120, -1.2, 1.2);
    positronStartPZP = new TH1D("positronStartPZP", "Positron Start Momentum Pz/P", 50, -1, 1);
    triggeredPositronStartZ = new TH1D("triggeredPositronStartZ", "Triggered Positron Start Position Z Coordinate", 20, -0.4, 0.6);
    triggeredPositronStartPZ = new TH1D("triggeredPositronStartPZ", "Triggered Positron Start Momentum Z Component", 75, -75, 75);
    triggeredPositronStartPZP = new TH1D("triggeredPositronStartPZP", "Triggered Positron Start Momentum Pz/P", 50, -1, 1);

    // All T1 Histograms require hits in B1, B2, Target, and T1 to fill
    positronT1R = new TH1D("positronT1R", "Positron T1 Position R Coordinate", 50, 0, 10);
    positronT1Z = new TH1D("positronT1Z", "Positron T1 Position Z Coordinate", 50, 0, 7);
    positronT1PR = new TH1D("positronT1PR", "Positron T1 Momentum R Component", 300, -75, 75);
    positronT1PZ = new TH1D("positronT1PZ", "Positron T1 Momentum Z Component", 300, -75, 75);
    positronT1PRP = new TH1D("positronT1PRP", "Positron T1 Momentum Pr/P", 100, -1, 1);
    positronT1PZP = new TH1D("positronT1PZP", "Positron T1 Momentum Pz/P", 100, -1, 1);
    positronT1ExpectedR = new TH1D("positronT1ExpectedR", "Expected Radius at center of T1 from positron start direction", 50, 0, 10);
    positronT1DifferenceR = new TH1D("positronT1DifferenceR", "Difference between expected and actual radial coordinates", 50,0,50);
    positronT1DifferencePZP = new TH1D("positronT1DifferencePZP", "Difference between start PZ/P and PZ/P at T1 center", 20, 0, 0.8);
    positronT1dEdX = new TH1D("positronT1dEdX", "dE/dX for positrons in T1", 100, 0, 0.4);
    positronT1dEdXvPathLength = new TH2D("positronT1dEdXvPathLength", "Edep vs. Path Length for positrons in T1", 100, 0, 10, 100, 0, 6);

    // All T3 Histograms require hits in B1, B2, Target, T1, and T3 to fill
    positronT3R = new TH1D("positronT3R", "Positron T3 Position R Coordinate", 50, 0, 20);
    positronT3Z = new TH1D("positronT3Z", "Positron T3 Position Z Coordinate", 50, 50, 70);
    positronT3PR = new TH1D("positronT3PR", "Positron T3 Momentum R Component", 300, -75, 75);
    positronT3PZ = new TH1D("positronT3PZ", "Positron T3 Momentum Z Component", 300, -75, 75);
    positronT3PRP = new TH1D("positronT3PRP", "Positron T3 Momentum Pr/P", 50, -1, 1);
    positronT3PZP = new TH1D("positronT3PZP", "Positron T3 Momentum Pz/P", 50, -1, 1);
    positronT3ExpectedR = new TH1D("positronT3ExpectedR", "Expected Radius at entrance to T3 from positron start direction", 50, 0, 50);
    positronT3DifferenceR = new TH1D("positronT3DifferenceR", "Difference between expected and actual radial coordinates at entrance to T3", 50,0,50);
    positronT3DifferencePZP = new TH1D("positronT3DifferencePZP", "Difference between start Pz/P and PZ/P at T3 center", 20, 0, 0.8);
    positronT3dEdX = new TH1D("positronT3dEdX", "dE/dX for positrons in T3", 100, 0, 0.4);
    positronT3dEdXvPathLength = new TH2D("positronT3dEdXvPathLength","Edep vs. Path Length for positrons in T3", 100, 0, 10, 100, 0, 6);
    triggeredPositronT3ExpectedR = new TH1D("triggeredPositronT3ExpectedR", "Expected Radius at entrance to T3 for triggered events", 50, 0, 50);
    triggeredPositronT3DifferenceR = new TH1D("triggeredPositronT3DifferenceR", "Difference between expected and actual radial coordinates at entrance to T3", 50, 0, 50);
    triggeredPositronT3DifferencePZP = new TH1D("triggeredPositronT3DifferencePZP", "Difference between start Pz/P and Pz/P at T3 center for triggered events", 20, 0, 0.8);

    hitTuple->SetAutoSave(32000);
    aTree->SetAutoSave(32000);
}

// Write any remaining ROOT objects and close the file.
//
void RunAction::CloseRoot() {

    hitTuple->AutoSave();

    // Particle Start Position/Momentum and Stop Position
    pionStartR->Write(); // Pion Start Position R Coordinate
    pionStartZ->Write(); // Pion Start Position Z Coordinate
    pionStartPR->Write(); // Pion Start Momentum R Component
    pionStartPZ->Write(); // Pion Start Momentum Z Component
    pionStartPRP->Write(); // Pion Start Momentum Pr/P
    pionStartPZP->Write(); // Pion Start Momentum Pz/P
    pionStopR->Write(); // Pion Stop Position R Coordinate
    pionStopZ->Write(); // Pion Stop Position Z Coordinate
    muonStartR->Write(); // Muon Start Position R Coordinate (Same As Pion Stop Coordinate - Redundant)
    muonStartZ->Write(); // Muon Start Position Z Coordinate (Same As Pion Stop Coordinate - Redundant)
    muonStartPR->Write(); // Muon Start Momentum R Component
    muonStartPZ->Write(); // Muon Start Momentum Z Component
    muonStartPRP->Write(); // Muon Start Momentum Pr/P
    muonStartPZP->Write(); // Muon Start Momentum Pz/P
    positronStartR->Write(); // Positron Start Position R Coordinate
    positronStartZ->Write(); // Positron Start Position Z Coordinate
    positronStartPR->Write(); // Positron Start Momentum R Component
    positronStartPZ->Write(); // Positron Start Momentum Z Component
    positronStartPRP->Write(); // Positron Start Momentum Pr/P
    positronStartPZP->Write(); // Positron Start Momentum Pz/P
    triggeredPositronStartZ->Write(); // Positron Start Position Z Coordinate
    triggeredPositronStartPZ->Write(); // Positron Start Momentum Z Component
    triggeredPositronStartPZP->Write(); // Positron Start Momentum Pz/P

    // All T1 Histograms require hits in B1, B2, Target, and T1 to fill
    positronT1R->Write(); // Positron T1 Position R Coordinate
    positronT1Z->Write(); // Positron T1 Position Z Coordinate
    positronT1PR->Write(); // Positron T1 Momentum R Component
    positronT1PZ->Write(); // Positron T1 Momentum Z Component
    positronT1PRP->Write(); // Positron T1 Momentum Pr/P
    positronT1PZP->Write(); // Positron T1 Momentum Pz/P
    positronT1ExpectedR->Write(); //Expected Radius at center of T1 from positron start direction
    positronT1DifferenceR->Write(); // Difference between expected and actual radial coordinates
    positronT1DifferencePZP->Write(); // Difference between start PZ/P and PZ/P at T1 center
    positronT1dEdX->Write(); // dE/dX for positrons in T1
    positronT1dEdXvPathLength->Write(); // dE/dX vs. Path Length for positrons in T1

    // All T3 Histograms require hits in B1, B2, Target, T1, and T3 to fill
    positronT3R->Write(); // Positron T3 Position R Coordinate
    positronT3Z->Write(); // Positron T3 Position Z Coordinate
    positronT3PR->Write(); // Positron T3 Momentum R Component
    positronT3PZ->Write(); // Positron T3 Momentum Z Component
    positronT3PRP->Write(); // Positron T3 Momentum Pr/P
    positronT3PZP->Write(); // Positron T3 Momentum Pz/P
    positronT3ExpectedR->Write(); // Expected Radius at center of T3 from positron start direction
    positronT3DifferenceR->Write(); // Difference between expected and actual radial coordinates
    positronT3DifferencePZP->Write(); // Difference between start PZ/P and PZ/P at T1 center
    positronT3dEdX->Write(); // dE/dX for positrons in T3
    positronT3dEdXvPathLength->Write(); // dE/dX vs. Path Length for positrons in T3
    triggeredPositronT3ExpectedR->Write();
    triggeredPositronT3DifferenceR->Write();
    triggeredPositronT3DifferencePZP->Write();

    aTree->Write();
    file->Close();
}

// Fill the hit Ntuple
//

void RunAction::FillTuple(G4double E1, G4double E2, G4double E3, G4double E4,
                          G4double E5, G4double E6, G4double E7, G4double E8,
                          G4double E9, G4double E10, G4double E11,
			  G4double E12, G4double E13, G4double E14,
                          G4double E15, G4double E16) {
  
  E16 = 0;
  hitTuple->Fill(E1,E2,E3,E4,E5,E6,E7,E8,E9,E10,E11,E12,E13,E14,E15);

}

void RunAction::FillTree(){aTree->Fill();}

void RunAction::SPionStart(G4ThreeVector position, 
                           G4double time, 
                           G4ThreeVector momentum, 
                           G4double Energy)
{
  PStartX[0]=position.x();
  PStartX[1]=position.y();
  PStartX[2]=position.z();
  PStartX[3]=time;
  PStartP[0]=momentum.x();
  PStartP[1]=momentum.y();
  PStartP[2]=momentum.z();
  PStartP[3]=Energy;
}

void RunAction::SPionDecay(G4ThreeVector position,
                           G4double time,
                           G4ThreeVector momentum,
                           G4double Energy)
{
  PDecayX[0]=position.x();
  PDecayX[1]=position.y();
  PDecayX[2]=position.z();
  PDecayX[3]=time;
  PDecayP[0]=momentum.x();
  PDecayP[1]=momentum.y();
  PDecayP[2]=momentum.z();
  PDecayP[3]=Energy;
}

void RunAction::SMuonStart(G4ThreeVector position,
                           G4double time,
                           G4ThreeVector momentum,
                           G4double Energy)
{
  MStartX[0]=position.x();
  MStartX[1]=position.y();
  MStartX[2]=position.z();
  MStartX[3]=time;
  MStartP[0]=momentum.x();
  MStartP[1]=momentum.y();
  MStartP[2]=momentum.z();
  MStartP[3]=Energy;
}

void RunAction::SMuonDecay(G4ThreeVector position,
                           G4double time,
                           G4ThreeVector momentum,
                           G4double Energy)
{
  MDecayX[0]=position.x();
  MDecayX[1]=position.y();
  MDecayX[2]=position.z();
  MDecayX[3]=time;
  MDecayP[0]=momentum.x();
  MDecayP[1]=momentum.y();
  MDecayP[2]=momentum.z();
  MDecayP[3]=Energy;
}

void RunAction::SMuonStop(G4ThreeVector position,
                           G4double time,
                           G4ThreeVector momentum,
                           G4double Energy) {
  MStopX[0]=position.x();
  MStopX[1]=position.y();
  MStopX[2]=position.z();
  MStopX[3]=time;
  MStopP[0]=momentum.x();
  MStopP[1]=momentum.y();
  MStopP[2]=momentum.z();
  MStopP[3]=Energy;
}

void RunAction::SPositronStart(G4ThreeVector position,
                               G4double time,
                               G4ThreeVector momentum,
                               G4double Energy)
{
  EStartX[0]=position.x();
  EStartX[1]=position.y();
  EStartX[2]=position.z();
  EStartX[3]=time;
  EStartP[0]=momentum.x();
  EStartP[1]=momentum.y();
  EStartP[2]=momentum.z();
  EStartP[3]=Energy;
}

void RunAction::SPositronStop(G4ThreeVector position,
                              G4double time,
                              G4ThreeVector momentum,
                              G4double Energy)
{
  EStopX[0]=position.x();
  EStopX[1]=position.y();
  EStopX[2]=position.z();
  EStopX[3]=time;
  EStopP[0]=momentum.x();
  EStopP[1]=momentum.y();
  EStopP[2]=momentum.z();
  EStopP[3]=Energy;
}

void RunAction::SPhotonuclear(G4ThreeVector position)
{
  EPhNX[0]=position.x();
  EPhNX[1]=position.y();
  EPhNX[2]=position.z();
}
    
void RunAction::SMuETG(G4double Energy)
{
  EMuTG=Energy;
}

void RunAction::SBhabha(G4double Energy)
{
  EBh=Energy;
 
}

void RunAction::SBhabhaEl(G4double EnergyEl)
{
  EBhEl=EnergyEl;
}


void RunAction::SAnnihilT1(G4int anhlT1)
{
  AnhlT1=anhlT1;
}

void RunAction::SAnnihilT2(G4int anhlT2)
{
  AnhlT2=anhlT2;
}

void RunAction::SAnnihilTg(G4int anhlTg)
{
  AnhlTg=anhlTg;
}


void RunAction::SInelasticB1(G4int inelB1)
{
  InelB1=inelB1;
}

void RunAction::SInelasticB1gamma(G4double EnergyGammaB1)
{
  InelB1gammaE=EnergyGammaB1;
}


void RunAction::SInelasticB2(G4int inelB2)
{
  InelB2=inelB2;
}

void RunAction::SInelasticB2gamma(G4double EnergyGammaB2)
{
  InelB2gammaE=EnergyGammaB2;
}

void RunAction::SInelasticTg(G4int inelTg)
{
  InelTg=inelTg;
}

void RunAction::SInelasticTggamma(G4double EnergyGammaTg)
{
  InelTggammaE=EnergyGammaTg;
}

void RunAction::SMscT2Angle(G4ThreeVector AngleT2)
{
  MscT2Angle=AngleT2;
}

void RunAction::SMscB1Angle(G4ThreeVector AngleB1)
{
  MscB1Angle=AngleB1;
}


void RunAction::SMuCap(G4int mucap)
{
  MuCap=mucap;
}


void RunAction::SRadCap(G4ThreeVector position,
			    G4double time,
			    G4ThreeVector momentum,
			    G4double Energy)
{
  RadCapX[0]=position.x();
  RadCapX[1]=position.y();
  RadCapX[2]=position.z();
  RadCapX[3]=time;
  RadCapP[0]=momentum.x();
  RadCapP[1]=momentum.y();
  RadCapP[2]=momentum.z();
  RadCapP[3]=Energy;
}      

 void RunAction::SGammaRadCap(G4ThreeVector position,
			    G4double time,
			    G4ThreeVector momentum,
			    G4double Energy)
{
  GammaRadCapX[0]=position.x();
  GammaRadCapX[1]=position.y();
  GammaRadCapX[2]=position.z();
  GammaRadCapX[3]=time;
  GammaRadCapP[0]=momentum.x();
  GammaRadCapP[1]=momentum.y();
  GammaRadCapP[2]=momentum.z();
  GammaRadCapP[3]=Energy;
}       

void RunAction::SGammaSignal(G4ThreeVector position,
				G4double time,
				G4ThreeVector momentum,
				G4double Energy)
{
  GammaSignalX[0]=position.x();
  GammaSignalX[1]=position.y();
  GammaSignalX[2]=position.z();
  GammaSignalX[3]=time;
  GammaSignalP[0]=momentum.x();
  GammaSignalP[1]=momentum.y();
  GammaSignalP[2]=momentum.z();
  GammaSignalP[3]=Energy;
}

void RunAction::SMuCapPID(G4int pid)
{
  MuCapPID=pid;
}
            
void RunAction::SZrPreP(G4double zrprep)
{
  MuZrPreP = zrprep;
}

void RunAction::SZrPostP(G4double zrpostp)
{
  MuZrPostP = zrpostp;
}

void RunAction::SZrFirstSlice(G4ThreeVector position,
			      G4ThreeVector momentum,
			      G4double Energy)
{
  ZrFirstSliceX[0]=position.x();
  ZrFirstSliceX[1]=position.y();
  ZrFirstSliceX[2]=position.z();
  ZrFirstSliceP[0]=momentum.x();
  ZrFirstSliceP[1]=momentum.y();
  ZrFirstSliceP[2]=momentum.z();
  ZrFirstSliceP[3]=Energy;
}
  
void RunAction::SZrSecondSlice(G4ThreeVector position,
			       G4ThreeVector momentum,
			       G4double Energy)
{
  ZrSecondSliceX[0]=position.x();
  ZrSecondSliceX[1]=position.y();
  ZrSecondSliceX[2]=position.z();
  ZrSecondSliceP[0]=momentum.x();
  ZrSecondSliceP[1]=momentum.y();
  ZrSecondSliceP[2]=momentum.z();
  ZrSecondSliceP[3]=Energy;
}

void RunAction::SZrThirdSlice(G4ThreeVector position,
			      G4ThreeVector momentum,
			      G4double Energy)
{
  ZrThirdSliceX[0]=position.x();
  ZrThirdSliceX[1]=position.y();
  ZrThirdSliceX[2]=position.z();
  ZrThirdSliceP[0]=momentum.x();
  ZrThirdSliceP[1]=momentum.y();
  ZrThirdSliceP[2]=momentum.z();
  ZrThirdSliceP[3]=Energy;
}

void RunAction::SZrFourthSlice(G4ThreeVector position,
			       G4ThreeVector momentum,
			       G4double Energy)
{
  ZrFourthSliceX[0]=position.x();
  ZrFourthSliceX[1]=position.y();
  ZrFourthSliceX[2]=position.z();
  ZrFourthSliceP[0]=momentum.x();
  ZrFourthSliceP[1]=momentum.y();
  ZrFourthSliceP[2]=momentum.z();
  ZrFourthSliceP[3]=Energy;
}

void RunAction::SZrFifthSlice(G4ThreeVector position,
			      G4ThreeVector momentum,
			      G4double Energy)
{
  ZrFifthSliceX[0]=position.x();
  ZrFifthSliceX[1]=position.y();
  ZrFifthSliceX[2]=position.z();
  ZrFifthSliceP[0]=momentum.x();
  ZrFifthSliceP[1]=momentum.y();
  ZrFifthSliceP[2]=momentum.z();
  ZrFifthSliceP[3]=Energy;
}

void RunAction::SMuonCapture(G4ThreeVector position,
			     G4double time,
			     G4ThreeVector momentum,
			     G4double Energy)
{
  MuCapX[0]=position.x();
  MuCapX[1]=position.y();
  MuCapX[2]=position.z();
  MuCapX[3]=time;
  MuCapP[0]=momentum.x();
  MuCapP[1]=momentum.y();
  MuCapP[2]=momentum.z();
  MuCapP[3]=Energy;
}

void RunAction::SBeamMuonBina(G4ThreeVector position,
			     G4double time,
			     G4ThreeVector momentum,
			     G4double Energy)
{
  BeamMuonBinaX[0]=position.x();
  BeamMuonBinaX[1]=position.y();
  BeamMuonBinaX[2]=position.z();
  BeamMuonBinaX[3]=time;
  BeamMuonBinaP[0]=momentum.x();
  BeamMuonBinaP[1]=momentum.y();
  BeamMuonBinaP[2]=momentum.z();
  BeamMuonBinaP[3]=Energy;
}

void RunAction::SMuStopZr(G4int stop)
{
  MuStopZr = stop;
}

void RunAction::SMuStopT1(G4int stop)
{
  MuStopT1 = stop;
}

void RunAction::SMuStopS3(G4int stop)
{
  MuStopS3 = stop;
}


void RunAction::SMuStopAfterZr(G4int stop)
{
  MuStopAfterZr = stop;
}

void RunAction::SMuonCaptureZr(G4ThreeVector position,
			       G4double time,
			       G4ThreeVector momentum,
			       G4double Energy)
{
  MuCapZrX[0]=position.x();
  MuCapZrX[1]=position.y();
  MuCapZrX[2]=position.z();
  MuCapZrX[3]=time;
  MuCapZrP[0]=momentum.x();
  MuCapZrP[1]=momentum.y();
  MuCapZrP[2]=momentum.z();
  MuCapZrP[3]=Energy;
}

void RunAction::SCaptureNeutron(G4ThreeVector position,
				G4double time,
				G4ThreeVector momentum,
				G4double Energy)
{
  CapNeutronX[0]=position.x();
  CapNeutronX[1]=position.y();
  CapNeutronX[2]=position.z();
  CapNeutronX[3]=time;
  CapNeutronP[0]=momentum.x();
  CapNeutronP[1]=momentum.y();
  CapNeutronP[2]=momentum.z();
  CapNeutronP[3]=Energy;
}

void RunAction::SCaptureProton(G4ThreeVector position,
				G4double time,
				G4ThreeVector momentum,
				G4double Energy)
{
  CapProtonX[0]=position.x();
  CapProtonX[1]=position.y();
  CapProtonX[2]=position.z();
  CapProtonX[3]=time;
  CapProtonP[0]=momentum.x();
  CapProtonP[1]=momentum.y();
  CapProtonP[2]=momentum.z();
  CapProtonP[3]=Energy;
}
void RunAction::SCaptureGamma(G4ThreeVector position,
				G4double time,
				G4ThreeVector momentum,
				G4double Energy)
{
  CapGammaX[0]=position.x();
  CapGammaX[1]=position.y();
  CapGammaX[2]=position.z();
  CapGammaX[3]=time;
  CapGammaP[0]=momentum.x();
  CapGammaP[1]=momentum.y();
  CapGammaP[2]=momentum.z();
  CapGammaP[3]=Energy;
}

void RunAction::SMuBrems(G4ThreeVector position,
				G4double time,
				G4ThreeVector momentum,
				G4double Energy)
{
  MuBremsX[0]=position.x();
  MuBremsX[1]=position.y();
  MuBremsX[2]=position.z();
  MuBremsX[3]=time;
  MuBremsP[0]=momentum.x();
  MuBremsP[1]=momentum.y();
  MuBremsP[2]=momentum.z();
  MuBremsP[3]=Energy;
}
           
void RunAction::SGammaStart(G4ThreeVector position,
			    G4double time,
			    G4ThreeVector momentum,
			    G4double Energy)
{
  GammaStartX[0]=position.x();
  GammaStartX[1]=position.y();
  GammaStartX[2]=position.z();
  GammaStartX[3]=time;
  GammaStartP[0]=momentum.x();
  GammaStartP[1]=momentum.y();
  GammaStartP[2]=momentum.z();
  GammaStartP[3]=Energy;
}

void RunAction::SGammaStartZr(G4ThreeVector position,
			      G4double time,
			      G4ThreeVector momentum,
			      G4double Energy)
{
  GammaStartZrX[0]=position.x();
  GammaStartZrX[1]=position.y();
  GammaStartZrX[2]=position.z();
  GammaStartZrX[3]=time;
  GammaStartZrP[0]=momentum.x();
  GammaStartZrP[1]=momentum.y();
  GammaStartZrP[2]=momentum.z();
  GammaStartZrP[3]=Energy;
}


void RunAction::SGammaBina(G4ThreeVector position,
			    G4double time,
			    G4ThreeVector momentum,
			    G4double Energy)
{
  GammaBinaX[0]=position.x();
  GammaBinaX[1]=position.y();
  GammaBinaX[2]=position.z();
  GammaBinaX[3]=time;
  GammaBinaP[0]=momentum.x();
  GammaBinaP[1]=momentum.y();
  GammaBinaP[2]=momentum.z();
  GammaBinaP[3]=Energy;
}        

 void RunAction::SElectronStart(G4ThreeVector position,
			    G4double time,
			    G4ThreeVector momentum,
			    G4double Energy)
{
  ElectronStartX[0]=position.x();
  ElectronStartX[1]=position.y();
  ElectronStartX[2]=position.z();
  ElectronStartX[3]=time;
  ElectronStartP[0]=momentum.x();
  ElectronStartP[1]=momentum.y();
  ElectronStartP[2]=momentum.z();
  ElectronStartP[3]=Energy;
}      

 void RunAction::SElectronStop(G4ThreeVector position,
			    G4double time,
			    G4ThreeVector momentum,
			    G4double Energy)
{
  ElectronStopX[0]=position.x();
  ElectronStopX[1]=position.y();
  ElectronStopX[2]=position.z();
  ElectronStopX[3]=time;
  ElectronStopP[0]=momentum.x();
  ElectronStopP[1]=momentum.y();
  ElectronStopP[2]=momentum.z();
  ElectronStopP[3]=Energy;
}      



void RunAction::SElectronBina(G4ThreeVector position,
			    G4double time,
			    G4ThreeVector momentum,
			    G4double Energy)
{
  ElectronBinaX[0]=position.x();
  ElectronBinaX[1]=position.y();
  ElectronBinaX[2]=position.z();
  ElectronBinaX[3]=time;
  ElectronBinaP[0]=momentum.x();
  ElectronBinaP[1]=momentum.y();
  ElectronBinaP[2]=momentum.z();
  ElectronBinaP[3]=Energy;
}        

void RunAction::SMuPolarization(G4ThreeVector mupolarization)
{
  //G4cout << "Pol = " << mupolarization << G4endl;
  MuPol[0] = mupolarization.x();
  MuPol[1] = mupolarization.y();
  MuPol[2] = mupolarization.z();
}
      
void RunAction::SEmomentum(G4ThreeVector emomentum)
{
  //G4cout << "Mom = " << emomentum << G4endl;
  EMom[0] = emomentum.x();
  EMom[1] = emomentum.y();
  EMom[2] = emomentum.z();

  double N1=0;
  double N2=0;
  double Q = 0;
  for (int i=0;i<3;i++){
    N1 += MuPol[i]*MuPol[i];
    N2 += EMom[i]*EMom[i];
    
    Q += MuPol[i]*EMom[i];
  }
  
  Q /= sqrt(N1*N2);
  
  //G4cout << Q << G4endl;
  
  EMom[0] = Q;
  EMom[1] = emomentum.y();
  EMom[2] = emomentum.z();

}    

void RunAction::ClearVariable(){
  G4ThreeVector tV(-10000,-10000,-10000);
  G4double temp=-10000;
  RunAction::SPionStart(tV,temp,tV,temp);
  RunAction::SPionDecay(tV,temp,tV,temp);
  RunAction::SMuonStart(tV,temp,tV,temp);
  RunAction::SMuonDecay(tV,temp,tV,temp);
  RunAction::SMuonStop(tV,temp,tV,temp);
  RunAction::SPositronStart(tV,temp,tV,temp);
  RunAction::SPositronStop(tV,temp,tV,temp);
  RunAction::SPhotonuclear(tV);
  RunAction::SBhabha(temp);
  RunAction::SBhabhaEl(temp);
  RunAction::SAnnihilT1(0);
  RunAction::SAnnihilT2(0);
  RunAction::SAnnihilTg(0);
  RunAction::SInelasticB1(0);
  RunAction::SInelasticB1gamma(temp);
  RunAction::SInelasticB2(0);
  RunAction::SInelasticB2gamma(temp);
  RunAction::SInelasticTg(0);
  RunAction::SInelasticTggamma(temp);
  RunAction::SMscT2Angle(tV);
  RunAction::SMscB1Angle(tV);
  RunAction::SMuCap(0);
  RunAction::SRadCap(tV,temp,tV,temp);
  RunAction::SGammaRadCap(tV,temp,tV,temp);
  RunAction::SGammaSignal(tV,temp,tV,temp);
  RunAction::SMuCapPID(0);
  RunAction::SMuStopZr(0);
  RunAction::SMuStopT1(0);
  RunAction::SMuStopS3(0);
  RunAction::SMuStopAfterZr(0);
  RunAction::SMuonCapture(tV,temp,tV,temp);
  RunAction::SBeamMuonBina(tV,temp,tV,temp);
  RunAction::SMuonCaptureZr(tV,temp,tV,temp);
  RunAction::SCaptureNeutron(tV,temp,tV,temp);
  RunAction::SCaptureProton(tV,temp,tV,temp);
  RunAction::SCaptureGamma(tV,temp,tV,temp);
  RunAction::SZrPreP(temp);
  RunAction::SZrPostP(temp);
  RunAction::SZrFirstSlice(tV,tV,temp);
  RunAction::SZrSecondSlice(tV,tV,temp);
  RunAction::SZrThirdSlice(tV,tV,temp);
  RunAction::SZrFourthSlice(tV,tV,temp);
  RunAction::SZrFifthSlice(tV,tV,temp);
  RunAction::SMuBrems(tV,temp,tV,temp);
  RunAction::SGammaStart(tV,temp,tV,temp);
  RunAction::SGammaStartZr(tV,temp,tV,temp);
  RunAction::SGammaBina(tV,temp,tV,temp);
  RunAction::SElectronBina(tV,temp,tV,temp);
  RunAction::SElectronStart(tV,temp,tV,temp);
  RunAction::SElectronStop(tV,temp,tV,temp);
}

void RunAction::SavePionStart(G4ThreeVector position, G4ThreeVector momentum)
{
  pionStartPosition = position;
  pionStartMomentum = momentum;
}

void RunAction::SavePionStop(G4ThreeVector position)
{
  pionStopPosition = position;
}

void RunAction::SaveMuonStart(G4ThreeVector position, G4ThreeVector momentum)
{
  muonStartPosition = position;
  muonStartMomentum = momentum;
}

void RunAction::SavePositronStart(G4ThreeVector position, 
                                  G4ThreeVector momentum)
{
  positronStartPosition = position;
  positronStartMomentum = momentum;
}

void RunAction::SaveT1(G4ThreeVector position, G4ThreeVector momentum)
{
  t1Position = position;
  t1Momentum = momentum;
}

void RunAction::SaveT1dEdX(G4double energyDeposit, G4double pathLength)
{
  t1energyDeposit = energyDeposit;
  t1pathLength = pathLength;
}

void RunAction::SaveT3(G4ThreeVector position, G4ThreeVector momentum)
{
  t3Position = position;
  t3Momentum = momentum;
}

void RunAction::SaveT3dEdX(G4double energyDeposit, G4double pathLength)
{
  t3energyDeposit = energyDeposit;
  t3pathLength = pathLength;
}

void RunAction::FillUntriggeredHistos() {
    // Particle Start Position/Momentum and Stop Position
    pionStartR->Fill( std::sqrt( std::pow(pionStartPosition.x(),2) +
                                 std::pow(pionStartPosition.y(),2) )/10. );
    pionStartZ->Fill(pionStartPosition.z()/10.);
    pionStartPR->Fill( std::sqrt( std::pow(pionStartMomentum.x(),2) +
                                  std::pow(pionStartMomentum.y(),2) ) );
    pionStartPZ->Fill(pionStartMomentum.z());
    pionStartPRP->Fill( std::sqrt( std::pow(pionStartMomentum.x(),2) +
                                   std::pow(pionStartMomentum.y(),2) )/
                                             pionStartMomentum.mag() );
    pionStartPZP->Fill(pionStartMomentum.z()/pionStartMomentum.mag());
}

void RunAction::FillTargetStopHistos() {

  if (pionStopPosition.x() != 0 && pionStopPosition.y() != 0 &&pionStopPosition.z() != 0 ) {// DvB: changed to match version 2.1.3.1 of CLHEP
       pionStopR->Fill(std::sqrt(std::pow(pionStopPosition.x(),2) +
                                 std::pow(pionStopPosition.y(),2))/10.);
       pionStopZ->Fill(pionStopPosition.z()/10.);
       muonStartR->Fill(std::sqrt(std::pow(muonStartPosition.x(),2) +
                                  std::pow(muonStartPosition.y(),2))/10.);
       muonStartZ->Fill(muonStartPosition.z()/10.);
       muonStartPR->Fill(std::sqrt(std::pow(muonStartMomentum.x(),2) +
                                   std::pow(muonStartMomentum.y(),2)));
       muonStartPZ->Fill(muonStartMomentum.z());
       muonStartPRP->Fill(std::sqrt(std::pow(muonStartMomentum.x(),2) +
                                    std::pow(muonStartMomentum.y(),2))/
                                                muonStartMomentum.mag());
       muonStartPZP->Fill(muonStartMomentum.z()/muonStartMomentum.mag());
       positronStartR->Fill(std::sqrt(std::pow(positronStartPosition.x(),2) +
                                      std::pow(positronStartPosition.y(),2))/
                                                                         10.);
       positronStartZ->Fill(positronStartPosition.z()/10.);
       positronStartPR->Fill(std::sqrt(std::pow(positronStartMomentum.x(),2) +
                                       std::pow(positronStartMomentum.y(),2)));
       positronStartPZ->Fill(positronStartMomentum.z());
       positronStartPRP->Fill(std::sqrt(std::pow(positronStartMomentum.x(),2) +
                                        std::pow(positronStartMomentum.y(),2))/
                                                  positronStartMomentum.mag());
       positronStartPZP->Fill(positronStartMomentum.z()/
                              positronStartMomentum.mag());

       G4double t3Scaling = (66.0 - positronStartPosition.z())/
                    (positronStartMomentum.z()/positronStartMomentum.mag());
       G4double t3ExpectedX = t3Scaling*positronStartMomentum.x()/
                    positronStartMomentum.mag() + positronStartPosition.x();
       G4double t3ExpectedY = t3Scaling*positronStartMomentum.y()/
                    positronStartMomentum.mag() + positronStartPosition.y();
       G4double t3ExpectedRadius = std::sqrt(std::pow(t3ExpectedX,2) + 
                                          std::pow(t3ExpectedY,2));

       positronT3ExpectedR->Fill(t3ExpectedRadius/10.);
       positronT3DifferenceR->Fill(std::abs(t3ExpectedRadius -
                                std::sqrt(std::pow(t3Position.x(),2) + 
                                          std::pow(t3Position.y(),2)))/10.);
       positronT3DifferencePZP->Fill(std::abs((positronStartMomentum.z()/
         positronStartMomentum.mag()) - (t3Momentum.z()/t3Momentum.mag())));

    }
  pionStopPosition.setX(0);   // DvB: changed to match version 2.1.3.1 of CLHEP
  pionStopPosition.setY(0);
  pionStopPosition.setZ(0);
}

void RunAction::FillSemiTriggeredHistos() {

    // All T1 Histograms require hits in B1, B2, Target, and T1 to fill
    positronT1R->Fill(std::sqrt(std::pow(t1Position.x(),2) + 
                                std::pow(t1Position.y(),2))/10.);
    positronT1Z->Fill(t1Position.z()/10.);
    positronT1PR->Fill(std::sqrt( std::pow(t1Momentum.x(),2) + 
                                  std::pow(t1Momentum.y(),2)));
    positronT1PZ->Fill(t1Momentum.z());
    positronT1PRP->Fill(std::sqrt( std::pow(t1Momentum.x(),2) + 
                                   std::pow(t1Momentum.y(),2))/
                                                t1Momentum.mag());
    positronT1PZP->Fill(t1Momentum.z()/t1Momentum.mag());

    G4double t1Scaling = (t1Position.z() - positronStartPosition.z())/
                      (positronStartMomentum.z()/positronStartMomentum.mag());
    G4double t1ExpectedX = t1Scaling*positronStartMomentum.x()/
                      positronStartMomentum.mag() + positronStartPosition.x();
    G4double t1ExpectedY = t1Scaling*positronStartMomentum.y()/
                      positronStartMomentum.mag() + positronStartPosition.y();
    G4double t1ExpectedRadius = std::sqrt(std::pow(t1ExpectedX,2) + 
                                          std::pow(t1ExpectedY,2));

    positronT1ExpectedR->Fill(t1ExpectedRadius/10.);
    positronT1DifferenceR->Fill(std::abs(t1ExpectedRadius - 
                                std::sqrt(std::pow(t1Position.x(),2) + 
                                          std::pow(t1Position.y(),2)))/10.);
    positronT1DifferencePZP->Fill(std::abs((positronStartMomentum.z()/
                                      positronStartMomentum.mag()) - 
                                      (t1Momentum.z()/t1Momentum.mag()) ) );

    positronT1dEdX->Fill(t1energyDeposit/t1pathLength);
    positronT1dEdXvPathLength->Fill(t1pathLength/3.0, t1energyDeposit);
}

void RunAction::FillTriggeredHistos() {

    triggeredPositronStartZ->Fill(positronStartPosition.z()/10.);
    triggeredPositronStartPZ->Fill(positronStartMomentum.z());
    triggeredPositronStartPZP->Fill(positronStartMomentum.z()/
                                    positronStartMomentum.mag());

    // All T3 Histograms require hits in B1, B2, Target, T1, and T3 to fill
    positronT3R->Fill(std::sqrt(std::pow(t3Position.x(),2) + 
                                std::pow(t3Position.y(),2))/10.);
    positronT3Z->Fill(t3Position.z()/10.);
    positronT3PR->Fill(std::sqrt(std::pow(t3Momentum.x(),2) + 
                                 std::pow(t3Momentum.y(),2)));
    positronT3PZ->Fill(t3Momentum.z());
    positronT3PRP->Fill(std::sqrt(std::pow(t3Momentum.x(),2) + 
                        std::pow(t3Momentum.y(),2))/t3Momentum.mag());
    positronT3PZP->Fill(pionStartPosition.z());

    G4double t3Scaling = (t3Position.z() - positronStartPosition.z())/
                     (positronStartMomentum.z()/positronStartMomentum.mag());
    G4double t3ExpectedX = t3Scaling*positronStartMomentum.x()/
                     positronStartMomentum.mag() + positronStartPosition.x();
    G4double t3ExpectedY = t3Scaling*positronStartMomentum.y()/
                     positronStartMomentum.mag() + positronStartPosition.y();
    G4double t3ExpectedRadius = std::sqrt(std::pow(t3ExpectedX,2) + 
                                          std::pow(t3ExpectedY,2));

    triggeredPositronT3ExpectedR->Fill(t3ExpectedRadius/10.);
    triggeredPositronT3DifferenceR->Fill(std::abs(t3ExpectedRadius - 
                                  std::sqrt(std::pow(t3Position.x(),2) +
                                            std::pow(t3Position.y(),2)))/10.);
    triggeredPositronT3DifferencePZP->Fill(std::abs((positronStartMomentum.z()/
          positronStartMomentum.mag()) - (t3Momentum.z()/t3Momentum.mag())));

    positronT3dEdX->Fill(t3energyDeposit/t3pathLength);
    positronT3dEdXvPathLength->Fill(t3pathLength/3.0, t3energyDeposit);
}


// Called by Geant4 at the beginning of a run
//
void RunAction::BeginOfRunAction(const G4Run* aRun) {

    G4cout << "### Run " << aRun->GetRunID() << " start ###" << G4endl;

    // Open Root
    OpenRoot();

    // save random number status
    //    G4RunManager::GetRunManager()->SetRandomNumberStore(true);
    G4RunManager::GetRunManager()->SetRandomNumberStoreDir("random/");

    if (autoSeed) {
       // automatic (time-based) random seeds for each run
       G4cout << "*******************" << G4endl;
       G4cout << "*** AUTOSEED ON ***" << G4endl;
       G4cout << "*******************" << G4endl;
       long seeds[2];
       time_t systime = time(NULL);
       seeds[0] = (long) systime;
       seeds[1] = (long) (systime*G4UniformRand());
       CLHEP::HepRandom::setTheSeeds(seeds);
    }

}

// Called by Geant4 at the end of a run
//
void RunAction::EndOfRunAction(const G4Run* aRun)
{

    G4cout << "### Run " << aRun->GetRunID() << " end ###" << G4endl;

    // Close ROOT
    CloseRoot();

    // save Rndm status
    CLHEP::HepRandom::saveEngineStatus("random/EndOfRun.rndm");

}

// Convert an integer to a character
//
char* RunAction::itoa(int val, int base) {

    static char buf[32] = {0};

    int i = 30;

    for(; val && i; --i, val /= base) {

        buf[i] = "0123456789abcdef"[val % base];

    }

    return &buf[i+1];

}
