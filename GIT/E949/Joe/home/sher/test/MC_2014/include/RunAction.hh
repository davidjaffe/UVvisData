#ifndef RunAction_h
#define RunAction_h 1

#include "G4ThreeVector.hh"

#include "G4UserRunAction.hh"

class TFile;
class TNtuple;
class TH1D;
class TH2D;
class TTree;

class G4Run;
class RunActionMessenger;

class RunAction : public G4UserRunAction {

  public:

    RunAction();
    ~RunAction();

    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);

    void SetAutoSeed (const G4bool val) { autoSeed = val; }

    void OpenRoot();
    void CloseRoot();

    void FillTuple(G4double E1,
                   G4double E2,
                   G4double E3,
                   G4double E4,
                   G4double E5,
                   G4double E6,
                   G4double E7,
                   G4double E8,
                   G4double E9,
                   G4double E10,
                   G4double E11,
                   G4double E12,
                   G4double E13,
                   G4double E14,
                   G4double E15,
                   G4double E16);
    void FillTree();

    void ClearVariable();

    void SetEventID(G4int eID) { eventID=eID; }

    void SPionStart(G4ThreeVector, G4double, G4ThreeVector, G4double);
    void SPionDecay(G4ThreeVector, G4double, G4ThreeVector, G4double);
    void SMuonStart(G4ThreeVector, G4double, G4ThreeVector, G4double);
    void SMuonDecay(G4ThreeVector, G4double, G4ThreeVector, G4double);
    void SMuonStop(G4ThreeVector, G4double, G4ThreeVector, G4double);
    void SPositronStart(G4ThreeVector, G4double, G4ThreeVector, G4double);
    void SPositronStop(G4ThreeVector, G4double, G4ThreeVector, G4double);
    void SPhotonuclear(G4ThreeVector);
    void SMuETG(G4double);
    void SBhabha(G4double);
    void SBhabhaEl(G4double);
    void SAnnihilT1(G4int);
    void SAnnihilT2(G4int);
    void SAnnihilTg(G4int);
    void SInelasticB1(G4int);
    void SInelasticB1gamma(G4double);
    void SInelasticB2(G4int);
    void SInelasticB2gamma(G4double);
    void SInelasticTg(G4int);
    void SInelasticTggamma(G4double);
    void SMscT2Angle(G4ThreeVector);
    void SMscB1Angle(G4ThreeVector);
    void SMuPolarization(G4ThreeVector);
    void SEmomentum (G4ThreeVector);
    void SMuCap(G4int);
    void SRadCap(G4ThreeVector, G4double, G4ThreeVector, G4double);
    void SGammaRadCap(G4ThreeVector, G4double, G4ThreeVector, G4double);
    void SGammaSignal(G4ThreeVector, G4double, G4ThreeVector, G4double);
    void SMuCapPID(G4int);
    void SMuStopZr(G4int);
    void SMuStopT1(G4int);
    void SMuStopS3(G4int);
    void SMuStopAfterZr(G4int);
    void SZrPreP (G4double);
    void SZrPostP (G4double);
    void SZrFirstSlice(G4ThreeVector, G4ThreeVector, G4double);
    void SZrSecondSlice(G4ThreeVector, G4ThreeVector, G4double);
    void SZrThirdSlice(G4ThreeVector, G4ThreeVector, G4double);
    void SZrFifthSlice(G4ThreeVector, G4ThreeVector, G4double);
    void SZrFourthSlice(G4ThreeVector, G4ThreeVector, G4double);
    void SMuonCapture(G4ThreeVector, G4double, G4ThreeVector, G4double);
    void SBeamMuonBina(G4ThreeVector, G4double, G4ThreeVector, G4double);
    void SMuonCaptureZr(G4ThreeVector, G4double, G4ThreeVector, G4double);
    void SCaptureNeutron(G4ThreeVector, G4double, G4ThreeVector, G4double);
    void SCaptureProton(G4ThreeVector, G4double, G4ThreeVector, G4double);
    void SCaptureGamma(G4ThreeVector, G4double, G4ThreeVector, G4double);
    void SMuBrems(G4ThreeVector, G4double, G4ThreeVector, G4double);
    void SGammaStart(G4ThreeVector, G4double, G4ThreeVector, G4double);
    void SGammaStartZr(G4ThreeVector, G4double, G4ThreeVector, G4double);
    void SGammaBina(G4ThreeVector, G4double, G4ThreeVector, G4double);
    void SElectronBina(G4ThreeVector, G4double, G4ThreeVector, G4double);
    void SElectronStart(G4ThreeVector, G4double, G4ThreeVector, G4double);
    void SElectronStop(G4ThreeVector, G4double, G4ThreeVector, G4double);
    
    void SavePionStart(G4ThreeVector, G4ThreeVector);
    void SavePionStop(G4ThreeVector);
    void SaveMuonStart(G4ThreeVector, G4ThreeVector);
    void SavePositronStart(G4ThreeVector, G4ThreeVector);
    void SaveT1(G4ThreeVector, G4ThreeVector);
    void SaveT3(G4ThreeVector, G4ThreeVector);
    void SaveT1dEdX(G4double, G4double);
    void SaveT3dEdX(G4double, G4double);

    void FillUntriggeredHistos();
    void FillTargetStopHistos();
    void FillSemiTriggeredHistos();
    void FillTriggeredHistos();

  private:

    RunActionMessenger* runActionMessenger;

    G4bool autoSeed;

    TFile*   file;
    TNtuple* hitTuple;
    TTree*   aTree;

    //Diagnostic Histograms

    // Particle Start Position/Momentum and Stop Position
    TH1D* pionStartR; // Pion Start Position R Coordinate
    TH1D* pionStartZ; // Pion Start Position Z Coordinate
    TH1D* pionStartPR; // Pion Start Momentum R Component
    TH1D* pionStartPZ; // Pion Start Momentum Z Component
    TH1D* pionStartPRP; // Pion Start Momentum Pr/P
    TH1D* pionStartPZP; // Pion Start Momentum Pz/P
    TH1D* pionStopR; // Pion Stop Position R Coordinate
    TH1D* pionStopZ; // Pion Stop Position Z Coordinate
    TH1D* muonStartR; // Muon Start Position R Coordinate (Same As Pion Stop Coordinate - Redundant)
    TH1D* muonStartZ; // Muon Start Position Z Coordinate (Same As Pion Stop Coordinate - Redundant)
    TH1D* muonStartPR; // Muon Start Momentum R Component
    TH1D* muonStartPZ; // Muon Start Momentum Z Component
    TH1D* muonStartPRP; // Muon Start Momentum Pr/P
    TH1D* muonStartPZP; // Muon Start Momentum Pz/P
    TH1D* positronStartR; // Positron Start Position R Coordinate
    TH1D* positronStartZ; // Positron Start Position Z Coordinate
    TH1D* positronStartPR; // Positron Start Momentum R Component
    TH1D* positronStartPZ; // Positron Start Momentum Z Component
    TH1D* positronStartPRP; // Positron Start Momentum Pr/P
    TH1D* positronStartPZP; // Positron Start Momentum Pz/P
    TH1D* triggeredPositronStartZ; // Triggered Positron Start Position Z Coordinate
    TH1D* triggeredPositronStartPZ; // Triggered Positron Momentum Z Component
    TH1D* triggeredPositronStartPZP; // Triggered Positron Momentum Pz/P

    // All T1 Histograms require hits in B1, B2, Target, and T1 to fill
    TH1D* positronT1R; // Positron T1 Position R Coordinate
    TH1D* positronT1Z; // Positron T1 Position Z Coordinate
    TH1D* positronT1PR; // Positron T1 Momentum R Component
    TH1D* positronT1PZ; // Positron T1 Momentum Z Component
    TH1D* positronT1PRP; // Positron T1 Momentum Pr/P
    TH1D* positronT1PZP; // Positron T1 Momentum Pz/P
    TH1D* positronT1ExpectedR; //Expected Radius at center of T1 from positron start direction
    TH1D* positronT1DifferenceR; // Difference between expected and actual radial coordinates
    TH1D* positronT1DifferencePZP; // Difference between start PZ/P and PZ/P at T1 entrance
    TH1D* positronT1dEdX; // dE/dX for positrons in T1
    TH2D* positronT1dEdXvPathLength; // dE/dX vs. Path Length for positrons in T1

    // All T3 Histograms require hits in B1, B2, Target, T1, and T3 to fill
    TH1D* positronT3R; // Positron T3 Position R Coordinate
    TH1D* positronT3Z; // Positron T3 Position Z Coordinate
    TH1D* positronT3PR; // Positron T3 Momentum R Component
    TH1D* positronT3PZ; // Positron T3 Momentum Z Component
    TH1D* positronT3PRP; // Positron T3 Momentum Pr/P
    TH1D* positronT3PZP; // Positron T3 Momentum Pz/P
    TH1D* positronT3ExpectedR; // Expected Radius at T3 Entrance from positron start direction
    TH1D* positronT3DifferenceR; // Difference between expected and actual radial coordinates
    TH1D* positronT3DifferencePZP; // Difference between start PZ/P and PZ/P at T3 entrance
    TH1D* positronT3dEdX; // dE/dX for positrons in T3
    TH2D* positronT3dEdXvPathLength; // dE/dX vs. Path Length for positrons in T3
    TH1D* triggeredPositronT3ExpectedR; // Triggered Positron Expected Radius At T3 Entrance
    TH1D* triggeredPositronT3DifferenceR; // Difference between expected and actual radial coordinates for triggered events at T3
    TH1D* triggeredPositronT3DifferencePZP; // Difference between start PZ/P and PZ/P at T3 entrance

    G4int eventID;
    G4int AnhlT1;
    G4int AnhlT2;
    G4int AnhlTg;
    G4int InelB1;
  //G4int InelB1gamma;
    G4int InelB2;
  //G4int InelB2gamma;
    G4int InelTg;
  //G4int InelTggamma;
    G4int MuCap; 
    G4int MuCapPID;
    G4int MuStopZr;
    G4int MuStopT1;
    G4int MuStopS3;
    G4int MuStopAfterZr;

    G4double PStartX[4],PStartP[4];
    G4double PDecayX[4],PDecayP[4];
    G4double MStartX[4],MStartP[4];
    G4double MDecayX[4],MDecayP[4];
    G4double MStopX[4],MStopP[4];
    G4double EStartX[4],EStartP[4];
    G4double EStopX[4],EStopP[4];
    G4double RadCapX[5],RadCapP[4];
    G4double GammaRadCapX[5],GammaRadCapP[4];
    G4double EPhNX[3];
    G4double EBh;
    G4double EBhEl;
    G4double EMuTG;
    G4double InelB1gammaE;
    G4double InelB2gammaE;
    G4double InelTggammaE;
  //G4double MscT2Angle;
    G4double MuZrPreP;
    G4double MuZrPostP;
    G4double ZrFirstSliceP[4], ZrFirstSliceX[4];
    G4double ZrSecondSliceP[4], ZrSecondSliceX[4];
    G4double ZrThirdSliceP[4], ZrThirdSliceX[4];
    G4double ZrFourthSliceP[4], ZrFourthSliceX[4];
    G4double ZrFifthSliceP[4], ZrFifthSliceX[4];
    G4double MuCapX[4],MuCapP[4];
    G4double BeamMuonBinaX[4],BeamMuonBinaP[4];
    G4double MuCapZrX[4],MuCapZrP[4];
    G4double CapNeutronX[4],CapNeutronP[4];
    G4double CapProtonX[4],CapProtonP[4];
    G4double GammaSignalX[4], GammaSignalP[4];
    G4double CapGammaX[4],CapGammaP[4];		
    G4double MuBremsX[4], MuBremsP[4];
    G4double GammaStartX[4], GammaStartP[4];
    G4double GammaBinaX[4], GammaBinaP[4];
    G4double GammaStartZrX[4], GammaStartZrP[4];
    G4double ElectronBinaX[4], ElectronBinaP[4];
    G4double ElectronStartX[4], ElectronStartP[4];
    G4double ElectronStopX[4], ElectronStopP[4];

    G4ThreeVector MscT2Angle;
    G4ThreeVector MscB1Angle;

    G4double MuPol[3];
    G4double EMom[3];

    G4ThreeVector pionStartPosition;
    G4ThreeVector pionStartMomentum;
    G4ThreeVector pionStopPosition;
    G4ThreeVector muonStartPosition;
    G4ThreeVector muonStartMomentum;
    G4ThreeVector positronStartPosition;
    G4ThreeVector positronStartMomentum;
    G4ThreeVector t1Position;
    G4ThreeVector t1Momentum;
    G4ThreeVector t3Position;
    G4ThreeVector t3Momentum;
    G4double t1energyDeposit;
    G4double t1pathLength;
    G4double t3energyDeposit;
    G4double t3pathLength;

    char* itoa(int val, int base);

};

#endif
