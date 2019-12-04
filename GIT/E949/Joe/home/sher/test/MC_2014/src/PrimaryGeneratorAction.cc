#include "PrimaryGeneratorAction.hh"
#include "ParticleGunMessenger.hh"

#include "Randomize.hh"

#include "G4Event.hh"
#include "G4PionPlus.hh"
#include "G4MuonPlus.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4Gamma.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "TRandom3.h"


#define min(X,Y) ((X) < (Y) ? (X) : (Y))

// Primary events are generated from a tab delineated file of the
// following format:
// X position (cm)
// theta (mr)
// Y position (cm)
// phi (mr)
// Momentum magnitude (GeV/c)
//
// The values are given for particles at the location of the target.
// For this simulation, they are propagated backward 1.5m.

// Default constructor
//
PrimaryGeneratorAction::PrimaryGeneratorAction()
{

  //Beam Histograms
  hE  = new TH1F("hE","",200,0.07,0.074);
  hx  = new TH1F("hx","",200,-3,3);
  hy  = new TH1F("hy","",200,-3,3);
  htx = new TH1F("htx","",200,-0.2,0.2);
  hty = new TH1F("hty","",200,-0.2,0.2);
  
  //Messenger for the Particle Gun
  fMessenger = new ParticleGunMessenger(this);
  
  // Create a particle gun to shoot 1 particle per event
  G4int n_particle = 1;
  particleGunPi      = new G4ParticleGun(n_particle);
  particleGunMu      = new G4ParticleGun(n_particle);
  particleGunE       = new G4ParticleGun(n_particle);
  particleGunEMin    = new G4ParticleGun(n_particle);
  particleGunMuMin   = new G4ParticleGun(n_particle);
  particleGunPiMin   = new G4ParticleGun(n_particle);
  particleGunGamma   = new G4ParticleGun(n_particle);
  particleGunENuBar  = new G4ParticleGun(n_particle);
  particleGunMuNu    = new G4ParticleGun(n_particle);
  particleGunNeutron = new G4ParticleGun(n_particle);
  particleGunProton  = new G4ParticleGun(n_particle);

   // Find the definition of the particles
  G4String particleName1  = "pi+"; // "pi+" for pions
  G4String particleName2  = "mu+";
  G4String particleName3  = "e+";
  G4String particleName4  = "mu-";
  G4String particleName5  = "pi-";
  G4String particleName6  = "e-";
  G4String particleName7  = "gamma";
  G4String particleName8  = "anti_nu_e";
  G4String particleName9  = "nu_mu";
  G4String particleName10 = "neutron";
  G4String particleName11 = "proton";
  
    
  // Get the Geant4 particle table
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  particleGunPi->SetParticleDefinition(particleTable->FindParticle(particleName1));
  particleGunMu->SetParticleDefinition(particleTable->FindParticle(particleName2));
  particleGunE->SetParticleDefinition(particleTable->FindParticle(particleName3));
  particleGunMuMin->SetParticleDefinition(particleTable->FindParticle(particleName4));
  particleGunPiMin->SetParticleDefinition(particleTable->FindParticle(particleName5));
  particleGunEMin->SetParticleDefinition(particleTable->FindParticle(particleName6));
  particleGunGamma->SetParticleDefinition(particleTable->FindParticle(particleName7));
  particleGunENuBar->SetParticleDefinition(particleTable->FindParticle(particleName8));
  particleGunMuNu->SetParticleDefinition(particleTable->FindParticle(particleName9));
  particleGunNeutron->SetParticleDefinition(particleTable->FindParticle(particleName10));
  particleGunProton->SetParticleDefinition(particleTable->FindParticle(particleName11));

  const char filename1[] = "/home/pienumgr/beam/alex_26sep.rays";  // original pienu code
  const char filename2[] = "/home/pienumgr/beam/BeamData_G4_AS.rays";  // change from old BeamData_G4.rays
  //const char filename2[] = ""; 

  //const char filename1[] = "/home/dorothea/MC/rays_at_F4_72MeV.txt"; // input from g4bl momentum and x-y distribution 
  //const char filename1[] = "/home/dorothea/MC/Gamma_rays.txt";      // change input here for gamma rays
  //const char filename1[] = "/home/dorothea/MC/Gamma_rays_F1.txt";   // inputs for 72 MeV/c slices
  //const char filename1[] = "/home/dorothea/MC/Gamma_rays_F2.txt"; 
  //const char filename1[] = "/home/dorothea/MC/Gamma_rays_F3.txt"; 
  //const char filename1[] = "/home/dorothea/MC/Gamma_rays_F4.txt"; 
  //const char filename1[] = "/home/dorothea/MC/Gamma_rays_F5.txt"; 

  //const char filename1[] = "/home/dorothea/MC/Gamma_rays_74MeV_F1.txt";   // inputs for 74 MeV/c slices
  //const char filename1[] = "/home/dorothea/MC/Gamma_rays_74MeV_F2.txt";
  //const char filename1[] = "/home/dorothea/MC/Gamma_rays_74MeV_F3.txt";
  //const char filename1[] = "/home/dorothea/MC/Gamma_rays_74MeV_F4.txt";
  //const char filename1[] = "/home/dorothea/MC/Gamma_rays_74MeV_F5.txt";

  infile1 = new std::ifstream(filename1);
  infile2 = new std::ifstream(filename2);

  if(!(*infile1)) {
    //G4cerr << "Error opening input file " << filename1 << G4endl;
     G4Exception("Error opening input file ",filename1,FatalErrorInArgument,"");
  }

  if(!(*infile2)) {
    //G4cerr << "Error opening input file " << filename2 << G4endl;
    G4Exception("Error opening input file ",filename2,FatalErrorInArgument,"");
  }

 
  // An extra check - probably not needed
  std::filebuf* inbuff = infile1->rdbuf();
  if(!inbuff->is_open()) {
    //G4cerr << "Input file could not be opened " << filename1 << G4endl;
    G4Exception("Input file could not be opened ",filename1,FatalErrorInArgument,"");
  }

  //Fill the histogram to sample
  //FillBeamHistograms();   // original pienu code
  //Definition of the Choleski Matrix

  L[0][0] = 1.0;  L[1][0] = 0.501133; L[2][0] = -0.130539;  L[3][0] = 0.0794952;  L[4][0] = -0.128491;
  L[0][1] = 0.0;  L[1][1] = 0.86537;  L[2][1] = 0.259344;   L[3][1] = 0.387393;   L[4][1] = 0.102995;
  L[0][2] = 0.0;  L[1][2] = 0.0;      L[2][2] = 0.956922;   L[3][2] = -0.0534969; L[4][2] = 0.594524;
  L[0][3] = 0.0;  L[1][3] = 0.0;      L[2][3] = 0.0;        L[3][3] = 0.916922;   L[4][3] = 0.0307296;
  L[0][4] = 0.0;  L[1][4] = 0.0;      L[2][4] = 0.0;        L[3][4] = 0.0;        L[4][4] = 0.786434;

  //Definition of the Choleski Matrix (for BeamData_G4.rays) - older file, created initially by Luca
/* 
 L[0][0] = 1.0;  L[1][0] = 0.339225;  L[2][0] = -0.0928713;  L[3][0] = 0.058468;   L[4][0] = -0.0904278;
  L[0][1] = 0.0;  L[1][1] = 0.940705;  L[2][1] = 0.191632;    L[3][1] = 0.384308;   L[4][1] = 0.0420166;
  L[0][2] = 0.0;  L[1][2] = 0.0;       L[2][2] = 0.977063;    L[3][2] = -0.0311674; L[4][2] = 0.620841;
  L[0][3] = 0.0;  L[1][3] = 0.0;       L[2][3] = 0.0;         L[3][3] = 0.920824;   L[4][3] = 0.0283415;
  L[0][4] = 0.0;  L[1][4] = 0.0;       L[2][4] = 0.0;         L[3][4] = 0.0;        L[4][4] = 0.777052;
*/

}
 
// Default destructor
// 
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    
  delete particleGunPi;
  delete particleGunMu;
  delete particleGunE;
  delete particleGunMuMin;
  delete particleGunPiMin;
  delete particleGunEMin;
  delete particleGunGamma;
  delete particleGunENuBar;
  delete particleGunMuNu;
  delete particleGunNeutron;
  delete particleGunProton;

  infile1->close();
  infile2->close();
}

// This method is called at the beginning of each event in order
// to generate a new primary vertex for tracking
//
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  // Variables to store the values read in from the file
  //G4double X,dX,Y,dY,P;
  G4double X,Y,P,Z,Px,Py,Pz,dX,dY; int pID;

  // Variables to store the final vertex position
  G4double xv,yv,zv,mdX,mdY,mdZ,kinEnergyPi,kinEnergyMu,kinEnergyE,kinEnergyGamma,kinEnergyENuBar,kinEnergyMuNu,kinEnergyNeutron,kinEnergyProton;

  // If we're at the end of the file, rewind.
  if (infile1->eof()) {
    infile1->clear();
    infile1->seekg(0);
  }

  if (infile2->eof()) {
    infile2->clear();
    infile2->seekg(0);
  }

  // Read in a line from the file
  if (!databeam){
    
       // Original Pienu code
       
    //Revmoc file
    *infile1 >> X >> dX >> Y >> dY >> P;
    // Convert the angles to radians and find the x and y length of a
    // vector if the total length is 1.
    mdX = std::sin(dX/1000.0);
    mdY = std::sin(dY/1000.0);
    mdZ = std::sqrt(1-mdX*mdX-mdY*mdY);

   //  //=================================================
//     // Saul's code to read file from G4Beamline
//     //=================================================

//     *infile1 >> X >> Y >> Z >> Px >> Py >> Pz >> pID;
    
//     // vector if the total length is 1.

//     // need length conversion only for input from G4Beamline, not for gamma ray input
//    //  //from mm to cm
//     X = X / 10.0 *cm;
//     Y = Y / 10.0 *cm;
//     Z = Z / 10.0 *cm;    

//     // from MeV/c to GeV/c
//     Px = Px / 1000.0; // * GeV;
//     Py = Py / 1000.0; // * GeV;
//     Pz = Pz / 1000.0; // * GeV;

//     // G4cout << "Px " << Px << " Py " << Py << " Pz " << Pz << G4endl;

//     P = sqrt(Px*Px + Py*Py + Pz*Pz);
    
//     if ( P > 0 )
//       {
// 	mdX = Px/P;
// 	mdY = Py/P;
// 	mdZ = Pz/P;
//     }
    
    //end of G4Beamline input change
    
  } else {
    //Data parameterization file
 
    //Sample the histograms and construct the correlated variables
    DoCorrelatedVariables();

    //Conversion not needed for the data beam
    P = r[0];    // 0.072 = 72 MeV
    //P = 0.072 * GeV;
    //G4cout << r[0] << G4endl;
    
    X = r[1];
    Y = r[2];
    mdX = r[3];
    mdY = r[4];
    if (mdX*mdX+mdY*mdY<=1) mdZ = std::sqrt(1-mdX*mdX-mdY*mdY);
    else mdZ = 1.0;
  }

  // Scale and shift the vector.  Everything is in centimeters.

  // Original Pienu code
  xv = -60*cm*(mdX/mdZ) + X*cm;
  yv = -60*cm*(mdY/mdZ) + Y*cm;
  zv = -60*cm;
  
  //Code when using input from G4Beamline
//    // Saul's z-shift:
//   //float shift_z = -12.5261*cm; //10*mm; // cm    ----  WC1_1/z = -125.261 (mm)  // for input from G4Beamline data
//   //float shift_z = -60*cm;
//    xv = X; // xv = shift_z*cm*(mdX/mdZ) + X*cm;
//    yv = Y; // yv = shift_z*cm*(mdY/mdZ) + Y*cm;
//    //zv = shift_z;  // default for g4bl input

//    zv = Z; // for gamma ray data, z position correct as input

   //G4cout << xv << " " << yv << " " << " " << zv << " " <<  mdX << " " << mdY << " " << mdZ << G4endl;

   // Momentum scaling factor to ensure Pion stops in the middle of the target
   G4double momentumFactor = 1.03828;//
   //G4double momentumFactor = 0.86; //For Muons stopped in target

   G4double NeutronMass = G4Neutron::NeutronDefinition()->GetPDGMass();
   G4double ProtonMass  = G4Proton::ProtonDefinition()->GetPDGMass();
   G4double PionMass    = G4PionPlus::PionPlusDefinition()->GetPDGMass();
   G4double EMass       = G4Positron::PositronDefinition()->GetPDGMass();
   G4double MuMass      = G4MuonPlus::MuonPlusDefinition()->GetPDGMass();
   G4double GammaMass   = G4Gamma::GammaDefinition()->GetPDGMass();

   //========================================
   // change particle mometum here
   //========================================

   //G4double MuonFactor = .9645;   // .9645 for 72 MeV/c, .99 for 74 MeV/c

   //G4double MuonFactor = .95;

   //G4double MuonFactor = .963;  // 72 MeV/c
   //G4double MuonFactor = 1.;

   G4double MuonFactor = .99;  // 74 MeV/c
   //G4double MuonFactor = 1.004; // 75 MeV/c

   G4double ElectronFactor74 = 1.02273;
   G4double ElectronFactor72 = 1.01094;

   G4double Momentum =   momentumFactor * MuonFactor * P * GeV;   // for muon simulation without g4bl input

   //G4double Momentum = momentumFactor * P * GeV;  // original code

   //G4double Momentum =  P * MeV;

   // calculate kinetic energy of particle
   kinEnergyPi      = std::sqrt(Momentum*Momentum + PionMass*PionMass) - PionMass;
   kinEnergyMu      = std::sqrt(Momentum*Momentum + MuMass*MuMass) - MuMass;
   kinEnergyE       = std::sqrt(Momentum*Momentum + EMass*EMass) - EMass;
   kinEnergyGamma   = std::sqrt(Momentum*Momentum + GammaMass*GammaMass) - GammaMass;
   kinEnergyENuBar  = Momentum;
   kinEnergyMuNu    = Momentum;
   kinEnergyNeutron = std::sqrt(Momentum*Momentum + NeutronMass*NeutronMass) - NeutronMass;
   kinEnergyProton  = std::sqrt(Momentum*Momentum + ProtonMass*ProtonMass) - ProtonMass;

   //cout << kinEnergyProton << endl;

   // Set the particle position, energy, and momentum direction
   particleGunPi->SetParticlePosition(G4ThreeVector(xv,yv,zv));
   particleGunMu->SetParticlePosition(G4ThreeVector(xv,yv,zv));
   particleGunE->SetParticlePosition(G4ThreeVector(xv,yv,zv));
   particleGunMuMin->SetParticlePosition(G4ThreeVector(xv,yv,zv));
   particleGunPiMin->SetParticlePosition(G4ThreeVector(xv,yv,zv));
   particleGunEMin->SetParticlePosition(G4ThreeVector(xv,yv,zv));
   particleGunGamma->SetParticlePosition(G4ThreeVector(xv,yv,zv));
   particleGunENuBar->SetParticlePosition(G4ThreeVector(xv,yv,zv));
   particleGunMuNu->SetParticlePosition(G4ThreeVector(xv,yv,zv));
   particleGunNeutron->SetParticlePosition(G4ThreeVector(xv,yv,zv));
   particleGunProton->SetParticlePosition(G4ThreeVector(xv,yv,zv));

   particleGunPi->SetParticleMomentumDirection(G4ThreeVector(0,0,1)); // pions in z-direction
   //particleGunPi->SetParticleMomentumDirection(G4ThreeVector(mdX,mdY,mdZ));
   particleGunMu->SetParticleMomentumDirection(G4ThreeVector(mdX,mdY,mdZ));
   particleGunE->SetParticleMomentumDirection(G4ThreeVector(mdX,mdY,mdZ));
   particleGunMuMin->SetParticleMomentumDirection(G4ThreeVector(mdX,mdY,mdZ));
   particleGunPiMin->SetParticleMomentumDirection(G4ThreeVector(mdX,mdY,mdZ));
   particleGunEMin->SetParticleMomentumDirection(G4ThreeVector(mdX,mdY,mdZ));
   particleGunGamma->SetParticleMomentumDirection(G4ThreeVector(mdX,mdY,mdZ));
   particleGunENuBar->SetParticleMomentumDirection(G4ThreeVector(mdX,mdY,mdZ));
   particleGunMuNu->SetParticleMomentumDirection(G4ThreeVector(mdX,mdY,mdZ));
   particleGunNeutron->SetParticleMomentumDirection(G4ThreeVector(mdX,mdY,mdZ));
   particleGunProton->SetParticleMomentumDirection(G4ThreeVector(mdX,mdY,mdZ));

   particleGunPi->SetParticleEnergy(kinEnergyPi);
   particleGunMu->SetParticleEnergy(kinEnergyMu);
   particleGunE->SetParticleEnergy(kinEnergyE);
   particleGunMuMin->SetParticleEnergy(kinEnergyMu);
   particleGunPiMin->SetParticleEnergy(kinEnergyPi);
   particleGunEMin->SetParticleEnergy(kinEnergyE);
   particleGunGamma->SetParticleEnergy(kinEnergyGamma);
   particleGunENuBar->SetParticleEnergy(kinEnergyENuBar);
   particleGunMuNu->SetParticleEnergy(kinEnergyMuNu);
   particleGunNeutron->SetParticleEnergy(kinEnergyNeutron);
   particleGunProton->SetParticleEnergy(kinEnergyProton);

   //POLARIZATION
   //particleGunMu->SetParticlePolarization(G4ThreeVector(mdX,mdY,mdZ));

   //NO-POLARIZATION
   //double rdmpol = G4UniformRand();
   //if (rdmpol>0.5) particleGunMu->SetParticlePolarization(G4ThreeVector(-mdX,-mdY,-mdZ));
   //else particleGunMu->SetParticlePolarization(G4ThreeVector(mdX,mdY,mdZ));

   //ANTI-POLARIZATION
   //particleGunMu->SetParticlePolarization(G4ThreeVector(-mdX,-mdY,-mdZ));

   //TRANSVERSE-POLARIZATION
   //particleGunMu->SetParticlePolarization(G4ThreeVector(0,0,1));

   //Begin Pileup code 
   if (pileup){

     G4double Rate = rate /s; //Beam rate (Hertz)

     G4int N = 0;
     G4int MAX=5;
     G4double tdecay;
     G4double tau = 2197.03*ns;//26.033 * ns;
     G4double t0[5]; //time
     G4double td[5]; //decay time

     //Time window
     G4double tmin = Tmin*ns;
     G4double tmax = Tmax*ns;

     for (int i=0;i<5;i++){t0[i]=0;td[i]=0;}

     //Generate the trigger particle
     t0[0] = 0;
     td[0] = -tau*log(G4UniformRand());

     G4double t=0;

     //For now, removed from the original algorithm 
     /*
     while (t>=-ttrig) {
       ttrig = -(1/Rate)*log(G4UniformRand());
       t     =  (1/Rate)*log(G4UniformRand());
     }


     //Try to place the first pileup particle before the trigger (->?)
     t = (1/Rate)*log(G4UniformRand()); 
     tdecay = (t - tau*log(G4UniformRand()))*ns;
     if ((t>=tmin)||(tdecay>=tmin && tdecay<=tmax)){    
       N=1;
       t0[N] = t*ns;
       td[N] = (-tau*log(G4UniformRand()))*ns;
     }
     */

     //Construct the other pileup particles before the trigger
     while ((t>=tmin)||(tdecay>=tmin && tdecay<=tmax)) {
       t = t - (-(1/Rate)*log(G4UniformRand()));
       tdecay = t - tau*log(G4UniformRand());

       if ((t>=tmin)||(tdecay>=tmin && tdecay<=tmax)){
	N++;
	 N = min(N,MAX);
	 t0[N] = t*ns;
	 td[N] = (-tau*log(G4UniformRand()))*ns;
       }
     }

     //Construct the other pileup particles after the trigger
     t = 0.0 * ns;
     while (t<=tmax) {
       t = (t + (-(1/Rate)*log(G4UniformRand())))*ns;

       if (t<tmax){
	 N++;
	 N = min(N,MAX);      
	 t0[N] = t*ns;
	 td[N] = (-tau*log(G4UniformRand()))*ns;
       }
     }

     //G4cout << "---" << G4endl;
     //for (int i=0;i<5;i++) G4cout << "N = " << i << " " << t0[i] << "  " << td[i] << G4endl;

     // for beam of mostly pions, some muons and electrons

     G4double part;
    //  for (int i=0;i<=N;i++) {
 //       part = G4UniformRand();
 //       if (part<0.85)                 {particleGunPi->SetParticleTime(t0[i]);  particleGunPi->GeneratePrimaryVertex(anEvent);} //pion+
 //       if (part>=0.85 && part<0.99)   {particleGunMu->SetParticleTime(t0[i]);  particleGunMu->GeneratePrimaryVertex(anEvent);} //muon+
 //       if (part>=0.99)                {particleGunE->SetParticleTime(t0[i]);   particleGunE->GeneratePrimaryVertex(anEvent);}  //positron
 //     }

     //====================================================
     // CHANGE PARTICLE TYPE HERE (twice!!)
     //====================================================
     // Original Pienu code
     for (int i=0;i<=N;i++) {
       part = G4UniformRand();
       //particleGunMu->SetParticleTime(t0[i]);  particleGunMu->GeneratePrimaryVertex(anEvent); //muon+
       particleGunMuMin->SetParticleTime(t0[i]); particleGunMuMin->GeneratePrimaryVertex(anEvent); //muon-
       //particleGunE->SetParticleTime(t0[i]);   particleGunE->GeneratePrimaryVertex(anEvent); //positrons
       //particleGunEMin->SetParticleTime(t0[i]);   particleGunEMin->GeneratePrimaryVertex(anEvent); //electrons
       //particleGunPiMin->SetParticleTime(t0[i]);  particleGunPiMin->GeneratePrimaryVertex(anEvent); //pi-
     }

  //    // Original Pienu code
 //  //  } else{  // no pileup

 //     //particleGunMu->GeneratePrimaryVertex(anEvent);    // muon+
 //     particleGunMuMin->GeneratePrimaryVertex(anEvent); // muon-
 //     //particleGunE->GeneratePrimaryVertex(anEvent);       // positrons
 //     //particleGunEMin->GeneratePrimaryVertex(anEvent);       // electrons
 //     //particleGunPiMin->GeneratePrimaryVertex(anEvent);   //pi-

 // //   }

     //   } else { // No pileup

     //     particleGunPi->GeneratePrimaryVertex(anEvent); //pion+
     //     //G4cout << particleGunMu->GetParticlePolarization() << G4endl;
     //     // particleGunMu->GeneratePrimaryVertex(anEvent); 

     //   }

     // for beam of only muons (positive or negative)


   } else { // No pileup

     // particleGunPi->GeneratePrimaryVertex(anEvent); //pion+
     
     //G4cout << "Particle Type " << pID << G4endl;

     //pID = 13;   // when not using G4Beamline, but simulating mu^- for muon capture experiment
     //pID = -13;  // for mu^+

     pID = 11;

     if (pID == 211)
       particleGunPi->GeneratePrimaryVertex(anEvent);    // positive pion

     if ( pID == -11 )
      particleGunE->GeneratePrimaryVertex(anEvent);      // positron
    
    if ( pID == 11 )
      particleGunEMin->GeneratePrimaryVertex(anEvent);   // electron
    
    if ( pID == 22 )
      particleGunGamma->GeneratePrimaryVertex(anEvent);  // gamma
    
    if ( pID == 13)
      particleGunMuMin->GeneratePrimaryVertex(anEvent);  // muon
    
    if ( pID == -13)
      particleGunMu->GeneratePrimaryVertex(anEvent);  //anti  muon
 
    if ( pID == -12)
      particleGunENuBar->GeneratePrimaryVertex(anEvent); // anti electon neutrino

    if ( pID == 14)
     particleGunMuNu->GeneratePrimaryVertex(anEvent);    // muon neutrino
    
    if ( pID == 2112)
      particleGunNeutron->GeneratePrimaryVertex(anEvent); // neutron
        if ( pID == 2212) 
      particleGunProton->GeneratePrimaryVertex(anEvent);  // proton
    

   }
    
 
}

void PrimaryGeneratorAction::FillBeamHistograms()
{
 
  G4double P, X, Y, dX, dY;
 
  while (1)
  {
    *infile2 >> X >> dX >> Y >> dY >> P;

    if (!infile2->good()) break;

    hE->Fill(P);
    hx->Fill(X);
    hy->Fill(Y);
    htx->Fill(dX);
    hty->Fill(dY); 

  }

  //Smooth the X and Y distributions
  hx->Smooth(200);
  hy->Smooth(200);

  //Fit the tx and ty distributions
  htx->Fit("gaus");
  hty->Fit("gaus");
  
  
  tx_mean=htx->GetFunction("gaus")->GetParameter(1);
  tx_sigma=htx->GetFunction("gaus")->GetParameter(2);
  
  ty_mean=hty->GetFunction("gaus")->GetParameter(1);
  ty_sigma=hty->GetFunction("gaus")->GetParameter(2);
  
}

void PrimaryGeneratorAction::DoCorrelatedVariables()
{

  TRandom3  *r3 = new TRandom3(0);
  gRandom = r3;

  //Sample the histograms
  r[0] = hE->GetRandom();
  r[1] = hx->GetRandom();
  r[2] = hy->GetRandom();

  //Generate tx and ty disributions using gaussian
  r[3] = r3->Gaus(tx_mean,tx_sigma);
  r[4] = r3->Gaus(ty_mean,ty_sigma);

  delete r3;

  //Normalize the variables
  r[0] = (r[0] - hE->GetMean()) / hE->GetRMS();
  r[1] = (r[1] - hx->GetMean()) / hx->GetRMS();
  r[2] = (r[2] - hy->GetMean()) / hy->GetRMS();
  r[3] = (r[3] - htx->GetMean()) / htx->GetRMS();
  r[4] = (r[4] - hty->GetMean()) / hty->GetRMS();

  //Apply the Choleski Matrix
  G4double rt[5];
  for (int i = 0; i < 5; i++)
  {
    rt[i] = 0;
    for (int j = 0; j < 5; j++)
    {
      rt[i] += L[i][j] * r[j];
    } 
  }

  //Transform back the (correlated) variables
  r[0] = rt[0] * hE->GetRMS()  + hE->GetMean();
  r[1] = rt[1] * hx->GetRMS()  + hx->GetMean();
  r[2] = rt[2] * hy->GetRMS()  + hy->GetMean();
  r[3] = rt[3] * htx->GetRMS() + htx->GetMean();
  r[4] = rt[4] * hty->GetRMS() + hty->GetMean();

}
