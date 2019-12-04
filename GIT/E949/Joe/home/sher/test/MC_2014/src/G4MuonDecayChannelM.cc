
///////////////////////////////////
//
// Muon Decay into e+/nu/nu_heavy
//
// L.Doria Dec.2011
//
///////////////////////////////////


#include "G4ParticleDefinition.hh"
#include "G4DecayProducts.hh"
#include "G4VDecayChannel.hh"
#include "G4MuonDecayChannelM.hh"
#include "Randomize.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"
#include "G4RotationMatrix.hh"

G4MuonDecayChannelM::G4MuonDecayChannelM(const G4String& theParentName, 
				       G4double        theBR)
                   :G4VDecayChannel("Muon Decay",1)
{
  // set names for daughter particles
  if (theParentName == "mu+") {
    SetBR(theBR);
    SetParent("mu+");
    SetNumberOfDaughters(3);
    SetDaughter(0, "e+");
    SetDaughter(1, "nu_h");
    SetDaughter(2, "nu_e");//anti_nu_mu
  } else if (theParentName == "mu-") {
    SetBR(theBR);
    SetParent("mu-");
    SetNumberOfDaughters(3);
    SetDaughter(0, "e-");
    SetDaughter(1, "anti_nu_e");
    SetDaughter(2, "nu_mu");
  } else {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4MuonDecayChannel:: constructor :";
      G4cout << " parent particle is not muon but ";
      G4cout << theParentName << G4endl;
    }
#endif
  }
}

G4MuonDecayChannelM::~G4MuonDecayChannelM()
{
}

//Mandelstam "t" variable
G4double G4MuonDecayChannelM::tvar(G4double M, G4double m, G4double E){
  return M*M - 2*M*E + m*m;
}

G4DecayProducts *G4MuonDecayChannelM::DecayIt(G4double) 
{

  // this version neglects muon polarization,and electron mass  
  //              assumes the pure V-A coupling
  //              the Neutrinos are correctly V-A. 
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) G4cout << "G4MuonDecayChannel::DecayIt ";
#endif

  if (parent == 0) FillParent();  
  if (daughters == 0) FillDaughters();
 
  //parent mass (muon)
  G4double M = parent->GetPDGMass();

  //daughters'mass
  G4double daughtermass[3]; 
  G4double sumofdaughtermass = 0.0;
  for (G4int index=0; index<3; index++){
    daughtermass[index] = daughters[index]->GetPDGMass();
    sumofdaughtermass += daughtermass[index];
  }

  G4double mh = daughtermass[1];

   //create parent G4DynamicParticle at rest
  G4ThreeVector dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle( parent, dummy, 0.0);
  //create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  //Daughters momentum
  G4double daughtermomentum[3];

  G4double Ee;
  G4double gam;
   
  G4double cthetae,cthetaeh,cthetah;
  G4double d2,d3,d4;

  //Maximum positron energy
  G4double Emax = (M*M - mh*mh)/M/2.0;
  G4double Eh,Ehmax,Ehmin;

  //Generate random positron energy Ee and positron angle cos(thetae) 
  do {
    Ee = G4UniformRand()*Emax;
    cthetae = 1 - 2.0*G4UniformRand();
    gam = G4UniformRand()*1e8;
    d2 = M*M*Ee*Ee*(1-mh*mh/tvar(M,0,Ee))*(1-mh*mh/tvar(M,0,Ee)) * 
      (  ((3-4*Ee/M) + (3-2*Ee/M)*mh*mh/tvar(M,0,Ee)) + ( (1-4*Ee/M) - (1+2*Ee/M)*mh*mh/tvar(M,0,Ee) )*cthetae );
    //G4cout << Emax/MeV << " " << Ee/MeV << " " << cthetae << " " << d2 << G4endl;
  }
  while (gam>d2);
  
  //Maximum nu_h energy
  Ehmax = (M-Ee)/2.0*(1+mh*mh/tvar(M,0,Ee)) + Ee/2.0*(1-mh*mh/tvar(M,0,Ee));
  Ehmin = (M-Ee)/2.0*(1+mh*mh/tvar(M,0,Ee)) - Ee/2.0*(1-mh*mh/tvar(M,0,Ee));

  //Generate random nu_h energy
  do {
    Eh = Ehmin + G4UniformRand()*(Ehmax-Ehmin);
    gam = G4UniformRand()*1e5;
    cthetaeh = (M*M+mh*mh - 2*M*(Eh+Ee) + 2*Eh*Ee) / (2*sqrt(Eh*Eh-mh*mh)*Ee);
    d3 = (Eh*Ee-sqrt(Eh*Eh-mh*mh)*Ee*cthetaeh)*(M-Eh-Ee-Ee*cthetae-sqrt(Eh*Eh-mh*mh)*cthetaeh*cthetae);
    //G4cout << sqrt(Eh*Eh-mh*mh) << " " << cthetaeh << " " << Eh << " " << d3 << G4endl;
  }
  while (gam>d3);

  //Generate ramdom cos(theta_h)
  G4double cmax = cos(acos(cthetae)-acos(cthetaeh));
  G4double cmin = cos(acos(cthetae)+acos(cthetaeh));
  G4double cd1,cd2;
  do {
    cthetah = cmin + G4UniformRand()*(cmax-cmin);
    gam = G4UniformRand()*6e6;

    cd1 = cos(acos(cthetah)-acos(cthetae));
    cd2 = cos(acos(cthetah)+acos(cthetae));

    d4 = ( (M-(Eh+sqrt(Eh*Eh-mh*mh)*cthetah) - (Ee+Ee*cthetae) ) * (Eh*Ee-sqrt(Eh*Eh-mh*mh)*Ee*cthetaeh) ) / 
      sqrt((cd1-cthetaeh)*(cthetaeh-cd2));

    //G4cout << cthetah << " " << d4 << G4endl;

  }
  while (gam>d4);

  //Random 3D-Rotation matrix
  G4double rphi=twopi*G4UniformRand()*rad;
  G4double rtheta=(std::acos(2.*G4UniformRand()-1.));
  G4double rpsi=twopi*G4UniformRand()*rad;

  G4RotationMatrix rot;
  rot.set(rphi,rtheta,rpsi);
  
  //Set the momenta of the daughter particles:

  //electron 0
  G4double sthetae = sqrt(1-cthetae*cthetae);
  daughtermomentum[0]= Ee;
  G4ThreeVector direction0(sthetae*cthetae,sthetae*sthetae,cthetae);
  direction0 *= rot;
  G4DynamicParticle * daughterparticle = new G4DynamicParticle (daughters[0], direction0 * daughtermomentum[0]);
  products->PushProducts(daughterparticle);

  //massive neutrino  1
  G4double sthetah = sqrt(1-cthetah*cthetah);
  daughtermomentum[1]=sqrt(Eh*Eh-mh*mh);
  G4ThreeVector direction1(sthetah*cthetah,sthetah*sthetah,cthetah);
  direction1 *= rot;
  G4DynamicParticle * daughterparticle1 = new G4DynamicParticle (daughters[1], direction1 * daughtermomentum[1]);
  products->PushProducts(daughterparticle1);

  //electronic neutrino 2
  daughtermomentum[2] = M - Ee - Eh;
  G4ThreeVector direction2 = -direction0 - direction1;
  direction2 *= rot;
  G4DynamicParticle * daughterparticle2 = new G4DynamicParticle (daughters[2], direction2 * daughtermomentum[2]);
  products->PushProducts(daughterparticle2);




  // output message
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cout << "G4MuonDecayChannel::DecayIt ";
    G4cout << "  create decay products in rest frame " <<G4endl;
    products->DumpInfo();
  }
#endif
  return products;
}






