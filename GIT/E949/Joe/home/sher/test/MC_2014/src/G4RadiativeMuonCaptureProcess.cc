
///////////////////////////////////
//
// Radiative Muon Capture
//
// D. vom Bruch Oct 2013
//
///////////////////////////////////


#include "G4ParticleDefinition.hh"
#include "G4DecayProducts.hh"
#include "G4VDecayChannel.hh"
#include "Randomize.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"
#include "G4RotationMatrix.hh"
#include "G4RadiativeMuonCaptureProcess.hh"
#include "G4Step.hh"
#include "G4VParticleChange.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ThreeVector.hh"
#include "RadiativeMuCapCrossSection.hh"
#include "G4CrossSectionDataStore.hh"

using namespace std;


G4RadiativeMuonCaptureProcess::G4RadiativeMuonCaptureProcess(const G4String& aName)
//:G4VDiscreteProcess(aName)
  :G4HadronicProcess(aName)
{

  G4CrossSectionDataStore* theStore = GetCrossSectionDataStore();
  theStore->AddDataSet(new RadiativeMuCapCrossSection );

  if (GetVerboseLevel()>0) {
    G4cout << GetProcessName() << " is created "<< G4endl;
  }
  
}


G4RadiativeMuonCaptureProcess::~G4RadiativeMuonCaptureProcess()
{
}


// inline void G4RadiativeMuonCaptureProcess::AddDataSet(G4VCrossSectionDataSet * aDataSet)  {
//   theCrossSectionDataStore->AddDataSet(aDataSet);
//   DataSet = aDataSet;
// }

// get cross section per element
G4double G4RadiativeMuonCaptureProcess::GetMicroscopicCrossSection(const G4DynamicParticle *aParticle, 
								   const G4Element *anElement, 
								   G4double aTemp ){
  RadiativeMuCapCrossSection *DataSet = new RadiativeMuCapCrossSection();

  G4double sigma;

  G4int Z = anElement->GetZ();
  G4int A = anElement->GetA();

  G4bool IsApplicable = DataSet->RadiativeMuCapCrossSection::IsZAApplicable(aParticle, Z, A);

  if (IsApplicable){
    
    sigma = DataSet->RadiativeMuCapCrossSection::GetCrossSection(aParticle, anElement, aTemp);

  }
  else {
    sigma = 0;
  }

  return sigma;

}




G4Element * G4RadiativeMuonCaptureProcess::ChooseAandZ(const G4DynamicParticle *aParticle, const G4Material *aMaterial )
{
  G4double currentZ = 0;
  G4double currentN = 0;
  const G4int numberOfElements = aMaterial->GetNumberOfElements();
  const G4ElementVector *theElementVector = aMaterial->GetElementVector();
  
  if( numberOfElements == 1 ) 
    {
      currentZ = G4double( ((*theElementVector)[0])->GetZ());
      currentN = (*theElementVector)[0]->GetN();
      targetNucleus.SetParameters(currentN, currentZ);
      return (*theElementVector)[0];
    }
  
  const G4double *theAtomicNumberDensity = aMaterial->GetAtomicNumDensityVector();
  G4double aTemp = aMaterial->GetTemperature();
  G4double crossSectionTotal = 0;
  G4int i;
  std::vector<G4double> runningSum;
  for( i=0; i < numberOfElements; ++i )
    {
      runningSum.push_back(theAtomicNumberDensity[i] *
			   dispatch->GetMicroscopicCrossSection( aParticle, (*theElementVector)[i], aTemp));
      crossSectionTotal+=runningSum[i];
    }
  
  G4double random = G4UniformRand();
  for( i=0; i < numberOfElements; ++i )
    { 
      if( random<=runningSum[i]/crossSectionTotal )
	{
	  currentZ = G4double( ((*theElementVector)[i])->GetZ());
	  currentN = ((*theElementVector)[i])->GetN();
	  targetNucleus.SetParameters(currentN, currentZ);
	  return (*theElementVector)[i];
	}
    }
  currentZ = G4double((*theElementVector)[numberOfElements-1]->GetZ());
  currentN = (*theElementVector)[numberOfElements-1]->GetN();
  targetNucleus.SetParameters(currentN, currentZ);
  return (*theElementVector)[numberOfElements-1];
}


G4double G4RadiativeMuonCaptureProcess::GetMeanFreePath(
   const G4Track& aTrack, G4double /*previousStepSize*/, G4ForceCondition* condition  ) {
    G4double mfp;

//   // G4StepPoint* pPrePoint = aStep.GetPreStepPoint();
// //   G4VPhysicalVolume* prePhysVol = pPrePoint->GetPhysicalVolume();
  
// //   G4double mfp;

// //   if ( prePhysVol->GetName() == "Zr") {

    G4double E     = aTrack.GetKineticEnergy();    // kinetic energy
    G4double N     = 6.022e23;                     // Avogadro's number (in 1/mole)
    G4double rho   = 6250*1000;                   // denisty of zirconium (in g/m^3)
    G4double A     = 91.224;                       // mass of one mole of zirconium (in u = g/mole)
    G4double sigma;  
    
    G4ThreeVector pVec = aTrack.GetMomentum();      // vector containing momentum
    G4double px = pVec.getX();
    G4double py = pVec.getY();
    G4double pz = pVec.getZ();
    G4double p = std::sqrt(px*px + py*py + pz*pz);
    
    
    //density of atoms per m^3
    G4double n = N*rho/A;                         //in 1/m^3
    
    // only calculate ARC cross section if momentum is between 20 and 50 MeV
    // otherwise set to something very small
    if (p > 20 && p < 50) {
      
      // cross section of ARC process;
      // parameters from fit to Kr cross section using 6th order polynomial
      G4double p0 =  9.043734e+00;
      G4double p1 = 1.812840e+00;
      G4double p2 = -2.418358e-01;
      G4double p3 = 1.081611e-02;
      G4double p4 = -2.335374e-04;
      G4double p5 = 2.483614e-06;
      G4double p6 = -1.045962e-08;
      
      G4double sigma_val = p0 + p1*p + p2*p*p + p3*p*p*p + p4*p*p*p*p + p5*p*p*p*p*p + p6*p*p*p*p*p*p;
      sigma = sigma_val*1e-28;
      
      G4cout << "momentum " << p <<  " E_kin = " << E << " ARC cross section " << sigma << G4endl;
      
    //G4double sigma = RadiativeMuCapCrossSection::GetCrossSection(const G4DynamicParticle* theProjectile, const G4Element* theTarget,  G4double theTemperature);
    
    mfp = 1/(sigma*n) * cm;
   
    }
    else {
      
      mfp = 10000 * cm;
    }
    
    
//     *condition = Forced;
    
    return mfp;
    
    //}
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4VParticleChange*
  G4RadiativeMuonCaptureProcess::PostStepDoIt( const G4Track& aTrack,
					       const G4Step& aStep)
{
  
 //  G4StepPoint* pPrePoint = aStep.GetPreStepPoint();
//   G4VPhysicalVolume* prePhysVol = pPrePoint->GetPhysicalVolume();
//   if ( prePhysVol->GetName() == "Zr") {

  G4double Ekin = aTrack.GetKineticEnergy();
  G4cout << "Ekin " << Ekin << G4endl;
  G4Track* sec;
   
  G4StepPoint* pPrePoint = aStep.GetPreStepPoint();
  G4VPhysicalVolume* prePhysVol = pPrePoint->GetPhysicalVolume();
  if ( prePhysVol->GetName() == "Zr") {

    //Initialize particle change
    aParticleChange.Initialize(aTrack);
    
    // Kill current track, then create scattered
    // phonon as a secondary
    aParticleChange.ProposeEnergy(0.);
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    aParticleChange.SetNumberOfSecondaries(1);
    
    
    G4double EBind_1S = 0.003643*GeV;
    //G4double EBind_2S = 0.001021*GeV;
    
    G4double postE = Ekin + EBind_1S;
    
    G4ThreeVector vgroup; // = new G4ThreeVector(x,y,z);
    const G4String name = "gamma";
    // search in particle table
    G4ParticleTable *pTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* photon = pTable->FindParticle(name);
    G4ThreeVector momentum = (G4ThreeVector)aTrack.GetMomentum();
    sec=new G4Track(new G4DynamicParticle(photon, momentum.unit(), postE), aTrack.GetGlobalTime(), aTrack.GetPosition());
    aParticleChange.AddSecondary(sec);
    G4cout << "PHOTON E_kin = " << postE <<  " in " << prePhysVol->GetName() <<  G4endl; 
  }
  else {

    //Initialize particle change
    aParticleChange.Initialize(aTrack);

    // Kill current track, then create scattered
    // phonon as a secondary
    aParticleChange.ProposeEnergy(0.);
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    aParticleChange.SetNumberOfSecondaries(1);

    G4ThreeVector vgroup; // = new G4ThreeVector(x,y,z);
    const G4String name = "mu-";
    // search in particle table
    G4ParticleTable *pTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* muon = pTable->FindParticle(name);
    G4ThreeVector momentum = (G4ThreeVector)aTrack.GetMomentum();
    sec=new G4Track(new G4DynamicParticle(muon, momentum.unit(), Ekin), aTrack.GetGlobalTime(), aTrack.GetPosition());
    aParticleChange.AddSecondary(sec);
    G4cout << "MUON E_kin = " << Ekin << " in " << prePhysVol->GetName() << G4endl; 

  }
 
  
  return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4RadiativeMuonCaptureProcess::IsApplicable(const G4ParticleDefinition& )
{
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* photon = pTable->FindParticle("gamma");
  //return (photon);
  return true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....




