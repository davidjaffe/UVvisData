//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// *                                                                  *
// * Parts of this code which have been  developed by QinetiQ Ltd     *
// * under contract to the European Space Agency (ESA) are the        *
// * intellectual property of ESA. Rights to use, copy, modify and    *
// * redistribute this software for general public use are granted    *
// * in compliance with any licensing, distribution and development   *
// * policy adopted by the Geant4 Collaboration. This code has been   *
// * written by QinetiQ Ltd for the European Space Agency, under ESA  *
// * contract 19770/06/NL/JD (Technology Research Programme).         *
// *                                                                  *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file hadronic/Hadr02/src/G4DPMJET2_5CrossSection.cc
/// \brief Implementation of the G4DPMJET2_5CrossSection class
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4DPMJET2_5CrossSection.cc
//
// Version:             0.A
// Date:                02/04/08
// Author:              P R Truscott
// Organisation:        QinetiQ Ltd, UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            19770/06/NL/JD
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//

#include "RadiativeMuCapCrossSection.hh"
#include "G4ParticleTable.hh"
#include "G4DynamicParticle.hh"
#include "G4IonTable.hh"

#include "G4HadronicException.hh"
#include "G4StableIsotopes.hh"
#include "G4HadTmpUtil.hh"

#include "globals.hh"

#include <iomanip>
#include <fstream>
#include <sstream>

#include "G4Element.hh"


using namespace std;

///////////////////////////////////////////////////////////////////////////////
//
RadiativeMuCapCrossSection::RadiativeMuCapCrossSection ():
  upperLimit ( 50 * MeV ), lowerLimit ( 20 * MeV ), maxA(240)
{
  //theCrossSectionIndex.clear();
  
//
//
// vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
// This next bit is provisional, stating that this cross-section estimator
// is applicable to hydrogen targets.  However, the cross-section will be
// set to zero.
//
  //ATmin = 1;
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//
}
///////////////////////////////////////////////////////////////////////////////
//
RadiativeMuCapCrossSection::~RadiativeMuCapCrossSection ()
{
//
// Go through the list of cross-section fit parameters and delete the arrays.
//
  G4cout << "RadiativeMuCapCrossSection::~RadiativeMuCapCrossSection" << G4endl;
  //G4cout << "Size: " << theCrossSectionIndex.size() << G4endl;
  /*  
  if(theCrossSectionIndex.size() > 0) {

    RadiativeMuCapCrossSectionIndex::iterator it;
    for (it=theCrossSectionIndex.begin(); it!=theCrossSectionIndex.end(); ++it)
      {
        RadiativeMuCapCrossSectionParamSet *ptr = it->second;
        for (RadiativeMuCapCrossSectionParamSet *ptr1=ptr; ptr1<ptr+maxA; ptr1++)
          { delete ptr1; }
      }
  }
  */
}
///////////////////////////////////////////////////////////////////////////////
//
G4bool RadiativeMuCapCrossSection::IsApplicable
  (const G4DynamicParticle* theProjectile, const G4Element* theTarget)
{
  
  G4bool result;

  G4int Z    = G4lrint(theTarget->GetZ());
  G4int A    = G4lrint(theTarget->GetA());
  
  result = IsZAApplicable(theProjectile, Z, A);
  
  G4cout << "RadiativeMuCapCrossSection::IsApplicable E(GeV)= "
	 << theProjectile->GetKineticEnergy()/GeV << " off "
	 << theTarget->GetName() << " - " << result << G4endl;

  return result;
}
///////////////////////////////////////////////////////////////////////////////
//
G4bool RadiativeMuCapCrossSection::IsZAApplicable
  (const G4DynamicParticle* theProjectile, G4int Z, G4int A)
{

  G4bool result;

  if (Z == 40 && Z == 91) {
    result = true;
  }
  else {
    result = false;
  }  

  return result;
}
///////////////////////////////////////////////////////////////////////////////
//
G4double RadiativeMuCapCrossSection::GetIsoZACrossSection
  (const G4DynamicParticle* theProjectile, G4int ZZ, G4int AA,
   G4double /*theTemperature*/)
{

  G4double E     = theProjectile->GetTotalEnergy();
  G4double N     = 6.022e23;                     // Avogadro's number (in 1/mole)
  G4double rho   = 6250*1000;                   // denisty of zirconium (in g/m^3)
  G4double A     = 91.224;                       // mass of one mole of zirconium (in u = g/mole)
  G4double m_mu  = 105.658;                      // muon mass in MeV
  G4double sigma;  
  
 //  G4ThreeVector pVec = aTrack.GetMomentum();      // vector containing momentum
//   G4double px = pVec.getX();
//   G4double py = pVec.getY();
//   G4double pz = pVec.getZ();
//   G4double p = std::sqrt(px*px + py*py + pz*pz);
  
  G4double p = std::sqrt(E*E - m_mu*m_mu);   // get momentum from total energy and muon mass
  
  // density of atoms per m^3
  G4double n = N*rho/A;                         // in 1/m^3
  
  // // only calculate ARC cross section if momentum is between 20 and 50 MeV
//   // otherwise set to something very small
//   if (p > 20 && p < 50) {
    
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
  //sigma = sigma / 10 * millibarn;
  
  return sigma;
}
///////////////////////////////////////////////////////////////////////////////
//
G4double RadiativeMuCapCrossSection::GetCrossSection
  (const G4DynamicParticle* theProjectile, const G4Element* theTarget,
  G4double theTemperature)
{
  G4double xsection;

  G4int Z = 40;
  G4int A = 91;

  xsection = GetIsoZACrossSection(theProjectile, Z, A, theTemperature);
  
  return xsection;
}
///////////////////////////////////////////////////////////////////////////////
//

///////////////////////////////////////////////////////////////////////////////
//
void RadiativeMuCapCrossSection::BuildPhysicsTable(const G4ParticleDefinition&)
{;}
///////////////////////////////////////////////////////////////////////////////
//
void RadiativeMuCapCrossSection::DumpPhysicsTable(const G4ParticleDefinition 
  &theProjectile)
{
  ;
}
