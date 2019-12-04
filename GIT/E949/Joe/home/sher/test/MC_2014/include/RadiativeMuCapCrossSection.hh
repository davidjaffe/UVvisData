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
/// \file hadronic/Hadr02/include/G4DPMJET2_5CrossSection.hh
/// \brief Definition of the G4DPMJET2_5CrossSection class
//
#ifndef RadiativeMuCapCrossSection_h
#define RadiativeMuCapCrossSection_h
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              RadiativeMuCapCrossSection
//
// Version:             0.A
// Date:                1/11/2013
// Author:              D. vom Bruch
// Organisation:        TRIUMF
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// Class Description
// Cross section for ARC process
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"

#include "G4VCrossSectionDataSet.hh"


#include <map>

///////////////////////////////////////////////////////////////////////////////
//
class RadiativeMuCapCrossSection : public G4VCrossSectionDataSet
{
public:
  RadiativeMuCapCrossSection();
  ~RadiativeMuCapCrossSection();
  virtual G4bool IsApplicable(const G4DynamicParticle* theProjectile,
			      const G4Element* theTarget);
  
  virtual G4bool IsZAApplicable(const G4DynamicParticle* theProjectile,
				G4int Z, G4int A);
  
  virtual G4double GetCrossSection(const G4DynamicParticle* theProjectile,
				   const G4Element* theTarget, G4double theTemperature);
  
  virtual 
  G4double GetIsoZACrossSection(const G4DynamicParticle* theProjectile,
				G4int ZZ, G4int AA, 
				G4double theTemperature);
  
  virtual void BuildPhysicsTable(const G4ParticleDefinition&);
  
  virtual void DumpPhysicsTable(const G4ParticleDefinition&);
  
private:
  


  //RadiativeMuCapCrossSectionIndex theCrossSectionIndex;
    
    const G4double upperLimit;
    const G4double lowerLimit;
    const G4int    maxA;
};
///////////////////////////////////////////////////////////////////////////////
//
#endif
