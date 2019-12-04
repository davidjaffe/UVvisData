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
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4MuonDecayChannel.hh,v 1.6 2006/06/29 19:23:35 gunter Exp $
// GEANT4 tag $Name: geant4-09-01-patch-01 $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      30 May 1997 H.Kurashige
// ------------------------------------------------------------
#ifndef G4RadiativeMuonCaptureProcess_h
#define G4RadiativeMuonCaptureProcess_h 1

#include "G4ios.hh"
#include "globals.hh"
//#include "G4VDiscreteProcess.hh"
#include "G4HadronicProcess.hh"
#include "G4CrossSectionDataStore.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "RadiativeMuCapCrossSection.hh"


//class G4RadiativeMuonCaptureProcess : public G4VDiscreteProcess 
class G4RadiativeMuonCaptureProcess : public G4HadronicProcess
{
  // Class Decription
  // This class simulates the photons produces from radiative muon capture into the 1S state
  // of zirconium according to the theoretical prediction of M. Pospelov and D. McKeen
  
public:  // With Description
  //Constructors 
  G4RadiativeMuonCaptureProcess(const G4String& processName = "RadiativeMuonCapture");
  //  Destructor
  virtual ~G4RadiativeMuonCaptureProcess();
  
  virtual G4VParticleChange* PostStepDoIt( const G4Track&, const G4Step& );
  
  virtual G4bool IsApplicable(const G4ParticleDefinition&);
  
protected:
  
  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition* );
  
  // access to the cross section data store
  inline G4CrossSectionDataStore* GetCrossSectionDataStore()
  { return theCrossSectionDataStore; }
  
  // get cross section per element
  virtual G4double GetMicroscopicCrossSection(const G4DynamicParticle *aParticle, 
					      const G4Element *anElement, 
					      G4double aTemp );

  //RadiativeMuCapCrossSection *DataSet;
  G4CrossSectionDataStore* theCrossSectionDataStore;

  G4Nucleus targetNucleus;
  G4RadiativeMuonCaptureProcess *dispatch;

  inline void AddDataSet(G4VCrossSectionDataSet * aDataSet);

  G4Element*  ChooseAandZ(const G4DynamicParticle* , const G4Material* );
  
public:  // With Description
  
  
};  


#endif
