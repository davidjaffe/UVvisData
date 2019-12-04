#include "G4Material.hh"
#include "G4NistManager.hh"

#include "Materials.hh"

Materials::Materials()
{
  nistMan = G4NistManager::Instance();

  nistMan->SetVerbose(1);

  CreateMaterials();
}

Materials::~Materials(){}

Materials* Materials::instance = 0;

Materials* Materials::GetInstance()
{
  if (instance == 0)
    {
      instance = new Materials();
    }
  return instance;
}

G4Material* Materials::GetMaterial(G4String material)  
{
  G4Material* mat =  nistMan->FindOrBuildMaterial(material);

  if (!mat) mat = G4Material::GetMaterial(material);
  if (!mat) {G4cout << material << "Not Found, Please Retry"<< G4endl;}

  return mat;
}

void Materials::CreateMaterials()
{
  G4double density;

  std::vector<G4int>  natoms;
  std::vector<G4double> fractionMass;
  std::vector<G4String> elements;

  // Materials Definitions
  // =====================

  nistMan->FindOrBuildMaterial("G4_Al");
  nistMan->FindOrBuildMaterial("G4_Ag");
  nistMan->FindOrBuildMaterial("G4_Zr"); //DvB
  nistMan->FindOrBuildMaterial("G4_Sn"); //DvB
  nistMan->FindOrBuildMaterial("G4_Mo"); //DvB
  nistMan->FindOrBuildMaterial("G4_Cu"); //DvB
  nistMan->FindOrBuildMaterial("G4_Galactic");
  nistMan->FindOrBuildMaterial("G4_AIR");

  nistMan->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  nistMan->FindOrBuildMaterial("G4_BERYLLIUM_OXIDE");
  nistMan->FindOrBuildMaterial("G4_SODIUM_IODIDE");
  nistMan->FindOrBuildMaterial("G4_CESIUM_IODIDE");

  nistMan->FindOrBuildMaterial("G4_MYLAR");

  nistMan->FindOrBuildMaterial("G4_Si");
  // Set the Birks Constant for NaI(Tl)
  GetMaterial("G4_SODIUM_IODIDE")->
                              GetIonisation()->SetBirksConstant(0.024*mm/MeV);
  // Set the Birks Constant for CsI(Tl)
  GetMaterial("G4_CESIUM_IODIDE")->
                              GetIonisation()->SetBirksConstant(0.024*mm/MeV);

  //Define Scintillator

  elements.push_back("C");     natoms.push_back(9);
  elements.push_back("H");     natoms.push_back(10);

  density = 1.032*g/cm3; 

  nistMan->ConstructNewMaterial("Scintillator", elements, natoms, density);

  // Set the Birks Constant for plastic scintillator
  GetMaterial("Scintillator")->
                              GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  elements.clear();
  natoms.clear();


//   // Zirconium 
//   elements.push_back("Zr"); fractionMass.push_back(100*perCent);
//   density = 6.501*g/cm3;
//   nistMan->ConstructNewMaterial("Zirconium", elements, fractionMass, density);

//   elements.clear();
//   natoms.clear();

  //Define Rohacell

  elements.push_back("C");     natoms.push_back(9);
  elements.push_back("H");     natoms.push_back(13);
  elements.push_back("Na");    natoms.push_back(1);
  elements.push_back("O");     natoms.push_back(2);

  density = 0.03*g/cm3;

  nistMan->ConstructNewMaterial("Rohacell", elements, natoms, density);

  elements.clear();
  natoms.clear();

  //Define G10

  elements.push_back("C");     fractionMass.push_back(36.8*perCent);
  elements.push_back("O");     fractionMass.push_back(30.7*perCent);
  elements.push_back("Si");    fractionMass.push_back(21.9*perCent);
  elements.push_back("H");     fractionMass.push_back(6.2*perCent);
  elements.push_back("B");     fractionMass.push_back(2.2*perCent);
  elements.push_back("Na");    fractionMass.push_back(2.2*perCent);

  density = 1.7*g/cm3;

  nistMan->ConstructNewMaterial("G10", elements, fractionMass, density);

  elements.clear();
  fractionMass.clear();

  //Define industrial Iron

  elements.push_back("Fe");    fractionMass.push_back(99*perCent);
  elements.push_back("C");     fractionMass.push_back(1*perCent);

  density = 7.87*g/cm3;

  nistMan->ConstructNewMaterial("Iron", elements, fractionMass, density);

  elements.clear();
  fractionMass.clear();

  //Define Isobutane

  elements.push_back("C");     natoms.push_back(4);
  elements.push_back("H");     natoms.push_back(10);

  density =  2.47*mg/cm3;

  G4double temperature = 293.15*kelvin;
  G4double pressure = 1*atmosphere;

  nistMan->ConstructNewMaterial("Isobutane", elements, natoms, density,
                                true, kStateGas, temperature, pressure);

  elements.clear();
  natoms.clear();

  //Define CarbonTetrafluoride (CF4)

  elements.push_back("C");     natoms.push_back(1);
  elements.push_back("F");     natoms.push_back(4);

  density     = 3.66*mg/cm3;

  nistMan->ConstructNewMaterial("CF4", elements, natoms, density,
                                true, kStateGas, temperature, pressure);

  elements.clear();
  natoms.clear();

  //Define WCGas

  density = 3.56*mg/cm3;

  G4Material* WCGas = new G4Material("WCGas",density,2,kStateGas,temperature,pressure);

  WCGas->AddMaterial(GetMaterial("Isobutane"), 20.0*perCent);
  WCGas->AddMaterial(GetMaterial("CF4"), 80.0*perCent);



  //Define Ceramic (Ceramic)
  elements.push_back("Al");   natoms.push_back(2);
  elements.push_back("O");   natoms.push_back(3);
  density  = 3.69*g/cm3;
  nistMan->ConstructNewMaterial("Ceramic", elements, natoms, density);
  elements.clear();
  natoms.clear(); 


  //Define Vacuum
  //   elements.push_back("H");   natoms.push_back(1);
  // elements.push_back("H");   natoms.push_back(1);
  density  = (1E-25)*g/cm3;
 elements.push_back("H");    fractionMass.push_back(99*perCent);
 elements.push_back("C");    fractionMass.push_back(1*perCent);

   nistMan->ConstructNewMaterial("Vac", elements, fractionMass, density);
   elements.clear();
   fractionMass.clear(); 


}
