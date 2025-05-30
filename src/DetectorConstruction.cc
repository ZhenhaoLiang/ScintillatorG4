/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include <G4Tubs.hh>
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include <G4RotationMatrix.hh>
#include <G4ThreeVector.hh>
#include <G4UnionSolid.hh>
#include <G4SubtractionSolid.hh>
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include <G4RotationMatrix.hh>
#include "G4SystemOfUnits.hh"
#include <G4VisAttributes.hh>

#define pi 3.14159265359

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction()
{
  DefineMaterial();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::DefineMaterial()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  Air = nist->FindOrBuildMaterial("G4_AIR");
  const G4int NUMENTRIE = 2;

  // 定义光子能量范围（单位：eV）
  G4double photonEnergyAir[NUMENTRIE] = {1.5*eV, 10.0*eV};  // ~200nm 到 ~800nm

  // 空气的折射率约为 1.0003
  G4double rIndexAir[NUMENTRIE] = {1.0003, 1.0003};

  // 创建光学属性表
  G4MaterialPropertiesTable* mptAir = new G4MaterialPropertiesTable();
  mptAir->AddProperty("RINDEX", photonEnergyAir, rIndexAir, NUMENTRIE);

  // 可选：加入吸收长度（这里假设非常大，近似透明）
  G4double absorptionAir[NUMENTRIE] = {1000.*m, 1000.*m};
  mptAir->AddProperty("ABSLENGTH", photonEnergyAir, absorptionAir, NUMENTRIE);

  // 将属性表附加到材料
  Air->SetMaterialPropertiesTable(mptAir);

  Water = nist->FindOrBuildMaterial("G4_WATER");
  Al = nist->FindOrBuildMaterial("G4_Al");
  //Co60
  G4Isotope* isoCo60 = new G4Isotope("Co60", 27, 60, 59.933817*g/mole); // Z=27, A=60
    
  // 创建富集Co-60的元素（100%富集）
  G4Element* elCo60 = new G4Element("EnrichedCo60", "Co*", 1);
  elCo60->AddIsotope(isoCo60, 1.0); // 100% Co-60
  
  // 创建金属钴材料（密度8.9 g/cm³），可选择使用天然钴或富集Co-60
  Co60 = new G4Material("Co60_Metal", 8.9*g/cm3, 1, kStateSolid);
  Co60->AddElement(elCo60, 1.0); // 使用富集Co-60

  //EJ-200
  // 基本材料组成 (C10H11)
  G4Element* elH = nist->FindOrBuildElement("H");
  G4Element* elC = nist->FindOrBuildElement("C");
  
  EJ200 = new G4Material("EJ200", 1.023*g/cm3, 2); // 密度1.023g/cm3
  EJ200->AddElement(elC, 10);  // 10个碳原子
  EJ200->AddElement(elH, 11);  // 11个氢原子
  
  // 设置光学特性
  const G4int NUMENTRIES = 32;
  G4double ppckov[NUMENTRIES] = 
      {2.00*eV, 2.25*eV, 2.50*eV, 2.75*eV, 3.00*eV, 3.25*eV, 
        3.50*eV, 3.75*eV, 4.00*eV, 4.25*eV, 4.50*eV, 4.75*eV,
        5.00*eV, 5.25*eV, 5.50*eV, 5.75*eV, 6.00*eV, 6.25*eV,
        6.50*eV, 6.75*eV, 7.00*eV, 7.25*eV, 7.50*eV, 7.75*eV,
        8.00*eV, 8.25*eV, 8.50*eV, 8.75*eV, 9.00*eV, 9.25*eV,
        9.50*eV, 9.75*eV};
  
  // 折射率 (典型值 ≈1.58)
  G4double rindexEJ200[NUMENTRIES];
  for(int i=0; i<NUMENTRIES; i++) rindexEJ200[i] = 1.58;
  
  // 吸收长度 (假设10m)
  G4double absorptionEJ200[NUMENTRIES];
  for(int i=0; i<NUMENTRIES; i++) absorptionEJ200[i] = 10.*m;
  
  // 光输出 (归一化的发射谱)
  G4double scintilEJ200[NUMENTRIES] = 
      {0.00, 0.01, 0.05, 0.12, 0.20, 0.28, 
        0.35, 0.42, 0.48, 0.55, 0.62, 0.68,
        0.75, 0.82, 0.88, 0.95, 1.00, 0.95,
        0.88, 0.80, 0.70, 0.60, 0.50, 0.40,
        0.30, 0.20, 0.10, 0.05, 0.01, 0.00,
        0.00, 0.00};
  
  // 添加到材料属性表
  G4MaterialPropertiesTable* mptEJ200 = new G4MaterialPropertiesTable();
  mptEJ200->AddProperty("RINDEX", ppckov, rindexEJ200, NUMENTRIES);
  mptEJ200->AddProperty("ABSLENGTH", ppckov, absorptionEJ200, NUMENTRIES);
  mptEJ200->AddProperty("FASTCOMPONENT", ppckov, scintilEJ200, NUMENTRIES);
  mptEJ200->AddConstProperty("SCINTILLATIONYIELD", 10000./MeV); // 光子产额
  mptEJ200->AddConstProperty("RESOLUTIONSCALE", 1.0);           // 分辨率
  mptEJ200->AddConstProperty("FASTTIMECONSTANT", 2.1*ns);       // 衰减时间
  mptEJ200->AddConstProperty("YIELDRATIO", 1.0);                // 快速成分比例
  
  EJ200->SetMaterialPropertiesTable(mptEJ200);

  // EJ-200 -> Air
  G4double ePhoton[] = {2.00*eV, 9.75*eV};
  const G4int num = sizeof(ePhoton)/sizeof(G4double);
  stickToAir = new G4OpticalSurface("StickAir");

  stickToAir->SetType(dielectric_dielectric);
  stickToAir->SetFinish(polished);
  stickToAir->SetModel(glisur);

  G4double reflectivityStickToAir[] = {0.5, 0.5};// 无反射
  assert(sizeof(reflectivityStickToAir) == sizeof(ePhoton));
  G4double efficiencyStickToAir[] = {0.5, 0.5};// 无吸收/再发射

  assert(sizeof(efficiencyStickToAir) == sizeof(ePhoton));
 
  G4MaterialPropertiesTable* stickToAirProperty =
  new G4MaterialPropertiesTable();

  stickToAirProperty->AddProperty("REFLECTIVITY", ePhoton,
                                  reflectivityStickToAir, num);

  stickToAirProperty->AddProperty("EFFICIENCY", ePhoton,
                                  efficiencyStickToAir, num);

  stickToAir->SetMaterialPropertiesTable(stickToAirProperty);

  G4cout << "EJ200 : density " <<  EJ200->GetDensity()/(g / cm3) << " , "
         << "NbOfKAtomsPerVolume " << EJ200->GetTotNbOfAtomsPerVolume()/(1. / cm3) /3.0<< G4endl;

  //==== Similarly here define EJ-276 ====

  //G4double EJ276_density = 1.099 * g / cm3;
  G4double H_density = 4.647e+22*(1. / cm3);
  G4double C_density = 4.944e+22*(1. / cm3);

  G4double H_MassDensity = (H_density/CLHEP::Avogadro)*(elH->GetAtomicMassAmu()*(g / mole));
  G4double C_MassDensity = (C_density/CLHEP::Avogadro)*(elC->GetAtomicMassAmu()*(g / mole));

  G4double HC_density = H_MassDensity + C_MassDensity;
  G4double H_Frac = H_MassDensity/HC_density;
  G4double C_Frac 	= C_MassDensity/HC_density;

  EJ276 = new G4Material("EJ276", HC_density, 2, kStateSolid);
  EJ276->AddElement(elH,H_Frac);
  EJ276->AddElement(elC,C_Frac);

  const G4int NUM = 15;
  G4double photonEnergy[NUM] = { 
      1.77*eV, 1.96*eV, 2.07*eV, 2.17*eV, 2.28*eV,  // 350-450nm
      2.38*eV, 2.48*eV, 2.58*eV, 2.76*eV, 2.88*eV,  // 450-550nm 
      3.00*eV, 3.10*eV, 3.26*eV, 3.44*eV, 3.54*eV   // 550-650nm
  };
  
  // 折射率 
  G4double rindex[NUM];
  for (int i=0; i<NUM; i++) rindex[i] = 1.58;  

  // 吸收长度
  G4double absorption[NUM];
  for (int i=0; i<NUM; i++) absorption[i] = 3.0*m;  

  // 发射光谱
  G4double fastComponent[NUM] = { 
      0.01, 0.15, 0.40, 0.75, 0.90,  // 快成分（中子响应）
      1.00, 0.95, 0.80, 0.60, 0.30,
      0.10, 0.05, 0.01, 0.00, 0.00
  };
  G4double slowComponent[NUM] = { 
      0.00, 0.05, 0.20, 0.45, 0.70,  // 慢成分（伽马响应）
      0.85, 0.95, 1.00, 0.90, 0.75,
      0.50, 0.30, 0.10, 0.05, 0.00
  };

  // 创建材料属性表
  G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("RINDEX",       photonEnergy, rindex,        NUM);
  mpt->AddProperty("ABSLENGTH",    photonEnergy, absorption,    NUM);
  mpt->AddProperty("FASTCOMPONENT",photonEnergy, fastComponent, NUM);
  mpt->AddProperty("SLOWCOMPONENT",photonEnergy, slowComponent, NUM);

  // 关键闪烁参数（根据EJ-276技术手册）
  mpt->AddConstProperty("SCINTILLATIONYIELD",  8000./MeV);  // 光产额
  mpt->AddConstProperty("RESOLUTIONSCALE",     1.0);       // 能量分辨率
  mpt->AddConstProperty("FASTTIMECONSTANT",    3.2*ns);    // 快成分衰减时间（中子）
  mpt->AddConstProperty("SLOWTIMECONSTANT",   42.0*ns);    // 慢成分衰减时间（伽马）
  mpt->AddConstProperty("YIELDRATIO",          0.5);       // 快/慢成分比例

  // 绑定到材料
  EJ276->SetMaterialPropertiesTable(mpt);

  G4cout << "EJ276 : density " <<  EJ276->GetDensity()/(g / cm3) << " , "
         << "NbOfAtomsPerVolume " << EJ276->GetTotNbOfAtomsPerVolume()/(1. / cm3) << G4endl;
}
G4VPhysicalVolume* DetectorConstruction::Construct()
{  
  
  // Scintillator Box
  //
  G4double  ScintillatorSize= 6.0*cm;
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;
 
  // World
  G4double world_sizeXYZ = 20.0*cm;
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXYZ, 0.5*world_sizeXYZ, 0.5*world_sizeXYZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        Al,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking   

  // Scintillator==============================================================
   
  G4Box* solidScintillator =    
    new G4Box("solidScintillator",                    //its name
        0.5*ScintillatorSize, 0.5*ScintillatorSize, 0.5*ScintillatorSize); //its size
      
  G4LogicalVolume* logicScintillator =                         
    new G4LogicalVolume(solidScintillator,            //its solid
                        EJ200,             //its material, Use EJ200 or EJ276
                        "logicScintillator");         //its name
               
  G4VPhysicalVolume* phyScint =  
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,0),    //at (0,0,0)
                    logicScintillator,       //its logical volume
                    "Scintillator",          //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking  
                    
  // Add surface
  auto EJ200WorldSurface = new G4LogicalBorderSurface("EJ200WorldSurface",phyScint, physWorld, stickToAir);
  
  //Source==============================================================
  G4RotationMatrix *CylinderRotate = new G4RotationMatrix();
  CylinderRotate->rotateX(90. * deg);
  G4double InnerSourceRadius = 2.5*0.5*cm;
  G4double OutterSourceRadius = 3.5*0.5*cm;
  G4double SourceHalfLength = 0.2*cm;

  G4Tubs* solidSourceCylinder =    
    new G4Tubs("SourceCylinder",                    //its name
        InnerSourceRadius,OutterSourceRadius, SourceHalfLength, 0.*deg, 360.*deg); //its size
      
  G4LogicalVolume* logicSourceCylinder =                         
    new G4LogicalVolume(solidSourceCylinder,            //its solid
                        Air,             //its material
                        "SourceCylinder");         //its name
               
  new G4PVPlacement(CylinderRotate,                       //no rotation
                    G4ThreeVector(0,SourceHalfLength+ScintillatorSize*0.5,0),         //at (0,0,0)
                    logicSourceCylinder,                //its logical volume
                    "SourceCylinder",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking      

  //Detector===============================================================  
  G4double PMTRadius = 2.54*cm;
  G4Tubs* solidDetector =    
    new G4Tubs("solidDetector",                    //its name
        0,PMTRadius, 1*cm,  0.*deg, 360.*deg); //its size
  G4LogicalVolume* logicDetector =                         
    new G4LogicalVolume(solidDetector,         //its solid
                        Air,          //its material
                        "Detector");           //its name

  G4VPhysicalVolume* phyDetector =
  new G4PVPlacement(CylinderRotate,                       //no rotation
                    G4ThreeVector(0,-ScintillatorSize*0.5-1.1*cm,0),   //at position
                    logicDetector,             //its logical volume
                    "Detector",                //its name
                    logicWorld,                //its mother volume  is contanier
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  new G4LogicalBorderSurface("EJ200WorldSurface",phyDetector, physWorld, stickToAir);

  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
