#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeneralParticleSource.hh"
#include "Randomize.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
    fGPS = new G4GeneralParticleSource();

    // 1. 设置粒子类型（例如伽马）
    fGPS->GetCurrentSource()->SetParticleDefinition(
        G4ParticleTable::GetParticleTable()->FindParticle("gamma"));

    // 2. 设置圆柱体源形状
    fGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Volume");
    fGPS->GetCurrentSource()->GetPosDist()->SetPosDisShape("Cylinder");
    fGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(0,5.5*cm,0));
    fGPS->GetCurrentSource()->GetPosDist()->SetRadius(1.5*cm);
    fGPS->GetCurrentSource()->GetPosDist()->SetHalfZ(2.0*cm); 

    // 3. 设置各向同性发射
    fGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("iso");
    //fGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("direction");
    //fGPS->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));


    // 4. === Here energy spectrum need to be improved ===
    fGPS->GetCurrentSource()->GetEneDist()->SetMonoEnergy(1.25*MeV);
}


PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fGPS;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    fGPS->GeneratePrimaryVertex(anEvent);
}