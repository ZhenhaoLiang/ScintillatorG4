/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
// #include "Run.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "g4root.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//G4String m_hDataFilename;
RunAction::RunAction()
{ 
  auto analysisManager = G4AnalysisManager::Instance();
 // G4AccumulableManager* analysisManager = G4AccumulableManager::Instance();
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);

  analysisManager->CreateNtuple("event", "Energy and Position");
  analysisManager->CreateNtupleDColumn("Energy");
  analysisManager->CreateNtupleDColumn("prex");
  analysisManager->CreateNtupleDColumn("prey");
  analysisManager->CreateNtupleDColumn("prez");
  analysisManager->CreateNtupleDColumn("postx");   
  analysisManager->CreateNtupleDColumn("posty");    //5
  analysisManager->CreateNtupleDColumn("postz");
  analysisManager->CreateNtupleSColumn("ptype");
  analysisManager->CreateNtupleDColumn("eventID");
  analysisManager->CreateNtupleDColumn("trackID");
  analysisManager->CreateNtupleDColumn("parentID");  //10
  analysisManager->CreateNtupleDColumn("dE"); 
  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

void RunAction::BeginOfRunAction(const G4Run*)
{
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  auto analysisManager = G4AnalysisManager::Instance();

  G4String filename = m_hDataFilename;//"event.root";
  analysisManager->OpenFile(filename);
  G4cout << "Using " << analysisManager->GetType() << G4endl;


}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

