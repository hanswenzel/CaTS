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
// ********************************************************************
//
//  CaTS (Calorimetry and Tracking Simulation)
//
//  Authors : Hans Wenzel
//            Soon Yung Jun
//            (Fermi National Accelerator Laboratory)
//
// History
//   February 9th, 2022 : first implementation
//
// ********************************************************************
//
/// \file SimEnergyDepositSD.cc
/// \brief Implementation of the CaTS::SimEnergyDepositSD class

// Geant4 headers
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
// project headers
#include "SimEnergyDeposit.hh"
#include "SimEnergyDepositSD.hh"
#include "ConfigurationManager.hh"
#define UNUSED(expr)                                                                               \
  do                                                                                               \
  {                                                                                                \
    (void) (expr);                                                                                 \
  } while(0)
using namespace std;
SimEnergyDepositSD::SimEnergyDepositSD(G4String name)
  : G4VSensitiveDetector(name)
{
  G4String HCname = name + "_HC";
  collectionName.insert(HCname);
  G4cout << collectionName.size() << "   TrackerSD name:  " << name
         << " collection Name: " << HCname << G4endl;
  fHCID   = -1;
  verbose = ConfigurationManager::getInstance()->isEnable_verbose();
  std::cout << "SimEnergyDepositSD: constructor" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SimEnergyDepositSD::Initialize(G4HCofThisEvent* hce)
{
  fSimEnergyDepositCollection =
    new SimEnergyDepositCollection(SensitiveDetectorName, collectionName[0]);
  if(fHCID < 0)
  {
    if(verbose)
      G4cout << "SimEnergyDepositSD::Initialize:  " << SensitiveDetectorName << "   "
             << collectionName[0] << G4endl;
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }
  hce->AddHitsCollection(fHCID, fSimEnergyDepositCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool SimEnergyDepositSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep == 0.)
    return false;
  G4Track* aTrack        = aStep->GetTrack();
  G4int TrackID          = aTrack->GetTrackID();
  SimEnergyDeposit* sdep = new SimEnergyDeposit(
    0, 0, TrackID, (float) aStep->GetPreStepPoint()->GetPosition().getX() / CLHEP::cm,
    (float) aStep->GetPreStepPoint()->GetPosition().getY() / CLHEP::cm,
    (float) aStep->GetPreStepPoint()->GetPosition().getZ() / CLHEP::cm,
    (float) aStep->GetPostStepPoint()->GetPosition().getX() / CLHEP::cm,
    (float) aStep->GetPostStepPoint()->GetPosition().getY() / CLHEP::cm,
    (float) aStep->GetPostStepPoint()->GetPosition().getZ() / CLHEP::cm,
    (float) aStep->GetPreStepPoint()->GetGlobalTime() / CLHEP::ns,
    (float) aStep->GetPostStepPoint()->GetGlobalTime() / CLHEP::ns, (float) edep / CLHEP::MeV);
  fSimEnergyDepositCollection->insert(sdep);
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SimEnergyDepositSD::EndOfEvent(G4HCofThisEvent*)
{
  G4int NbHits = fSimEnergyDepositCollection->entries();
  if(verbose)
    G4cout << " Number of SimEnergyDeposit:  " << NbHits << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
