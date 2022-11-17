//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: myTrackerSD.cc,v 1.3.4.1 2001/06/28 19:07:31 gunter Exp $
// GEANT4 tag $Name:  $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#include "mySensitiveDetector.hh"
// #include "myDetectorConstruction.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include <new>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

mySensitiveDetector::mySensitiveDetector(G4String name, myDetectorConstruction* detector)
  : G4VSensitiveDetector(name)
  , Detector(detector)
  , SDName(name)
{
  G4String HCname;
  collectionName.insert("hitsCollection");
  HCID = -1;
  // first = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

mySensitiveDetector::~mySensitiveDetector() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void mySensitiveDetector::Initialize(G4HCofThisEvent* HCE)
{
  mySDHitsCollection = new myHitsCollection(SensitiveDetectorName, collectionName[0]);
  if(HCID < 0)
  {
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(mySDHitsCollection);
  }
  HCE->AddHitsCollection(HCID, mySDHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4bool mySensitiveDetector::ProcessHits(G4Step* aStep, G4TouchableHistory*) { return false; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void mySensitiveDetector::EndOfEvent(G4HCofThisEvent* HCE)
{
  //  if(numHits > 0) {
  //	G4cout << "\t"  << SDName << " " << numHits << G4endl;
  //	numHits = 0;
  //  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
