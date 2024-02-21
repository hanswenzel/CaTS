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
// ********************************************************************
//
//  CaTS (Calorimetry and Tracking Simulation)
//
//  Authors : Hans Wenzel
//            Soon Yung Jun
//            (Fermi National Accelerator Laboratory)
//
// History
//   October 18th, 2021 : first implementation
//
// ********************************************************************
//
/// \file PhotonSD.cc
/// \brief Implementation of the CaTS::PhotonSD class

// Geant4 headers
#include "ConfigurationManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4OpticalPhoton.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
// project headers
#include "PhotonSD.hh"
#ifdef WITH_G4CXOPTICKS
#  include "scuda.h"
#  include "SEvt.hh"
#  include "G4CXOpticks.hh"
#  include "NP.hh"
#endif

PhotonSD::PhotonSD(G4String name)
  : G4VSensitiveDetector(name)
{
  G4String HCname = name + "_HC";
  collectionName.insert(HCname);
  G4cout << collectionName.size() << "   PhotonSD name:  " << name << " collection Name: " << HCname
         << G4endl;
  fHCID = -1;
}

void PhotonSD::Initialize(G4HCofThisEvent* hce)
{
  fPhotonHitsCollection = new PhotonHitsCollection(SensitiveDetectorName, collectionName[0]);
  if(fHCID < 0)
  {
    if(ConfigurationManager::getInstance()->isEnable_verbose())
      G4cout << "PhotonSD::Initialize:  " << SensitiveDetectorName << "   " << collectionName[0]
             << G4endl;
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }
  hce->AddHitsCollection(fHCID, fPhotonHitsCollection);
}

G4bool PhotonSD::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist)
{
  G4TouchableHistory* theTouchable =
    (G4TouchableHistory*) (aStep->GetPreStepPoint()->GetTouchable());
  G4VPhysicalVolume* physVol = theTouchable->GetVolume();
  G4int CopyNr               = DetectorID(physVol);
  G4Track* theTrack          = aStep->GetTrack();
  // we only deal with optical Photons:
  if(theTrack->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
  {
    return false;
  }
  G4double theEdep              = theTrack->GetTotalEnergy() / CLHEP::eV;
  const G4VProcess* thisProcess = theTrack->GetCreatorProcess();
  G4String processname;
  if(thisProcess != NULL)
    processname = thisProcess->GetProcessName();
  else
    processname = "No Process";
  int theCreationProcessid;
  if(processname == "Cerenkov")
  {
    theCreationProcessid = 0;
  }
  else if(processname == "Scintillation")
  {
    theCreationProcessid = 1;
  }
  else
  {
    theCreationProcessid = -1;
  }
  PhotonHit* newHit = new PhotonHit(
    CopyNr, theCreationProcessid, etolambda(theEdep), theTrack->GetGlobalTime(),
    aStep->GetPostStepPoint()->GetPosition(), aStep->GetPostStepPoint()->GetMomentumDirection(),
    aStep->GetPostStepPoint()->GetPolarization());
  fPhotonHitsCollection->insert(newHit);
  theTrack->SetTrackStatus(fStopAndKill);
  return true;
}

void PhotonSD::EndOfEvent(G4HCofThisEvent*)
{
   if(ConfigurationManager::getInstance()->isEnable_verbose())
  {
  G4int NbHits = fPhotonHitsCollection->entries();
  G4cout << " PhotonSD::EndOfEvent Number of PhotonHits:  " << NbHits << G4endl;
  }
}
#ifdef WITH_G4CXOPTICKS
void PhotonSD::AddOpticksHits()
{
  SEvt* sev             = SEvt::Get_EGPU();
  unsigned int num_hits = sev->GetNumHit(0);

  for(int idx = 0; idx < int(num_hits); idx++)
  {
    sphoton hit;
    sev->getHit(hit, idx);
    G4ThreeVector position     = G4ThreeVector(hit.pos.x, hit.pos.y, hit.pos.z);
    G4ThreeVector direction    = G4ThreeVector(hit.mom.x, hit.mom.y, hit.mom.z);
    G4ThreeVector polarization = G4ThreeVector(hit.pol.x, hit.pol.y, hit.pol.z);
    int theCreationProcessid;
    if(OpticksPhoton::HasCerenkovFlag(hit.flagmask))
    {
      theCreationProcessid = 0;
    }
    else if(OpticksPhoton::HasScintillationFlag(hit.flagmask))
    {
      theCreationProcessid = 1;
    }
    else
    {
      theCreationProcessid = -1;
    }
    PhotonHit* newHit = new PhotonHit(DetectorID(position), theCreationProcessid, hit.wavelength,
                                      hit.time, position, direction, polarization);
    fPhotonHitsCollection->insert(newHit);
    if(ConfigurationManager::getInstance()->isEnable_verbose())
    {
      G4cout << " Process ID: " << theCreationProcessid << " PhotonSD  pos.:" << hit.pos.x << "  "
             << hit.pos.y << "  "
             << "  " << hit.pos.z << "  mom.:  " << hit.mom.x << "  " << hit.mom.y << "  "
             << hit.mom.z << "  pol.:  " << hit.pol.x << "  "
             << "  " << hit.pol.y << "  " << hit.pol.z << " iiindex: " << hit.iindex << "  "
             << "  wavel.:  " << hit.wavelength << "  time:  " << hit.time
             << "  boundary flag:  " << hit.boundary_flag << "  identy:  " << hit.identity
             << "  orient_idx: " << hit.orient_idx << "  flagmask:  " << hit.flagmask << G4endl;
    }
  }
}
#endif
