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
/// \file MCEventAction.cc
/// \brief Implementation of the CaTS::MCEventAction  class
// Geant4 headers
#include <G4UserEventAction.hh>
#include <G4ios.hh>
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include <G4String.hh>
#include <G4Types.hh>
#include <G4VHit.hh>
#include <G4VHitsCollection.hh>
#include "G4SDManager.hh"
// project headers:
#include "MCEventAction.hh"
#include "ConfigurationManager.hh"
#include "Event.hh"
#include "PhotonSD.hh"
#include "PhotonHit.hh"
#include "InteractionHit.hh"
#include "lArTPCHit.hh"
#include "TrackerHit.hh"
#include "MscHit.hh"
#include "CalorimeterHit.hh"
#include "DRCalorimeterHit.hh"
#include "SimEnergyDepositHit.hh"
#ifdef WITH_ROOT
#  include "RootIO.hh"
#endif
// stl headers
#include <map>
#include <utility>
#include <algorithm>
#include <istream>
#include <string>

#ifdef WITH_G4CXOPTICKS
#  include "SEvt.hh"
#  include "NP.hh"
#  include "G4CXOpticks.hh"
#endif

namespace
{
  G4Mutex opticks_mutex = G4MUTEX_INITIALIZER;
}

MCEventAction ::MCEventAction()
  : G4UserEventAction()
{
#ifdef WITH_ROOT
  if(ConfigurationManager::getInstance()->isWriteHits())
  {
    RootIO::GetInstance();
  }
#endif
}  // namespace MCEventAction::MCEventAction()

void MCEventAction ::BeginOfEventAction(const G4Event* event)
{  // if(verbose)
  G4cout << "MCEventAction BeginOfEventAction Event:   " << event->GetEventID() << G4endl;
}

void MCEventAction ::EndOfEventAction(const G4Event* event)
{
  G4bool verbose       = ConfigurationManager::getInstance()->isEnable_verbose();
  G4HCofThisEvent* HCE = event->GetHCofThisEvent();
  if(HCE == nullptr)
    return;
#ifdef WITH_ROOT
  Event* CaTSEvt = new Event();
  CaTSEvt->SetEventNr(event->GetEventID());
  std::map<G4String, std::vector<G4VHit*>>* hcmap = CaTSEvt->GetHCMap();
#endif  // end WITH_ROOT

#ifdef WITH_G4CXOPTICKS
  if(ConfigurationManager::getInstance()->isEnable_opticks())
  {
    if(ConfigurationManager::getInstance()->isEnable_verbose())
    {
      G4int inum_photon          = SEvt::GetNumPhotonFromGenstep(0);
      G4int inum_genstep         = SEvt::GetNumGenstepFromGenstep(0);
      G4int num_PhotonCollected  = SEvt::GetNumPhotonCollected(0);
      G4int num_PhotonGenstepMax = SEvt::GetNumPhotonGenstepMax(0);
      G4int num_Hit              = SEvt::GetNumHit(0);
      std::cout << "------------------------------------------------------------------"
                << std::endl;
      std::cout << "MCEventAction: GetNumPhotonFromGenstep: " << inum_photon << std::endl;
      std::cout << "MCEventAction: GetNumGenstepFromGenstep: " << inum_genstep << std::endl;
      std::cout << "MCEventAction: GetNumPhotonCollected:  " << num_PhotonCollected << std::endl;
      std::cout << "MCEventAction: GetNumPhotonGenstepMax: " << num_PhotonGenstepMax << std::endl;
      std::cout << "MCEventAction: GetNumHit:            " << num_Hit << std::endl;
      std::cout << "------------------------------------------------------------------"
                << std::endl;
    }

    cudaDeviceSynchronize();
    //        G4AutoLock lock(&opticks_mutex);
    // G4RunManager* rm     = G4RunManager::GetRunManager();
    // const G4Event* event = rm->GetCurrentEvent();
    G4int eventid = event->GetEventID();
    G4CXOpticks::Get()->simulate(eventid, true);

    /*

      //hjw   G4CXOpticks* g4cxok = G4CXOpticks::Get();
      G4int eventid       = event->GetEventID();
      //static void SEvt::SetIndex(eventid);
      G4int num_genstep = SEvt::GetNumGenstepFromGenstep(0);
      G4int num_photon  = SEvt::GetNumPhotonCollected(0);
      if(verbose)
      {
        G4cout << "MCEndOfEventAction: num_photon: " << num_photon << G4endl;
        G4cout << "MCEndOfEventAction: num_genstep: " << num_genstep << G4endl;
      }
      if(num_photon > 0)
      {
         G4CXOpticks::Get()->simulate(eventid);
        //      cudaDeviceSynchronize();
      }
  */
    // unsigned int num_hits = SEvt::GetNumHit(0);
    //   SEvt* sev             = SEvt::Get();
    /*
        if(verbose)
        {
            G4int inum_photon          = SEvt::GetNumPhotonFromGenstep(0);
            G4int inum_genstep         = SEvt::GetNumGenstepFromGenstep(0);
      G4int num_PhotonCollected  = SEvt::GetNumPhotonCollected(0);
      G4int num_PhotonGenstepMax = SEvt::GetNumPhotonGenstepMax(0);
      G4int num_Hit              = SEvt::GetNumHit(0);
          G4cout << "MCEndOfEventAction: num_hits: " << num_hits << G4endl;
          G4cout << "MCEndOfEventAction: num_photon: " << num_photon << G4endl;
          //      G4cout << "McEndOfEventAction: num_genstep: " << num_genstep << G4endl;
          // std::cout << sev->descFull();
        }
    */

    /*
        if(num_hits > 0)
        {
          G4HCtable* hctable = G4SDManager::GetSDMpointer()->GetHCtable();
          for(G4int i = 0; i < hctable->entries(); ++i)
          {
            std::string sdn   = hctable->GetSDname(i);
            std::size_t found = sdn.find("Photondetector");
            if(found != std::string::npos)
            {
              //          std::cout << "Photondetector: " << sdn << std::endl;
              PhotonSD* aSD = (PhotonSD*) G4SDManager::GetSDMpointer()->FindSensitiveDetector(sdn);
              aSD->AddOpticksHits();
            }
          }
        }
    */
    // U4HitExtra* hit_extra_ptr = way_enabled ? &hit_extra : nullptr;

    //     const std::string EvDescr = SEvt::desc();
    //  const NP* getPhoton() const ;
    //     std::cout << SEvt::descFull();
  }
#endif  //  WITH_G4CXOPTICKS

  //
  // Now we deal with the Geant4 Hit collections.
  //
  if(verbose)
    G4cout << "Number of collections:  " << HCE->GetNumberOfCollections() << G4endl;
#ifdef WITH_ROOT
  if(ConfigurationManager::getInstance()->isWriteHits())
  {
    for(int i = 0; i < HCE->GetNumberOfCollections(); i++)
    {
      G4VHitsCollection* hc      = HCE->GetHC(i);
      G4String hcname            = hc->GetName();
      std::vector<std::string> y = split(hcname, '_');
      std::string Classname      = y[1];
      if(verbose)
        G4cout << "Classname: " << Classname << G4endl;
      if(Classname == "lArTPC")
      {
        std::vector<G4VHit*> hitsVector;
        G4int NbHits = hc->GetSize();
        for(G4int ii = 0; ii < NbHits; ii++)
        {
          G4VHit* hit       = hc->GetHit(ii);
          lArTPCHit* tpcHit = dynamic_cast<lArTPCHit*>(hit);
          hitsVector.push_back(tpcHit);
        }
        hcmap->insert(std::make_pair(hcname, hitsVector));
      }
      else if(Classname == "Photondetector")
      {
        std::vector<G4VHit*> hitsVector;
        G4int NbHits = hc->GetSize();
        if(verbose)
          G4cout << "Photondetector size: " << hc->GetSize() << G4endl;
        for(G4int ii = 0; ii < NbHits; ii++)
        {
          G4VHit* hit     = hc->GetHit(ii);
          PhotonHit* pHit = dynamic_cast<PhotonHit*>(hit);
          hitsVector.push_back(pHit);
        }
        hcmap->insert(std::make_pair(hcname, hitsVector));
      }
      else if(Classname == "Target")
      {
        std::vector<G4VHit*> hitsVector;
        G4int NbHits = hc->GetSize();
        if(verbose)
          G4cout << "Interaction size: " << hc->GetSize() << G4endl;
        for(G4int ii = 0; ii < NbHits; ii++)
        {
          G4VHit* hit           = hc->GetHit(ii);
          InteractionHit* iaHit = dynamic_cast<InteractionHit*>(hit);
          hitsVector.push_back(iaHit);
        }
        hcmap->insert(std::make_pair(hcname, hitsVector));
      }
      else if(Classname == "Tracker")
      {
        std::vector<G4VHit*> hitsVector;
        G4int NbHits = hc->GetSize();
        if(verbose)
          G4cout << "Tracker size: " << hc->GetSize() << G4endl;
        for(G4int ii = 0; ii < NbHits; ii++)
        {
          G4VHit* hit      = hc->GetHit(ii);
          TrackerHit* tHit = dynamic_cast<TrackerHit*>(hit);
          hitsVector.push_back(tHit);
        }
        hcmap->insert(std::make_pair(hcname, hitsVector));
      }
      else if(Classname == "Msc")
      {
        std::vector<G4VHit*> hitsVector;
        G4int NbHits = hc->GetSize();
        if(verbose)
          G4cout << "Msc size: " << hc->GetSize() << G4endl;
        for(G4int ii = 0; ii < NbHits; ii++)
        {
          G4VHit* hit    = hc->GetHit(ii);
          MscHit* mscHit = dynamic_cast<MscHit*>(hit);
          hitsVector.push_back(mscHit);
        }
        hcmap->insert(std::make_pair(hcname, hitsVector));
      }
      else if(Classname == "Calorimeter")
      {
        std::vector<G4VHit*> hitsVector;
        G4int NbHits = hc->GetSize();
        if(verbose)
          G4cout << "Calorimeter size: " << hc->GetSize() << G4endl;
        for(G4int ii = 0; ii < NbHits; ii++)
        {
          G4VHit* hit          = hc->GetHit(ii);
          CalorimeterHit* cHit = dynamic_cast<CalorimeterHit*>(hit);
          hitsVector.push_back(cHit);
        }
        hcmap->insert(std::make_pair(hcname, hitsVector));
      }
      else if(Classname == "DRCalorimeter")
      {
        std::vector<G4VHit*> hitsVector;
        G4int NbHits = hc->GetSize();
        if(verbose)
          G4cout << "DRCalorimeter size: " << hc->GetSize() << G4endl;
        for(G4int ii = 0; ii < NbHits; ii++)
        {
          G4VHit* hit             = hc->GetHit(ii);
          DRCalorimeterHit* drHit = dynamic_cast<DRCalorimeterHit*>(hit);
          hitsVector.push_back(drHit);
        }
        hcmap->insert(std::make_pair(hcname, hitsVector));
      }
      else if(Classname == "SimEnergyDeposit")
      {
        std::vector<G4VHit*> hitsVector;
        G4int NbHits = hc->GetSize();
        if(verbose)
          G4cout << "SimEnergyDeposit size: " << hc->GetSize() << G4endl;
        for(G4int ii = 0; ii < NbHits; ii++)
        {
          G4VHit* hit                = hc->GetHit(ii);
          SimEnergyDepositHit* seHit = dynamic_cast<SimEnergyDepositHit*>(hit);
          hitsVector.push_back(seHit);
        }
        hcmap->insert(std::make_pair(hcname, hitsVector));
      }
      else
      {
        G4cout << "SD type: " << Classname << " unknown" << G4endl;
      }
    }
    G4AutoLock lock(&opticks_mutex);
    RootIO::GetInstance()->Write(CaTSEvt);
    CaTSEvt->Reset();
    delete CaTSEvt;
  }  // end enableio
#endif
  // if(verbose)
  G4cout << "MCEventAction EndOfEventAction Event:   " << event->GetEventID() << G4endl;
}

std::vector<std::string>& MCEventAction::split(const std::string& s, char delim,
                                               std::vector<std::string>& elems)
{
  std::stringstream ss(s);
  std::string item;
  while(std::getline(ss, item, delim))
  {
    elems.push_back(item);
  }
  return elems;
}

std::vector<std::string> MCEventAction::split(const std::string& s, char delim)
{
  std::vector<std::string> elems;
  return split(s, delim, elems);
}
