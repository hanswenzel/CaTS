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
/// \file DetectorConstruction.cc
/// \brief Implementation of the CaTS::DetectorConstruction class

// Geant4 headers
#include "G4GDMLParser.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"
// project headers
#include "CalorimeterSD.hh"
#include "ColorReader.hh"
#include "ConfigurationManager.hh"
#include "DRCalorimeterSD.hh"
#include "DetectorConstruction.hh"
#include "InteractionSD.hh"
#include "MscSD.hh"
#include "PhotonSD.hh"
#include "RadiatorSD.hh"
#include "TrackerSD.hh"
#include "SimTrajectorySD.hh"
#include "lArTPCSD.hh"
// c++ headers
#include <iostream>
DetectorConstruction::DetectorConstruction(G4String fname)
  : G4VUserDetectorConstruction()
  , gdmlFile(fname)
{}

DetectorConstruction::~DetectorConstruction() {}
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  verbose = ConfigurationManager::getInstance()->isEnable_verbose();
  ReadGDML();
  const G4GDMLAuxMapType* auxmap = parser->GetAuxMap();
  if(verbose)
  {
    G4cout << "Found " << auxmap->size() << " volume(s) with auxiliary information." << G4endl
           << G4endl;
  }
  for(G4GDMLAuxMapType::const_iterator iter = auxmap->begin(); iter != auxmap->end(); iter++)
  {
    if(verbose)
    {
      G4cout << "Volume " << ((*iter).first)->GetName()
             << " has the following list of auxiliary information: " << G4endl;
    }
    for(G4GDMLAuxListType::const_iterator vit = (*iter).second.begin(); vit != (*iter).second.end();
        vit++)
    {
      if(verbose)
      {
        G4cout << "--> Type: " << (*vit).type << " Value: " << (*vit).value << G4endl;
      }
      if((*vit).type == "StepLimit")
      {
        G4UserLimits* fStepLimit = new G4UserLimits(atof((*vit).value));
        ((*iter).first)->SetUserLimits(fStepLimit);
      }
    }
  }
  G4VPhysicalVolume* worldPhysVol = parser->GetWorldVolume();
  if(ConfigurationManager::getInstance()->isDumpgdml())
  {
    std::ifstream ifile;
    ifile.open(ConfigurationManager::getInstance()->getGDMLFileName());
    if(ifile)
    {
      G4cout << "****************************************************" << G4endl;
      G4cout << ConfigurationManager::getInstance()->getGDMLFileName() << " already exists!!!"
             << G4endl;
      G4cout << "No new gdml dump created!!!" << G4endl;
      G4cout << "****************************************************" << G4endl;
    }
    else
    {
      G4cout << "Writing: " << ConfigurationManager::getInstance()->getGDMLFileName() << G4endl;
      parser->Write(ConfigurationManager::getInstance()->getGDMLFileName(), worldPhysVol);
    }
  }
  return worldPhysVol;
}
void DetectorConstruction::ConstructSDandField()
{
  G4SDManager* SDman             = G4SDManager::GetSDMpointer();
  const G4GDMLAuxMapType* auxmap = parser->GetAuxMap();
  if(verbose)
  {
    G4cout << "Found " << auxmap->size() << " volume(s) with auxiliary information." << G4endl
           << G4endl;
  }
  std::map<std::string, int> mapofSensedets = {
    { "PhotonDetector", 0 }, { "Target", 1 },      { "Tracker", 2 },
    { "SimTrajectory", 3 },  { "Msc", 4 },         { "lArTPC", 5 },
    { "Radiator", 6 },       { "Calorimeter", 7 }, { "DRCalorimeter", 8 }
  };
  enum SensDet
  {
    PhotonDetector,
    Target,
    Tracker,
    SimTrajectory,
    Msc,
    lArTPC,
    Radiator,
    Calorimeter,
    DRCalorimeter
  };
  for(auto const& [logVol, listType] : *auxmap)
  {
    for(auto const& auxtype : listType)
    {
      if(auxtype.type == "SensDet")
      {
        if(verbose)
        {
          G4cout << "Found sensitive Detector: " << auxtype.value << G4endl;
        }
        if(mapofSensedets.find(auxtype.value) == mapofSensedets.end())
        {
          G4cout << "Unknown type of sensitive Detector: " << auxtype.value << G4endl;
        }
        else
        {
          G4String name;
          switch(mapofSensedets[auxtype.value])
          {
            case PhotonDetector: {
              name                = logVol->GetName() + "_Photondetector";
              PhotonSD* aPhotonSD = new PhotonSD(name);
              SDman->AddNewDetector(aPhotonSD);
              logVol->SetSensitiveDetector(aPhotonSD);
              break;
            }
            case Target: {
              name                          = logVol->GetName() + "_Target";
              InteractionSD* aInteractionSD = new InteractionSD(name);
              SDman->AddNewDetector(aInteractionSD);
              logVol->SetSensitiveDetector(aInteractionSD);
              break;
            }
            case Tracker: {
              name                  = logVol->GetName() + "_Tracker";
              TrackerSD* aTrackerSD = new TrackerSD(name);
              SDman->AddNewDetector(aTrackerSD);
              logVol->SetSensitiveDetector(aTrackerSD);
              break;
            }
            case SimTrajectory: {
              name                              = logVol->GetName() + "_SimTrajectory";
              SimTrajectorySD* aSimTrajectorySD = new SimTrajectorySD(name);
              SDman->AddNewDetector(aSimTrajectorySD);
              logVol->SetSensitiveDetector(aSimTrajectorySD);
              break;
            }
            case Msc: {
              name          = logVol->GetName() + "_Msc";
              MscSD* aMscSD = new MscSD(name);
              SDman->AddNewDetector(aMscSD);
              logVol->SetSensitiveDetector(aMscSD);
              break;
            }
            case lArTPC: {
              name                = logVol->GetName() + "_lArTPC";
              lArTPCSD* alArTPCSD = new lArTPCSD(name);
              SDman->AddNewDetector(alArTPCSD);
              logVol->SetSensitiveDetector(alArTPCSD);
              break;
            }
            case Radiator: {
              name                    = logVol->GetName() + "_Radiator";
              RadiatorSD* aRadiatorSD = new RadiatorSD(name);
              SDman->AddNewDetector(aRadiatorSD);
              logVol->SetSensitiveDetector(aRadiatorSD);
              break;
            }
            case Calorimeter: {
              name                          = logVol->GetName() + "_Calorimeter";
              CalorimeterSD* aCalorimeterSD = new CalorimeterSD(name);
              SDman->AddNewDetector(aCalorimeterSD);
              logVol->SetSensitiveDetector(aCalorimeterSD);
              break;
            }
            case DRCalorimeter: {
              name                              = logVol->GetName() + "_DRCalorimeter";
              DRCalorimeterSD* aDRCalorimeterSD = new DRCalorimeterSD(name);
              SDman->AddNewDetector(aDRCalorimeterSD);
              logVol->SetSensitiveDetector(aDRCalorimeterSD);
              break;
            }
          }
          if(verbose)
          {
            G4cout << "Attaching sensitive Detector: " << auxtype.value
                   << " to Volume:  " << logVol->GetName() << G4endl;
          }
        }
      }
      else if(auxtype.type == "Solid")
      {
        if(auxtype.value == "True")
        {
          G4VisAttributes* visibility = new G4VisAttributes();
          visibility->SetForceSolid(true);
          G4VisAttributes* visatt = new G4VisAttributes(logVol->GetVisAttributes()->GetColour());
          visatt->SetVisibility(true);
          visatt->SetForceSolid(true);
          visatt->SetForceAuxEdgeVisible(true);
          logVol->SetVisAttributes(visatt);
        }
      }
    }
  }
}
void DetectorConstruction::ReadGDML()
{
  fReader = new ColorReader;
  parser  = new G4GDMLParser(fReader);
  parser->Read(gdmlFile, false);
  G4VPhysicalVolume* World = parser->GetWorldVolume();
  //----- GDML parser makes world invisible, this is a hack to make it
  // visible again...
  G4LogicalVolume* pWorldLogical = World->GetLogicalVolume();
  pWorldLogical->SetVisAttributes(0);
  if(verbose)
  {
    G4cout << "Found World:  " << World->GetName() << G4endl;
    G4cout << "World LV:  " << World->GetLogicalVolume()->GetName() << G4endl;
    G4LogicalVolumeStore* pLVStore = G4LogicalVolumeStore::GetInstance();
    G4cout << "Found " << pLVStore->size() << " logical volumes." << G4endl << G4endl;
    G4PhysicalVolumeStore* pPVStore = G4PhysicalVolumeStore::GetInstance();
    G4cout << "Found " << pPVStore->size() << " physical volumes." << G4endl << G4endl;
  }
}
void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}
