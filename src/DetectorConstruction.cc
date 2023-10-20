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
#include "G4AutoDelete.hh"
#include "G4GDMLParser.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include <G4Region.hh>
#include <G4RegionStore.hh>
#include <G4ProductionCuts.hh>
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4UserLimits.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

// project headers
#include "F01FieldSetup.hh"
#include "EMPHATICMagneticField.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4PropagatorInField.hh"
#include "G4TransportationManager.hh"
#include "G4ClassicalRK4.hh"

#include "CalorimeterSD.hh"
#include "ColorReader.hh"
#include "ConfigurationManager.hh"
#include "DRCalorimeterSD.hh"
#include "DetectorConstruction.hh"
#include "InteractionSD.hh"
#include "MscSD.hh"
#include "PhotonSD.hh"
#ifdef WITH_G4OPTICKS
#  include "RadiatorSD.hh"
#endif
#ifdef WITH_G4CXOPTICKS
// #  include "OPTICKS_LOG.hh"
#  include "G4CXOpticks.hh"
// #  include "OpticksMode.hh"
#  include <cuda_runtime.h>
#  include "SEventConfig.hh"
#endif

#include "TrackerSD.hh"
#include "lArTPCSD.hh"
#include "SimTrajectorySD.hh"
#include "SimEnergyDepositSD.hh"
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

  // if(fMagFieldMap != "")
  //{
  /*
    G4String fMagFieldMap           = "magField.root";
    EMPHATICMagneticField* magField = new EMPHATICMagneticField(fMagFieldMap);
    G4FieldManager* fieldMgr =
    G4TransportationManager::GetTransportationManager()->GetFieldManager();
    fieldMgr->SetDetectorField(magField);
    G4Mag_UsualEqRhs* fEquation = new G4Mag_UsualEqRhs(magField);

    G4MagIntegratorStepper* pStepper = new G4ClassicalRK4(fEquation);

    G4ChordFinder* pChordFinder = new G4ChordFinder(magField, 1.e-1 * mm, pStepper);
    pChordFinder->SetDeltaChord(1.0e-3 * mm);
    fieldMgr->SetChordFinder(pChordFinder);
    fieldMgr->SetDeltaOneStep(1.0e-3 * mm);
    fieldMgr->SetDeltaIntersection(1.0e-4 * mm);
    G4PropagatorInField* fieldPropagator =
      G4TransportationManager::GetTransportationManager()->GetPropagatorInField();
    fieldPropagator->SetMinimumEpsilonStep(1.e-5 * mm);
    fieldPropagator->SetMaximumEpsilonStep(1.e-2 * mm);
    */
  //}
  /*
    fUseFSALstepper = true;
    F01FieldSetup* fieldSetup =
      new F01FieldSetup(G4ThreeVector(0.0, 3.3 * tesla, 0.0), fUseFSALstepper);
    G4AutoDelete::Register(fieldSetup);  // Kernel will delete the F01FieldSetup
    fEmFieldSetup.Put(fieldSetup);
  */

  ReadGDML();
  const G4GDMLAuxMapType* auxmap = parser->GetAuxMap();
  if(verbose)
  {
    G4cout << "Found " << auxmap->size() << " volume(s) with auxiliary information." << G4endl
           << G4endl;
  }
  for(auto const& [logVol, listType] : *auxmap)
  {
    if(verbose)
    {
      G4cout << "Volume " << logVol->GetName()
             << " has the following list of auxiliary information: " << G4endl;
    }
    for(auto const& auxtype : listType)
    {
      G4double value             = atof(auxtype.value);
      G4double val_unit          = 1;  //--no unit
      G4String provided_category = "NONE";
      G4cout << auxtype.type << G4endl;
      if((auxtype.unit) && (auxtype.unit != ""))
      {  // -- if provided and non-NULL
        val_unit          = G4UnitDefinition::GetValueOf(auxtype.unit);
        provided_category = G4UnitDefinition::GetCategory(auxtype.unit);
        if(verbose)
        {
          G4cout << " Unit parsed = " << auxtype.unit
                 << " from unit category: " << provided_category.c_str();
        }
        value *=
          val_unit;  //-- Now do something with the value, making sure that the unit is appropriate
      }
      if(auxtype.type == "StepLimit")
      {
        //-- check that steplimit has valid length unit category
        G4String steplimit_category = "Length";
        if(provided_category == steplimit_category)
        {
          G4cout << "Valid steplimit unit category obtained: " << provided_category.c_str()
                 << G4endl;
          // -- convert length to mm
          value                    = (value / CLHEP::mm) * CLHEP::mm;
          G4UserLimits* fStepLimit = new G4UserLimits(value);
          G4AutoDelete::Register(fStepLimit);
          logVol->SetUserLimits(fStepLimit);
          G4cout << "StepLimit for log Volume: " << logVol->GetName() << "  " << value << "  "
                 << value / CLHEP::cm << " cm" << G4endl;
        }
        else if(provided_category == "NONE")
        {  //--no unit category provided, use the default CLHEP::mm
          G4cout << "StepLimit in geometry file does not have a unit!"
                 << " Defaulting to mm..." << G4endl;
          value *= CLHEP::mm;
          G4UserLimits* fStepLimit = new G4UserLimits(value);
          G4AutoDelete::Register(fStepLimit);
          logVol->SetUserLimits(fStepLimit);
          G4cout << "StepLimit for log Volume: " << logVol->GetName() << "  " << value << "  "
                 << value / CLHEP::cm << " cm" << G4endl;
        }
        else
        {  //--wrong unit category provided
          G4cout << "StepLimit does not have a valid length unit!" << G4endl;
          G4cout << "Category of unit provided = " << provided_category << G4endl;
          exit(EXIT_FAILURE);
        }
      }
      if(auxtype.type == "Solid")
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
  /*
//
// dump material properties:
//

const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
G4int nMaterials                     = G4Material::GetNumberOfMaterials();
for(G4int m = 0; m < nMaterials; m++)
{
  const G4Material* aMaterial = (*materialTable)[m];
  G4cout << "Material Name:  " << aMaterial->GetName() << G4endl;
  G4MaterialPropertiesTable* aMaterialPropertiesTable =
aMaterial->GetMaterialPropertiesTable(); if(aMaterialPropertiesTable != nullptr)
    aMaterialPropertiesTable->DumpTable();
}
*/
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
  std::map<std::string, int> mapofSensedets = { { "PhotonDetector", 0 },   { "Target", 1 },
                                                { "Tracker", 2 },          { "SimTrajectory", 3 },
                                                { "SimEnergyDeposit", 4 }, { "Msc", 5 },
                                                { "lArTPC", 6 },           { "Radiator", 7 },
                                                { "Calorimeter", 8 },      { "DRCalorimeter", 9 } };
  enum SensDet
  {
    PhotonDetector,
    Target,
    Tracker,
    SimTrajectory,
    SimEnergyDeposit,
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
            case SimEnergyDeposit: {
              name                                    = logVol->GetName() + "_SimEnergyDeposit";
              SimEnergyDepositSD* aSimEnergyDepositSD = new SimEnergyDepositSD(name);
              SDman->AddNewDetector(aSimEnergyDepositSD);
              logVol->SetSensitiveDetector(aSimEnergyDepositSD);
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
#ifdef WITH_G4OPTICKS
            case Radiator: {
              name                    = logVol->GetName() + "_Radiator";
              RadiatorSD* aRadiatorSD = new RadiatorSD(name);
              SDman->AddNewDetector(aRadiatorSD);
              logVol->SetSensitiveDetector(aRadiatorSD);
              break;
            }
#endif
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
      }  // end if sensdet
    }
  }
  /*
  G4MagneticField* magField;
  magField                 = new G4UniformMagField(G4ThreeVector(0., 3.0 * tesla, 0.0));
  G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldMgr->SetDetectorField(magField);
  fieldMgr->CreateChordFinder(magField);
*/
  // if(fMagFieldMap != "")
  //{
  /*
    G4String fMagFieldMap           = "magField.root";
    EMPHATICMagneticField* magField = new EMPHATICMagneticField(fMagFieldMap);
    G4FieldManager* fieldMgr =
    G4TransportationManager::GetTransportationManager()->GetFieldManager();
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField);
    */
  //  G4Mag_UsualEqRhs* fEquation = new G4Mag_UsualEqRhs(magField);

  // G4MagIntegratorStepper* pStepper = new G4ClassicalRK4(fEquation);

  // G4ChordFinder* pChordFinder = new G4ChordFinder(magField, 1.e-1 * mm, pStepper);
  /*
  pChordFinder->SetDeltaChord(1.0e-3 * mm);
  fieldMgr->SetChordFinder(pChordFinder);
  fieldMgr->SetDeltaOneStep(1.0e-3 * mm);
  fieldMgr->SetDeltaIntersection(1.0e-4 * mm);
  G4PropagatorInField* fieldPropagator =
    G4TransportationManager::GetTransportationManager()->GetPropagatorInField();
  fieldPropagator->SetMinimumEpsilonStep(1.e-5 * mm);
  fieldPropagator->SetMaximumEpsilonStep(1.e-2 * mm);
*/
  //}
  /*
    fUseFSALstepper = true;
    F01FieldSetup* fieldSetup =
      new F01FieldSetup(G4ThreeVector(0.0, 3.3 * tesla, 0.0), fUseFSALstepper);
    G4AutoDelete::Register(fieldSetup);  // Kernel will delete the F01FieldSetup
    fEmFieldSetup.Put(fieldSetup);
    */
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
#ifdef WITH_G4CXOPTICKS
  if(ConfigurationManager::getInstance()->isEnable_opticks())
  {
    // G4CXOpticks gx;  // Simulate is the default RGMode
    // if(opticksMode != 0)
    std::cout << "DetectorConstruction setGeometry" << std::endl;
    cudaDeviceReset();
    // SEventConfig::SetMaxPhoton(100000000);
    // G4CXOpticks* g4cx =
    //if (opticksMode != 0)
      G4CXOpticks::SetGeometry(World);
    // gx.setGeometry(World);
    // SEventConfig::SetMaxPhoton(100000000);
    std::cout << SEventConfig::Desc() << std::endl;
  }
#endif
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
