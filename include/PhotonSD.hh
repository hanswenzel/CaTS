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
//   October 18th, 2021 : first implementation
//
// ********************************************************************
//
/// \file PhotonSD.hh
/// \brief Definition of the CaTS::PhotonSD class

#pragma once

#include "G4VSensitiveDetector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "PhotonHit.hh"
#include <G4String.hh>
#include <G4Types.hh>
class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;
class PhotonSD : public G4VSensitiveDetector
{
 public:
  PhotonSD(G4String);
  ~PhotonSD() = default;
  void Initialize(G4HCofThisEvent*) final;
  G4bool ProcessHits(G4Step*, G4TouchableHistory* ROhist) final;
#ifdef WITH_G4CXOPTICKS
  void AddOpticksHits();
#endif
  void EndOfEvent(G4HCofThisEvent* hitCollection) final;

 private:
  PhotonHitsCollection* fPhotonHitsCollection{ 0 };
  G4int fHCID{ 0 };
  static constexpr G4double hc = (h_Planck * c_light) / (CLHEP::eV * CLHEP::nm);
  inline G4double etolambda(G4double E)
  {
    // input photon energy in eV
    // return wavelength in nm:
    // lambda = h c/e
    //
    return hc / E;
  }
  inline G4int DetectorID(G4VPhysicalVolume* pV) { return pV->GetCopyNo(); }
  // #ifdef WITH_G4OPTICKS
  inline G4int DetectorID(G4ThreeVector pos)
  {
    const G4double rmax       = 31.;
    const G4double rmaxsquare = rmax * rmax;
    const G4double zpos[5]    = { -950., 950., 0.0, -950., 950. };
    const G4double ypos[5]    = { 450., 450., 0.0, -450., -450. };
    if((pos.y() - ypos[0]) * (pos.y() - ypos[0]) + (pos.z() - zpos[0]) * (pos.z() - zpos[0]) <
       rmaxsquare)
    {
      return 0;
    }
    else if((pos.y() - ypos[1]) * (pos.y() - ypos[1]) + (pos.z() - zpos[1]) * (pos.z() - zpos[1]) <
            rmaxsquare)
    {
      return 1;
    }
    else if((pos.y() - ypos[2]) * (pos.y() - ypos[2]) + (pos.z() - zpos[2]) * (pos.z() - zpos[2]) <
            rmaxsquare)
    {
      return 2;
    }
    else if((pos.y() - ypos[3]) * (pos.y() - ypos[3]) + (pos.z() - zpos[3]) * (pos.z() - zpos[3]) <
            rmaxsquare)
    {
      return 3;
    }
    else if((pos.y() - ypos[4]) * (pos.y() - ypos[4]) + (pos.z() - zpos[4]) * (pos.z() - zpos[4]) <
            rmaxsquare)
    {
      return 4;
    }
    else
    {
      G4cout << "x: " << pos.x() << "  y: " << pos.y() << "z: " << pos.z() << G4endl;
      return -1;
    }
  }
  // #endif
};
