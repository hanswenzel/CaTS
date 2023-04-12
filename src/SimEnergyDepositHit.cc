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
//   February 9th, 2022 : first implementation
//
// ********************************************************************
//
/// \file SimEnergyDepositHit.cc
/// \brief Implementation of the CaTS::SimEnergyDepositHit class
// Geant4 headers
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "SimEnergyDepositHit.hh"
G4ThreadLocal G4Allocator<SimEnergyDepositHit>* SimEnergyDepositHitAllocator = nullptr;
SimEnergyDepositHit::SimEnergyDepositHit()
  : G4VHit()
{}

SimEnergyDepositHit::SimEnergyDepositHit(unsigned int znph, unsigned int znelec, unsigned int ztid,
                                   float zx, float zy, float zz, float zxe, float zye, float zze,
                                   double zt, double zte, float zedep)
  : nph(znph)
  , nelec(znelec)
  , tid(ztid)
  , x(zx)
  , y(zy)
  , z(zz)
  , xe(zxe)
  , ye(zye)
  , ze(zze)
  , t(zt)
  , te(zte)
  , edep(zedep)
{}

SimEnergyDepositHit::SimEnergyDepositHit(const SimEnergyDepositHit& right)
{
  nph   = right.nph;
  nelec = right.nelec;
  tid   = right.tid;
  x     = right.x;
  y     = right.y;
  z     = right.z;
  xe    = right.xe;
  ye    = right.ye;
  ze    = right.ze;
  t     = right.t;
  te    = right.te;
  edep  = right.edep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const SimEnergyDepositHit& SimEnergyDepositHit::operator=(const SimEnergyDepositHit& right)
{
  nph   = right.nph;
  nelec = right.nelec;
  tid   = right.tid;
  x     = right.x;
  y     = right.y;
  z     = right.z;
  xe    = right.xe;
  ye    = right.ye;
  ze    = right.ze;
  t     = right.t;
  te    = right.te;
  edep  = right.edep;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool SimEnergyDepositHit::operator==(const SimEnergyDepositHit& right) const
{
  return (this == &right) ? true : false;
}
void SimEnergyDepositHit::Draw() {}