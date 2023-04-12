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
/// \file SimEnergyDepositHit.hh
/// \brief Definition of the CaTS::SimEnergyDepositHit class

#pragma once
#include "globals.hh"
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
class SimEnergyDepositHit : public G4VHit
{
 public:
  const SimEnergyDepositHit& operator=(const SimEnergyDepositHit&);
  bool operator==(const SimEnergyDepositHit&) const;
  SimEnergyDepositHit();
  SimEnergyDepositHit(unsigned int znph, unsigned int znelec, unsigned int ztid, float zx, float zy,
                   float zz, float zxe, float zye, float zze, double zt, double zte, float zedep);
  SimEnergyDepositHit(const SimEnergyDepositHit&);
  ~SimEnergyDepositHit() = default;
  inline void* operator new(size_t);
  inline void operator delete(void*);
  void Draw() final;
  void SetEdep(float edep);
  inline float GetEdep() const { return edep; };
  void SetT(float t);
  inline float GetT() const { return t; };
  void SetZ(float z);
  inline float GetZ() const { return z; };
  void SetY(float y);
  inline float GetY() const { return y; };
  void SetX(float x);
  inline float GetX() const { return x; }
  double GetTe() const { return te; };
  float GetZe() const { return ze; };
  float GetYe() const { return ye; };
  float GetXe() const { return xe; };

 private:
  unsigned int nph{ 0 };
  unsigned int nelec{ 0 };
  unsigned int tid{ 0 };
  float x{ 0.0 };
  float y{ 0.0 };
  float z{ 0.0 };
  float xe{ 0.0 };
  float ye{ 0.0 };
  float ze{ 0.0 };
  double t{ 0.0 };
  double te{ 0.0 };
  float edep{ 0.0 };
};
using SimEnergyDepositHitCollection = G4THitsCollection<SimEnergyDepositHit>;
extern G4ThreadLocal G4Allocator<SimEnergyDepositHit>* SimEnergyDepositHitAllocator;

inline void* SimEnergyDepositHit::operator new(size_t)
{
  if(!SimEnergyDepositHitAllocator)
  {
    SimEnergyDepositHitAllocator = new G4Allocator<SimEnergyDepositHit>;
  }
  return (void*) SimEnergyDepositHitAllocator->MallocSingle();
}

inline void SimEnergyDepositHit::operator delete(void* aHit)
{
  SimEnergyDepositHitAllocator->FreeSingle((SimEnergyDepositHit*) aHit);
}