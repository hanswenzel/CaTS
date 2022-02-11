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
/// \file SimTrajectory.hh
/// \brief Definition of the CaTS::SimTrajectory class

#pragma once
#include <vector>
#include "globals.hh"
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
class SimStep;

class SimTrajectory : public G4VHit
{
 public:
  SimTrajectory();
  SimTrajectory(G4int id, G4int pdg, G4int parent);
  ~SimTrajectory();
  SimTrajectory(const SimTrajectory&);
  const SimTrajectory& operator=(const SimTrajectory&);
  G4bool operator==(const SimTrajectory&) const;
  inline void* operator new(size_t);
  inline void operator delete(void*);
  void Draw() final;
  G4int getTrackID() const { return fTrackID; }
  void setTrackID(const G4int& fTrackID_) { fTrackID = fTrackID_; }

  G4int getPDGcode() const { return fPDGcode; }
  void setPDGcode(const G4int& fPDGcode_) { fPDGcode = fPDGcode_; }

  G4int getParentID() const { return fParentID; }
  void setParentID(const G4int& fParentID_) { fParentID = fParentID_; }

  std::vector<SimStep*>* getTrajectory() const { return fTrajectory; }
  void setTrajectory(std::vector<SimStep*>* fTrajectory_) { fTrajectory = fTrajectory_; }

  std::vector<G4int>* getDaughters() const { return fDaughters; }
  void setDaughters(std::vector<G4int>* fDaughters_) { fDaughters = fDaughters_; }

  void AddSimStep(SimStep*);

 private:
  G4int fTrackID{ 0 };
  G4int fPDGcode{ 0 };
  G4int fParentID{ 0 };
  std::vector<G4int>* fDaughters{ nullptr };
  std::vector<SimStep*>* fTrajectory{ nullptr };
};
using SimTrajectoryCollection = G4THitsCollection<SimTrajectory>;
extern G4ThreadLocal G4Allocator<SimTrajectory>* SimTrajectoryAllocator;

inline void* SimTrajectory::operator new(size_t)
{
  if(!SimTrajectoryAllocator)
  {
    SimTrajectoryAllocator = new G4Allocator<SimTrajectory>;
  }
  return (void*) SimTrajectoryAllocator->MallocSingle();
}

inline void SimTrajectory::operator delete(void* aHit)
{
  SimTrajectoryAllocator->FreeSingle((SimTrajectory*) aHit);
}