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
/// \file SimTrajectory.cc
/// \brief Implementation of the CaTS::SimTrajectory class

#include "SimTrajectory.hh"
#include "SimStep.hh"
G4ThreadLocal G4Allocator<SimTrajectory>* SimTrajectoryAllocator = nullptr;
SimTrajectory::SimTrajectory()
  : G4VHit()
  , fTrackID(0)
  , fTrajectory(0)
{
  fTrajectory = new std::vector<SimStep*>();
}

SimTrajectory::SimTrajectory(G4int id)
{
  fTrackID    = id;
  fTrajectory = new std::vector<SimStep*>();
}

SimTrajectory::SimTrajectory(const SimTrajectory& orig) {}

SimTrajectory::~SimTrajectory()
{
  for(auto step = fTrajectory->begin(); step != fTrajectory->end(); ++step)
  {
    delete *step;
  }
  delete fTrajectory;
}
void SimTrajectory::Draw() {}