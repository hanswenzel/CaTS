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
/// \file SimStep.cc
/// \brief Implementation of the CaTS::SimStep class

#include "SimStep.hh"

SimStep::SimStep()
  : x(0.0)
  , y(0.0)
  , z(0.0)
  , len(0.0)
  , t(0.0)
  , edep(0.0)
{}

SimStep::SimStep(float xx = 0.0, float yy = 0.0, float zz = 0.0, float ll = 0.0, float tt = 0.0,
                 float ed = 0.0)
  : x(xx)
  , y(yy)
  , z(zz)
  , len(ll)
  , t(tt)
  , edep(ed)
{}

SimStep::SimStep(const SimStep& right)
{
  x    = right.x;
  y    = right.y;
  z    = right.z;
  len  = right.len;
  t    = right.t;
  edep = right.edep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const SimStep& SimStep::operator=(const SimStep& right)
{
  x    = right.x;
  y    = right.y;
  z    = right.z;
  len  = right.len;
  t    = right.t;
  edep = right.edep;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool SimStep::operator==(const SimStep& right) const { return (this == &right) ? true : false; }
