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
/// \file SimStep.hh
/// \brief Definition of the CaTS::SimStep class

#pragma once

class SimStep
{
 public:
  const SimStep& operator=(const SimStep&);
  bool operator==(const SimStep&) const;
  SimStep();
  SimStep(float xx, float yy, float zz, float ll, float tt, float ed);
  SimStep(const SimStep& orig);
  ~SimStep() = default;
  float getX() const { return x; }
  void setX(float x_) { x = x_; }
  float getY() const { return y; }
  void setY(float y_) { y = y_; }
  float getZ() const { return z; }
  void setZ(float z_) { z = z_; }
  float getLen() const { return len; }
  void setLen(float len_) { len = len_; }
  float getT() const { return t; }
  void setT(float t_) { t = t_; }
  float getEdep() const { return edep; }
  void setEdep(float edep_) { edep = edep_; }

 private:
  float x;
  float y;
  float z;
  float len;
  float t;
  float edep;
};