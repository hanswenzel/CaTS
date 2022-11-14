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
///////////////////////////////////////////////////////////////////////////////
// File: CCalMagneticField.hh
// Description: A class for control of the Magnetic Field of the detector.
//              The field is assumed to be uniform.
///////////////////////////////////////////////////////////////////////////////
#ifndef EMPHATICMagneticField_H
#define EMPHATICMagneticField_H

#include "G4UniformMagField.hh"
#include "G4ThreeVector.hh"
#include <map>
#include <vector>

class G4FieldManager;

class EMPHATICMagneticField : public G4MagneticField
{
 public:
  EMPHATICMagneticField(const G4String& name);
  ~EMPHATICMagneticField();

  // Access functions
  void MagneticField(const double Point[3], double Bfield[3]) const;
  CLHEP::Hep3Vector MagneticField(const CLHEP::Hep3Vector Point) const;
  virtual void GetFieldValue(const double Point[3], double* Bfield) const;
  // G4double GetConstantFieldvalue() const {return fval;}

 protected:
  // Find the global Field Manager
  G4FieldManager* GetGlobalFieldManager();

 private:
  std::map<int, std::map<int, std::map<int, std::vector<double>>>> field;
  double step;
  double start[3];
  G4bool fVerbosity;
};

#endif
