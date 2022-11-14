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
// File: CCalMagneticField.cc
// Description: User Field class implementation.
///////////////////////////////////////////////////////////////////////////////
#include <fstream>

#include "EMPHATICMagneticField.hh"

#include "G4SystemOfUnits.hh"
#include "G4FieldManager.hh"

#include "TFile.h"
#include "TTree.h"
#include "math.h"
// #define ddebug
#define debug

// Constructor and destructor:

EMPHATICMagneticField::EMPHATICMagneticField(const G4String& filename)
  : start{ -50, -50, -50 }
  , step(0)
{
#ifdef debug
  fVerbosity = true;
#else
  fVerbosity = false;
#endif

  /*for(int i = 0; i < 250; i++){
    for(int j = 0; j < 250; j++){
      for(int k = 0; k < 250; k++){
        std::vector<double> temp(3, 0);
        field[i][j][k] = temp;
      }
    }
  }*/

  TFile mfFile(filename.c_str(), "READ");
  G4cout << " ==> Opening file " << filename << " to read magnetic field..." << G4endl;

  if(mfFile.IsZombie())
  {
    G4cerr << "Error opening file" << G4endl;
    exit(-1);
  }
  TTree* tree = (TTree*) mfFile.Get("magField");

  double x;
  double y;
  double z;
  double Bx;
  double By;
  double Bz;
  tree->SetBranchAddress("x", &x);
  tree->SetBranchAddress("y", &y);
  tree->SetBranchAddress("z", &z);
  tree->SetBranchAddress("Bx", &Bx);
  tree->SetBranchAddress("By", &By);
  tree->SetBranchAddress("Bz", &Bz);
  Bx           = Bx * 100000;
  By           = By * 100000;
  Bz           = Bz * 100000;
  int nEntries = tree->GetEntries();

  tree->GetEntry(0);
  double xVal = x;
  double yVal = y;
  double zVal = z;

  // step = 0;
  tree->GetEntry(1);
  if(abs(xVal - x) > step)
    step = abs(xVal - x);
  else if(abs(yVal - y) > step)
    step = abs(yVal - y);
  else
    step = abs(zVal - z);

  // start = {-60, -60, -60};

  for(int i = 0; i < nEntries; i++)
  {
    tree->GetEntry(i);
    int indX   = (int) (x - start[0]) / step;
    int indY   = (int) (y - start[1]) / step;
    int indZ   = (int) (z - start[2]) / step;
    fVerbosity = true;
    if(fVerbosity)
    {
      G4cout << "(x, y, z) = (" << x << ", " << y << ", " << z << ") cm,    (ix, iy, iz) = ("
             << indX << ", " << indY << ", " << indZ << "),    (Bx, By, Bz) = (" << Bx << ", " << By
             << ", " << Bz << ") kG" << G4endl;
    }
    std::vector<double> temp;
    temp.push_back(Bx);
    temp.push_back(By);
    temp.push_back(Bz);

    field[indX][indY][indZ] = temp;
  }

  ///////////////////////////////////////////////////////////////
  // Close the file
  G4cout << " ==> Closing file " << filename << G4endl;
  mfFile.Close();
}

EMPHATICMagneticField::~EMPHATICMagneticField() {}

// Member functions

void EMPHATICMagneticField::MagneticField(const double x[3], double B[3]) const
{
  double indX = (x[0] / 10 - start[0]) / step;
  double indY = (x[1] / 10 - start[1]) / step;
  double indZ = (x[2] / 10 - start[2]) / step;

  int ix[2] = { int(floor(indX)), int(ceil(indX)) };
  int iy[2] = { int(floor(indY)), int(ceil(indY)) };
  int iz[2] = { int(floor(indZ)), int(ceil(indZ)) };

  bool skip = false;

  if(field.find(ix[0]) == field.end())
    skip = true;
  else if(field.find(ix[1]) == field.end())
    skip = true;
  else
  {
    if(field.at(ix[0]).find(iy[0]) == field.at(ix[0]).end())
      skip = true;
    else if(field.at(ix[0]).find(iy[1]) == field.at(ix[0]).end())
      skip = true;
    else if(field.at(ix[1]).find(iy[0]) == field.at(ix[1]).end())
      skip = true;
    else if(field.at(ix[1]).find(iy[1]) == field.at(ix[1]).end())
      skip = true;
    else
    {
      if(field.at(ix[0]).at(iy[0]).find(iz[0]) == field.at(ix[0]).at(iy[0]).end())
        skip = true;
      else if(field.at(ix[0]).at(iy[0]).find(iz[1]) == field.at(ix[0]).at(iy[0]).end())
        skip = true;
      else if(field.at(ix[0]).at(iy[1]).find(iz[0]) == field.at(ix[0]).at(iy[1]).end())
        skip = true;
      else if(field.at(ix[0]).at(iy[1]).find(iz[1]) == field.at(ix[0]).at(iy[1]).end())
        skip = true;
      else if(field.at(ix[1]).at(iy[0]).find(iz[0]) == field.at(ix[1]).at(iy[0]).end())
        skip = true;
      else if(field.at(ix[1]).at(iy[0]).find(iz[1]) == field.at(ix[1]).at(iy[0]).end())
        skip = true;
      else if(field.at(ix[1]).at(iy[1]).find(iz[0]) == field.at(ix[1]).at(iy[1]).end())
        skip = true;
      else if(field.at(ix[1]).at(iy[1]).find(iz[1]) == field.at(ix[1]).at(iy[1]).end())
        skip = true;
    }
  }

  if(skip)
  {
    B[0] = 0;
    B[1] = 0;
    B[2] = 0;
    return;
  }

  double sumx = 0;
  double sumy = 0;
  double sumz = 0;
  double norm = 0;

  for(int i = 0; i < 2; i++)
  {
    for(int j = 0; j < 2; j++)
    {
      for(int k = 0; k < 2; k++)
      {
        double dist = sqrt((indX - ix[i]) * (indX - ix[i]) + (indY - iy[j]) * (indY - iy[j]) +
                           (indZ - iz[k]) * (indZ - iz[k]));
        sumx += field.at(ix[i]).at(iy[j]).at(iz[k]).at(0) * dist;
        sumy += field.at(ix[i]).at(iy[j]).at(iz[k]).at(1) * dist;
        sumz += field.at(ix[i]).at(iy[j]).at(iz[k]).at(2) * dist;
        norm += dist;
      }
    }
  }

  B[0] = (sumx / norm) * kilogauss;
  B[1] = (sumy / norm) * kilogauss;
  B[2] = (sumz / norm) * kilogauss;

  if(fVerbosity)
  {
    G4cout << "(x, y, z) = (" << x[0] << ", " << x[1] << ", " << x[2] << ") cm,    (Bx, By, Bz) = ("
           << B[0] << ", " << B[1] << ", " << B[2] << ") kG" << G4endl;
  }
}

CLHEP::Hep3Vector EMPHATICMagneticField::MagneticField(const CLHEP::Hep3Vector point) const
{
  G4double x[3], B[3];
  CLHEP::Hep3Vector v;

  x[0] = point.x();
  x[1] = point.y();
  x[2] = point.z();
  EMPHATICMagneticField::MagneticField(x, B);
  v.setX(B[0]);
  v.setY(B[1]);
  v.setZ(B[2]);
  return v;
}

void EMPHATICMagneticField::GetFieldValue(const double x[3], double* B) const
{
  EMPHATICMagneticField::MagneticField(x, B);
}
