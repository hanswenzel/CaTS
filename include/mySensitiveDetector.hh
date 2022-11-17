//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: myTrackerSD.hh,v 1.3.4.1 2001/06/28 19:07:29 gunter Exp $
// GEANT4 tag $Name:  $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#ifndef mySensitiveDetector_h
#define mySensitiveDetector_h 1

#include "G4VSensitiveDetector.hh"
#include "myHit.hh"

#include "G4VSensitiveDetector.hh"
#include "G4ScintillationTrackInformation.hh"
#include <G4MaterialPropertyVector.hh>
#include <G4String.hh>
#include <G4Types.hh>

class G4Step;
class G4HCofThisEvent;
class myDetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

class mySensitiveDetector : public G4VSensitiveDetector
{

public:
  mySensitiveDetector(G4String name, myDetectorConstruction *);
  ~mySensitiveDetector();
  
  virtual void Initialize(G4HCofThisEvent*);
  virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  //#ifdef WITH_G4CXOPTICKS
   //void ProcessOpticksHits(G4Step*, G4TouchableHistory*);
  //#endif
  virtual void EndOfEvent(G4HCofThisEvent*);

private:

  G4int materialIndex;
  const G4Material* aMaterial;
  G4MaterialPropertiesTable* aMaterialPropertiesTable;

  myDetectorConstruction *Detector;
  G4String SDName;
  G4int numHits=0;
  myHitsCollection* mySDHitsCollection;
  G4int HCID;

};

#endif

