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
//   October 18th, 2021 : first implementation
//
// ********************************************************************
//
/// \file MCTrackingAction.cc
/// \brief Implementation of the CaTS::MCTrackingAction class

// Geant4 headers
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "G4ClassificationOfNewTrack.hh"
#include "G4GenericMessenger.hh"
// project headers
#include "MCTrackingAction.hh"

MCTrackingAction::MCTrackingAction()
  : G4UserTrackingAction()
{
  // DefineCommands();
}

MCTrackingAction::~MCTrackingAction()
{
  // delete fMessenger;
}

void MCTrackingAction::Print()
{
  /*
  G4cout << "===================================================" << G4endl;
  G4cout << " MCTrackingAction configuration:      " << G4endl;
  G4cout << " Kill Pi0s :                        " << fkillPi0 << G4endl;
  G4cout << " Kill etas :                        " << fkilleta << G4endl;
  G4cout << " Kill Gammas from neutron Capture:  " << fkillGammafromnCapture
         << G4endl;
  G4cout << "===================================================" << G4endl;
*/
}
/*
void MCTrackingAction::DefineCommands()
{
  fMessenger = new G4GenericMessenger(this, "/CaTS/MCTrackingAction/", "select particles to kill");
auto& killPi0Cmd = fMessenger->DeclareProperty("killPi0", fkillPi0);
killPi0Cmd.SetGuidance("kill Pi0 (true/false)");
killPi0Cmd.SetStates(G4State_PreInit, G4State_Init, G4State_Idle);
killPi0Cmd.SetDefaultValue("false");
auto& killetaCmd = fMessenger->DeclareProperty("killeta", fkilleta);
killetaCmd.SetGuidance("kill eta (true/false)");
killetaCmd.SetStates(G4State_PreInit, G4State_Init, G4State_Idle);
killetaCmd.SetDefaultValue("false");
auto& killGammafromnCaptureCmd = fMessenger->DeclareProperty(
"killGammafromnCapture", fkillGammafromnCapture);
killGammafromnCaptureCmd.SetGuidance("kill GammafromnCapture (true/false)");
killGammafromnCaptureCmd.SetStates(G4State_PreInit, G4State_Init,
 G4State_Idle);
killGammafromnCaptureCmd.SetDefaultValue("false");
fMessenger->DeclareMethod("Print", &MCTrackingAction::Print)
.SetGuidance("Print MCTrackingAction configuration")
.SetStates(G4State_PreInit, G4State_Idle);

}
*/