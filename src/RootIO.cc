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
/// \file RootIO.cc
/// \brief Implementation of the CaTS::RootIO class

#ifdef WITH_ROOT
#include "RootIO.hh"
#include <G4String.hh>
#include <G4ios.hh>
#include <string>
#include "ConfigurationManager.hh"
#include "Event.hh"
#include "TBranch.h"
#include "TFile.h"
#include "TObject.h"
#include "TSystem.h"
#include "TTree.h"
#include "TROOT.h"
#include "ROOT/TBufferMerger.hxx"
//-----------------------
#include <random>
#include <thread>
//-----------------------
static RootIO* instance = 0;
//using ROOT::TBufferMerger;

RootIO::RootIO()
{
  TSystem ts;
  gSystem->Load("libCaTSClassesDict");
  G4String FileName = ConfigurationManager::getInstance()->getfname();
  G4cout << "Opening File: " << FileName << G4endl;
  ROOT::EnableThreadSafety();
  merger = new ROOT::TBufferMerger(std::unique_ptr<TMemFile>(new TMemFile(FileName.c_str(), "RECREATE")));
  fFile = merger->GetFile();
  TTree::SetMaxTreeSize(1000 * Long64_t(2000000000));
  // Create a ROOT Tree and one superbranch
  ftree = new TTree("Events", "ROOT tree containing Hit collections");
  ftree->ResetBit(kMustCleanup);
  int flush = 32;
  ftree->SetAutoFlush(-(flush) * 1024 * 1024); // Flush at exceeding 32MB
  Int_t branchStyle = 1;
  TTree::SetBranchStyle(branchStyle);
}

RootIO::~RootIO() {
  G4cout << "RootIO destructor !!!!!!"<< G4endl;
  // delete merger;
}

RootIO* RootIO::GetInstance()
{
  if(instance == 0)
  {
    instance = new RootIO();
  }
  return instance;
}

void RootIO::Write(Event* fevent)
{
  G4cout << "writing Event: " << fevent->GetEventNumber() << G4endl;
  if(ConfigurationManager::getInstance()->isEnable_verbose())
    G4cout << "writing Event: " << fevent->GetEventNumber() << G4endl;
  if((fevent->GetEventNumber()) % 1000 == 0)
    G4cout << "writing Event: " << fevent->GetEventNumber() << G4endl;
  if(!fevtinitialized)
  {
    Int_t bufsize = 64000;
    fevtbranch    = ftree->Branch("event.", &fevent, bufsize, 0);
    fevtbranch->SetAutoDelete(kFALSE);
    fevtinitialized = true;
  }
  fnb += ftree->Fill();
  fFile->Write("", TObject::kOverwrite);
}
void RootIO::Close() {
  G4cout << "closing !!!!!"<<G4endl;
}
#endif
