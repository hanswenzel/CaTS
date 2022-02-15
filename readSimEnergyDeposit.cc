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
//---------------------------------------------------------------------
//*            |\___/|                                                *
//*            )     (                                                *
//*           =\     /=                                               *
//*             )===(                                                 *
//*            /     \     CaTS: Calorimeter and Tracker Simulation   *
//*            |     |     is a flexible and extend-able framework    *
//*           /       \    for the simulation of various detector     *
//*	          \       /    systems                                    *
//*            \__  _/     https://github.com/hanswenzel/CaTS         *
//*	             ( (                                                  *
//*	              ) )                                                 *
//*              (_(                                                  *
//* CaTS also serves as an example that demonstrates how to use       *
//* opticks from within Geant4 for the creation and propagation of    *
//* optical photons.                                                  *
//* see https://bitbucket.org/simoncblyth/opticks.git).               *
//* Ascii Art by Joan Stark: https://www.asciiworld.com/-Cats-2-.html *
//---------------------------------------------------------------------
//
/// \file readSimEnergyDeposit.cc
/// \brief example how to read the  CaTS::SimEnergyDeposit
//
// Root headers
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include "TTree.h"
// Project headers
#include "Event.hh"
#include "SimEnergyDeposit.hh"

int main(int argc, char** argv)
{
  // initialize ROOT
  TSystem ts;
  gSystem->Load("libCaTSClassesDict");
  if(argc < 4)
  {
    G4cout << "Program requires 3 arguments: name of input file, name of "
              "output file, Volume that sensitive detector is attached to"
           << G4endl;
    exit(1);
  }
  TFile* outfile = new TFile(argv[2], "RECREATE");
  outfile->cd();
  TH1F* energy = new TH1F("energy", "total energy", 100, 0., 2000.);

  TFile fo(argv[1]);
  fo.GetListOfKeys()->Print();
  Event* event = new Event();
  TTree* Tevt  = (TTree*) fo.Get("Events");
  Tevt->SetBranchAddress("event.", &event);
  TBranch* fevtbranch = Tevt->GetBranch("event.");
  Int_t nevent        = fevtbranch->GetEntries();
  G4cout << "Nr. of Events:  " << nevent << G4endl;
  std::string CollectionName = argv[3];
  CollectionName             = CollectionName + "_SimEnergyDeposit_HC";
  for(Int_t i = 0; i < nevent; i++)
  {
    fevtbranch->GetEntry(i);
    auto* hcmap = event->GetHCMap();
    for(const auto& ele : *hcmap)
    {
      auto hits = ele.second;
      if(ele.first.compare(CollectionName) == 0)
      {
        auto hits    = ele.second;
        G4int NbHits = hits.size();
        G4cout << "Event: " << i << "  Number of Hits:  " << NbHits << G4endl;
        // np->Fill(NbHits);
        double tote = 0.0;
        for(G4int ii = 0; ii < NbHits; ii++)
        {
          SimEnergyDeposit* simEnergyDeposit = dynamic_cast<SimEnergyDeposit*>(hits.at(ii));
          tote                               = tote + simEnergyDeposit->GetEdep();
        }
        energy->Fill(tote);
      }
    }
  }
  outfile->cd();
  outfile->Write();
}