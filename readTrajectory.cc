
// Include files
#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TKey.h"
//
#include <iostream>
#include <map>
#include <string>
#include "include/SimStep.hh"
#include "include/SimTrajectory.hh"

int main(int argc, char** argv) {
    // initialize ROOT
    TSystem ts;
    gSystem->Load("liblArTestClassesDict");
    if (argc < 2) std::cout << "Missing name of the file to read!" << std::endl;

    TFile fo(argv[1]);
    std::cout << " reading: "<<argv[1]<<std::endl;
    std::map<int, SimTrajectory*>* tmap;
    std::vector<SimStep*>* trajectory;
    //fo.GetListOfKeys()->Print();

    TIter next(fo.GetListOfKeys());
    TKey *key;
    //double tot_en;
    while ((key = (TKey*) next())) {
        //        fo.GetObject(key->GetName(), hits);

        //tot_en = 0;
        //      std::std::cout << "Collection: " << key->GetName() << std::endl;
        // std::std::cout << "Number of hits: " << hits->size() << std::endl;
        // for (size_t i=0;i!=hits->size();i++)
        // {
        //   (*hits)[i]->Print();
        // }   

        fo.GetObject(key->GetName(), tmap);

        std::cout << "Collection: " << key->GetName() << std::endl;
        key->Print();

        const char* Collectionname = key->GetName();
        std::string s(Collectionname);
        std::string s2("Simtrajectory");
        std::size_t found = s.find(s2);
        //s.compare(0, 13, s)
        if (found == 0) {
            std::cout << "Number of trajectories: " << tmap->size() << std::endl;
            for (auto itr = tmap->begin(); itr != tmap->end(); itr++) {
                std::cout << "ID:  " << itr->first << "  steps:  " << (itr->second)->GetTrajectory()->size() << '\n';
                trajectory = (itr->second)->GetTrajectory();
                for (auto vitr = trajectory->begin(); vitr != trajectory->end(); vitr++) {
                    SimStep* st = *vitr;
                    std::cout << " Edep:  " << st->GetEdep()
                            << "  Len:  " << st->GetLen()
                            << "  Time: " << st->GetT()
                            << "  X:    " << st->GetX()
                            << "  Y:    " << st->GetY()
                            << "  Z:    " << st->GetZ()
                            << std::endl;
                }
            }
        }



    }
}

