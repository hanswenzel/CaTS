# CaTS
![Cats logo](images/medium.CaTSlogo.png)

# CaTS: Calorimeter and Tracker Simulation

CaTS is a flexible and extend-able framework (based on geant4 and ROOT)
for the simulation of calorimeter and tracking detectors. 
It also serves as an Example that demonstrates how to use opticks
from within Geant4 for the creation and propagation of optical photons.
https://bitbucket.org/simoncblyth/opticks.git.

The components of CaTS are:


Detector Description:      described in gdml input file (e.g. crystalcal.gdml)
(Geometry, Materials,
 optical properties,
 sensitive detector)

Input modules:                

    GPS
    Particle Gun



Physics Lists:                  choice of all Reference Physics Lists
                                          optical physics processes (Cerenkov, Rayleigh,
                                          Scintillation etc.) are added (talk to)  
                
Sensitive Detectors:        (+ corresponding Hit classes)        

    TrackerSD(Hit)  registers the step points of charged particles.
    CalorimeterSD(Hit) registering energy deposit
    DRCalorimeterSD(Hit) besides registering energy deposit counts produced Cerenkov photons
    MsCSD(Hit) used to study multiple scattering on a thin layer
    PhotonSD(Hit) sensitive detector that registers optical photons.
    lArTPCSD(Hit) sensitive detector that registers ionization and collects Gensteps (Scintillation and Cerenkov) to be processed by Opticks.


This requires the opticks environment to be set up properly see: [Instructions](Instructions.md)

 source setup_opticks.sh

To get started : 

```bash
git clone https://github.com/hanswenzel/CaTS.git
cd CaTS/

# to setup the opticks environment source the setup_opticks.sh file described in:
# https://github.com/hanswenzel/CaTS/blob/master/Instructions.md


source (path to opticks WORK_DIR)/setup_opticks.sh 
cd ../
mkdir CaTS-build
cd CaTS-build

cmake -GNinja -DCMAKE_BUILD_TYPE=Release \
  -DWITH_G4CXOPTICKS=ON \
  -DCMAKE_PREFIX_PATH="${LOCAL_BASE}/opticks/externals;${LOCAL_BASE}/opticks" \
  -DOPTICKS_PREFIX=${LOCAL_BASE}/opticks \
  -DCMAKE_MODULE_PATH=${OPTICKS_HOME}/cmake/Modules \
  -DCMAKE_INSTALL_PREFIX=../CaTS-install \
  ../CaTS

ninja install
cd ../CaTS-install/bin
time ./CaTS -g  simpleLArTPC.gdml -pl 'FTFP_BERT+OPTICAL+STEPLIMIT'  -t 1 -m time.mac >& time.log
The command line variables are 
-g 'name of gdml file defining the geometry'
-pl 'name of the Physics list and physiscs constructors'
-m 'name of the Geant4 macro'
-t 'n number of Geant4 threads'

Note! Only the -g command line variable is mandatory! If you don't specify the macro file interactive mode is assumed:
Note! For the moment one can use only 1 Geant4 thread when using G4CXOpticks'

./CaTS -g simpleLArTPC.gdml -pl 'FTFP_BERT+OPTICAL+STEPLIMIT'

```

to compile CaTS without Opticks do:

```bash
cmake  -GNinja -DCMAKE_BUILD_TYPE=Release  -DWITH_G4CXOPTICKS=OFF      -DCMAKE_MODULE_PATH="../CaTS/cmake/Modules"   -DCMAKE_INSTALL_PREFIX=../CaTS-install   ../CaTS
```
if you don't provide the -pl argument the default physics list configuration:
'FTFP_BERT+OPTICAL+STEPLIMIT'
is used


![example1](images/display.png)
to look at the hit collection and make a few histograms:

    ./readPhotonHits NewHits_point_Run0.root PhotonHistos.root Det

The 3 arguments here are: name of input file (hits), name of output file for the histograms and the logical Volume that sensitive detector (PhotonSD)is attached to.
One can then use root to look at the plots:

```bash
root histos.root
   ------------------------------------------------------------------
  | Welcome to ROOT 6.22/06                        https://root.cern |
  | (c) 1995-2020, The ROOT Team; conception: R. Brun, F. Rademakers |
  | Built for linuxx8664gcc on Dec 13 2020, 13:28:00                 |
  | From tags/v6-22-06@v6-22-06                                      |
  | Try '.help', '.demo', '.license', '.credits', '.quit'/'.q'       |
   ------------------------------------------------------------------

root [0] 
Attaching file histos.root as _file0...
(TFile *) 0x559578f76f60
root [1] TBrowser b
(TBrowser &) Name: Browser Title: ROOT Object Browser
```

![example2](images/position.png)

For comparison one might want to disable Opticks and use Geant4 to generate and propagate optical photons:
    
    time ./CaTS -g  simpleLArTPC.gdml -pl 'FTFP_BERT+OPTICAL+STEPLIMIT'  -m time_G4.mac
