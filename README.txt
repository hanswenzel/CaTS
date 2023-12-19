
using:
Optix 7.5 (june 2022)       latest as of Aug 1 2023 7.7
cuda 11.7 (may 2022)        latest as of Aug 1 2023 12.2
Geant4 11.1.p02 (June 2023) latest as of Aug 1 2023 11.1.p02
CLHEP 2.4.6.2               latest as of Aug 1 2023 2.4.6.4
Root 6_28_04                latest stable  as of Aug 1 2023  6_28_04


First of all make sure that all the necessary system tools and development libraries are available on the System.
For Ubuntu ( LTS 22.04) we provide the script in CaTS
https://github.com/hanswenzel/CaTS/blob/main/scripts/checkpr.sh


CLHEP:
wget https://proj-clhep.web.cern.ch/proj-clhep/dist1/clhep-2.4.6.2.tgz 
tar xzvf clhep-2.4.6.2.tgz 
cd 2.4.6.2/
mkdir CLHEP-build
cd  CLHEP-build
cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DCLHEP_BUILD_CXXSTD=-std=c++17 ../CLHEP
ninja
sudo ninja install


geant4:
wget https://geant4-data.web.cern.ch/releases/geant4-v11.1.2.tar.gz
tar xzvf geant4-v11.1.2.tar.gz
mkdir  geant4-v11.1.2-build


cd   geant4-v11.1.2-build
cmake  -GNinja -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../geant4-v11.1.2-install -DGEANT4_BUILD_BUILTIN_BACKTRACE=OFF -DGEANT4_BUILD_VERBOSE_CODE=OFF -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_SYSTEM_CLHEP=ON -DGEANT4_USE_GDML=ON -DGEANT4_USE_SYSTEM_EXPAT=ON -DGEANT4_USE_SYSTEM_ZLIB=ON  -DGEANT4_USE_QT=ON -DGEANT4_BUILD_MULTITHREADED=OFF -DGEANT4_USE_OPENGL_X11=ON ../geant4-v11.1.2
ninja
ninja install
. ../geant4-v11.1.2-install/bin/geant4.sh

Note we specifically set -DGEANT4_BUILD_MULTITHREADED=OFF since there seem to be some incompabilties opticks 

git clone --branch latest-stable https://github.com/root-project/root.git root_6_28_04_src
mkdir root_6_28_04_build
cd root_6_28_04_build
cmake -GNinja -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=../root_6_28_04-install  -Droot7=ON   -Dxrootd=OFF ../root_6_28_04_src/
cmake --build . --target install
source ../root_6_28_04-install/bin/thisroot.sh
root


cuda:
https://developer.nvidia.com/cuda-downloads
legacy versions can be found in 
https://developer.nvidia.com/cuda-toolkit-archive
e.g. 11.7
https://developer.nvidia.com/cuda-11-7-0-download-archive

cuda-samples:
git clone https://github.com/nvidia/cuda-samples
cd cuda-samples
make --ignore-errors
try
bin/x86_64/linux/release/deviceQueryDrv
After the cuda-samples are build copy them to /usr/local/cuda or create a symbolic link /usr/local/cuda/samples that
points to the installation of cuda-samples. Opticks uses some of the sample headers and therefore building opticks
will not suceed without them.

optix:
https://developer.nvidia.com/designworks/optix/download
legacy  versions can be found in. ../
https://developer.nvidia.com/designworks/optix/downloads/legacy
e.g. 7.5
https://developer.nvidia.com/optix/downloads/7.5.0/linux64-x86_64



opticks:
git clone https://bitbucket.org/simoncblyth/opticks.git

now with the tag:

git clone https://github.com/simoncblyth/opticks.git
cd opticks
git checkout v0.2.2
git fetch --all --tags



change opticks/optickscore/OpticksSwitches.h

so that:

#define WITH_SKIPAHEAD 1



mkdir -p ${WORK_DIR}/local/opticks/externals/
cd ${WORK_DIR}/local/opticks/externals/
ln -s ${OptiX_INSTALL_DIR} OptiX
cd ${WORK_DIR}
opticks-externals-install >& install_ext.log &
tail -f install_ext.log

cd ${WORK_DIR}
opticks-full  >& install_full.log &
tail -f install_full.log




For CaTS the installation should be the one below but it looks like Simon made
changes to opticks and CaTS needs to be modified for that. 



git clone https://github.com/hanswenzel/CaTS.git
cd CaTS/
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
