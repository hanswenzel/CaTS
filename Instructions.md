# Opticks configurations
## irago (Hans home machine)

<table>
  <tr>
    <th>Product </th>
    <th>Versio</th>
  </tr>
  <tr>
    <td>$${\color{blue} cuda}$$ </td>
    <td>$${\color{blue} 11.7}$$</td>
  </tr>
  <tr>
    <td>$${\color{blue} Optix }$$ </td>
    <td>7.5</td>
  </tr>
  <tr>
    <td>$${\color{blue} Geant4}$$ </td>
    <td>$${\color{blue} 11.1.p02}$$</td>
  </tr>
  <tr>
    <td>$${\color{blue} CLHEP}$$ </td>
    <td>$${\color{blue} 2.4.6.2}$$</td>
  </tr>
    <tr>
    <td>$${\color{blue} Root}$$ </td>
    <td>$${\color{blue} \ 6 \_ 28 \_ 04} $$</td>
  </tr>
    <tr>
    <td>$${\color{blue} Opticks}$$ </td>
    <td>$${\color{blue}  v0.2.7}$$</td>
  </tr>
    </tr>
    <tr>
    <td>$${\color{blue} CaTS}$$ </td>
    <td>$${\color{blue}  v2.0.4}$$</td>
  </tr>
    </tr>
    <tr>
    <td>$${\color{blue} Opticks}$$ </td>
    <td>$${\color{blue} v0.2.7}$$ </td>
  </tr>
    </tr>
 <tr>
    <td>$${\color{blue} GPU}$$ </td>
    <td>$${\color{blue} NVIDIA \ GeForce \  RTX\   2070}$$</td>
  </tr>
</table>


<p align="left"> $${\color{blue} Optix \ 7.5 }$$ </p>

$${\color{blue} cuda \ 11.7}$$ 

$${\color{blue} Geant4 \  11.1.p02}$$

$${\color{blue} CLHEP \ 2.4.6.2 }$$

$${\color{blue} Root \ 6 \_ 28 \_ 04  }$$

$${\color{blue} Opticks \  v0.2.7}$$ 

$${\color{blue} CaTS  \ v2.0.4}$$  

$${\color{blue} NVIDIA \ GeForce \  RTX\   2070.  }$$ 

$${\color{blue} NVIDIA \ driver \ version \  515.43.04  }$$ 

$${\color{blue} OS: \  Ubuntu \ 22.04.3 \ LTS  }$$ 
to check version use the following command:  
ubuntu: lsb_release -a   
cuda: nvcc --version  

# Prerequisites
First of all make sure that all the necessary system tools and development libraries are available on the System. For Ubuntu we provide the script: [checkpr.sh](scripts/checkpr.sh) that ensures the system is ready. Opticks requires Geant4, nvidia cuda and nvidia Optix among other libraries. CaTS in addition will require ROOT. If all these libraries and development headers are available on your machine skip directly to  (**Building opticks vs. existing libraries**). On a 'blank' computing system it makes sense to build CLHEP, then Geant4 and finally ROOT in that order assuring that all the necessary development libraries and headers are installed.   

# Building CLHEP
Check the release note for the required clhep version for the Geant4 release you are using. For example for Geant4 11.2 this can be found at: https://geant4.web.cern.ch/download/release-notes/notes-v11.2.0.html

CLHEP can be found at:
https://proj-clhep.web.cern.ch/proj-clhep/

to build it from scratch using cmake (used cmake version > 3.22.0) 

    cd to the directory where you want to build clhep (replace version)
    wget https://proj-clhep.web.cern.ch/proj-clhep/dist1/clhep- $${\color{blue} 2.4.5.1}$$ .tgz 
    tar xzvf clhep-2.4.5.1.tgz
    cd 2.4.5.1/
    mkdir CLHEP-build
    cd  CLHEP-build
    cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../CLHEP-install -DCLHEP_BUILD_CXXSTD=-std=c++17 ../CLHEP
    ninja
    ninja install

**Note** the default install directory is /usr/local but one needs root privileges to install it there:

    cd to the directory where you want to build clhep
    wget https://proj-clhep.web.cern.ch/proj-clhep/dist1/clhep-2.4.5.1.tgz 
    tar xzvf clhep-2.4.5.1.tgz
    cd 2.4.5.1/
    mkdir CLHEP-build
    cd  CLHEP-build
    cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DCLHEP_BUILD_CXXSTD=-std=c++17 ../CLHEP
    ninja
    sudo ninja install

# Building Geant4

Geant4 versions are available at:
https://geant4.web.cern.ch/support/download


    cd to the directory where you want to install Geant4
    wget https://geant4-data.web.cern.ch/releases/geant4-v11.0.2.tar.gz
    tar xzvf geant4-v11.0.2.tar.gz
    mkdir geant4-v11.0.2-build
    cd  geant4-v11.0.2-build
    
    cmake  -GNinja -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../geant4-v11.0.2-install -DGEANT4_BUILD_BUILTIN_BACKTRACE=OFF -DGEANT4_BUILD_VERBOSE_CODE=OFF -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_SYSTEM_CLHEP=ON -DGEANT4_USE_GDML=ON -DGEANT4_USE_SYSTEM_EXPAT=ON -DGEANT4_USE_SYSTEM_ZLIB=ON  -DGEANT4_USE_QT=ON -DGEANT4_BUILD_MULTITHREADED=ON -DGEANT4_USE_OPENGL_X11=ON ../geant4-v11.0.2 

    ninja
    ninja install
    . ../geant4-v11.0.2-install/bin/geant4.sh


check the output for any error. 




# Building ROOT 
Instructions how to build ROOT from ssource can be found at:
https://root.cern/install/build_from_source/

    # cd  to the diretory where you want to install root
    git clone --branch latest-stable https://github.com/root-project/root.git root_src
    mkdir root-build
    cd root-build
    
    cmake -GNinja -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=../root-install  -Droot7=ON   -Dxrootd=OFF ../root_src/ 
    cmake --build . --target install
    
    or if you use make instead of ninja:
    cmake  -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=../root-install  -Droot7=ON   -Dxrootd=OFF ../root_src/ 
    # speed up the make process
    new=" -j$(($(grep -c ^processor /proc/cpuinfo) - 1))" 
    case ":${MAKEFLAGS:=$new}:" in
        *:"$new":*)  ;;
        *) MAKEFLAGS="$MAKEFLAGS:$new"  ;;
    esac
    cmake --build . --target install
    
check the output for any error, install any development packages that might be necessary. 


    . ../root-install/bin/thisroot.sh
    root
    

# Installing CUDA

cuda (11.5) is available at the NVIDIA web site just follow the instruction depending on the system you are using. 

https://developer.nvidia.com/cuda-downloads

**Note** this will also install the corresponding NVIDIA graphics driver you might have to reboot.  

A good way to check that things are working properly is to build the cuda samples and execute them

    # cd to the directory where you want to build the cuda samples. E.g. the commands deviceQueryDrv and deviceQuery provide useful information. 
    mkdir cuda-test
    cd cuda-test
    cp -r /usr/local/cuda-11.3/samples .
    cd  samples/
    which nvcc
    make 
    bin/x86_64/linux/release/deviceQuery
    bin/x86_64/linux/release/deviceQueryDrv
    
**Note** later versions of cuda don't provide the samples anymore. Instead the samples are provided in a git repository 
    
    git clone https://github.com/nvidia/cuda-samples
    cd cuda-samples
    make --ignore-errors 
    
After the cuda-samples are build copy them to /usr/local/cuda or create a symbolic link /usr/local/cuda/sampes that points to the installation of cuda-samples. Opticks uses some of the sample headers and therefore building opticks will not suceed without them. 


Two tools for monitoring Nvidia GPUs On Linux can be found here:
https://www.linuxuprising.com/2019/06/2-tools-for-monitoring-nvidia-gpus-on.html

# Installing Optix (6.5)

https://developer.nvidia.com/designworks/optix/download

Optix comes with precompiled samples and one might want to try them:

    # cd to the Optix installation directory
    cd SDK-precompiled-samples
    export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH
    # execute e.g.:
    ./optixMDLSphere
    ./optixSphere
    ./optixTutorial
    # etc.





# Building opticks vs. existing libraries

This are instructions how to build opticks making use of preinstalled libraries available on the system. These libraries include CLHEP, xerces-c, boost and  Geant4.
For geant 4 we use the current version at the time of writing which is Geant4.10.7.p2. We make use of the fact that the om-cmake function of om.bash is sensitive to the CMAKE_PREFIX_PATH envvar so that we can point to the directories where the libraries are installed and avoid having to rebuild them.  In principle just cut and paste the following line to a file change the envars of the different directories to match your system and source the resulting script.

    cd to the directory where you want to install Opticks (the WORK_DIR environmental variable will point to this directory). 
    
    
Here we are using branch of Opticks which can be found in github. It contains some adjustments that we had to do to make opticks work with Geant4 version 11 and up which introduced some changes to the Geant4 API. Also there has been some rearrangement of the cuda header files in later cuda version that needed to be accounted for.

git clone https://github.com/hanswenzel/opticks opticks.v0.1.7
git checkout  v0.1.7
git status

Tagged versions of Simon's Blyth opticks can be found in:

    git clone https://github.com/simoncblyth/opticks.git
    cd opticks
    git checkout tags/v0.1.6 -b v0.1.6-branch
    git status
    
The development version (a. k. a. the latest and greatest) can be found in the following repository: 

    git clone https://bitbucket.org/simoncblyth/opticks.git
    
    
But again here we do:
git clone https://github.com/hanswenzel/opticks opticks.v0.1.7
git checkout  v0.1.7
git status
  
change opticks/optickscore/OpticksSwitches.h

so that:

    #define WITH_SKIPAHEAD 1

is set. To create a setup_opticks.sh you can use the cat statement below. But you have to edit the created file so that the environmental variables WORK_DIR, OptiX_INSTALL_DIR, OPTICKS_COMPUTE_CAPABILITy, CUDA_INSTALL_DIR, CUDA_SAMPLES, G4INSTALL correspond to the correct values of your system. 

    cat > setup_opticks.sh << +EOF
      # ----------------------------------------------------------------------------------------------------------------------
    # --- you need to modify the following environmental variables so that point to the specific directories on your system
    # --- 
    export WORK_DIR=/data/software/gpu4112
    export OptiX_INSTALL_DIR=/data/software/NVIDIA-OptiX-SDK-6.5.0-linux64
    export OPTICKS_COMPUTE_CAPABILITY=86
    export CUDA_INSTALL_DIR=/usr/local/cuda
    export CUDA_SAMPLES=${CUDA_INSTALL_DIR}/samples
    export G4INSTALL=/data/software/geant4-v11.0.2-install
    export LOCAL_BASE=${WORK_DIR}/local
    # ----------------------------------------------------------------------------------------------------------------------
    export CMAKE_PREFIX_PATH=${G4INSTALL}:${LOCAL_BASE}/opticks/externals:${OptiX_INSTALL_DIR}:${WORK_DIR}/opticks/cmake/Modules/:${WORK_DIR}/local /opticks:${WORK_DIR}/local/opticks:${WORK_DIR}/local/opticks/externals/
    export PYTHONPATH=$WORK_DIR
    export OPTICKS_HOME=${WORK_DIR}/opticks
    export OPTICKS_PREFIX=${WORK_DIR}/local/opticks                            
    export OPTICKS_INSTALL_PREFIX=$LOCAL_BASE/opticks
    export OPTICKS_OPTIX_PREFIX=/data/software/NVIDIA-OptiX-SDK-6.5.0-linux64/
    export OPTICKS_CUDA_PREFIX=${CUDA_INSTALL_DIR}
    export OPTICKS_EMBEDDED_COMMANDLINE_EXTRA="--rngmax 10 --rtx 1 --skipaheadstep 10000"
    #   
    # setup Geant4 and root
    #
    . ${G4INSTALL}/bin/geant4.sh
    . /data/software/root-install/bin/thisroot.sh
    opticks-(){ . ${OPTICKS_HOME}/opticks.bash && opticks-env $* ; }
    op(){ op.sh $* ; }
    o(){ cd $(opticks-home) ; hg st ; }
    _path_prepend() {
        if [ -n "$2" ]; then
            case ":$(eval "echo \$$1"):" in
                    *":$2:"*) :;;
                *) eval "export $1=$2\${$1:+\":\$$1\"}" ;;
         esac
        else
            case ":$PATH:" in
                *":$1:"*) :;;
                *) export PATH="$1${PATH:+":$PATH"}" ;;
            esac    
        fi
    }   

    _path_append() {
        if [ -n "$2" ]; then
            case ":$(eval "echo \$$1"):" in
                *":$2:"*) :;;
                *) eval "export $1=\${$1:+\"\$$1:\"}$2" ;;
            esac
        else
            case ":$PATH:" in
             *":$1:"*) :;;
                *) export PATH="${PATH:+"$PATH:"}$1" ;;
            esac
        fi
    }
    _path_prepend "${LOCAL_BASE}/bin"
    # make sure to add the compiler options
    new=" -fPIC" 
    case ":${CXXFLAGS:=$new}:" in
        *:"$new":*)  ;;
        *) CXXFLAGS="$CXXFLAGS:$new"  ;;
    esac
    new=" -fPIC" 
    case ":${CFLAGS:=$new}:" in
        *:"$new":*)  ;;
        *) CFLAGS="$CFLAGS:$new"  ;;
    esac
    # speed up the make process
    new=" -j$(($(grep -c ^processor /proc/cpuinfo) - 1))" 
    case ":${MAKEFLAGS:=$new}:" in
        *:"$new":*)  ;;
        *) MAKEFLAGS="$MAKEFLAGS:$new"  ;;
    esac
    # deal with the $LD_LIBRARYPATH
    new=${OptiX_INSTALL_DIR}/lib64/
    case ":${LD_LIBRARY_PATH:=$new}:" in
     *:"$new":*)  ;;
        *) LD_LIBRARY_PATH="$new:$LD_LIBRARY_PATH"  ;;
    esac
    new=${OPTICKS_HOME}/externals/lib
    case ":${LD_LIBRARY_PATH:=$new}:" in
        *:"$new":*)  ;;
        *) LD_LIBRARY_PATH="$new:$LD_LIBRARY_PATH"  ;;
    esac
    new=${CUDA_INSTALL_DIR}/lib64/
    case ":${LD_LIBRARY_PATH:=$new}:" in
         *:"$new":*)  ;;
         *) LD_LIBRARY_PATH="$new:$LD_LIBRARY_PATH"  ;;
    esac
    new=${LOCAL_BASE}/opticks/lib/
    case ":${LD_LIBRARY_PATH:=$new}:" in
        *:"$new":*)  ;;
        *) LD_LIBRARY_PATH="$new:$LD_LIBRARY_PATH"  ;;
    esac

    opticks-
    new=${CUDA_INSTALL_DIR}/bin
    case ":${PATH:=$new}:" in
        *:"$new":*)  ;;
        *) PATH="$new:$PATH"  ;;
    esac
    new=${OPTICKS_HOME}/bin/
    case ":${PATH:=$new}:" in
     *:"$new":*)  ;;
        *) PATH="$new:$PATH"  ;;
    esac
    new=${OPTICKS_HOME}/ana/
    case ":${PATH:=$new}:" in
        *:"$new":*)  ;;
        *) PATH="$new:$PATH"  ;;
    esac
    new=${LOCAL_BASE}/opticks/lib/
    case ":${PATH:=$new}:" in
        *:"$new":*)  ;;
        *) PATH="$new:$PATH"  ;;
    esac
    new=${CUDA_SAMPLES}/bin/x86_64/linux/release/
    case ":${PATH:=$new}:" in
     *:"$new":*)  ;;
        *) PATH="$new:$PATH"  ;;
    esac
    oinfo-(){
        echo 'LD_LIBRARY_PATH:';
        echo '================';
        echo  ${LD_LIBRARY_PATH}| tr : \\n;
        echo;
        echo 'PATH:';
        echo '=====';
        echo  ${PATH}| tr : \\n;
        echo;
        echo 'CMAKE_PREFIX_PATH:';
        echo '==================';
        echo  ${CMAKE_PREFIX_PATH}| tr : \\n;
        }
    dinfo-(){    
        nvidia-smi;
        ${CUDA_SAMPLES}/bin/x86_64/linux/release/deviceQuery
    }
    +EOF

Change setup_opticks.sh so that the environmental variables are set to correspond to your installation as instructed in the script. Then get some information  about your system:
    source setup_opticks.sh
    oinfo-
    dinfo-
 To build the opticks external packages do:
 
    mkdir -p ${WORK_DIR}/local/opticks/externals/
    cd ${WORK_DIR}/local/opticks/externals/
    ln -s ${OptiX_INSTALL_DIR} OptiX
    cd ${WORK_DIR}
    opticks-externals-install >& install_ext.log &

scan the log file or any errors and correct them. Now build opticks:

    cd ${WORK_DIR}
    opticks-full  >& install_full.log &

scan the log file or any errors and correct them. Before you run 
    opticks-t 

you want to to create the geocache using e.g. one of the gdml files provided by CaTS
(https://github.com/hanswenzel/CaTS)

    geocache-
    geocache-create- --gdmlpath  ${WORK_DIR}/local/opticks/opticksdata/export/juno1808/g4_00.gdml
    export OPTICKS_KEY=`output from above command"
    opticks-t
    
if the geocache-create- command doesn't work try:


    ${WORK_DIR}/local/opticks/lib/OKX4Test --okx4test --g4codegen --deletegeocache --gdmlpath  ${WORK_DIR}/CaTS/gdml/simpleLArTPC.gdml
    export OPTICKS_KEY=`output from above command"
    opticks-t

Also CaTS ( see instructions how to build and run CaTS below) rebuilds the geocache each time you run it and it will print out what to set the OPTICKS_KEY to:
    e.g.  export OPTICKS_KEY=CaTS.X4PhysicalVolume.World_PV.6a511c07e6d72b5e4d71b74bd548e8fd
    
# Building CaTS

Instructions for building and running CaTS can be found here: [README](README.md)
