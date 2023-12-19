#------------------------------------------------------------------------------------
# optionally enable logging in relevant classes
#export SEvt=INFO
#export G4CXOpticks=INFO
#export SEventConfig=INFO
#export CSGOptiX=INFO
#export SSim__stree_level=2
#------------------------------------------------------------------------------------
# Setting envvar GEOM to an identifier for your geometry
#will be saved in:
#/tmp/${USER}/opticks/GEOM/$GEOM/CaTS/ALL${VERSION}/...
export GEOM=simpleLArTPC
#
# You have to change the following environmental variables so that they point
# to the correct directories in your installation
#
#------------------------------------------------------------------------------------
export WORK_DIR=/data3/${user}/alexei5
export OptiX_INSTALL_DIR=/home/wenzel/NVIDIA-OptiX-SDK-7.5.0-linux64-x86_64
export OPTICKS_COMPUTE_CAPABILITY=75
export CUDA_INSTALL_DIR=/usr/local/cuda
export CUDA_SAMPLES=/home/wenzel/cuda-samples
export G4INSTALL=/data3/wenzel/geant4-v11.1.2_nomt-install
export ROOTINSTALL=/data3/wenzel/root_6_28_04-install
#------------------------------------------------------------------------------------
export LOCAL_BASE=${WORK_DIR}/local
export CMAKE_PREFIX_PATH=${G4INSTALL}:${LOCAL_BASE}/opticks/externals:${OptiX_INSTALL_DIR}:${WORK_DIR}/opticks/cmake/Modules/:${WORK_DIR}/local/opticks:${WORK_DIR}/local/opticks:${WORK_DIR}/local/opticks/externals/:/usr/local/lib/CLHEP-2.4.6.2
export PYTHONPATH=$WORK_DIR
export OPTICKS_HOME=${WORK_DIR}/opticks
export OPTICKS_PREFIX=${WORK_DIR}/local/opticks                            
export OPTICKS_INSTALL_PREFIX=$LOCAL_BASE/opticks
export OPTICKS_OPTIX_PREFIX=${OptiX_INSTALL_DIR}
export OPTICKS_CUDA_PREFIX=${CUDA_INSTALL_DIR}
export OPTICKS_EMBEDDED_COMMANDLINE_EXTRA="--rngmax 100 --rtx 1 --skipaheadstep 10000"
export OPTICKS_MAX_PHOTON=10000000
#
# setup Geant4 and root
#
. ${G4INSTALL}/bin/geant4.sh
. ${ROOTINSTALL}/bin/thisroot.sh
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
#export PATH=${LOCAL_BASE}/bin:${PATH}
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
export OPTICKS_KEY=CaTS.X4PhysicalVolume.World_PV.6a511c07e6d72b5e4d71b74bd548e8fd
