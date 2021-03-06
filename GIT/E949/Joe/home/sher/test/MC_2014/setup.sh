# change to build directory of geant4 version 9.6

unset  CLHEP_BASE_DIR
unset  CLHEP_INCLUDE_DIR
unset  CLHEP_LIB
unset  CLHEP_LIB_DIR

unset  G4ANALYSIS_USE
unset  G4DEBUG
unset  G4INCLUDE
unset  G4INSTALL

unset  G4LEDATA
unset  G4LEVELGAMMADATA
unset  G4NEUTRONHPDATA
unset  G4RADIOACTIVEDATA
unset  G4ABLADATA

unset  G4LIB
unset  G4LIB_BUILD_G3TOG4
unset  G4LIB_BUILD_SHARED
unset  G4LIB_BUILD_STATIC
unset  G4LIB_BUILD_ZLIB
unset  G4LIB_BUILD_GDML
unset  G4LIB_USE_G3TOG4
unset  G4LIB_USE_GRANULAR
unset  G4LIB_USE_ZLIB

unset  G4SYSTEM

unset  G4UI_BUILD_WIN32_SESSION
unset  G4UI_BUILD_XAW_SESSION
unset  G4UI_BUILD_XM_SESSION
unset  G4UI_USE_TCSH
unset  G4UI_USE_WIN32
unset  G4UI_USE_XAW
unset  G4UI_USE_XM
unset  G4UI_USE_QT

unset  G4VIS_BUILD_DAWN_DRIVER
unset  G4VIS_BUILD_OIWIN32_DRIVER
unset  G4VIS_BUILD_OIX_DRIVER
unset  G4VIS_BUILD_OPENGLWIN32_DRIVER
unset  G4VIS_BUILD_OPENGLXM_DRIVER
unset  G4VIS_BUILD_OPENGLX_DRIVER
unset  G4VIS_BUILD_RAYTRACERX_DRIVER
unset  G4VIS_BUILD_VRML_DRIVER
unset  G4VIS_BUILD_OPENGLQT_DRIVER

unset  G4VIS_USE_DAWN
unset  G4VIS_USE_OIWIN32
unset  G4VIS_USE_OIX
unset  G4VIS_USE_OPENGLWIN32
unset  G4VIS_USE_OPENGLX
unset  G4VIS_USE_OPENGLXM
unset  G4VIS_USE_RAYTRACERX
unset  G4VIS_USE_VRML
unset  G4VIS_USE_OPENGLQT


cd /home/dorothea/apps/geant4.9.6/geant4.9.6-build/
# source environment script
source geant4make.sh
# go back to work directory
#cd /home/dorothea/MC/
cd /home/sher/test/MC_2014/

export PATH=$PATH:$G4WORKDIR/bin/$G4SYSTEM

export G4WORKDIR=$PWD
# MCOUTPUT:  Redefine this variable to direct MC output, but be careful 
# with  inclusion of the "/" at the end. The code take the env variable
# MCOUTPUT and appends the number to it, if no "/" is present at the end
# it could produce a funny filename, for example /home/sher/mctest51.root
# instead of /home/sher/mctest5/1.root

export MCOUTPUT=$PWD/data/

export PYTHONPATH=$ROOTSYS/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/pienu/packages/log4cpp/lib

export XERCESCROOT=/home/pienumgr/apps64/xerces-c
export ROOTSYS=/home/pienu/packages/root/
#export CLHEP_BASE_DIR=/home/pienumgr/apps64/2011/clhep
#export CLHEP_INCLUDE_DIR=/home/pienumgr/apps64/2011/clhep/include

#export CLHEP_BASE_DIR=/home/dorothea/apps/CLHEP/install 
#export CLHEP_INCLUDE_DIR=/home/dorothea/apps/CLHEP/install/include
#export CLHEP_LIB_DIR=/home/dorothea/apps/CLHEP/install/libcase 

export PATH=$PATH:/home/pienumgr/apps/wired
export PATH=$PATH:$G4WORKDIR/bin/$G4SYSTEM
export PATH=$PATH:$CLHEP_LIB_DIR
#export PATH=$PATH:/home/pienumgr/apps64/2011/clhep/lib
export CRYHOME=/home/luca/cry_v1.5
export CRYDATAPATH=/home/luca/cry_v1.5/data


export PATH=$PATH:/home/pienumgr/apps/wired
export PATH=$PATH:/home/pienumgr/apps64/2011/clhep/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/pienumgr/apps64/2011/geant4/lib/geant4/Linux-g++
export CRYHOME=/home/luca/cry_v1.5
export CRYDATAPATH=/home/luca/cry_v1.5/data
export PATH=$PATH:$G4WORKDIR/bin/$G4SYSTEM
export PYTHONPATH=$ROOTSYS/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/pienu/packages/log4cpp/lib

