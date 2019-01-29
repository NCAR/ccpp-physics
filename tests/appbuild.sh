#!/bin/bash

set -xeu

SECONDS=0

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 PATHTR MACHINE_ID [ MAKE_OPT [ BUILD_NR ] ]"
  exit 1
fi

readonly PATHTR=$1
readonly APP=$2
readonly BUILD_NAME=fv3${3:+_$3}

hostname

echo "Compiling app $APP into $BUILD_NAME.exe"

cd ${PATHTR}/..

rm -f $PATHTR/../NEMS/exe/NEMS.x
rm -f $PATHTR/../NEMS/src/conf/modules.nems

if [[ $APP = CCPP || $APP = CCPP_static_trans ]]; then
  if [[ $MACHINE_ID = theia.* ]]; then
    echo "Move original modulefile modulefiles/theia.intel/fv3 aside and replace with modulefiles/theia.intel/fv3.intel-15.1.133"
    cd modulefiles/theia.intel
    mv -v fv3 fv3.original
    ln -svf fv3.intel-15.1.133 fv3
    cd ../..
  else
    echo "ERROR, appbuild.sh with APP CCPP not configured for build target ${MACHINE_ID}"
    exit 1
  fi
fi
set +e
./NEMS/NEMSAppBuilder app="$APP"
RC=$?
set -e
cd ${PATHTR}/..
if [[ $APP = CCPP || $APP = CCPP_static_trans ]]; then
  if [[ $MACHINE_ID = theia.* ]]; then
    echo "Reinstantiate original modulefile modulefiles/theia.intel/fv3"
    cd modulefiles/theia.intel
    rm -v fv3
    mv -v fv3.original fv3
    cd ../..
  else
    echo "ERROR, appbuild.sh with APP CCPP not configured for build target ${MACHINE_ID}"
    exit 1
  fi
fi
if [[ $RC -ne 0 ]]; then
  echo "ERROR in './NEMS/NEMSAppBuilder app=$APP'"
  exit $RC
fi

cp $PATHTR/../NEMS/exe/NEMS.x ${PATHTR}/../tests/$BUILD_NAME.exe
cp $PATHTR/../NEMS/src/conf/modules.nems ${PATHTR}/../tests/modules.$BUILD_NAME

cd $PATHTR/../NEMS/src
gmake clean
cd $PATHTR
gmake cleanall

# A few things that "cleanall" doesn't clean:
rm -rf FV3_INSTALL
rm -rf nems_dir

elapsed=$SECONDS
echo "Elapsed time $elapsed seconds. Compiling app ${APP} finished"
