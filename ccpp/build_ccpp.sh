#!/bin/bash

set +x
set -eu

# List of valid/tested machines
VALID_MACHINES=( theia.intel theia.gnu theia.pgi cheyenne.intel cheyenne.gnu cheyenne.pgi macosx.gnu linux.gnu )

###################################################################################################

function usage   {
  echo "Usage: "
  echo "build_ccpp.sh MACHINE_ID [ MAKE_OPT ] [ clean_before ] [ clean_after ]"
  echo "    Where: MACHINE      [required] can be : ${VALID_MACHINES[@]}"
  echo "           CCPP_DIR     [required] is the target installation directory for CCPP"
  echo "           MAKE_OPT     [optional] can be any of the NEMSfv3gfs MAKE_OPT options; used:"
  echo "                                   SION=Y/N   (default N)"
  echo "                                   DEBUG=Y/N  (default N)"
  echo "                                   OPENMP=Y/N (default Y)"
  echo "                                   HYBRID=Y/N (default Y)"
  echo "                                   32BIT=Y/N (default N)"
  echo "                                   STATIC=Y/N (default N, STATIC=Y requires HYBRID=N)"
  echo "                                   SUITE=name_of_sdf_without_path.xml (only if STATIC=Y)"
  echo "           clean_before [optional] can be 'YES' (default) or 'NO'"
  echo "           clean_after  [optional] can be 'YES' (default) or 'NO'"
  exit 1
}

function checkvalid {
# Ensure value ($2) of variable ($1) is contained in list of validvalues ($3)
  if [ $# -lt 3 ]; then
    echo $FUNCNAME requires at least 3 arguments: stopping
    exit 1
  fi

  var_name=$1 && shift
  input_val=$1 && shift
  valid_vars=($@)

  for x in ${valid_vars[@]}; do
    if [ "$input_val" == "$x" ]; then
      echo "${var_name}=${input_val} is valid."
      return
    fi
  done
  echo "ERROR: ${var_name}=${input_val} is invalid."
  usage
  exit 1
}

function trim {
    local var="$1"
    # remove leading whitespace characters
    var="${var#"${var%%[![:space:]]*}"}"
    # remove trailing whitespace characters
    var="${var%"${var##*[![:space:]]}"}"   
    echo -n "$var"
}

###################################################################################################

# Check and process command line arguments

if [[ $# -lt 2 ]]; then usage; fi

readonly MACHINE_ID=$1
readonly CCPP_DIR=$2
readonly MAKE_OPT=${3:-}
readonly clean_before=${4:-YES}
readonly clean_after=${5:-YES}

checkvalid MACHINE_ID $MACHINE_ID ${VALID_MACHINES[@]}

# Run ccpp_prebuild.py from the top-level directory before building the CCPP framework and physics
cd ..
if [[ "${MAKE_OPT}" == *"STATIC=Y"* ]]; then
  if [[ "${MAKE_OPT}" == *"HYBRID=N"* ]]; then
    if [[ "${MAKE_OPT}" == *"SUITE="* ]]; then
      # Extract name of suite definition file
      TMP=${MAKE_OPT#*SUITE=}
      SUITE=${TMP% *}
      STATICFLAGS="--static --suite=ccpp/suites/${SUITE}"
    else
      echo "Error, option STATIC=Y requires suite definition file as SUITE=xyz.xml (w/o path)"
      exit 1
    fi
  else
    echo "Error, option STATIC=Y requires HYBRID=N"
    exit 1
  fi
else
  STATICFLAGS=""
fi
if [[ "${MAKE_OPT}" == *"DEBUG=Y"* ]]; then
  DEBUGFLAG="--debug"
else
  DEBUGFLAG=""
fi
./ccpp/framework/scripts/ccpp_prebuild.py --config=./ccpp/config/ccpp_prebuild_config.py ${STATICFLAGS} ${DEBUGFLAG}
cd -

# Generate CCPP cmake flags from MAKE_OPT
CCPP_CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=${CCPP_DIR} -DNCEPLIBS_DIR=${NCEPLIBS_DIR} -DMPI=ON"
CCPP_MAKE_FLAGS=""
if [[ "${MAKE_OPT}" == *"SION=Y"* ]]; then
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DSIONLIB=${SIONLIB}"
fi
if [[ "${MAKE_OPT}" == *"DEBUG=Y"* ]]; then
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DCMAKE_BUILD_TYPE=Debug"
  CCPP_MAKE_FLAGS="${CCPP_MAKE_FLAGS} VERBOSE=1"
elif [[ "${MAKE_OPT}" == *"REPRO=Y"* ]]; then
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DCMAKE_BUILD_TYPE=Bitforbit"
  CCPP_MAKE_FLAGS="${CCPP_MAKE_FLAGS} VERBOSE=1"
else
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DCMAKE_BUILD_TYPE=Release"
  CCPP_MAKE_FLAGS="${CCPP_MAKE_FLAGS} VERBOSE=1"
fi
if [[ "${MAKE_OPT}" == *"OPENMP=N"* ]]; then
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DOPENMP=OFF"
else
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DOPENMP=ON"
fi
if [[ "${MAKE_OPT}" == *"HYBRID=N"* ]]; then
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DTEMPLOG=ON"
else
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DTEMPLOG=OFF"
fi
if [[ "${MAKE_OPT}" == *"32BIT=Y"* ]]; then
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DDYN32=ON"
else
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DDYN32=OFF"
fi
if [[ "${MAKE_OPT}" == *"STATIC=Y"* ]]; then
  if [[ "${MAKE_OPT}" == *"HYBRID=N"* ]]; then
    CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DSTATIC=ON"
  else
    echo "Error, option STATIC=Y requires HYBRID=N"
    exit 1
  fi
fi

CCPP_CMAKE_FLAGS=$(trim "${CCPP_CMAKE_FLAGS}")
CCPP_MAKE_FLAGS=$(trim "${CCPP_MAKE_FLAGS}")

# Generate additional CCPP cmake flags depending on machine / compiler
if [[ "${MACHINE_ID}" == "macosx.gnu" ]]; then
  # Intel MKL (for FFTW)
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DMKL_DIR=${MKL_DIR}"
  # ESMF (for DGEMM) - replace with MKL version in the future?
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DESMF_LIB_DIR=${ESMF_LIB}"
  # netCDF (needed when linking ESMF)
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DNETCDF_DIR=${NETCDF}"
elif [[ "${MACHINE_ID}" == "linux.gnu" ]]; then
  # ESMF (for DGEMM) - replace with MKL version in the future?
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DESMF_LIB_DIR=${ESMF_LIB}"
  # netCDF (needed when linking ESMF)
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DNETCDF_DIR=${NETCDF}"
fi

# Build and install CCPP

echo "Building CCPP with options '${CCPP_CMAKE_FLAGS}' ..."
PATH_CCPP=${PWD}
PATH_CCPP_BUILD=${PWD}/build
PATH_CCPP_INC=${PWD}/include
PATH_CCPP_LIB=${PWD}/lib

if [ $clean_before = YES ]; then
    rm -fr ${PATH_CCPP_BUILD}
    rm -fr ${PATH_CCPP_INC}
    rm -fr ${PATH_CCPP_LIB}
fi
mkdir -p ${PATH_CCPP_BUILD}
cd ${PATH_CCPP_BUILD}
cmake ${CCPP_CMAKE_FLAGS} ${PATH_CCPP}
make ${CCPP_MAKE_FLAGS}
make ${CCPP_MAKE_FLAGS} install

if [ $clean_after = YES ]; then
    rm -fr ${PATH_CCPP_BUILD}
fi
