#!/bin/bash

set +x
set -eu

# List of valid/tested machines
VALID_MACHINES=( gaea.intel jet.intel theia.intel theia.gnu theia.pgi cheyenne.intel cheyenne.gnu cheyenne.pgi endeavor.intel macosx.gnu linux.gnu )

###################################################################################################

function usage   {
  echo "Usage: "
  echo "build_ccpp.sh MACHINE_ID CCPP_DIR ESMF_MK [ 'MAKE_OPT' ] [ clean_before ] [ clean_after ]"
  echo "    Where: MACHINE      [required] can be : ${VALID_MACHINES[@]}"
  echo "           CCPP_DIR     [required] is the target installation directory for CCPP"
  echo "           ESMF_MK      [required] is the location of the ESMF makefile fragement"
  echo "           MAKE_OPT     [optional] can be any of the NEMSfv3gfs MAKE_OPT options,"
  echo "                                   enclosed in a single string; used:"
  echo "                                   SION=Y/N   (default N)"
  echo "                                   DEBUG=Y/N  (default N)"
  echo "                                   REPRO=Y/N  (default N)"
  echo "                                   TRANSITION=Y/N (default N)"
  echo "                                   OPENMP=Y/N (default Y)"
  echo "                                   HYBRID=Y/N (default Y)"
  echo "                                   32BIT=Y/N  (default N, affects dynamics/fast physics only)"
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
readonly ESMF_MK=$3
readonly MAKE_OPT=${4:-}
readonly clean_before=${5:-YES}
readonly clean_after=${6:-YES}

checkvalid MACHINE_ID $MACHINE_ID ${VALID_MACHINES[@]}

# Generate CCPP cmake flags from MAKE_OPT
CCPP_CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=${CCPP_DIR} -DNCEPLIBS_DIR=${NCEPLIBS_DIR} -DNETCDF_DIR=${NETCDF} -DMPI=ON"
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
  if [[ "${MACHINE_ID}" == "jet.intel" ]]; then
    CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DSIMDMULTIARCH=ON"
  fi
  CCPP_MAKE_FLAGS="${CCPP_MAKE_FLAGS} VERBOSE=1"
fi
if [[ "${MAKE_OPT}" == *"TRANSITION=Y"* ]]; then
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DTRANSITION=ON"
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
if [[ "${MAKE_OPT}" == *"HYBRID=N"* ]]; then
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DHYBRID=OFF"
else
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DHYBRID=ON"
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
elif [[ "${MACHINE_ID}" == "gaea.intel" ]]; then
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DESMF_LIB_DIR=${ESMF_LIB}"
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DLIBXML2_LIB_DIR=${LIBXML2_LIB_DIR} -DLIBXML2_INCLUDE_DIR=${LIBXML2_INCLUDE_DIR}"
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DCRAY=ON"
  if [[ "${MAKE_OPT}" == *"STATIC=Y"* ]]; then
    # DH* At this time, it is not possible to use the dynamic CCPP
    # build on gaea. While compiling/linking works, the model crashes
    # immediately. This may be related to 64bit/32bit mismatches
    # in the MPI libraries (missing "-fPIC" flags when the MPI libraries
    # were compiled on the system?) - to be investigated.
    ## FOR DYNAMIC BUILD, SET ENVIRONMENT VARIABLE
    #export CRAYPE_LINK_TYPE=dynamic
    echo "Dynamic CCPP build not supported on gaea at this time."
    exit 1
    # *DH
  fi
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
    rm -f ${ESMF_MK}
fi
mkdir -p ${PATH_CCPP_BUILD}
cd ${PATH_CCPP_BUILD}
cmake ${CCPP_CMAKE_FLAGS} ${PATH_CCPP}
make ${CCPP_MAKE_FLAGS}
make ${CCPP_MAKE_FLAGS} install

if [ $clean_after = YES ]; then
    rm -fr ${PATH_CCPP_BUILD}
fi

# Generate ESMF makefile fragment

# Explicitly append libxml2, with or without path
CCPP_XML2_LIB="${LIBXML2_LIB_DIR:+-L${LIBXML2_LIB_DIR} }-lxml2"
set -u
if ( echo "${MAKE_OPT}" | grep STATIC=Y ) ; then
  # Set linker flags for static build
  CCPP_LINK_OBJS="-L${PATH_CCPP_LIB} -lccpp -lccppphys ${CCPP_XML2_LIB}"
else
  # Set link objects
  if ( echo "$MACHINE_ID" | grep gaea ) ; then
    CCPP_LINK_OBJS="-dynamic -L${PATH_CCPP_LIB} -lccpp ${CCPP_XML2_LIB} ${CRAY_PMI_POST_LINK_OPTS} -lpmi"
  else
    CCPP_LINK_OBJS="-L${PATH_CCPP_LIB} -lccpp ${CCPP_XML2_LIB}"
  fi
fi
echo "ESMF_DEP_INCPATH=-I${PATH_CCPP_INC}" > ${ESMF_MK}
echo "ESMF_DEP_LINK_OBJS=${CCPP_LINK_OBJS}" >> ${ESMF_MK}
