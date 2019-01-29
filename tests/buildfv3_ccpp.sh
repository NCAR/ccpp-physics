#!/bin/bash

SITE="stampede.intel"
FV3GFS_DIR=$(pwd)/..

./compile.sh ${FV3GFS_DIR}/FV3 ${SITE} "CCPP=Y REPRO=Y HYBRID=N" 64bit YES YES |& tee make.out.32bit
