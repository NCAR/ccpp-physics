#!/bin/bash

SITE="stampede.intel"
FV3GFS_DIR=$(pwd)/..

./compile.sh ${FV3GFS_DIR}/FV3 ${SITE} "DEBUG=N 32BIT=Y" 32bit YES NO |& tee make.out.32bit
