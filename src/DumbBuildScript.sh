#!/bin/bash

if [ ${1} ] ; then
  if [ ${1} == "-NVCC" ] || [ ${1} == "-nvcc" ] ; then
    COMPILE_NVCC=1
  else
    COMPILE_NVCC=0
  fi
  if [ ${1} == "-GCC" ] || [ ${1} == "-gcc" ] ; then
    COMPILE_GCC=1
  else
    COMPILE_GCC=0
  fi
else
  COMPILE_NVCC=1
  COMPILE_GCC=1
fi

module load gcc/6.4.0
module load cuda/11.2.2
module load valgrind/3.12.0

ALL_DIR="Constants Cuda DataTypes Math Parsing Random Reporting Topology UnitTesting"
ALL_FILES=""
for DIR in ${ALL_DIR} ; do
    for CPPFI in `ls ${DIR}/*.cpp` ; do
      ALL_FILES="${ALL_FILES} ${CPPFI}"
    done
done

if [ ${COMPILE_NVCC} -eq 1 ] ; then
  EXE_FLAGS="-std=c++11 -g"
  LINK="-L. -L${CUDA_HOME}/lib64 -L${CUDA_HOME}/lib64/stubs"
  LINK="${LINK} -lcurand -lcublas -lcusolver -lcudart -lcudadevrt -lnvidia-ml"
  INCL="-I. -I${CUDA_HOME}/include"
  ALL_OBJ_FILES=""
  for CPPFI in ${ALL_FILES} ; do
    FBASE=${CPPFI%%.cpp}
    FBASE=`cut -d "/" -f2 <<< ${FBASE}`
    #nvcc ${EXE_FLAGS} -c lib/${FBASE}.o ${CPPFI} ${LINK} ${INCL}
    ALL_OBJ_FILES="${ALL_OBJ_FILES} lib/${FBASE}.o"
  done
  nvcc ${EXE_FLAGS} -o omni.cuda omni.cpp -DOMNI_USE_CUDA ${ALL_FILES} ${LINK} ${INCL}
  cd ../test/
  for DIR in ${ALL_DIR} ; do
    if [ ! -d ${DIR} ] ; then
      mkdir ${DIR}
    fi
    for CPPFI in `ls ${DIR}/*.cpp` ; do
      FBASE=${CPPFI%%.cpp}
      FBASE=`cut -d "/" -f2 <<< ${FBASE}`
    done
  done
  cd ../src/
fi
if [ ${COMPILE_GCC} -eq 1 ] ; then
  EXE_FLAGS="-std=c++11 -Wall -g -o"
  LINK="-L."
  INCL="-I."
  g++ ${EXE_FLAGS} omni omni.cpp ${ALL_FILES} ${LINK} ${INCL}
  cd ../test/
  for DIR in ${ALL_DIR} ; do
    for CPPFI in `ls ${DIR}/*.cpp` ; do
      FBASE=${CPPFI%%.cpp}
      FBASE=`cut -d "/" -f2 <<< ${FBASE}`
    done
  done
  cd ../src/
fi
