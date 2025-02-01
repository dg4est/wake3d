#!/bin/bash -e
# Folder structure:
# ================
# project_root
#   install            #install directory
#   builds             #build directory
#   src                #source codes

# ======================= #
# project directory paths
# ======================= #
CURRENT_PATH="$(pwd)"
MY_PATH="$( cd "$( dirname "$0" )" && pwd )"
PROJECT_ROOT=${MY_PATH}

# ====================== #
# folder directory paths
# ====================== #
INSTALL_DIRECTORY=${PROJECT_ROOT}/install
INSTALL_WAKE_DIRECTORY=${INSTALL_DIRECTORY}/wake3d
BUILD_DIRECTORY=${PROJECT_ROOT}/builds
BUILD_WAKE_DIRECTORY=${BUILD_DIRECTORY}/wake3d

# ============== #
# wake3d sources
# ============== #
SOURCES_DIRECTORY=${PROJECT_ROOT}/src
WAKE_DIRECTORY=${SOURCES_DIRECTORY}
WAKE_SOURCES_DIRECTORY=${WAKE_DIRECTORY}

# ================== #
# compiling defaults
# ================== #
BUILD_WAKE=1
BUILD_TEST=0
BUILD_CLEAN=0
COMPILE_FAIL=0

# ================= #
# compiler defaults
# ================= #
CC=mpicc
CXX=mpicxx

C_FLAGS=
CXX_FLAGS=

# ======================== #
# compiler option defaults
# ======================== #
BUILD_SUFFIX="_release"
BUILD_TYPE="Release"

# ======================== #
# make and install command
# ======================== #
MAKE_CMD="make -j4 install"

# ============== #
# print strings
# ============== #
opt_str="[OPTION] "

eC="\x1B[0m"
GC="\x1B[1;32m"
rC="\x1B[0;41m"
gC="\x1B[0;42m"
yC="\x1B[0;33m"
oC="\x1B[3;93m"
aC="\x1B[3;92m"
mC="\x1B[0;43m"

help() {
    echo -e "Usage: ${yC}$0 <OPTIONS> <COMPILER OPTIONS>${eC}"
    echo " "
    echo "  [OPTION]:"
    echo "    --help     -h      displays this help message"
    echo "    --clean    -c      removes build directory: builds/wake3d_{}"
    echo "    --release  -opt    compile the project in optimized mode"
    echo "    --debug    -deb    compile the project in debug mode"
    echo " "
    echo "  [COMPILER OPTIONS]:"
    echo "     CC=<arg>     cc=<arg>    sets the C compiler"
    echo "    CXX=<arg>    cxx=<arg>    sets the C++ compiler"
    echo " "
    echo "      C_FLAGS=<arg>    c_flags=<arg>    sets the C compiler flags"
    echo "    CXX_FLAGS=<arg>  cxx_flags=<arg>    sets the C++ compiler flags"
    echo " "
    echo -e "  ${aC}Recommended Options:${eC}"
    echo -e "    Default (-go): ${yC}./configure.sh -go${eC}"
    echo -e "     ${oC}sets${eC} CC=mpicc CXX=mpicxx"
    echo "  "
    echo -e "    Intel MPI (-impi): ${yC}./configure.sh -impi${eC}"
    echo -e "     ${oC}sets${eC} CC=mpiicc CXX=mpiicpc"
}

# ----------------------------- #
# Start the compilation process #
# ----------------------------- #
cd $PROJECT_ROOT

# ============ #
# parse inputs
# ============ #
for var in "$@"
do
  if [ "$var" == "--help" -o "$var" == "-help" -o "$var" == "-h" ]; then
    help
    exit 0

  elif [ "$var" == "--clean" -o "$var" == "-clean" -o "$var" == "-c" -o "$var" == "-dc" ]; then
    echo ${opt_str} "Clean build directories."
    BUILD_CLEAN=1
    BUILD_WAKE=0

  elif [ "$var" == "--release" -o "$var" == "-release" -o "$var" == "-opt" ]; then
    echo ${opt_str} "Compiling in optimized mode"
    BUILD_SUFFIX="_release"
    BUILD_TYPE="Release"

  elif [ "$var" == "--debug" -o "$var" == "-debug" -o "$var" == "-deb" ]; then
    echo ${opt_str} "Compiling in debug mode"
    BUILD_SUFFIX="_debug"
    BUILD_TYPE="Debug"

  elif [ "$var" == "--tests" -o "$var" == "-tests" ]; then
    echo ${opt_str} "Turning on unit testing"
    BUILD_TEST=1

  elif [ "${var:0:3}" == "CC=" -o "${var:0:3}" == "cc=" ]; then
    CC=${var:3}
    echo -e "[OPTION]       C Compiler: $yC$CC$eC"

  elif [ "${var:0:4}" == "CXX=" -o "${var:0:4}" == "cxx=" ]; then
    CXX=${var:4}
    echo -e "[OPTION]     CXX Compiler: $yC$CXX$eC"

  elif [ "${var:0:8}" == "C_FLAGS=" -o "${var:0:8}" == "c_flags=" ]; then
    C_FLAGS=${var:8}
    echo -e "[OPTION]       C Compiler Flags: $yC$C_FLAGS$eC"

  elif [ "${var:0:10}" == "CXX_FLAGS=" -o "${var:0:10}" == "cxx_flags=" ]; then
    CXX_FLAGS=${var:10}
    echo -e "[OPTION]     CXX Compiler Flags: $yC$CXX_FLAGS$eC"

  elif [ "$var" == "-go" ]; then
    CC=mpicc
    CXX=mpicxx
    FC=mpif90

  elif [ "$var" == "-impi" ]; then
    CC=mpiicc
    CXX=mpiicpc
    FC=mpiifort

  fi
done

# ========================= #
# display command line args
# ========================= #
echo " "
echo "$0 $@"

# ----------------------------------------------------- #
# After reading in cmd arg options, set remaining paths #
# ----------------------------------------------------- #
# ====================================== #
# install/build location compiled source
# ====================================== #
COMPILE_INSTALL_WAKE_DIRECTORY="${INSTALL_WAKE_DIRECTORY}${BUILD_SUFFIX}"
COMPILE_BUILD_WAKE_DIRECTORY="${BUILD_WAKE_DIRECTORY}${BUILD_SUFFIX}"

# ============== #
# compiler paths
# ============== #
CC_PATH="`which $CC`"
CXX_PATH="`which $CXX`"
LD_PATH="`which ld`"

# ====================== #
# check source directory
# ====================== #
if [ ! -d "${WAKE_SOURCES_DIRECTORY}" ]; then
  echo " "
  echo "Error:"
  echo "${WAKE_SOURCES_DIRECTORY} does not exist."
  exit 1
fi

# ======================= #
# check install directory
# ======================= #
if [ ! -d "${INSTALL_DIRECTORY}" ]; then
  echo  "${INSTALL_DIRECTORY} does not exist. Making it..."
  mkdir "${INSTALL_DIRECTORY}"
fi

# ====================== #
# check builds directory
# ====================== #
if [ ! -d "${BUILD_DIRECTORY}" ]; then
  echo  "${BUILD_DIRECTORY} does not exist. Making it..."
  mkdir "${BUILD_DIRECTORY}"
fi

# =================================================================== #
if [ $BUILD_CLEAN == 1 ]; then
  echo " "
  echo "Clean: removing ${COMPILE_BUILD_WAKE_DIRECTORY}..."
  echo "Clean: removing ${COMPILE_INSTALL_WAKE_DIRECTORY}..."
  echo " "
  rm -rf $COMPILE_BUILD_WAKE_DIRECTORY
  rm -rf $COMPILE_INSTALL_WAKE_DIRECTORY
  rm -rf $BUILD_DIRECTORY
  rm -rf $INSTALL_DIRECTORY
  unlink bin
  exit 0
fi

if [ $BUILD_WAKE == 1 ]; then
  echo " "
  echo -e "${mC} ===================== Building WAKE3D ==================== ${eC}"
  echo         "   Compiling Options:"
  echo         "          Build Type: ${BUILD_TYPE}"
  echo         " "
  echo         "                  CC: ${CC}"
  echo         "                 CXX: ${CXX}"
  echo         " "
  echo         "            CC Flags: ${C_FLAGS}"
  echo         "           CXX Flags: ${CXX_FLAGS}"
  echo         " "
  echo         "          Build Type: ${BUILD_TYPE}"
  echo         "      Build Location: ${COMPILE_BUILD_WAKE_DIRECTORY}"
  echo         "    Install Location: ${COMPILE_INSTALL_WAKE_DIRECTORY}"
  echo         " Executable Location: ${COMPILE_BUILD_WAKE_DIRECTORY}/bin"
  echo -e "${mC} ========================================================== ${eC}"
  echo " "

  # move to the build directory
  cd $BUILD_DIRECTORY

  if [ ! -d $COMPILE_BUILD_WAKE_DIRECTORY ]; then
    mkdir $COMPILE_BUILD_WAKE_DIRECTORY
  fi
  cd $COMPILE_BUILD_WAKE_DIRECTORY

  cmake -D CMAKE_C_COMPILER=${CC_PATH}                              \
        -D CMAKE_CXX_COMPILER=${CXX_PATH}                           \
        -D CMAKE_C_FLAGS=${C_FLAGS}                                 \
        -D CMAKE_CXX_FLAGS=${CXX_FLAGS}                             \
        -D CMAKE_LINKER=${LD_PATH}                                  \
        -D CMAKE_INSTALL_PREFIX=${COMPILE_INSTALL_WAKE_DIRECTORY}   \
        -D CMAKE_BUILD_TYPE=${BUILD_TYPE}                           \
        -G "Unix Makefiles" ${PROJECT_ROOT} | tee cmake_config.out

  ${MAKE_CMD}
  cd ${CURRENT_PATH}

  if [ ! -d "${COMPILE_INSTALL_WAKE_DIRECTORY}" ]; then
    echo "ERROR:"
    echo "${COMPILE_INSTALL_WAKE_DIRECTORY} does not exist."
    COMPILE_FAIL=1
  fi
fi

if [ ${COMPILE_FAIL} == 0 ]; then
  echo " "
  echo -e " ========================================================== "
  echo -e " ${gC}WAKE3D build successful! ${eC}"
  echo    "   Compiling Options:"
  echo    "          Build Type: ${BUILD_TYPE}"
  echo    "          Unit Tests: ${UNIT_TEST}"
  echo    " "
  echo    "                  CC: ${CC}"
  echo    "                 CXX: ${CXX}"
  echo    " "
  echo    "             C Flags: ${C_FLAGS}"
  echo    "           CXX Flags: ${CXX_FLAGS}"
  echo    " "
  echo    "          Build Type: ${BUILD_TYPE}"
  echo    "      Build Location: ${COMPILE_BUILD_WAKE_DIRECTORY}"
  echo -e "                    : ${GC}make clean; make -j install${eC} in this directory"
  echo    "    Install Location: ${COMPILE_INSTALL_WAKE_DIRECTORY}"
  echo    " Executable Location: ${COMPILE_BUILD_WAKE_DIRECTORY}/bin"
  echo -e " ========================================================== "
  echo    " "
else
  echo " "
  echo         "======================"
  echo -e "${rC} WAKE3D build FAILED! ${eC}"
  echo         "======================"
  echo " "
  exit 1
fi
# =================================================================== #

# Create hyperlink to bin directory
ln -sf ${COMPILE_BUILD_WAKE_DIRECTORY}/bin bin

# Copy sample input files to the executable directory
#cp ${PROJECT_ROOT}/example/input.* ${COMPILE_BUILD_WAKE_DIRECTORY}/bin/.

echo "================================================"
echo -e "${gC} Finished Successfully...${eC}"
echo -e " Executable: ${GC}${COMPILE_BUILD_WAKE_DIRECTORY}/bin/wake3d.mpi${eC}"
echo "================================================"
