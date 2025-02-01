```shell
__/\\\______________/\\\_____/\\\\\\\\\_____/\\\________/\\\__/\\\\\\\\\\\\\\\_____/\\\\\\\\\\___/\\\\\\\\\\\\____
 _\/\\\_____________\/\\\___/\\\\\\\\\\\\\__\/\\\_____/\\\//__\/\\\///////////____/\\\///////\\\_\/\\\////////\\\__
  _\/\\\_____________\/\\\__/\\\/////////\\\_\/\\\__/\\\//_____\/\\\______________\///______/\\\__\/\\\______\//\\\_
   _\//\\\____/\\\____/\\\__\/\\\_______\/\\\_\/\\\\\\//\\\_____\/\\\\\\\\\\\_____________/\\\//___\/\\\_______\/\\\_
    __\//\\\__/\\\\\__/\\\___\/\\\\\\\\\\\\\\\_\/\\\//_\//\\\____\/\\\///////_____________\////\\\__\/\\\_______\/\\\_
     ___\//\\\/\\\/\\\/\\\____\/\\\/////////\\\_\/\\\____\//\\\___\/\\\_______________________\//\\\_\/\\\_______\/\\\_
      ____\//\\\\\\//\\\\\_____\/\\\_______\/\\\_\/\\\_____\//\\\__\/\\\______________/\\\______/\\\__\/\\\_______/\\\__
       _____\//\\\__\//\\\______\/\\\_______\/\\\_\/\\\______\//\\\_\/\\\\\\\\\\\\\\\_\///\\\\\\\\\/___\/\\\\\\\\\\\\/__
        ______\///____\///_______\///________\///__\///________\///__\///////////////____\/////////_____\////////////___
```

## Getting Started

Configure and build WAKE3D
```shell
Usage: ./configure.sh <OPTIONS> <COMPILER OPTIONS>
 
  [OPTION]:
    --help     -h      displays this help message
    --clean    -c      removes build directory: builds/wake3d_{}
    --release  -opt    compile the project in optimized mode
    --debug    -deb    compile the project in debug mode
 
  [COMPILER OPTIONS]:
     CC=<arg>     cc=<arg>    sets the C compiler
    CXX=<arg>    cxx=<arg>    sets the C++ compiler
 
      C_FLAGS=<arg>    c_flags=<arg>    sets the C compiler flags
    CXX_FLAGS=<arg>  cxx_flags=<arg>    sets the C++ compiler flags
 
  Recommended Options:
    Default (-go): ./configure.sh -go
     sets CC=mpicc CXX=mpicxx
  
    Intel MPI (-impi): ./configure.sh -impi
     sets CC=mpiicc CXX=mpiicpc
```
