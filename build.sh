# Reset
Reset='\033[0m'           # Text Reset
Red='\033[0;31m'          # Red
Green='\033[0;32m'        # Green
Yellow='\033[0;33m'       # Yellow
Blue='\033[0;34m'         # Blue
Magenta='\033[0;35m'      # Purple
Cyan='\033[0;36m'         # Cyan
White='\033[0;37m'        # White

# Bold
BoldRed='\033[1;31m'         # Red
BoldGreen='\033[1;32m'       # Green
BoldYellow='\033[1;33m'      # Yellow
BoldBlue='\033[1;34m'        # Blue
BoldMagenta='\033[1;35m'     # Purple
BoldCyan='\033[1;36m'        # Cyan
BoldWhite='\033[1;37m'       # White

clear
rm -r build
mkdir build
cd build

echo -e "${BoldMagenta}-- BUILDING NAVTOOLS${Reset}";

# -DCMAKE_C_COMPILER=/usr/local/gcc-14.1.0/bin/gcc-14.1.0 \
# -DCMAKE_CXX_COMPILER=/usr/local/gcc-14.1.0/bin/g++-14.1.0 \

case "$OSTYPE" in
  linux*)
    echo -e "${BoldMagenta}-- OS: linux${Reset}";
    cmake .. \
        -DCMAKE_C_COMPILER=/usr/local/gcc-14.1.0/bin/gcc-14.1.0 \
        -DCMAKE_CXX_COMPILER=/usr/local/gcc-14.1.0/bin/g++-14.1.0 \
        -DINSTALL_NAVTOOLS_EXAMPLES=True \
        -DCMAKE_INSTALL_PREFIX=../build \
        -DCMAKE_BUILD_TYPE=Debug
        ;;
  darwin*)
    echo -e "${BoldMagenta}-- OS: mac${Reset}"; 
    cmake .. \
        -DCMAKE_C_COMPILER=clang-17 \
        -DCMAKE_CXX_COMPILER=clang++-17 \
        -DINSTALL_NAVTOOLS_EXAMPLES=True \
        -DCMAKE_INSTALL_PREFIX=../build \
        -DCMAKE_BUILD_TYPE=Debug
        ;;
  msys*)
    echo -e "${BoldMagenta}-- OS: windows${Reset}";
    cmake .. \
        -G "MinGW Makefiles" \
        -DCMAKE_CXX_COMPILER=C:/MinGW/bin/g++.exe \
        -DCMAKE_C_COMPILER=C:/MinGW/bin/gcc.exe \
        -DINSTALL_NAVTOOLS_EXAMPLES=True \
        -DCMAKE_INSTALL_PREFIX=../build \
        -DCMAKE_BUILD_TYPE=Debug
        ;;
  solaris*)
    echo -e "${BoldMagenta}-- OS: solaris${Reset}";;
  bsd*)
    echo -e "${BoldMagenta}-- OS: bsd${Reset}";;
  *)
    echo -e "${BoldMagenta}-- OS: unknown${Reset}";;
esac

cmake --build . -- -j4
# make
# make install
cd ..
