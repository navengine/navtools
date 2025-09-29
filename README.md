# Navtools
A foundational C++ template library for basic navigation utilities.

## Prerequisites
For simplicity, you can run:
```sh
./scripts/check_dependencies.sh
```

Otherwise:
1) First a valid C++20 compiler should be installed. I recommend `clang`, which is easy on Ubuntu and MacOS. For example, using MacOS homebrew to get clang++20:
```zsh
brew install llvm@20
```
 - For ubuntu follow this tutorial [Install Clang 20](https://ubuntuhandbook.org/index.php/2023/09/how-to-install-clang-17-or-16-in-ubuntu-22-04-20-04/). At minimum, run the following commands (*it is recommended to follow the tutorial to update install to follow new Debian policies*):
```sh
wget https://apt.llvm.org/llvm.sh -O /tmp/llvm.sh
chmod u+x /tmp/llvm.sh
sudo ./tmp/llvm.sh 20
sudo apt update -y
sudo apt install -y clang-tidy-20 clang-format-20 clang-tools-20 llvm-20-dev lld-20 lldb-20 llvm-20-tools libomp-20-dev libc++-20-dev libc++abi-20-dev libclang-common-20-dev libclang-20-dev libclang-cpp20-dev liblldb-20-dev libunwind-20-dev
```

2) Next an updated version of `cmake` is also necessary. For MacOS:
```zsh
brew install cmake
```
 - On Ubuntu, this can be installed from the shell script at the cmake downloads page, for example using cmake 4.1.1:
```sh
wget https://github.com/Kitware/CMake/releases/download/v4.1.1/cmake-4.1.1-linux-x86_64.sh -O /tmp/cmake-4.1.1-linux-x86_64.sh
chmod +x /tmp/cmake-4.1.1-linux-x86_64.sh
sudo /tmp/cmake-4.1.1-linux-x86_64.sh --skip-license --prefix=/opt/cmake
export PATH=/opt/cmake/bin:$PATH
```

3) Make sure Eigen is installed on your computer.
```sh
(ubuntu) sudo apt install -y libeigen3-dev
(macos)  brew install eigen
```

## Building
To build the the standalone C++ project use the build script:
```sh
./scripts/build.sh
```

# TODO:
1) Finish the "wgs84" header.
2) Test "wgs84" models.
3) Test "utils"