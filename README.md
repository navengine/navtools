# navtools
Common functions and definitions used across navengine.

## Prerequisites
First a valid C++20 compiler should be installed. I recommend clang++-18 which is easy on ubuntu 24.04.
```sh
sudo apt install clang-18
```

Make sure Eigen and Pybind11 are installed on your computer.
```sh
sudo apt install libeigen3-dev
sudo apt install python3-pybind11
```

## Building
To build the the standalone C++ project use the build script:
```sh
./build.sh
```
You can turn off the python installation by setting the `build_python` option to `False` in the `build.sh` file.

To build the python project, first create a virtual environment, and then pip install the project (replace the "." with the path to the navtools folder):
```sh
python3 -m venv .venv
. .venv/bin/activate
pip install numpy
pip install .
```

## Python Linting (Optional)
This should already be applied, but in case it is not, you can simply do it as follows:
```sh
pip install pybind11-stubgen
pybind11-stubgen navtools -o src
```