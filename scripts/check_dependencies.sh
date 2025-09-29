#!/bin/bash

OS=$(uname -s)

if [ "$OS" == "Darwin" ]; then
  PACKAGES=("git" "llvm@20" "cmake" "eigen")

  for pkg in "${PACKAGES[@]}"; do
    # Check if the package is in the list of installed Homebrew packages
    if brew list --formula | grep -q "^${pkg}\$"; then
      echo "✅ ${pkg}: Installed (Homebrew)"
    else
      echo "❌ ${pkg}: NOT Installed. Attempting to install with Homebrew..."
      if brew install "${pkg}"; then
        echo "✅ ${pkg}: Successfully Installed (Homebrew)"
      else
        echo "❌ ${pkg}: Homebrew installation FAILED."
      fi
    fi
  done

elif [ "$OS" == "Linux" ]; then
  PACKAGES=("build-essential" "git" "libeigen3-dev") 
  sudo apt update -y &> /dev/null

  for pkg in "${PACKAGES[@]}"; do
    if dpkg -l | grep -q "^ii.*${pkg}"; then
      echo "✅ ${pkg}: Installed (apt)"
    else
      echo "❌ ${pkg}: NOT Installed. Attempting to install with apt..."
      if sudo apt install -y "${pkg}"; then
        echo "✅ ${pkg}: Successfully Installed (apt)"
      else
        echo "❌ ${pkg}: apt installation FAILED."
      fi
    fi
  done

  CURRENT_DIR="$(pwd)"

  # 1) download a new version of cmake
  CMAKE_INSTALL_PREFIX="/opt/cmake"
  CMAKE_FILE_BASE="cmake-4.1.1-linux-x86_64"
  CMAKE_INSTALLER="${CMAKE_FILE_BASE}.sh"
  CMAKE_DOWNLOAD_URL="https://github.com/Kitware/CMake/releases/download/v4.1.1/${CMAKE_INSTALLER}"
  if [ ! -d "$CMAKE_INSTALL_PREFIX/$CMAKE_FILE_BASE" ]; then
    echo "--- Installing CMake 4.1.1 ---"
    sudo apt install -y wget &> /dev/null
    if wget "${CMAKE_DOWNLOAD_URL}" -O /tmp/"${CMAKE_INSTALLER}"; then
      chmod +x /tmp/"${CMAKE_INSTALLER}"
      if sudo ./tmp/"${CMAKE_INSTALLER}" --skip-license --prefix="${CMAKE_INSTALL_PREFIX}"; then
        export PATH="${CMAKE_INSTALL_PREFIX}/bin:$PATH"
        echo "✅ CMake 4.1.1 installation script executed successfully."
        rm -f /tmp/"${CMAKE_INSTALLER}"
      else
        echo "❌ CMake installation script FAILED to execute."
        rm -f /tmp/"${CMAKE_INSTALLER}"
        exit 1
      fi
    else
      echo "❌ Failed to download CMake installer from ${CMAKE_DOWNLOAD_URL}."
      exit 1
    fi
  else 
    # Set PATH even if already installed, in case it was not set previously
    export PATH="${CMAKE_INSTALL_PREFIX}/bin:$PATH"
    echo "✅ CMake 4.1.1 already installed at ${CMAKE_INSTALL_PREFIX}"
  fi

  # 2) Add apt repo for clang++20 and install clang++20
  LLVM_INSTALLER="llvm.sh"
  NEW_LLVM_POLICY_TEXT="[arch=amd64 signed-by=/etc/apt/keyrings/apt.llvm.org.asc] "
  LLVM_LIST_FILE_PATTERN="/etc/apt/sources.list.d/archive_uri-http_apt_llvm_org_*.list"
  if ! [ -f "/usr/bin/clang++-20" ]; then
    echo "--- Setting up LLVM PPA for Clang++20 and dependencies ---"
    sudo apt install -y wget &> /dev/null
    if wget https://apt.llvm.org/llvm.sh -O /tmp/"${LLVM_INSTALLER}"; then
      chmod u+x /tmp/"${LLVM_INSTALLER}"
      if sudo ./tmp/"${LLVM_INSTALLER}" 20; then

        # update debian policies
        sudo mkdir -p /etc/apt/keyrings
        sudo mv /etc/apt/trusted.gpg.d/apt.llvm.org.asc /etc/apt/keyrings/
        for FILE in $LLVM_LIST_FILE_PATTERN; do
          if [[ -f "$FILE" ]]; then
            sed -i "s|\(deb \)\(http\)|\1${NEW_LLVM_POLICY_TEXT}\2|" "$FILE"
          else
            echo "No LLVM list file found matching the pattern: $LLVM_LIST_FILE_PATTERN"
          fi
        done

        # finish install clang-20
        sudo apt update -y &> /dev/null
        if sudo apt install -y ${CLANG_PACKAGES}; then
          echo "✅ Clang-20 environment installed."
          rm -f /tmp/"${CMAKE_INSTALLER}"
        else
          echo "❌ Clang-20 environment installation FAILED."
          rm -f /tmp/"${CMAKE_INSTALLER}"
          exit 1
        fi
      else
        echo "❌ LLVM install script for clang-20 failed."
        rm -f /tmp/"${CMAKE_INSTALLER}"
        exit 1
      fi
      
    else
      echo "❌ Failed to download llvm install script."
      exit 1
    fi
  else 
    echo "✅ Clang-20 Installed"
  fi

  # 3) cd back to starting dir
  cd "$CURRENT_DIR"

else
  echo "❌ Unsupported operating system: $OS"
  exit 1
fi