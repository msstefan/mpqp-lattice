FROM mcr.microsoft.com/devcontainers/cpp:0-ubuntu-22.04

ARG REINSTALL_CMAKE_VERSION_FROM_SOURCE="none"

# Optionally install the cmake for vcpkg
COPY ./reinstall-cmake.sh /tmp/

RUN if [ "${REINSTALL_CMAKE_VERSION_FROM_SOURCE}" != "none" ]; then \
        chmod +x /tmp/reinstall-cmake.sh && /tmp/reinstall-cmake.sh ${REINSTALL_CMAKE_VERSION_FROM_SOURCE}; \
    fi \
    && rm -f /tmp/reinstall-cmake.sh

# [Optional] Uncomment this section to install additional vcpkg ports.
# RUN su vscode -c "${VCPKG_ROOT}/vcpkg install <your-port-name-here>"

# Install Eigen3, CDD-lib (Fukuda), and libigl dependencies
 RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
     && apt-get -y install --no-install-recommends libeigen3-dev libcdd-dev mesa-common-dev libgl1-mesa-dev libglu1-mesa-dev libxrandr-dev libxinerama-dev libxi-dev libxmu-dev libblas-dev libxcursor-dev libboost-dev

