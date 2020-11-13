#!/bin/bash

make DEBUG=1 -j 4
make sdk_install_pkg DEBUG=1
export DEB_BUILD_OPTIONS="nostrip"
make deb_sgx_enclave_common_pkg DEBUG=1
./linux/installer/bin/sgx_linux_x64_sdk_*.bin
sudo dpkg -i ./linux/installer/deb/libsgx-enclave-common/libsgx-enclave-common_*.deb

