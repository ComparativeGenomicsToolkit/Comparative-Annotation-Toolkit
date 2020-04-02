#!/bin/bash -ex

sudo sed -i -e 's/^Defaults\tsecure_path.*$//' /etc/sudoers

echo "GO Version:"
go version

echo "Python Version:"
python --version
pip install --user sregistry[all]

echo "sregistry Version:"
sregistry version

# Install Singularity

export PATH="${GOPATH}/bin:${PATH}"

mkdir -p "${GOPATH}/src/github.com/sylabs"
cd "${GOPATH}/src/github.com/sylabs"

git clone -b vault/release-3.3 https://github.com/sylabs/singularity
cd singularity
./mconfig -v -p /usr/local
make -j `nproc 2>/dev/null || echo 1` -C ./builddir all
sudo make -C ./builddir install