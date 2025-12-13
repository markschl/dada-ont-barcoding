#!/bin/bash
# This script installs Dorado

set -euo pipefail


out=bin
version=1.3.0  # default version
if [[ "$#" >= 1 ]]; then
    version=$1
fi

if dorado --version > /dev/null 2>&1; then
    echo "Dorado already installed systemwide"
elif bin/dorado/bin/dorado --version > /dev/null 2>&1; then
    echo "Dorado already installed in bin/dorado"
else
    os="$(uname -s | tr '[:upper:]' '[:lower:]')"
    arch=$(arch)
    if [[ "$os" != "linux" && "$os" != "darwin" ]]; then
        echo "Linux (Windows WSL) or Mac expected" >&2
        exit 1
    fi
    echo "Installing dorado v$v to bin/dorado"
    if [[ "$os" == "darwin" ]]; then
        os=osx
    fi
    if [[ "$arch" == "x86_64" ]]; then
        arch=x64
    fi
    name=dorado-$version-$os-$arch
    wget -O - "https://cdn.oxfordnanoportal.com/software/analysis/$name.tar.gz" | tar -xzf -
    mv $name "$out"/dorado
fi
