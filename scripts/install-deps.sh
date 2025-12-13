#!/bin/bash
# This script installs all dependencies (except Dorado) into the "bin" directory.
# Alternatively, consider system-wide installation of at least part of the tools
# (Ubuntu: sudo apt install ...)
# or Conda (but does currently not work with standard RStudio)

set -e

out=bin
cores=1

exists() {
    bash -c "command -v $1" > /dev/null
}

arch=$(arch)
os="$(uname -s | tr '[:upper:]' '[:lower:]')"
if [[ "$os" != "linux" && "$os" != "darwin" ]]; then
    echo "Linux (Windows WSL) or Mac expected" >&2
    exit 1
fi

# seqtool

if ! exists st; then
    v=0.4.0-beta.4
    echo "Installing seqtool v$v to bin/st"
    case "${os}" in
        linux*)  name=unknown-linux-gnu;;
        darwin*) name=apple-darwin;;
    esac
    url=https://github.com/markschl/seqtool/releases/download/v$v/seqtool-$arch-$name.tar.xz
    wget -qO - "$url" | tar -xJf -
    mv seqtool-$arch-$name/st "$out"
    rm -R seqtool-$arch-$name
fi


# VSEARCH

if ! exists vsearch; then
    v=2.30.1
    echo "Installing VSEARCH v$v to bin/vsearch"
    platform=$os
    if [[ "$os" == "darwin" ]]; then
        platform=mac
    fi
    url=https://github.com/torognes/vsearch/releases/download/v$v/vsearch-$v-$platform-$arch.tar.gz
    wget -qO - "$url" | tar -xzf -
    mv vsearch-$v-$platform-$arch/bin/vsearch "$out"
    rm -R vsearch-$v-$platform-$arch
fi

# Minimap2

if ! exists minimap2; then
    echo "Installing Minimap2 to bin/minimap2"
    v=2.30
    if [[ "$os-$arch" == "linux-x86_64" ]]; then
        url=https://github.com/lh3/minimap2/releases/download/v"$v"/minimap2-"$v"_x64-linux.tar.bz2
        wget -qO - "$url" | tar -xjf -
        mv minimap2-"$v"_x64-linux/minimap2 "$out"
        rm -R minimap2-"$v"_x64-linux
    elif exists make; then
        echo "No $os / $arch binary of Minimap2, building from source..." >&2
        url=https://github.com/lh3/minimap2/releases/download/v"$v"/minimap2-$v.tar.bz2
        wget -qO - "$url" | tar -xjf -
        (cd minimap2-"$v" && make -s -j $cores && mv minimap2 ../bin)
        rm -R minimap2-"$v"
    else
        echo "No $os / $arch binary of Minimap2, and 'make' not available for building from source" >&2
    fi
fi

# Samtools

if ! exists samtools; then
    if ! exists make; then
        echo "No 'make' available for building Samtools from source" >&2
        exit 1
    fi
    echo "Installing Samtools to bin/samtools"
    v=1.22.1
    url=https://github.com/samtools/samtools/releases/download/$v/samtools-$v.tar.bz2
    wget -qO - "$url" | tar -xjf -
    (cd samtools-$v && ./configure && make -j $cores && mv samtools ../bin)
    rm -R samtools-$v
fi
