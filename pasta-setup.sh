#!/bin/bash
set -euo pipefail

# 1. System deps
sudo apt install -y python3 python3-pip openjdk-21-jdk gcc cmake mafft fasttree hmmer

# 2. Python deps
pip3 install --break-system-packages dendropy
pip3 install --break-system-packages git+https://github.com/uym2/TreeShrink.git

# 3. Build tronko-build
cd ~/tronko/tronko-build && make

# 4. sate-tools-linux (PASTA runtime deps — opal.jar, hmmbuild, etc.)
git clone https://github.com/smirarab/sate-tools-linux.git ~/sate-tools-linux
mkdir -p ~/tronko/pasta/bin
cd ~/tronko/pasta/bin && ln -sf ~/sate-tools-linux/* .

# Override bundled mafft (v7.149b/2014) with system mafft
ln -sf "$(which mafft)" ~/tronko/pasta/bin/mafft

# 5. PASTA scripts not in sate-tools-linux
cp ~/tronko/pasta/resources/scripts/hmmeralign ~/tronko/pasta/bin/hmmeralign
sed -i '1s|#! /usr/bin/env python|#!/usr/bin/env python3|' ~/tronko/pasta/bin/hmmeralign
chmod +x ~/tronko/pasta/bin/hmmeralign
ln -sf "$(which run_treeshrink.py)" ~/tronko/pasta/bin/treeshrink

# 6. Build VeryFastTree
cd /tmp && git clone https://github.com/citiususc/veryfasttree.git
cd veryfasttree && mkdir -p build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local && make -j$(nproc)
sudo make install

# 7. Run
cd ~/tronko
TREE_BACKEND=veryfasttree bash pasta-3builds.sh 2>&1 | tee pasta-build.log
