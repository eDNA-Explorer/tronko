# Ubuntu 24.04: the release with clang 19, which BWA-MEM3 (the aligner tronko-assign runs) needs
FROM ubuntu:24.04

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    make \
    zlib1g-dev \
    libzstd-dev \
    libc6-dev \
    git \
    ca-certificates \
    clang-19 \
    libomp-19-dev \
    libdeflate-dev \
    autoconf \
    automake \
    libtool \
    cmake \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy source code
COPY . /app/

# Build tronko-assign and the BWA-MEM3 it runs (pinned v0.14.0, built with clang 19)
WORKDIR /app/tronko-assign
RUN make clean && make && make bwa-mem3

# Build tronko-build
WORKDIR /app/tronko-build
RUN make clean && make

# Set the working directory back to root
WORKDIR /app

# Create a simple test script
RUN echo '#!/bin/bash\n\
echo "Testing tronko-assign with logging..."\n\
echo "Available options:"\n\
./tronko-assign/tronko-assign -h\n\
echo ""\n\
echo "Testing verbose logging (should show help and logging info):"\n\
./tronko-assign/tronko-assign -V 2 -h\n\
' > test_logging.sh && chmod +x test_logging.sh

CMD ["/bin/bash"]