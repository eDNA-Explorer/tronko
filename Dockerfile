FROM ubuntu:22.04

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    make \
    zlib1g-dev \
    libc6-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy source code
COPY . /app/

# Build tronko-assign
WORKDIR /app/tronko-assign
RUN make clean && make

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