# Use NVIDIA CUDA base image
FROM nvidia/cuda:12.8.1-devel-ubuntu22.04

# Set non-interactive mode for apt-get to avoid prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# Install necessary packages
RUN apt-get update && \
    apt-get clean && \
    apt-get install -y \
    wget \
    python3-pip \
    tzdata \
    gcovr \
    git \
    libeigen3-dev \
    libboost-all-dev \
    gcc \
    g++ \
    cmake \
    && apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install CUDA 12.x
RUN echo "Checking and installing CUDA 12.x..." && \
    wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-keyring_1.1-1_all.deb && \
    dpkg -i cuda-keyring_1.1-1_all.deb && \
    apt-get update -y && \
    apt-get -y install cuda-toolkit-12-8 && \
    export PATH=/usr/local/cuda/bin:$PATH && \
    export CUDADIR=/usr/local/cuda && \
    export CXXFLAGS="-I/usr/local/cuda/include $CXXFLAGS" && \
    export LDFLAGS="-L/usr/local/cuda/lib64 $LDFLAGS"

# Set the working directory
WORKDIR /app

# Copy the code from the CI/CD pipeline context into the Docker image
RUN git clone https://github.com/psivant/stormm.git

# Build STORMM using CMake with CUDA enabled
RUN cmake -S stormm -B stormmbuild \
    -DCMAKE_BUILD_TYPE=RELEASE \
    -DSTORMM_ENABLE_CUDA=YES \
    -DSTORMM_ENABLE_RDKIT=NO \
    -DCUSTOM_NVCC_THREADS=8

RUN export STORMM_HOME=/app/stormm
WORKDIR /app/stormmbuild
RUN make -j

# Set the default command to run bash
CMD ["/bin/bash"]