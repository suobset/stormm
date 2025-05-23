# Use NVIDIA CUDA base image
FROM nvidia/cuda:12.4.1-devel-ubuntu22.04

# Set non-interactive mode for apt-get to avoid prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# Install necessary packages
RUN apt-get update && \
    apt-get clean && \
    apt-get install -y \
    wget \
    python3-pip \
    git \
    libeigen3-dev \
    libboost-all-dev \
    gcc \
    g++ \
    cmake 

# Install CUDA 12.x
RUN echo "Checking and installing CUDA 12.x..."
RUN wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-keyring_1.1-1_all.deb
RUN dpkg -i cuda-keyring_1.1-1_all.deb
RUN apt-get update -y
RUN apt-get -y install cuda-toolkit-12-4
ENV PATH=/usr/local/cuda/bin:$PATH
ENV CUDADIR=/usr/local/cuda
ENV CXXFLAGS="-I/usr/local/cuda/include $CXXFLAGS"
ENV LDFLAGS="-L/usr/local/cuda/lib64 $LDFLAGS"

# Set the working directory
WORKDIR /app

# Copy the repository from the build context to the container
RUN mkdir -p /app/stormm
COPY . /app/stormm

# Build STORMM using CMake with CUDA enabled
RUN cmake -S stormm -B stormmbuild \
    -DSTORMM_ENABLE_CUDA=YES \
    -DSTORMM_ENABLE_RDKIT=NO \
    -DCUSTOM_GPU_ARCH=89

ENV STORMM_HOME=/app/stormm
ENV STORMM_SOURCE=/app/stormm
ENV STORMM_BUILD=/app/stormmbuild
WORKDIR /app/stormmbuild
RUN make -j

# After building all apps
COPY entrypoint /usr/local/bin/entrypoint
RUN chmod +x /usr/local/bin/entrypoint

# Set the default command to run bash
RUN echo "To run this container with GPU support, use the --gpus flag with docker run (e.g. docker run --gpus all stormm-config)"
ENTRYPOINT ["/usr/local/bin/entrypoint"]
