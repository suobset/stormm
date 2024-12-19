# Use the official Debian base image
FROM debian:latest

# Set non-interactive mode for apt-get to avoid prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install necessary packages
RUN apt-get update && \
    apt-get install -y \
    cmake \
    gcc \
    g++ \
    make \
    git \
    wget \
    curl \
    vim \
    emacs \
    neovim \
    ca-certificates \
    gnupg2 \
    lsb-release \
    software-properties-common \
    && apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Add contrib repository and update package list
RUN add-apt-repository contrib && \
    apt-get update

# Install CUDA Toolkit
RUN apt-get -y install cuda-toolkit-12-6

# Set up environment variables for CUDA
ENV PATH=/usr/local/cuda-12.6/bin:${PATH}
ENV LD_LIBRARY_PATH=/usr/local/cuda-12.6/lib64:${LD_LIBRARY_PATH}

# Create a directory for the repository
WORKDIR /app

# Clone the Git repository
RUN git clone https://github.com/psivant/stormm.git

# Build STORMM using CMake with CUDA enabled
RUN cmake -S stormm -B stormmbuild -DSTORMM_ENABLE_CUDA=OFF -DSTORMM_ENABLE_RDKIT=OFF
RUN export STORMM_HOME=/app/stormm
WORKDIR /app/stormmbuild
RUN make 

# Set the default command to run bash
CMD ["/bin/bash"]
