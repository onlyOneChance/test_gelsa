#!/bin/bash

# # Run nvidia-smi and capture its output
# nvidia_smi_output=$(nvidia-smi)

# # Search for the line containing CUDA Version and extract CUDA version
# cuda_version=$(echo "$nvidia_smi_output" | sed -n '/CUDA Version/ s/.*: \([0-9.]*\).*/\1/p')

# # Print CUDA version
# echo "CUDA version: $cuda_version"

# Run nvidia-smi and capture its output
nvidia_smi_output=$(nvidia-smi)

# Search for the line containing CUDA Version and extract CUDA version
cuda_version=$(echo "$nvidia_smi_output" | sed -n '/CUDA Version/ s/.*: \([0-9.]*\).*/\1/p')

# Check if CUDA version is captured
if [ -z "$cuda_version" ]; then
    echo "Please install an NVIDIA GPU the corresponding version NVIDIA driver and the corresponding version of CUDA.\
    Now building by cpu"
    # If CUDA version is not detected, exit or perform other actions
else
    # Print CUDA version
    echo "CUDA version: $cuda_version"
fi


#Compiler Dependencies:
apt-get update && apt-get install -y wget sudo python3.10 python3.10-dev python3-pip
sudo update-alternatives --install /usr/bin/python python /usr/bin/python3.10 1
pip install numpy scipy pandas statsmodels 
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
sudo apt-get update
sudo apt-get -y install cuda-toolkit-12-3

# #Get Latest Release:
# wget https://bitbucket.org/charade/lsa/get/release.tar.gz

#Python Dependencies:
cd ./lsa/ &&
rm -rf ./compcore.so &&
cd ../pybind_compcore/ &&
/usr/local/cuda-$cuda_version/bin/nvcc -Xcompiler -fPIC \
-ccbin /usr/bin/gcc \
-std=c++14 \
-c ./compcore.cu \
-o ./libcompcore.o && 
g++ -std=c++14 -fPIC -shared \
./*.cpp \
./libcompcore.o \
-I /usr/include/python3.10 \
-L /usr/lib/python3.10 \
-lpython3.10 \
-I../pybind11/include \
-I/usr/local/cuda-$cuda_version/include \
-L/usr/local/cuda-$cuda_version/lib64 \
-lcudart \
-O3 -o ../lsa/compcore.so &&
cd .. &&
python setup.py install

cd lsa

m


echo "CUDA version: $cuda_version"