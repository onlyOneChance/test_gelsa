/home/liy/packages/cuda-12.2/bin/nvcc -Xcompiler -fPIC \
-ccbin /usr/bin/gcc \
-std=c++14 \
-c ./compcore.cu \
-o ./libcompcore.o

g++ -std=c++17 -O3 -o compcore \
./*.cpp \
./libcompcore.o \
-I/home/liy/packages/cuda-12.2/include \
-L/home/liy/packages/cuda-12.2/lib64 \
-lcudart

./compcore

export LD_LIBRARY_PATH=/home/liy/packages/cuda-12.2/lib64:$LD_LIBRARY_PATH

#稳定设置
echo 'export LD_LIBRARY_PATH=/home/liy/packages/cuda-12.2/lib64:$LD_LIBRARY_PATH' >> ~/.bashrc
source ~/.bashrc

nvidia-smi pmon
watch -n 1 nvidia-smi
clear&&make clean&&make&&./compcore
