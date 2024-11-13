CC = /home/liy/packages/cuda-12.2/bin/nvcc
CXX = g++
CXXFLAGS = -std=c++14
CUDALIBDIR = /home/liy/packages/cuda-12.2/lib64

all: compcore

libcompcore.o: compcore.cu
	$(CC) -Xcompiler -fPIC -ccbin /usr/bin/gcc $(CXXFLAGS) -c $< -o $@

compcore: libcompcore.o
	$(CXX) -std=c++17 -O3 -o $@ ./*.cpp $< -I/home/liy/packages/cuda-12.2/include -L$(CUDALIBDIR) -lcudart

clean:
	rm -f compcore libcompcore.o