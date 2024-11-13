#include"device_launch_parameters.h"
#include <thrust/extrema.h>
#include <stdio.h>

using namespace std;

__device__ void kernel_2(float* Gpu_ptr_X, float* Gpu_ptr_Y, int ROW, int COL, int x, int y, int* r, float* score, int* x_0, int* y_0, int* x_1, int* y_1)              //动态规划判断
{
	int shift = ROW / 2;
	int block_size = gridDim.y;
	int thread_size = blockDim.x;

	extern __shared__ float ptr_s[];

	// 使用单线程为共享内存赋初值=======性能损失0.5s/w
	// if (threadIdx.x == 0)
	// 	for(int i=0;i<ROW;++i)
	// 	{
	// 		for(int j=0;j<COL;++j)
	// 		{
	// 			ptr_s[i*COL+j]=0;
	// 		}		
	// 	}	
    //  使用单线程计算矩阵==============性能损失1s/w
	// if (threadIdx.x == 0)
	// for (int i = 0; i < ROW; ++i)
    // {
    // 	if (i <= shift)
    //     	for (int j = 0; j < COL - shift + i; ++j)
    //         	ptr_s[i * COL + j] = Gpu_ptr_X[shift - i + j] * Gpu_ptr_Y[j];
    //     else
    //         for (int j = 0; j < COL + shift - i; ++j)
    //         ptr_s[i * COL + j] = Gpu_ptr_X[j] * Gpu_ptr_Y[i - shift + j];
    // }

	int tid = threadIdx.x;
	while (tid < ROW * COL)
	{
		ptr_s[tid]=0;
		tid += thread_size;
	}
	
	int id = threadIdx.x;
	int cid = id;
	int h_cid = id;
	int rx = x * COL + id + shift;
	int ry = y * COL + id;
	int cid_0 = cid;
	int rx_0 = rx;
	int ry_0 = ry;
	int d_r = ROW - 1;

	//动态规划矩阵计算
	while (cid < (shift + 1) * COL)
	{
		while (rx < (x + 1) * COL)
		{
			ptr_s[cid] = Gpu_ptr_X[rx] * Gpu_ptr_Y[ry];
			ptr_s[h_cid + d_r * COL] = Gpu_ptr_X[ry + (x - y) * COL] * Gpu_ptr_Y[rx - (x - y) * COL];
			h_cid += thread_size;
			cid += thread_size;
			rx += thread_size;
			ry += thread_size;
		}
		cid_0 += COL;
		rx_0 -= 1;
		cid = cid_0;
		rx = rx_0;
		ry = ry_0;
		d_r -= 1;
		h_cid = id;
	}
	__syncthreads();

	//对某一block中矩阵进行动态规划
	if (threadIdx.x == 0)
	{
		int I = 0;
		int P = 0;
		int J = 0;
		for (int i = 0; i < ROW; ++i)
		{
			int p = 0;                                         //正相关起始位置
			int n = 0;                                         //负相关起始位置
			float p_s = 0;
			float n_s = 0;
			float max_p = 0;
			float max_n = 0;

			for (int j = 0; j < COL; ++j)
			{
				p_s = max(p_s + ptr_s[i * COL + j], (float)0.0);            //状态转移方程
				n_s = max(n_s - ptr_s[i * COL + j], (float)0.0);            //状态转移方程
	
				if (p_s == 0) p = j + 1;
				else  if (p_s > max_p) max_p = p_s;
				if (n_s == 0) n = j + 1;
				else  if (n_s > max_n) max_n = n_s;

				if (max_p > max_n) { if (max_p > score[block_size * x + y]) { score[block_size * x + y] = max_p; r[block_size * x + y] = 1; I = i, P = p, J = j; } }
				else 			   { if (max_n > score[block_size * x + y]) { score[block_size * x + y] = max_n; r[block_size * x + y] =-1; I = i, P = n, J = j; } }
			}
		}

		if (I <= shift)
		{
			x_0[block_size * x + y] = shift - I + P + 1; y_0[block_size * x + y] = P + 1; x_1[block_size * x + y] = shift - I + J + 1;  y_1[block_size * x + y] = J + 1;
		}
		else
		{
			x_0[block_size * x + y] = P + 1; y_0[block_size * x + y] = I - shift + P + 1; y_1[block_size * x + y] = I - shift + J + 1; x_1[block_size * x + y] = J + 1;
		}
		score[block_size * x + y] *= r[block_size * x + y];
		score[block_size * x + y] /= COL;
	}
};

__global__ void kernel_1(float* Gpu_ptr_X, float* Gpu_ptr_Y,int num_01,int num_02,int shift,int COL,int* r, float* score,int* x_0, int* y_0, int* x_1, int* y_1)
{
	int ROW = 2 * shift + 1;
	int x = blockIdx.x;
	int y = blockIdx.y;

	kernel_2(Gpu_ptr_X,Gpu_ptr_X,ROW, COL, x, y, r, score, x_0, y_0, x_1, y_1);
};

extern "C"
void gpu_compcore(float* x,float* y,long long num_01,long long num_02,int shift,long long COL,float* score_data,int* x_0,int* x_1,int* y_0,int* y_1)
{
	int gpuCount = -1;
	cudaGetDeviceCount(&gpuCount);//获取显卡设备数量
	// cout << "gpuGount:\n" << gpuCount << endl;
	// if (gpuCount < 0)
	// 	cout << "np device!\n" << endl;
	cudaSetDevice(0);//选取显卡设备 

	cudaDeviceProp prop;
	int count;
	cudaGetDeviceCount(&count);
	// printf("gpu num %d\n", count);
	cudaGetDeviceProperties(&prop, 0);
	// printf("max thread num: %d\n", prop.maxThreadsPerBlock);
	// printf("max grid dimensions: %d, %d, %d)\n",prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);

	int ROW = 2 * shift + 1;
	float* Gpu_ptr_X = 0;
	float* Gpu_ptr_Y = 0;
	float* Gpu_score = 0;

	int* Gpu_r = 0;
	int* Gpu_x_0 = 0;
	int* Gpu_y_0 = 0;
	int* Gpu_x_1 = 0;
	int* Gpu_y_1 = 0;

	cudaMalloc((void**)&Gpu_ptr_X, num_01*COL * sizeof(float));
	cudaMalloc((void**)&Gpu_ptr_Y, num_02*COL * sizeof(float));
	cudaMalloc((void**)&Gpu_score, num_01*num_02 * sizeof(float));
	cudaMalloc((void**)&Gpu_r, num_01 * num_02 * sizeof(int));
	cudaMalloc((void**)&Gpu_x_0, num_01 * num_02 * sizeof(int));                                                                                               
	cudaMalloc((void**)&Gpu_y_0, num_01 * num_02 * sizeof(int));                                                                                               
	cudaMalloc((void**)&Gpu_x_1, num_01 * num_02 * sizeof(int));                                                                                               
	cudaMalloc((void**)&Gpu_y_1, num_01 * num_02 * sizeof(int));

	cudaMemcpy(Gpu_ptr_X, x, num_01*COL * sizeof(float), cudaMemcpyHostToDevice);//============================================设置流操作
	cudaMemcpy(Gpu_ptr_Y, y, num_02*COL * sizeof(float), cudaMemcpyHostToDevice);//============================================设置流操作

	dim3 grid(num_01, num_02);

	kernel_1 << <grid, 32, ROW* COL * sizeof(float) >> > (Gpu_ptr_X,Gpu_ptr_Y,num_01,num_02,shift,COL,Gpu_r,Gpu_score,Gpu_x_0,Gpu_y_0,Gpu_x_1,Gpu_y_1);

	cudaMemcpy(score_data,Gpu_score,num_01 * num_02 * sizeof(float), cudaMemcpyDeviceToHost);

	cudaMemcpy(x_0,Gpu_x_0,num_01 * num_02 * sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(y_0,Gpu_y_0,num_01 * num_02 * sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(x_1,Gpu_x_1,num_01 * num_02 * sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(y_1,Gpu_y_1,num_01 * num_02 * sizeof(int), cudaMemcpyDeviceToHost);

	cudaFree(Gpu_ptr_X);
	cudaFree(Gpu_ptr_Y);
	cudaFree(Gpu_r);
	cudaFree(Gpu_score);
	cudaFree(Gpu_x_0);
	cudaFree(Gpu_y_0);
	cudaFree(Gpu_x_1);
	cudaFree(Gpu_y_1);
}
