#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<stdio.h>
#include<cuda.h>
#include<math.h>
#include<fstream>
#include<cuda_runtime.h>
#include<cooperative_groups.h>
#include<cuda_runtime_api.h>

#define BLOCK_SIZE 2
#define GRID_SIZE 2

__global__ void test(int A[BLOCK_SIZE][BLOCK_SIZE], int B[BLOCK_SIZE][BLOCK_SIZE],int C[BLOCK_SIZE][BLOCK_SIZE])
{

    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int j = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < BLOCK_SIZE && j < BLOCK_SIZE)
        C[i][j] = A[i][j] + B[i][j];
    printf("%d \n");

}

int main()
{

    int d_A[BLOCK_SIZE][BLOCK_SIZE];
    int d_B[BLOCK_SIZE][BLOCK_SIZE];
    int d_C[BLOCK_SIZE][BLOCK_SIZE];

    int C[BLOCK_SIZE][BLOCK_SIZE];

    for(int i=0;i<BLOCK_SIZE;i++)
      for(int j=0;j<BLOCK_SIZE;j++)
      {
        d_A[i][j]=i+j;
        d_B[i][j]=i+j;
      }


    // dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE); 
    // dim3 dimGrid(GRID_SIZE, GRID_SIZE); 

    // cudaMemcpy(C,d_A,BLOCK_SIZE*BLOCK_SIZE , cudaMemcpyDeviceToHost);
    // cudaMemcpy(C,d_B,BLOCK_SIZE*BLOCK_SIZE , cudaMemcpyDeviceToHost);
    // cudaMemcpy(C,d_C,BLOCK_SIZE*BLOCK_SIZE , cudaMemcpyDeviceToHost);
    
    test<<<dimGrid, dimBlock>>>(d_A,d_B,d_C); 



    // for(int i=0;i<BLOCK_SIZE;i++)
    //     {
    //       for(int j=0;j<BLOCK_SIZE;j++)
    //     {
    //       printf("%d\n",C[i][j]);

    //     }
    // }
}