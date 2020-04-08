#include<stdio.h>
#include<stdlib.h>
#include<cuda.h>
#include<math.h>
#include<cuda_runtime.h>
#include<cooperative_groups.h>
#include<cuda_runtime_api.h>
#include<bits/stdc++.h>

__global__ void cuda_hello(){
    printf("Hello World from GPU!\n");
}

int main() {
    cuda_hello<<<1,1>>>(); 

    cudaDeviceSynchronize();
    return 0;
}