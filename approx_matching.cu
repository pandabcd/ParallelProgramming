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


using namespace std;

#define num_threads 50
#define num_edges 4
#define num_vertices1 2
#define num_vertices2 2


__device__ int d_edges_u[num_edges];
__device__ int d_edges_v[num_edges];

__device__ unsigned int d_visited_1[num_vertices1]={0};
__device__ unsigned int d_visited_2[num_vertices2]={0};
__device__ unsigned int d_matched[num_edges]={0};

__device__ unsigned int d_first_edge[num_vertices1];

__global__ 
void get_approx_matching(){
	int tid = blockIdx.x*1024 + threadIdx.x;
	if(tid<num_edges){

		int visited1 = atomicExch(&d_visited_1[d_edges_u[tid]], 1);   
		if(!visited1){
			printf("[%d]%d is unvisitedd \n", tid,d_edges_u[tid]);

			int visited2 = atomicExch(&d_visited_2[d_edges_v[tid]], 1);
			if(!visited2)
			{
				   
				printf("[%d]%d is unvisited::::: %d \n",tid, d_edges_v[tid], visited2);
				printf("Pairing %d with %d \n", d_edges_u[tid], d_edges_v[tid]);
			}
			else{
				printf("[%d]%d is visitedd \n",tid, d_edges_v[tid]);

			}
		}
	}	
}


int main(){
	int fc = num_vertices1;
	int set1[num_vertices1], set2[num_vertices2];
	int first_edge[num_vertices1];
	int edges_u[num_edges], edges_v[num_edges];

	ifstream fin;
    fin.open("FC_" + to_string(fc) + "_" + to_string(fc) + ".txt", ios::in);
    int u, v;

    for(int i=0;i<fc;i++){
    	fin >> first_edge[i];
    }

    cout << "Printing all the edges: \n";

    // Vertices with 0 edges are implicitly ignored while reading the file itself
    for(int i=0;i<num_edges;i++){
            fin >> u >> v;
            // Check if not a repeat, then add
            // set1.push_back(u);
            // set2.push_back(v);
            edges_u[i] = u;
            edges_v[i] = v;
            cout << u << " " << v <<endl;
    }
    
    cudaMemcpyToSymbol(d_first_edge, first_edge, num_vertices1*sizeof(int),0,cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_edges_u, edges_u, num_edges*sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_edges_v, edges_v, num_edges*sizeof(int), 0, cudaMemcpyHostToDevice);


	get_approx_matching<<<1, num_threads>>>();

	cudaDeviceSynchronize();
	return 0;
}