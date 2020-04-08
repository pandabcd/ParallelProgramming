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
#define num_edges 25
#define num_vertices1 5
#define num_vertices2 5


// Stores edges
__device__ int d_edges_u[num_edges];
__device__ int d_edges_v[num_edges];

// __device__ unsigned int d_visited_1[num_vertices1]={0};

__device__ unsigned int d_visited_2[num_vertices2]={0};  //visited or not
__device__ unsigned int d_matched[num_edges]={0};		// matched or not

__device__ unsigned int d_first_edge[num_vertices1+1];  // making vertex disjoint


__global__ 
void get_approx_matching(){
	int tid = blockIdx.x*1024 + threadIdx.x;
	if(tid<num_edges){
		// printf("[%d] %d \n", tid, d_first_edge[tid]);
		for(int i=d_first_edge[tid]; i<d_first_edge[tid+1];i++)
		{
			int visited2 = atomicExch(&d_visited_2[d_edges_v[i]], 1);
			if(!visited2)
			{
				printf("Pairing %d with %d (edge number %d)\n", d_edges_u[i], d_edges_v[i], i);
				d_matched[i] = 1;
				return;
			}
		}
	}
	
}


int main(){
	int fc = num_vertices1;
	int set1[num_vertices1], set2[num_vertices2];
	int first_edge[num_vertices1+1]; // one artifial index at last for easier coding
	int edges_u[num_edges], edges_v[num_edges];
	int matched[num_edges]={0};

	ifstream fin;
    fin.open("FC_" + to_string(fc) + "_" + to_string(fc) + ".txt", ios::in);
    int u, v;

    for(int i=0;i<num_vertices1;i++){
    	fin >> first_edge[i];
    }
    first_edge[num_vertices1] = num_edges;

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
    
    cudaMemcpyToSymbol(d_first_edge, first_edge, (num_vertices1+1)*sizeof(int),0,cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_edges_u, edges_u, num_edges*sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_edges_v, edges_v, num_edges*sizeof(int), 0, cudaMemcpyHostToDevice);


	get_approx_matching<<<1, num_threads>>>();

	cudaMemcpyFromSymbol(matched, d_matched, num_edges*sizeof(int), 0, cudaMemcpyDeviceToHost);


	cudaDeviceSynchronize();
	cout << "Printing matched edges"<<endl;
	for(int i=0;i<num_edges;i++){
		if(matched[i]){
			cout << edges_u[i] << " " << edges_v[i] << endl;
		}
	}

	
	return 0;
}