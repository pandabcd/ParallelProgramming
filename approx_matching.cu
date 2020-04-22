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

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

using namespace std;

#define num_threads 50
#define num_edges 25
#define num_vertices1 5
#define num_vertices2 5

__device__ unsigned int d_degree[num_vertices1+num_vertices2+1];    //Is this required?
__device__ unsigned int d_flat_adj_list[2*num_edges];
__device__ unsigned int d_list_ptr[num_vertices1+num_vertices2+2];

__device__ unsigned int d_matched_vertices[num_vertices1+num_vertices2+1]={0};
__device__ unsigned int d_matched_edges[2*num_edges]={0};
__device__ unsigned int d_visited[num_vertices1+num_vertices2+1]={0};
		


// Every vertex gets a node
__global__ 
void get_approx_matching(){
	int tid = blockIdx.x*1024 + threadIdx.x;
	int vertex = tid + 1;	// The world is 1-indexed
	if(vertex<=num_vertices1){

		printf("[%d]Looking from %d to %d \n" ,tid, d_list_ptr[vertex], d_list_ptr[vertex+1]);
		for(int i=d_list_ptr[vertex];i<d_list_ptr[vertex+1];i++){


			// Problem in here.... You can do it :)
			printf("[%d]working %d \n",tid, d_list_ptr[vertex]);
			int visited = atomicExch(&d_visited[d_list_ptr[vertex]], 1);    // Index of connected vertex
			printf("inside %d \n", visited);
			if(!visited)
			{
				printf("Pairing %d with %d \n", vertex, d_flat_adj_list[d_list_ptr[vertex]]);
				// d_matched[i] = 1;
				return;
			}
		}

	}
}


__global__
void vertex_disjoint_bfs(){

}


int main(){
	int fc = num_vertices1;
	
	int degree[num_vertices1+num_vertices2+1]={0};      //store degree of each vertex
	int flat_adj_list[2*num_edges];
	int list_ptr[num_vertices1+num_vertices2+2];        //1-indexed and extra element at the end for easy size access  // Pointer to the start of adjacency list
	int list_ptr_copy[num_vertices1+num_vertices2+2];    // Temporrary stuff, gotta sleep
	// Only required for results
	int matched_vertices[num_vertices1+num_vertices2+1]={0};
	int matched_edges[2*num_edges]={0};

	// to and from of edges
	int edges_u[num_edges], edges_v[num_edges];			// Make this dynamic memory and free it once we have our 2 pass initialisation phase
	

	
	ifstream fin;
    fin.open("FC_" + to_string(fc) + "_" + to_string(fc) + ".txt", ios::in);
    int u, v;

    cout << "Printing all the edges: \n";

    // Vertices with 0 edges are implicitly ignored while reading the file itself
    for(int i=0;i<num_edges;i++){
            fin >> u >> v;
            cout << u << " " << v <<endl;
            edges_u[i] = u;
            edges_v[i] = v;
            degree[u]++;
            degree[v]++;
    }

    // Get pointer to adjacency list using prefix sum (no opti here since other parts are more complex anyway)
    // Index 0 will never be used.... the last elem
    list_ptr[1] = 0;
    list_ptr_copy[1] = list_ptr[1];
    for(int i=2;i<=num_vertices1+num_vertices2;i++){
    	list_ptr[i] = list_ptr[i-1] + degree[i-1];
    	list_ptr_copy[i] = list_ptr[i];
    }
    list_ptr[num_vertices1+num_vertices2+1] = 2*num_edges;       //For easy coding
    list_ptr_copy[num_vertices1+num_vertices2+1] = 2*num_edges;



    for(int i=0;i<num_edges;i++){
    	flat_adj_list[list_ptr_copy[edges_u[i]]] = edges_v[i];
    	flat_adj_list[list_ptr_copy[edges_v[i]]] = edges_u[i];
    	list_ptr_copy[edges_u[i]]++;
    	list_ptr_copy[edges_v[i]]++;
    }

    cout << "Printing flat adjacency list for 4: " << endl;
    // for(int i=0;i<2*num_edges;i++){
    // 	cout << flat_adj_list[i] << endl;
    // }

    for(int i=list_ptr[4];i<list_ptr[5];i++){
    	cout << flat_adj_list[i] << endl;
    }

   
    cudaMemcpyToSymbol(d_degree, degree, (num_vertices1+num_vertices2+1)*sizeof(int),0,cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(d_flat_adj_list, flat_adj_list, (2*num_edges)*sizeof(int),0,cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(d_list_ptr, list_ptr, (num_vertices1+num_vertices2+2)*sizeof(int),0,cudaMemcpyHostToDevice);
	
    // cout<< list_ptr[0];
    cout<<endl<<endl;
	get_approx_matching<<<1, num_threads>>>();

	// cudaMemcpyFromSymbol(matched, d_matched, num_edges*sizeof(int), 0, cudaMemcpyDeviceToHost);


	cudaDeviceSynchronize();
	// cout << "Printing matched edges"<<endl;
	// for(int i=0;i<num_edges;i++){
	// 	if(matched[i]){
	// 		cout << edges_u[i] << " " << edges_v[i] << endl;
	// 	}
	// }

	
	return 0;
}