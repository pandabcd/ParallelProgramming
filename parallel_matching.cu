#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<fstream>
#include<vector>

using namespace std;

#define num_threads 10
// #define num_edges 700000
// #define num_vertices1 10000
// #define num_vertices2 10000

// #define num_edges 1000000
// #define num_vertices1 1000
// #define num_vertices2 1000


#define lli long long int

// #define num_edges 2998468
// #define num_vertices1 100000
// #define num_vertices2 100000



const lli num_edges = 4;
const lli num_vertices1 = 2;
const lli num_vertices2 = 2;

__device__ int d_flat_adj_list[2*num_edges];
__device__ int d_degree[num_vertices1+num_vertices2+1]={0};      //store degree of each vertex
__device__ int d_list_ptr[num_vertices1+num_vertices2+2];        //1-indexed and extra element at the end for easy size access  // Pointer to the start of adjacency list
__device__ int d_list_ptr_copy[num_vertices1+num_vertices2+2];    // Temporrary stuff, gotta sleep

__device__ bool d_is_matched_edge[(num_vertices1+ num_vertices2 + 1)*(num_vertices1 + num_vertices2+1)] = {0} ;     // Adjacency matrix (1-indexed)
__device__ bool d_is_matched_vertex[num_vertices1 + num_vertices2 + 1] = {0};	//is the vertex matched
__device__ int d_partner_vertex[num_vertices1 + num_vertices2 + 1];
__device__ int d_visited[num_vertices1 + num_vertices2 + 1] = {0};
__device__ int d_bfs_parent[num_vertices1 +  num_vertices2 + 1];
__device__ bool d_is_parent_change[num_vertices1 + num_vertices2 + 1] = {0};

__device__ int d_frontier[num_vertices1 + num_vertices2+1] = {0};
__device__ int d_next_frontier[num_vertices1+num_vertices2+1] = {0};



int *h_flat_adj_list;
int *h_degree;
int * h_list_ptr;
int *h_list_ptr_copy;

bool *h_is_matched_edge;
bool *h_is_matched_vertex;
int *h_partner_vertex;
int *h_visited;
int *h_bfs_parent;
bool *h_is_parent_change;

int fc = num_vertices1;
// int num_aug_paths = 0;

int *h_frontier;
int *h_next_frontier;


__device__ 
int get_is_matched_edge(int i, int j){
	return d_is_matched_edge[i*(num_vertices1 + num_vertices2+1) + j ];
}

__device__ 
void set_is_matched_edge(int i, int j, int value){
	d_is_matched_edge[i*(num_vertices1 + num_vertices2+1) + j ] = value;
}
__device__
void clear_visited(){
	int tid = blockIdx.x*1024 + threadIdx.x;
	int vertex1 = tid + 1;

	if(vertex1<=num_vertices1 + num_vertices2){
		d_visited[vertex1] = 0;
	}
}


__device__ 
void clear_bfs_parent(){
	int tid = blockIdx.x*1024 + threadIdx.x;
	int vertex1 = tid + 1;

	if(vertex1<=num_vertices1 + num_vertices2){
		d_bfs_parent[vertex1] = vertex1;
	}
}

__device__ 
void initialise_partner_vertex(){
	int tid = blockIdx.x*1024 + threadIdx.x;
	int vertex1 = tid + 1;

	if(vertex1<=num_vertices1 + num_vertices2){
		d_partner_vertex[vertex1] = -1;
	}
}

__device__ 
void clear_is_parent_change(){
	int tid = blockIdx.x*1024 + threadIdx.x;
	int vertex1 = tid + 1;

	if(vertex1<=num_vertices1 + num_vertices2){
		d_is_parent_change[vertex1] = -1;
	}
}

__device__ 
void copy_frontier(int *my_frontier, int *my_next_frontier){
	
	for (int i=1;i<=num_vertices1+num_vertices2;i++){
		my_frontier[i] = my_next_frontier[i];
	}
}

__device__ 
void clear_frontier(int *my_frontier, int *my_next_frontier ){
	for (int i=1;i<=num_vertices1+num_vertices2;i++){
			my_frontier[i] = 0;
			my_next_frontier[i] = 0;
	}
}
__device__
void vertex_disjoint_bfs(int binary_level, int vertex){
	// int frontier_element = vertex;
	// printf("Frontier element: %d \n", frontier_element );
	if(!d_frontier[vertex]){
		return;
	}

	// my_frontier[vertex] = 1;
	d_frontier[vertex] = 1;


	// Iterate all frontier elements
	while(frontier_element!=-1){
	
		int vertex = frontier_element;
		
		// Make this atomic
		d_visited[vertex] = true;
		
		// cout << "Frontier: " << frontier_element << endl;
		// cout << "Continuining for vertex: " << vertex << endl;
		
		bool found_path = false;
		int start_edge = d_list_ptr[vertex];
		int end_edge = d_list_ptr[vertex + 1]; 
		
		// cout << "Start-End edge " <<  start_edge << " " << end_edge  << endl;
		printf ("Start-End edge %d %d \n", start_edge, end_edge);
		for(int j=start_edge;j<end_edge;j++){
			if(found_path)
				break;


			int neighbor = d_flat_adj_list[j];

			

			int visited = atomicExch(&d_visited[neighbor], 1);

			if(!visited){
				// We want to alternate between unmatched and matched edges, otherwise we ignore
				d_visited[neighbor] = true;
				// cout << "Processing: " << vertex << " " << neighbor << endl;
				// exit(0);
				d_bfs_parent[neighbor] = vertex;

				if( binary_level==0 && get_is_matched_edge(vertex, neighbor)==0 && d_is_matched_vertex[neighbor]==1 ){
					// next_frontier.push_back(neighbor);
					d_next_frontier[neighbor] = 1;
				}

				// is_matched_vertex is implicitly true since the edge is matched
				// In level 1, we are only interested in matched edges
				else if( binary_level==1 && get_is_matched_edge(vertex, neighbor)==1 ){
					// next_frontier.push_back(neighbor);
					d_next_frontier[neighbor] = 1;
					// If I have found a path to the next level; I have to break
					// found_path = 1;
					return;
				}

				// Changing parent change only for this node
				else if(binary_level==0 && get_is_matched_edge(vertex, neighbor)==0 && d_is_matched_vertex[neighbor]==0){
					// cout << "Found a aug. path with " << neighbor << " with parent: " << vertex << endl;
					printf("Found a aug. path with %d with parent: %d \n", neighbor, vertex);
					d_is_parent_change[neighbor] = 1;
					// num_aug_paths++ ;
					// remove this return so that multiple paths can be found 
					return;
				}
			}
		}

		// frontier_element = get_frontier_element(vertex);
		// Getting next frontier element
		for(int x = vertex+1; x <=num_vertices1 + num_vertices2; x++){
			if(d_frontier[x]==1){
				frontier_element = x;
				break;
			}
		}
		break;
	}
	copy_frontier(d_frontier, d_next_frontier);
	// Only working for one level for now
	// bfs(binary_level = !binary_level);

}


__global__
void vertex_disjoint_bfs_util(){
	
	clear_visited();
	clear_bfs_parent();
	clear_is_parent_change();

	// int *my_frontier = new int[num_vertices1+num_vertices2+1];
	// int *my_next_frontier = new int[num_vertices1+num_vertices2+1];
	// clear_frontier(my_frontier, my_next_frontier );

	initialise_partner_vertex();
	//Can add fairness here

	
	// int num_aug_paths = 1000;

	int tid = blockIdx.x*1024 + threadIdx.x;
	int vertex = tid+1;
	if(vertex > num_vertices1)
		return;

	if(!d_visited[vertex] && !d_is_matched_vertex[vertex]){
		d_frontier[vertex] = 1;
		vertex_disjoint_bfs(0, vertex);
		__syncthreads();
	}

	printf("Working \n" );

	// if(num_aug_paths > 0){
	// 	update_matchings();
	// }

	// return num_aug_paths;
}



int main(){
	h_is_matched_edge = (bool *)calloc( (num_vertices1+ num_vertices2 + 1)*(num_vertices1 + num_vertices2+1), sizeof(bool));

	h_flat_adj_list = (int *)malloc(2*num_edges*sizeof(int));
	h_degree = (int *)malloc((num_vertices1+num_vertices2+1)*sizeof(int));
	h_list_ptr = (int *)malloc((num_vertices1+num_vertices2+2)*sizeof(int));
	h_list_ptr_copy = (int *)malloc((num_vertices1+num_vertices2+2)*sizeof(int));
	h_is_matched_vertex = (bool *)malloc((num_vertices1+num_vertices2+1)*sizeof(bool));
	h_partner_vertex = (int *)malloc((num_vertices1+num_vertices2+1)*sizeof(int));
	h_visited = (int *)malloc((num_vertices1+num_vertices2+1)*sizeof(int));
	h_bfs_parent = (int *)malloc((num_vertices1+num_vertices2+1)*sizeof(int));
	h_is_parent_change = (bool *)malloc((num_vertices1+num_vertices2+1)*sizeof(bool));
	h_frontier = (int *)malloc((num_vertices1+num_vertices2+1)*sizeof(int));
	h_next_frontier = (int *)malloc((num_vertices1+num_vertices2+1)*sizeof(int));

	// Add a check for null memory

	memset(h_degree, 0, num_vertices1 + num_vertices2 +1);
	// memset(h_is_matched_edge, 0, (num_vertices1 + num_vertices2 +1)*(num_vertices1+num_vertices2+1));
	memset(h_is_matched_vertex, 0, num_vertices1 + num_vertices2 +1);
	memset(h_visited, 0, num_vertices1 + num_vertices2 +1);
	memset(h_is_parent_change, 0, num_vertices1 + num_vertices2 +1);
	memset(h_frontier, 0, num_vertices1 + num_vertices2 +1);
	memset(h_next_frontier, 0, num_vertices1 + num_vertices2 +1);



	// to and from of edges
	// int h_edges_u[num_edges], h_edges_v[num_edges];			// Make this dynamic memory and free it once we have our 2 pass initialisation phase
	int *h_edges_u, *h_edges_v;
	h_edges_u = (int *)malloc((num_edges)*sizeof(int));
	h_edges_v = (int *)malloc((num_edges)*sizeof(int));


	ifstream fin;
    fin.open("FC_" + to_string(fc) + "_" + to_string(fc) + ".txt", ios::in);
    // fin.open("random_" + to_string(num_vertices1) + "_" + to_string(num_vertices2) + ".txt", ios::in);
    int u, v;

    // cout << "Printing all the edges: \n";

    // Vertices with 0 edges are implicitly ignored while reading the file itself
    for(int i=0;i<num_edges;i++){
    		// cout << i << endl;
            fin >> u >> v;
            h_edges_u[i] = u;
            h_edges_v[i] = v;
            h_degree[u]++;
            h_degree[v]++;
    }

    cout << "Done reading edges" << endl;

    // Get pointer to adjacency list using prefix sum (no opti here since other parts are more complex anyway)
    // Index 0 will never be used.... the last elem
    h_list_ptr[1] = 0;
    h_list_ptr_copy[1] = h_list_ptr[1];
    for(int i=2;i<=num_vertices1+num_vertices2;i++){
    	h_list_ptr[i] = h_list_ptr[i-1] + h_degree[i-1];
    	h_list_ptr_copy[i] = h_list_ptr[i];
    }
    h_list_ptr[num_vertices1+num_vertices2+1] = 2*num_edges;       //For easy coding
    h_list_ptr_copy[num_vertices1+num_vertices2+1] = 2*num_edges;  // list_ptr has the start of the adj list ; list_ptr_copy has the current position

    
    for(int i=0;i<num_edges;i++){
    	h_flat_adj_list[h_list_ptr_copy[h_edges_u[i]]] = h_edges_v[i];
    	h_flat_adj_list[h_list_ptr_copy[h_edges_v[i]]] = h_edges_u[i];
    	h_list_ptr_copy[h_edges_u[i]]++;
    	h_list_ptr_copy[h_edges_v[i]]++;
    }
    

	cudaMemcpyToSymbol(d_is_matched_edge, h_is_matched_edge, (num_vertices1+ num_vertices2 + 1)*(num_vertices1 + num_vertices2+1)*sizeof(int),0,cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_flat_adj_list, h_flat_adj_list, 2*num_edges*sizeof(int)*sizeof(int),0,cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_degree, h_degree, (num_vertices1+num_vertices2+1)*sizeof(int),0,cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_list_ptr, h_list_ptr, (num_vertices1+num_vertices2+2)*sizeof(int),0,cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_is_matched_vertex, h_is_matched_vertex, (num_vertices1+num_vertices2+1)*sizeof(int),0,cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_visited, h_visited, (num_vertices1+num_vertices2+1)*sizeof(int),0,cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_frontier, h_frontier, (num_vertices1+num_vertices2+2)*sizeof(int),0,cudaMemcpyHostToDevice);
	



    // for(int i=1;i<=num_vertices1+num_vertices2;i++){
    // 	for(int j=1;j<=num_vertices1+num_vertices2+1;j++){
    // 		h_is_matched_edge[j*num_vertices2 + i] = 0;
    // 	}
    // }
    // sleep(20000);

    // initialise_partner_vertex();
    // cout << "Partner vertex initialized " << endl;
  	

  	// for(int i=1;i<=num_vertices1+num_vertices2;i++){
  	// 	cout << h_degree[i] << " ";
  	// }

  	// for(int i=0;i<2*num_edges;i++){
  	// 	cout << h_flat_adj_list[i] << " ";
  	// }
  	// cout << endl;
   //  for(int i=0;i<=num_vertices1+num_vertices2;i++){
  	// 	cout << h_list_ptr[i] << " ";
  	// // }

    // cout << " ------------------------" <<endl;
    // for(int i=1;i<=num_vertices1;i++){
    // 	for(int j=1;j<=num_vertices2;j++){
    // 		get_is_matched_edge(i,j);
    // 	}
    // }


    // cout << get_frontier_element(9265);
  	// int x = check_matching();

  	cout << "Matching checked " << endl;

  	vertex_disjoint_bfs_util<<<1, num_threads>>>();
  	cudaDeviceSynchronize();
  


}