#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fstream>
#include<time.h>
#include<sys/time.h>
#include<string>

using namespace std;

// #define num_threads 1000

#define lli long long int

int size[5] = {100, 500, 1000, 5000, 10000};

// int edges_2[5] = {200, 447, 1969, 4991, 200001};
int edges_8[5] = {801, 20000, 79580, 1999218, 8000000};


const lli num_edges = 8000000;
const lli num_vertices1 = 10000;
const lli num_vertices2 = 10000;


__device__ int d_flat_adj_list[2*num_edges];
__device__ int d_degree[num_vertices1+num_vertices2+1]={0};      //store degree of each vertex
__device__ int d_list_ptr[num_vertices1+num_vertices2+2];        //1-indexed and extra element at the end for easy size access  // Pointer to the start of adjacency list
__device__ int d_list_ptr_copy[num_vertices1+num_vertices2+2];    // 


__device__ bool d_matched_edge[2*num_edges];					// Tells for every edge in the list if the edge is matched or not
__device__ bool d_is_matched_vertex[num_vertices1 + num_vertices2 + 1] = {0};	//is the vertex matched
__device__ int d_partner_vertex[num_vertices1 + num_vertices2 + 1];
__device__ int d_visited[num_vertices1 + num_vertices2 + 1] = {0};
__device__ int d_bfs_parent[num_vertices1 +  num_vertices2 + 1];
__device__ bool d_is_parent_change[num_vertices1 + num_vertices2 + 1] = {0};

__device__ int d_frontier[num_vertices1 + num_vertices2+1] = {0};
__device__ int d_next_frontier[num_vertices1+num_vertices2+1] = {0};

__device__ int num_aug_paths = 10000000;						//Any number not equal to 0 works



int *h_flat_adj_list;
int *h_degree;
int * h_list_ptr;
int *h_list_ptr_copy;

bool *h_matched_edge;
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
bool get_matched_edge(int x, int y){
	int vertex = x;
	int start_edge = d_list_ptr[vertex];
	int end_edge = d_list_ptr[vertex + 1]; 
	for(int i = start_edge; i<end_edge;i++){
		if(d_flat_adj_list[i]==y){
			return d_matched_edge[i];
		}
	}
	printf("Error! Querying for an edge which is not present \n");
	return -1;
}

__device__
void set_matched_edge(int x, int y, int value){
	bool edge_present = false;
	int vertex = x;
	int start_edge = d_list_ptr[vertex];
	int end_edge = d_list_ptr[vertex + 1]; 
	for(int i = start_edge; i<end_edge;i++){
		if(d_flat_adj_list[i] == y){
			d_matched_edge[i] = value;
			edge_present = true;
			break;
		}
	}
	
	vertex = y;
	start_edge = d_list_ptr[vertex];
	end_edge = d_list_ptr[vertex + 1]; 
	for(int i = start_edge; i<end_edge;i++){
		if(d_flat_adj_list[i] == x){
			d_matched_edge[i] = value;
			edge_present = true;
			break;
		}
	}

	if(!edge_present){
		printf("Error! Querying for an edge which is not present \n");
	}
}

void print_matchings(){
	cout << "Matchings: " << endl;
    for(int i=1;i<=num_vertices1+num_vertices2; i++){
    	cout<< i << " " << h_partner_vertex[i] << endl;
    }
}

int get_matched_edge_h(int x, int y){
	int vertex = x;
	int start_edge = h_list_ptr[vertex];
	int end_edge = h_list_ptr[vertex + 1]; 
	for(int i = start_edge; i<end_edge;i++){
		if(h_flat_adj_list[i] == y){
			return h_matched_edge[i];
		}
	}
	cout << "Error! Querying for an edge which is not present";
	exit(0);
}


__device__
void match_edges(int u, int v){
	set_matched_edge(u,v,1);
	set_matched_edge(v,u,1);
	d_is_matched_vertex[u] = 1;
	d_is_matched_vertex[v] = 1;
	d_partner_vertex[u] = v;
	d_partner_vertex[v] = u;

}

// Unmatching edges also unmatches the vertices since the graph is a matching
__device__
void unmatch_edges(int u, int v){
	set_matched_edge(u,v,0);
	set_matched_edge(v,u,0);
	if(d_partner_vertex[u]==v){
		d_is_matched_vertex[u] = 0;
		d_partner_vertex[u] = -1;
	}
	if(d_partner_vertex[v]==u){
		d_is_matched_vertex[v] = 0;
		d_partner_vertex[v] = -1;
	}
}

// Make this parallel
__global__
void update_matchings(){
	int tid = blockIdx.x*1024 + threadIdx.x;
	for(int i=tid; i<=num_vertices1+num_vertices2; i+=num_vertices1){
		int vertex = i;
		if(d_is_parent_change[vertex] == true){
			
			// cout << "Found aug. path till " << vertex << endl;
			// There should always be odd number of vertices in aug. path
			int path_length = 1;
			int parent = d_bfs_parent[vertex];
			while(parent!=vertex){
				// cout << vertex << " " <<parent << endl;
				if(path_length%2==1){
					match_edges(vertex, parent);
					// printf("Matching %d and %d \n", vertex, parent);
				}
				else{
					unmatch_edges(vertex, parent);
					// printf("Unmatching %d and %d \n", vertex, parent);
				}
				vertex =  d_bfs_parent[vertex];
				parent = d_bfs_parent[vertex];
				path_length++;
			}
		}
	}
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
void vertex_disjoint_bfs(int binary_level, int vertex, int tid){
	
		int visited_self = atomicExch(&d_visited[vertex], 1);
		if(visited_self && binary_level==0){
			return;
		}
		d_visited[vertex] = true;
		
		
		bool found_path = false;
		int start_edge = d_list_ptr[vertex];
		int end_edge = d_list_ptr[vertex + 1]; 
		
		for(int j=start_edge;j<end_edge;j++){
			if(found_path)
				break;


			int neighbor = d_flat_adj_list[j];


			if(neighbor > num_vertices1 + num_vertices2){
				printf("[%d]Error(neighbor out of range: vertex, neighbor : %d, %d \n", tid, vertex, neighbor);
			}
			
			int visited = atomicExch(&d_visited[neighbor], 1);

			if(!visited){
				// We want to alternate between unmatched and matched edges, otherwise we ignore
				d_visited[neighbor] = true;
				d_bfs_parent[neighbor] = vertex;

				if( binary_level==0 && get_matched_edge(vertex, neighbor)==0 && d_is_matched_vertex[neighbor]==1 ){
						d_next_frontier[neighbor] = 1;
						if(binary_level==1)
							printf("Going odd %d \n", vertex);
						vertex_disjoint_bfs(!binary_level, neighbor, tid);
				}

				// In level 1, we are only interested in matched edges
				else if( binary_level==1 && get_matched_edge(vertex, neighbor)==1 ){
					d_next_frontier[neighbor] = 1;
					vertex_disjoint_bfs(!binary_level, neighbor, tid);
					return;
				}

				// Changing parent change only for this node
				else if(binary_level==0 && get_matched_edge(vertex, neighbor)==0 && d_is_matched_vertex[neighbor]==0){
					d_is_parent_change[neighbor] = 1;
					// atomicAdd(&num_aug_paths, 1);
					num_aug_paths++;
					return;
				}
			}
	}
}


__global__
void vertex_disjoint_bfs_util(){

	// parallelise these functions
	clear_visited();
	clear_bfs_parent();
	clear_is_parent_change();

	// clear_frontier(my_frontier, my_next_frontier );
	initialise_partner_vertex();
	int tid = blockIdx.x*1024 + threadIdx.x;
	int vertex = tid+1;
	if(vertex > num_vertices1)
		return;

	if(vertex >=  num_vertices1+num_vertices2+1)
		printf("[%d] Error \n");

	if(!d_visited[vertex] && !d_is_matched_vertex[vertex]){
		d_frontier[vertex] = 1;
		vertex_disjoint_bfs(0, vertex, tid);
	}

}

int check_matching(){
	int total_matched = 0;
	for(int i=1;i<=num_vertices1+num_vertices2;i++){
		int vertex = i;
		int num_matched = 0;


		for(int j=h_list_ptr[i];j<h_list_ptr[i+1];j++){
			int neighbor = h_flat_adj_list[j];
			// cout << vertex << " " << neighbor << endl;
			if(get_matched_edge_h(vertex, neighbor)){
				num_matched++;
			}
		}


		if(num_matched==1){
			// cout << "Hi" << endl;
			total_matched++;
		}
		if(num_matched>1){
			cout << vertex << endl;
			cout << "Error! Not a matching!";
			exit(0);
		}
	}
	return total_matched/2;
}


int main(){

	struct timespec start, end;

	// h_is_matched_edge = (bool *)calloc( (num_vertices1+ num_vertices2 + 1)*(num_vertices1 + num_vertices2+1), sizeof(bool));

	h_matched_edge = (bool *)calloc(2*num_edges, sizeof(bool));
	
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
	fin.open("random_10000_10000_high.txt", ios::in);
    int u, v;

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
    

    clock_gettime( CLOCK_REALTIME,&start);

	cudaMemcpyToSymbol(d_matched_edge, h_matched_edge, (2*num_edges)*sizeof(int),0,cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_flat_adj_list, h_flat_adj_list, 2*num_edges*sizeof(int),0,cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_degree, h_degree, (num_vertices1+num_vertices2+1)*sizeof(int),0,cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_list_ptr, h_list_ptr, (num_vertices1+num_vertices2+2)*sizeof(int),0,cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_is_matched_vertex, h_is_matched_vertex, (num_vertices1+num_vertices2+1)*sizeof(int),0,cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_visited, h_visited, (num_vertices1+num_vertices2+1)*sizeof(int),0,cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_frontier, h_frontier, (num_vertices1+num_vertices2+2)*sizeof(int),0,cudaMemcpyHostToDevice);

	int h_num_aug_paths = 1000;
	
  	cudaDeviceSynchronize();

  	while(h_num_aug_paths>0){
  		h_num_aug_paths = 0;
  		cudaMemcpyToSymbol(num_aug_paths, &h_num_aug_paths, (1)*sizeof(int),0,cudaMemcpyHostToDevice);
	  	
	  	vertex_disjoint_bfs_util<<<10, 1024>>>();
	  	update_matchings<<<10, 1024>>>();
	  	cudaDeviceSynchronize();	  	
	  	cudaMemcpyFromSymbol(&h_num_aug_paths, num_aug_paths, sizeof(num_aug_paths),0,cudaMemcpyDeviceToHost);
	  	
	  	break;
	}
  	clock_gettime( CLOCK_REALTIME,&end);
  	cudaMemcpyFromSymbol(h_matched_edge, d_matched_edge, sizeof(d_matched_edge),0,cudaMemcpyDeviceToHost);
  	cudaMemcpyFromSymbol(h_partner_vertex, d_partner_vertex, sizeof(d_partner_vertex),0,cudaMemcpyDeviceToHost);
  	
  printf("Number of augmenting paths(actual number may be higher): %d \n", h_num_aug_paths);

  	int num_matches = check_matching();
  	
  	printf("Number of matchings: %d \n", num_matches);

  	double elapsed = (end.tv_sec-start.tv_sec)*1000000000 + end.tv_nsec-start.tv_nsec;
  	printf("Time elapsed %lf\n", elapsed/1e6);
	

  	cudaDeviceSynchronize();
  
}