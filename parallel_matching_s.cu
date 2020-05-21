#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<fstream>
#include<vector>

using namespace std;

// #define num_threads 50
// #define num_edges 700000
// #define num_vertices1 10000
// #define num_vertices2 10000

// #define num_edges 1000000
// #define num_vertices1 1000
// #define num_vertices2 1000


#define num_edges 2998468
// #define num_vertices1 100000
// #define num_vertices2 100000

const long long num_vertices1 = 100000;
const long long num_vertices2 = 100000;

#define long long int lli
// vector<int> adj_list[num_vertices1 + num_vertices2 + 1];			// Do we need this? YES
	// vector<bool> is_matched_edge[num_vertices1 + num_vertices2 + 1];    // Adjacency matrix with boolean indicators
	// bool is_matched_vertex[num_vertices1 + num_vertices2 + 1] = {0};	// Is the vertex matched
	// int partner_vertex[num_vertices1 + num_vertices2 + 1];				// Get the vertex with which this vertex is matched. Initialised as -1  

	// int visited[num_vertices1+num_vertices2+1] = {0} ;			// Visited array for each vertex
	// int bfs_parent[num_vertices1+num_vertices2+1] ;				// Parent of the vertex. Required to find the augmenting path
	// int is_parent_change[num_vertices1+num_vertices2+1] = {0};	// Denotes if the parent changed in the last round
	// int num_aug_paths = 0;										// Counts number of augmenting paths found
int h_fc = num_vertices1;

int h_flat_adj_list[2*num_edges];
int h_degree[num_vertices1+num_vertices2+1]={0};      //store degree of each vertex
int h_list_ptr[num_vertices1+num_vertices2+2];        //1-indexed and extra element at the end for easy size access  // Pointer to the start of adjacency list
int h_list_ptr_copy[num_vertices1+num_vertices2+2];    // Temporrary stuff, gotta sleep

// bool h_is_matched_edge[(num_vertices1+ num_vertices2 + 1)*(num_vertices1 + num_vertices2+1)] = {0} ;     // Adjacency matrix (1-indexed)

bool *h_is_matched_edge;
bool h_is_matched_vertex[num_vertices1 + num_vertices2 + 1] = {0};	//is the vertex matched
int h_partner_vertex[num_vertices1 + num_vertices2 + 1];
int h_visited[num_vertices1 + num_vertices2 + 1] = {0};
int h_bfs_parent[num_vertices1 +  num_vertices2 + 1];
bool h_is_parent_change[num_vertices1 + num_vertices2 + 1] = {0};

int fc = num_vertices1;
int num_aug_paths = 0;


// Only required for results
// int matched_vertices[num_vertices1+num_vertices2+1]={0};
// int matched_edges[2*num_edges]={0};

// vector<int> frontier;
// int aug_path_end = -1;
int frontier[num_vertices1 + num_vertices2] = {0};
int next_frontier[num_vertices1+num_vertices2] = {0};


int get_is_matched_edge(int i, int j){
	// cout << i << " " << j << endl;
	return h_is_matched_edge[i*(num_vertices1 + num_vertices2+1) + j ];
}

void set_is_matched_edge(int i, int j, int value){
	h_is_matched_edge[i*(num_vertices1 + num_vertices2+1) + j ] = value;
}
// Checks if the matching is correct and also returns the total number of vertices matched
int check_matching(){
	int total_matched = 0;
	for(int i=1;i<=num_vertices1+num_vertices2;i++){
		int vertex = i;
		int num_matched = 0;
		int start_edge = h_list_ptr[vertex];
		int end_edge = h_list_ptr[vertex+1];

		for(int j=start_edge;j<end_edge;j++){

			int neighbor = h_flat_adj_list[j];
			// cout << "vertex-neighbor " << vertex << " " <<neighbor <<endl;
			if(get_is_matched_edge(vertex, neighbor)){
				// cout << "Matched" << endl;
				// cout << vertex << " " << neighbor <<endl;
				num_matched++;
			}
		}
		if(num_matched==1){
			total_matched++;
		}
		if(num_matched>1){
			// cout << vertex << endl;
			// cout << "Error! Not a matching!";
			// exit(0);
		}
	}
	cout << "Matching is correct! " << endl;
	return total_matched/2;
}

void clear_visited(){
	for(int i=1;i<=num_vertices1+num_vertices2;i++){
		h_visited[i] = 0;
	}
}

void clear_bfs_parent(){
	for(int i=1;i<=num_vertices1+num_vertices2;i++){
		h_bfs_parent[i] = i;
	}
}

void initialise_partner_vertex(){
	for(int i=1;i<=num_vertices1+num_vertices2;i++){
		h_partner_vertex[i] = -1;
	}
}

void clear_is_parent_change(){
	for(int i=1;i<=num_vertices1+num_vertices2;i++){
		h_is_parent_change[i] = 0;
	}
}





void print_matchings(){
	cout << "Matchings: " << endl;
    for(int i=1;i<=num_vertices1+num_vertices2; i++){
    	cout<< i << " " << h_partner_vertex[i] << endl;
    }
}

void match_edges(int u, int v){
	// h_is_matched_edge[u][v] = 1;
	// h_is_matched_edge[v][u] = 1;
	// cout << "Matching " << u << " " << v << endl;	

	set_is_matched_edge(u,v,1);
	set_is_matched_edge(v,u,1);
	h_is_matched_vertex[u] = 1;
	h_is_matched_vertex[v] = 1;
	h_partner_vertex[u] = v;
	h_partner_vertex[v] = u;
}

// Unmatching edges also unmatches the vertices since the graph is a matching
void unmatch_edges(int u, int v){
	// h_is_matched_edge[u][v] = 0;
	// h_is_matched_edge[v][u] = 0;
	

	// cout << "UnMatching " << u << " " << v << endl;
	set_is_matched_edge(u,v,0);
	set_is_matched_edge(v,u,0);
	h_is_matched_vertex[u] = 0;
	h_is_matched_vertex[v] = 0;
	h_partner_vertex[u] = -1;
	h_partner_vertex[v] = -1;
}


void update_matchings(){
	for(int i=1;i<=num_vertices1+num_vertices2;i++){
		int vertex = i;
		if(h_is_parent_change[vertex] == true){
			
			// cout << "Found aug. path till " << vertex << endl;
			// There should always be odd number of vertices in aug. path
			int path_length = 1;
			int parent = h_bfs_parent[vertex];
			while(parent!=vertex){
				// cout << vertex << " " <<parent << endl;
				if(path_length%2==1){
					match_edges(vertex, parent);
					// cout << "Matching " << vertex <<  " and " << parent << endl; 
				}
				else{
					unmatch_edges(vertex, parent);
					// cout << "Unmatching " << vertex <<  " and " << parent << endl;
				}
				vertex =  h_bfs_parent[vertex];
				parent = h_bfs_parent[vertex];
				path_length++;
				// cout << vertex << " " << parent << endl;
			}
			
			// cout << ". The path length is: " << path_length << endl;
			// break;
		}

		// return here to stop after updating only one path : Important for experiments
	}
}


int get_frontier_element(int ele){
	for(int i=ele+1;i<=num_vertices1+num_vertices2+1;i++){
		if(frontier[i]){
			return i;
		}
	}
	return -1;
}

void copy_frontier(){
	for(int i=0;i<=num_vertices1+num_vertices2;i++){
		frontier[i] = next_frontier[i];
		next_frontier[i] = 0;
	}
}

void bfs(bool binary_level){

	// vector<int> next_frontier;
	
	int frontier_element = get_frontier_element(0);
	// int frontier_element = 9265;

	// cout << "Frontier: " << frontier_element << endl;
	// if(not frontier.empty()){
	if(frontier_element!=-1){
		// for(int i=0;i<frontier.size();i++){
		// Iterate all frontier elements
		while(frontier_element!=-1){

			// int vertex = frontier[i];
			int vertex = frontier_element;
			h_visited[vertex] = true;
			
			// cout << "Frontier: " << frontier_element << endl;
			// cout << "Continuining for vertex: " << vertex << endl;
			bool found_path = false;
			int start_edge = h_list_ptr[vertex];
			int end_edge = h_list_ptr[vertex + 1]; 
			
			// cout << "Start-End edge " <<  start_edge << " " << end_edge  << endl;
			for(int j=start_edge;j<end_edge;j++){
				if(found_path)
					break;


				int neighbor = h_flat_adj_list[j];

				// cout << "Vertex- neighbor " << vertex << " " << neighbor <<endl; 

				if(!h_visited[neighbor]){
					// We want to alternate between unmatched and matched edges, otherwise we ignore
					h_visited[neighbor] = true;
					// cout << "Processing: " << vertex << " " << neighbor << endl;
					// exit(0);
					h_bfs_parent[neighbor] = vertex;

					if( binary_level==0 && get_is_matched_edge(vertex, neighbor)==0 && h_is_matched_vertex[neighbor]==1 ){
						// next_frontier.push_back(neighbor);
						next_frontier[neighbor] = 1;
					}

					// is_matched_vertex is implicitly true since the edge is matched
					// In level 1, we are only interested in matched edges
					else if( binary_level==1 && get_is_matched_edge(vertex, neighbor)==1 ){
						// next_frontier.push_back(neighbor);
						next_frontier[neighbor] = 1;
						// If I have found a path to the next level; I have to break
						// found_path = 1;
						return;
					}

					// Changing parent change only for this node
					else if(binary_level==0 && get_is_matched_edge(vertex, neighbor)==0 && h_is_matched_vertex[neighbor]==0){
						// cout << "Found a aug. path with " << neighbor << " with parent: " << vertex << endl;
						h_is_parent_change[neighbor] = 1;
						num_aug_paths++ ;
						// remove this return so that multiple paths can be found 
						return;
					}
				}
			}

			frontier_element = get_frontier_element(vertex);
		}
		// frontier.clear();
		// frontier.assign(next_frontier.begin(), next_frontier.end());
		copy_frontier();
		bfs(binary_level = !binary_level);
	}
	
}

void clear_frontier(){
	for(int i=0;i<num_vertices1+num_vertices2+1;i++){
		frontier[i] = 0;
	}
}

int bfs_util(){
	clear_visited();
	clear_bfs_parent();
	clear_is_parent_change();
	// frontier.clear();
	clear_frontier();

	//Can add fairness here

	num_aug_paths = 0;

	// Special style bfs
	for(int i=1;i<=num_vertices1;i++){
		if(!h_visited[i] && !h_is_matched_vertex[i]){
			// frontier.clear();
			clear_frontier();
			// frontier.push_back(i);
			frontier[i] = 1;
			bfs(0);
			// cout << "Loop";
		}
		// break;	
	}

	// cout << "Printing parents: " << endl;
	// for(int i=1;i<=num_vertices2+num_vertices1;i++){
	// 	cout << i << " " << h_bfs_parent[i] <<endl;
	// }

	if(num_aug_paths > 0){
		update_matchings();
	}

	return num_aug_paths;

}


int main(){
	cout << "hi " << endl;
	h_is_matched_edge = (bool *)malloc((num_vertices1+ num_vertices2 + 1)*(num_vertices1 + num_vertices2+1)*sizeof(int));
	h_flat_adj_list = (int *)malloc(2*num_edges*sizeof(int));
	h_list_ptr = (int *)malloc((num_vertices1+num_vertices2+2)*sizeof(int));
	h_list_ptr_copy = (int *)malloc((num_vertices1+num_vertices2+2)*sizeof(int));
	h_is_matched_vertex = (int *)malloc((num_vertices1+num_vertices2+1)*sizeof(int));
	h_partner_vertex = (int *)malloc((num_vertices1+num_vertices2+1)*sizeof(int));
	h_visited = (int *)malloc((num_vertices1+num_vertices2+1)*sizeof(int));
	h_bfs_parent = (int *)malloc((num_vertices1+num_vertices2+1)*sizeof(int));
	h_is_parent_change = (int *)malloc((num_vertices1+num_vertices2+1)*sizeof(int));


// int h_flat_adj_list[2*num_edges];
// int h_degree[num_vertices1+num_vertices2+1]={0};      //store degree of each vertex
// int h_list_ptr[num_vertices1+num_vertices2+2];        //1-indexed and extra element at the end for easy size access  // Pointer to the start of adjacency list
// int h_list_ptr_copy[num_vertices1+num_vertices2+2];    // Temporrary stuff, gotta sleep

// // bool h_is_matched_edge[(num_vertices1+ num_vertices2 + 1)*(num_vertices1 + num_vertices2+1)] = {0} ;     // Adjacency matrix (1-indexed)

// bool h_is_matched_vertex[num_vertices1 + num_vertices2 + 1] = {0};	//is the vertex matched
// int h_partner_vertex[num_vertices1 + num_vertices2 + 1];
// int h_visited[num_vertices1 + num_vertices2 + 1] = {0};
// int h_bfs_parent[num_vertices1 +  num_vertices2 + 1];
// bool h_is_parent_change[num_vertices1 + num_vertices2 + 1] = {0};

	

	// to and from of edges
	int h_edges_u[num_edges], h_edges_v[num_edges];			// Make this dynamic memory and free it once we have our 2 pass initialisation phase
	

	
	ifstream fin;
    // fin.open("FC_" + to_string(fc) + "_" + to_string(fc) + ".txt", ios::in);
    fin.open("random_" + to_string(num_vertices1) + "_" + to_string(num_vertices2) + ".txt", ios::in);
    int u, v;

    // cout << "Printing all the edges: \n";

    // Vertices with 0 edges are implicitly ignored while reading the file itself
    for(int i=0;i<num_edges;i++){
            fin >> u >> v;
            // cout << u << " " << v <<endl;
            h_edges_u[i] = u;
            h_edges_v[i] = v;
            h_degree[u]++;
            h_degree[v]++;
    }

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


    // for(int i=1;i<=num_vertices1+num_vertices2;i++){
    // 	for(int j=1;j<=num_vertices1+num_vertices2+1;j++){
    // 		h_is_matched_edge[j*num_vertices2 + i] = 0;
    // 	}
    // }


    initialise_partner_vertex();
  	

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
  	int x = check_matching();
    bfs_util();
    print_matchings();


    x = check_matching();
    cout << "Number of matchings: " << x << endl;








    // int x = check_matching();
    // cout << "Total matches before running code: " << x << endl;
    
    
    // int aug_paths = bfs_util();
    // cout << "Main : Number of augmenting paths " << aug_paths << endl;
    // // print_matchings();

    // while(aug_paths>0)
    // {	
    // 	aug_paths = bfs_util();
    // 	cout << "Main : Number of augmenting paths " << aug_paths << endl;
    // 	// print_matchings();
    // 	break;
    // }

    // x = check_matching();
    // cout << "Total matches " << x/2 << endl;


}