#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<fstream>
#include<vector>

using namespace std;

// #define num_threads 50
// #define num_edges 25
// #define num_vertices1 5
// #define num_vertices2 5

// #define num_edges 700000
// #define num_vertices1 10000
// #define num_vertices2 10000

// #define int long long

// #define num_edges 2998468
// #define num_vertices1 100000
// #define num_vertices2 100000

#define num_edges 291
#define num_vertices1 100
#define num_vertices2 100

int fc = num_vertices1;
	

vector<int> adj_list[num_vertices1 + num_vertices2 + 1];			// Do we need this?
vector<bool> is_matched_edge[num_vertices1 + num_vertices2 + 1];    // Adjacency matrix with boolean indicators
bool is_matched_vertex[num_vertices1 + num_vertices2 + 1] = {0};	// Is the vertex matched
int partner_vertex[num_vertices1 + num_vertices2 + 1];				// Get the vertex with which this vertex is matched. Initialised as -1  

int visited[num_vertices1+num_vertices2+1] = {0} ;			// Visited array for each vertex
int bfs_parent[num_vertices1+num_vertices2+1] ;				// Parent of the vertex. Required to find the augmenting path
int is_parent_change[num_vertices1+num_vertices2+1] = {0};	// Denotes if the parent changed in the last round
int num_aug_paths = 0;										// Counts number of augmenting paths found

// Only required for results
// int matched_vertices[num_vertices1+num_vertices2+1]={0};
// int matched_edges[2*num_edges]={0};

vector<int> frontier;
// int aug_path_end = -1;

// Checks if the matching is correct and also returns the total number of vertices matched
int check_matching(){
	int total_matched = 0;
	for(int i=1;i<=num_vertices1+num_vertices2;i++){
		int vertex = i;
		int num_matched = 0;
		for(int j=0;j<adj_list[vertex].size();j++){
			int neighbor = adj_list[vertex][j];
			if(is_matched_edge[vertex][neighbor]){
				// cout << vertex << " " << neighbor <<endl;
				num_matched++;
			}
		}
		if(num_matched==1){
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

void clear_visited(){
	for(int i=1;i<=num_vertices1+num_vertices2;i++){
		visited[i] = 0;
	}
}

void clear_bfs_parent(){
	for(int i=1;i<=num_vertices1+num_vertices2;i++){
		bfs_parent[i] = i;
	}
}

void initialise_partner_vertex(){
	for(int i=1;i<=num_vertices1+num_vertices2;i++){
		partner_vertex[i] = -1;
	}
}

void clear_is_parent_change(){
	for(int i=1;i<=num_vertices1+num_vertices2;i++){
		is_parent_change[i] = 0;
	}
}

void print_matchings(){
	cout << "Matchings: " << endl;
    for(int i=1;i<=num_vertices1+num_vertices2; i++){
    	cout<< i << " " << partner_vertex[i] << endl;
    }
}

void match_edges(int u, int v){
	is_matched_edge[u][v] = 1;
	is_matched_edge[v][u] = 1;
	is_matched_vertex[u] = 1;
	is_matched_vertex[v] = 1;
	partner_vertex[u] = v;
	partner_vertex[v] = u;
}

// Unmatching edges also unmatches the vertices since the graph is a matching
void unmatch_edges(int u, int v){
	is_matched_edge[u][v] = 0;
	is_matched_edge[v][u] = 0;
	if(partner_vertex[u]==v){
		is_matched_vertex[u] = 0;
		partner_vertex[u] = -1;
	}
	if(partner_vertex[v]==u){
		is_matched_vertex[v] = 0;
		partner_vertex[v] = -1;
	}
}

void debug(){
	match_edges(3,4);
	// match_edges(3,5);
	// match_edges(3,8);
	// // match_edges(3,9);

	// match_edges(10,5);
	// // match_edges(5,11);
	// match_edges(11,6);
}


void update_matchings(){
	for(int i=1;i<=num_vertices1+num_vertices2;i++){
		int vertex = i;
		if(is_parent_change[vertex] == true){
			
			// cout << "Found aug. path till " << vertex << endl;
			// There should always be odd number of vertices in aug. path
			int path_length = 1;
			int parent = bfs_parent[vertex];
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
				vertex =  bfs_parent[vertex];
				parent = bfs_parent[vertex];
				path_length++;
				// cout << vertex << " " << parent << endl;
			}
			
			// cout << ". The path length is: " << path_length << endl;
			// break;
		}

		// return here to stop after updating only one path : Important for experiments
	}
}


// Wanna alternate between unmatched and matched EDGES
// All VERTICES except the first and the last should be matched
void bfs(bool binary_level){

	// cout << "Size of frontier: " << frontier.size() << endl;
	vector<int> next_frontier;
	if(not frontier.empty()){
		// cout << "Frontier elements: " ;
		// for(int i=0;i<frontier.size();i++){
		// 	cout << frontier[i] << " ";
		// }
		// cout << endl;
		for(int i=0;i<frontier.size();i++){

			int vertex = frontier[i];
			visited[vertex] = true;
		
			// cout << "Continuining for vertex: " << vertex << endl;
			bool found_path = false;
			for(int j=0;j<adj_list[vertex].size();j++){
				if(found_path)
					break;

				int neighbor = adj_list[vertex][j];

				if(!visited[neighbor]){
					// We want to alternate between unmatched and matched edges, otherwise we ignore
					visited[neighbor] = true;
					// cout << vertex << " " << neighbor << endl;

					bfs_parent[neighbor] = vertex;

					if( binary_level==0 && is_matched_edge[vertex][neighbor]==0 && is_matched_vertex[neighbor]==1 ){
						next_frontier.push_back(neighbor);
					}

					// is_matched_vertex is implicitly true since the edge is matched
					// In level 1, we are only interested in matched edges
					else if( binary_level==1 && is_matched_edge[vertex][neighbor]==1 ){
						next_frontier.push_back(neighbor);
						// If I have found a path to the next level; I have to break
						// found_path = 1;
						// return;
					}

					// Changing parent change only for this node
					else if(binary_level==0 && is_matched_edge[vertex][neighbor]==0 && is_matched_vertex[neighbor]==0){
						// cout << "Found a aug. path with " << neighbor << " with parent: " << vertex << endl;
						is_parent_change[neighbor] = 1;
						num_aug_paths++ ;
						// remove this return so that multiple paths can be found 
						return;
					}
				}
			}
		}
		frontier.clear();
		frontier.assign(next_frontier.begin(), next_frontier.end());
		// cout << endl;
		bfs(binary_level = !binary_level);
	}
	
}

int bfs_util(){
	clear_visited();
	clear_bfs_parent();
	clear_is_parent_change();
	frontier.clear();

	//Can add fairness here
	
	// Debug
	

	num_aug_paths = 0;
	// Special style bfs
	for(int i=1;i<=num_vertices1;i++){
		if(!visited[i] && !is_matched_vertex[i]){
			frontier.clear();
			frontier.push_back(i);
			bfs(0);
		}
		// break;	
	}

	// cout << "Printing parents: " << endl;
	// for(int i=1;i<=num_vertices2+num_vertices1;i++){
	// 	cout << i << " " << bfs_parent[i] <<endl;
	// }

	if(num_aug_paths > 0){
		update_matchings();
	}

	return num_aug_paths;

}


int main(){
	ifstream fin;
    // fin.open("FC_" + to_string(fc) + "_" + to_string(fc) + ".txt", ios::in);

    fin.open("random_" + to_string(num_vertices1) + "_" + to_string(num_vertices2) + ".txt", ios::in);


    cout << "Check the file being read";
    int u, v;

    

    // Vertices with 0 edges are implicitly ignored while reading the file itself
    // Reading edges ; Format is assumed to be u-v where u is from set 1 and v is from set 2.
    for(int i=0;i<num_edges;i++){
            fin >> u >> v;
            // cout << u << " " << v <<endl;
            adj_list[u].push_back(v);
            adj_list[v].push_back(u);
    }
    cout << "Done reading edges \n";

   

    for(int i=1;i<=num_vertices1+num_vertices2;i++){
    	for(int j=1;j<=num_vertices1+num_vertices2+1;j++){
    		is_matched_edge[i].push_back(0);
    	}
    }


    initialise_partner_vertex();
   
    // debug();

    int x = check_matching();
    cout << "Total matches before running code: " << x/2 << endl;
    
    // exit(0);

    int aug_paths = bfs_util();
    cout << "Main : Number of augmenting paths " << aug_paths << endl;
    // print_matchings();

    while(aug_paths>0)
    {	
    	aug_paths = bfs_util();
    	cout << "Main : Number of augmenting paths " << aug_paths << endl;
    	// print_matchings();
    	break;
    }

    // print_matchings();
    x = check_matching();
    cout << "Total matches " << x << endl;

    // cout << "Matchings: " << endl;
    // for(int i=1;i<=num_vertices1+num_vertices2; i++){
    // 	cout<< i << " " << partner_vertex[i] << endl;
    // }

    

    return 0;
}