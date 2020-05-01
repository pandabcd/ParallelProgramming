#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<fstream>
#include<vector>

using namespace std;

// #define num_threads 50
#define num_edges 9
#define num_vertices1 3
#define num_vertices2 3


int fc = num_vertices1;
	

vector<int> adj_list[num_vertices1 + num_vertices2 + 1];			// Do we need this?
vector<bool> is_matched_edge[num_vertices1 + num_vertices2 + 1];    // Adjacency matrix with boolean indicators
bool is_matched_vertex[num_vertices1 + num_vertices2 + 1] = {0};	// Is the vertex matched

int visited[num_vertices1+num_vertices2+1] = {0} ;			// Visited array for each vertex
int bfs_parent[num_vertices1+num_vertices2+1] ;				// Parent of the vertex. Required to find the augmenting path
int is_parent_change[num_vertices1+num_vertices2+1] = {0};	// Denotes if the parent changed in the last round
int num_aug_paths = 0;										// Counts number of augmenting paths found

// Only required for results
// int matched_vertices[num_vertices1+num_vertices2+1]={0};
// int matched_edges[2*num_edges]={0};

vector<int> frontier;
int aug_path_end = -1;

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

void clear_is_parent_change(){
	for(int i=1;i<=num_vertices1+num_vertices2;i++){
		is_parent_change[i] = 0;
	}
}

void match_edges(int u, int v){
	is_matched_edge[u][v] = 1;
	is_matched_edge[v][u] = 1;
	is_matched_vertex[u] = 1;
	is_matched_vertex[v] = 1;
}

// Unmatching edges also unmatches the vertices since the graph is a matching
void unmatch_edges(int u, int v){
	is_matched_edge[u][v] = 0;
	is_matched_edge[v][u] = 0;
	is_matched_vertex[u] = 0;
	is_matched_vertex[v] = 0;
}

void debug(){
	match_edges(2,4);
	match_edges(3,5);
	// match_edges(3,5);
}


void update_matchings(){
	for(int i=1;i<=num_vertices1+num_vertices2;i++){
		int vertex = i;
		if(is_parent_change[vertex] == true){
			
			cout << "Found aug. path till " << vertex << endl;
			// There should always be odd number of vertices in aug. path
			int path_length = 1;
			int parent = bfs_parent[vertex];
			while(parent!=vertex){
				if(path_length%2==1){
					match_edges(vertex, parent);
					cout << "Matching " << vertex <<  " and " << parent << endl; 
				}
				else{
					unmatch_edges(vertex, parent);
					
					cout << "Unmatching " << vertex <<  " and " << parent << endl;
				}
				vertex =  bfs_parent[vertex];
				parent = bfs_parent[vertex];
				path_length++;
				cout << vertex << endl;
			}
			
			cout << ". The path length is: " << path_length << endl;
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
			for(int j=0;j<adj_list[vertex].size();j++){
				int neighbor = adj_list[vertex][j];


				if(!visited[neighbor]){
					// We want to alternate between unmatched and matched edges, otherwise we ignore
					visited[neighbor] = true;
					bfs_parent[neighbor] = vertex;

					if( binary_level==0 && is_matched_edge[vertex][neighbor]==0 && is_matched_vertex[neighbor]==1 ){
						next_frontier.push_back(neighbor);
					}

					// is_matched_vertex is implicitly true since the edge is matched
					// In level 1, we are only interested in matched edges
					else if( binary_level==1 && is_matched_edge[vertex][neighbor]==1 ){
						next_frontier.push_back(neighbor);
					}

					// Changing parent change only for this node
					else if(binary_level==0 && is_matched_edge[vertex][neighbor]==0 && is_matched_vertex[neighbor]==0){
						cout << "Found a aug. path with " << neighbor << " with parent: " << vertex << endl;
						is_parent_change[neighbor] = 1;
						num_aug_paths++ ;
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
	debug();

	num_aug_paths = 0;
	// Special style bfs
	for(int i=1;i<=num_vertices1;i++){
		if(!visited[i]){
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
    fin.open("FC_" + to_string(fc) + "_" + to_string(fc) + ".txt", ios::in);
    int u, v;

    cout << "Printing all the edges: \n";

    // Vertices with 0 edges are implicitly ignored while reading the file itself
    for(int i=0;i<num_edges;i++){
            fin >> u >> v;
            cout << u << " " << v <<endl;
            adj_list[u].push_back(v);
            adj_list[v].push_back(u);

            is_matched_edge[u].push_back(0);
            is_matched_edge[v].push_back(0);
    }

    int aug_paths = bfs_util();
    cout << "Number of augmenting paths " << aug_paths << endl;
    // while(aug_paths>0)
    // {	
    // 	aug_paths = bfs_util();
    // 	cout << "Number of augmenting paths " << aug_paths << endl;
    // 	break;
    // }

    cout << "BFS parents: " << endl;
    for(int i=1;i<=num_vertices1+num_vertices2; i++){
    	cout<< i << bfs_parent[i] << endl;
    }


}