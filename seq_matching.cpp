#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<fstream>
#include<vector>

using namespace std;

#define num_threads 50
#define num_edges 25
#define num_vertices1 5
#define num_vertices2 5


int fc = num_vertices1;
	

vector<int> adj_list[num_vertices1 + num_vertices2 + 1];
int visited[num_vertices1+num_vertices2+1];
int bfs_parent[num_vertices1+num_vertices2+1];

// Only required for results
// int matched_vertices[num_vertices1+num_vertices2+1]={0};
// int matched_edges[2*num_edges]={0};

vector<int> frontier;

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

void bfs(){
	cout << "Size of frontier: " << frontier.size() << endl;
	vector<int> next_frontier;
	if(not frontier.empty()){
		for(int i=0;i<frontier.size();i++){

			int vertex = frontier[i];
			cout << vertex <<endl;
			visited[vertex] = true;

			for(int j=0;j<adj_list[vertex].size();j++){
				int neighbor = adj_list[vertex][j];

				if(!visited[neighbor]){
					visited[neighbor] = true;
					bfs_parent[neighbor] = vertex;
					next_frontier.push_back(neighbor);
				}
			}
		}
		frontier.clear();
		frontier.assign(next_frontier.begin(), next_frontier.end());
		bfs();
	}
	
}

void bfs_util(){
	clear_visited();
	clear_bfs_parent();
	frontier.clear();
	//Can add fairness here
	frontier.push_back(1);
	// cout << "good" << endl;
	bfs();
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
    }

    bfs_util();

    for(int i=0;i<num_vertices1+num_vertices2; i++){
    	cout<< i << bfs_parent[i] << endl;
    }


}