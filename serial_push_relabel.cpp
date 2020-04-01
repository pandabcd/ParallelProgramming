
// C++ program to implement push-relabel algorithm for 
// getting maximum flow of graph 
#include <bits/stdc++.h> 
using namespace std; 

  
struct Edge 
{ 
    // To store current flow and capacity of edge 
    int flow, capacity; 
  
    // An edge u--->v has start vertex as u and end 
    // vertex as v. 
    int u, v; 
  
    Edge(int flow, int capacity, int u, int v) 
    { 
        this->flow = flow; 
        this->capacity = capacity; 
        this->u = u; 
        this->v = v; 
    } 
}; 
  
// Represent a Vertex 
struct Vertex 
{ 
    int h, e_flow; 
  
    Vertex(int h, int e_flow) 
    { 
        this->h = h; 
        this->e_flow = e_flow; 
    } 
}; 
  
// To represent a flow network 
class Graph 
{ 
    int V;    // No. of vertices 
    vector<Vertex> ver; 
    vector<Edge> edge; 
  
    // Function to push excess flow from u 
    bool push(int u); 
  
    // Function to relabel a vertex u 
    void relabel(int u); 
  
    // This function is called to initialize 
    // preflow 
    void preflow(int s); 
  
    // Function to reverse edge 
    void updateReverseEdgeFlow(int i, int flow); 
  
public: 
    Graph(int V);  // Constructor 
  
    // function to add an edge to graph 
    void addEdge(int u, int v, int w); 
  
    // returns maximum flow from s to t 
    int getMaxFlow(int s, int t); 

    // Modification
    bool check_bipartite(int start, int end, int source, int sink);
}; 
  
Graph::Graph(int V) 
{ 
    this->V = V; 
  
    // all vertices are initialized with 0 height 
    // and 0 excess flow 
    for (int i = 0; i < V; i++) 
        ver.push_back(Vertex(0, 0)); 
} 
  
void Graph::addEdge(int u, int v, int capacity) 
{ 
    // flow is initialized with 0 for all edge 
    edge.push_back(Edge(0, capacity, u, v)); 
} 
  
void Graph::preflow(int s) 
{ 
    // Making h of source Vertex equal to no. of vertices 
    // Height of other vertices is 0. 
    ver[s].h = ver.size(); 
  
    // 
    for (int i = 0; i < edge.size(); i++) 
    { 
        // If current edge goes from source 
        if (edge[i].u == s) 
        { 
            // Flow is equal to capacity 
            edge[i].flow = edge[i].capacity; 
  
            // Initialize excess flow for adjacent v 
            ver[edge[i].v].e_flow += edge[i].flow; 
  
            // Add an edge from v to s in residual graph with 
            // capacity equal to 0 
            edge.push_back(Edge(-edge[i].flow, 0, edge[i].v, s)); 
        } 
    } 
} 
  
// returns index of overflowing Vertex 
int overFlowVertex(vector<Vertex>& ver, int s, int t) 
{ 
    for (int i = 0; i < ver.size(); i++) 
       if (i!=s && i!=t && ver[i].e_flow > 0) 
            return i; 
  
    // -1 if no overflowing Vertex 
    return -1; 
} 
  
// Update reverse flow for flow added on ith Edge 
void Graph::updateReverseEdgeFlow(int i, int flow) 
{ 
    int u = edge[i].v, v = edge[i].u; 
  
    for (int j = 0; j < edge.size(); j++) 
    { 
        if (edge[j].v == v && edge[j].u == u) 
        { 
            edge[j].flow -= flow; 
            return; 
        } 
    } 
  
    // adding reverse Edge in residual graph 
    Edge e = Edge(0, flow, u, v); 
    edge.push_back(e); 
} 
  
// To push flow from overflowing vertex u 
bool Graph::push(int u) 
{ 
    // Traverse through all edges to find an adjacent (of u) 
    // to which flow can be pushed 
    for (int i = 0; i < edge.size(); i++) 
    { 
        // Checks u of current edge is same as given 
        // overflowing vertex 
        if (edge[i].u == u) 
        { 
            // if flow is equal to capacity then no push 
            // is possible 
            if (edge[i].flow == edge[i].capacity) 
                continue; 
  
            // Push is only possible if height of adjacent 
            // is smaller than height of overflowing vertex 
            if (ver[u].h > ver[edge[i].v].h) 
            {   
                // Flow to be pushed is equal to minimum of 
                // remaining flow on edge and excess flow. 
                int flow = min(edge[i].capacity - edge[i].flow, 
                               ver[u].e_flow); 
  
                // Reduce excess flow for overflowing vertex 
                ver[u].e_flow -= flow; 
  
                // Increase excess flow for adjacent 
                ver[edge[i].v].e_flow += flow; 
  
                // Add residual flow (With capacity 0 and negative 
                // flow) 
                edge[i].flow += flow; 
  
                updateReverseEdgeFlow(i, flow); 

                cout << "Pushing " << flow << " from " << edge[i].u << " to " << edge[i].v <<endl;
                
  
                return true; 
            } 
        } 
    } 
    return false; 
} 
  
// function to relabel vertex u 
void Graph::relabel(int u) 
{ 
    // Initialize minimum height of an adjacent 
    int mh = INT_MAX; 
  
    // Find the adjacent with minimum height 
    for (int i = 0; i < edge.size(); i++) 
    { 
        if (edge[i].u == u) 
        { 
            // if flow is equal to capacity then no 
            // relabeling 
            if (edge[i].flow == edge[i].capacity) 
                continue; 
  
            // Update minimum height 
            if (ver[edge[i].v].h < mh) 
            { 
                mh = ver[edge[i].v].h; 
  
                // updating height of u 
                ver[u].h = mh + 1; 

                cout << "Relabling vertex: " <<u<<endl;
            } 
        } 
    } 
} 
  
// main function for printing maximum flow of graph 
int Graph::getMaxFlow(int s, int t) 
{ 
    preflow(s); 
  
    // loop untill none of the Vertex is in overflow 
    while (overFlowVertex(ver,s,t) != -1) 
    { 
        int u = overFlowVertex(ver,s,t); 
        if (!push(u)) 
            relabel(u); 
    } 
  
    // ver.back() returns last Vertex, whose 
    // e_flow will be final maximum flow 
    // return ver.back().e_flow; 
    return ver[t].e_flow;
} 

bool check_set(int start, int end, int vertex){
    return(vertex>=start && vertex<=end);
    //     return true;
    // }
    // return false;
}




bool Graph::check_bipartite(int start, int end, int source, int sink){
    //start and end are indices of vertex of a set(and are not supposed to have edges)
    for (int i = 0; i < edge.size(); i++) {
        if(edge[i].u==source || edge[i].v==sink){
            continue;
        }
        // Probably not needed
        if(edge[i].u==sink || edge[i].v==source){
            continue;
        }

        if( check_set(start, end, edge[i].u) && check_set(start, end, edge[i].v)){
            cout << "Edges " << edge[i].u << " and " <<  edge[i].v << " belong to the same set and has an edge. Invalid graph input!";
            return false;
        }
    }

    return true;

}

// Driver program to test above functions 
int main() 
{   
    int fc = 20;
    int V = 2*fc + 2; 
    Graph g(V); 
        

    int pseudo_source = 0;
    int pseudo_sink = 2*fc+1;

    ifstream fin;
    fin.open("FC_" + to_string(fc) + "_" + to_string(fc) + ".txt", ios::in);
    int u, v;

    cout << "Printing all the edges: ";

    for(int i=0;i<fc;i++){
        for(int j=0;j<fc;j++){
            fin >> u >> v;
            g.addEdge(u,v,1);
            cout << u+1 << " " << v+1 <<endl;
        }
    }
   
    for(int i=1;i<=fc;i++){
        g.addEdge(pseudo_source, i, 1);   
        cout << pseudo_source << " " << i <<endl;             
    }

    for(int i=fc+1;i<=2*fc;i++){
     cout << i << " " << pseudo_sink <<endl;
        g.addEdge(i, pseudo_sink, 1);                
    }


    if(!g.check_bipartite(1,fc,pseudo_source, pseudo_sink) || !g.check_bipartite(fc+1,2*fc,pseudo_source,pseudo_sink)){
        cout << "Error! Input graph is not bipartite. Aborting!";
        return(0);
    }

    cout<<"Finding maximum flow from " << pseudo_source << " to " << pseudo_sink <<endl;
    cout << "Maximum flow is " << g.getMaxFlow(pseudo_source, pseudo_sink); 
    // cout << "Max flow: " << g.getMaxFlow(0,5);
    return 0; 
} 