#ifndef ANALYSER_HPP_
#define ANALYSER_HPP_

#include "graph.hpp"
#include <stack>
#include <list>
#include <limits.h>
#include <fstream>
#include <algorithm>
#include "edge.hpp"
#define INF INT_MAX
#define INF2 0x3f3f3f3f

using namespace std;

typedef pair<int, int> pi;

class analyser {
    public:
        graph a;
        int num_nodes;
        analyser() {

        }

        analyser(graph A) {
            a = A;
            num_nodes = a.get_number_nodes();
        }

        int countEdges() {
            return a.get_edges();
        }

        int findInDegree(int x) { // x = node number
            int degree = 0;          
            vector<vector<point>>* lis = a.get_lis2();
            for (auto i = (lis -> at(x)).begin(); i != (lis -> at(x)).end(); i++) {
                ++degree;              
            }
            return degree;  
        }

        int findOutDegree(int x) {
            int degree = 0;          
            vector<vector<point>>* lis = a.get_lis();
            for (auto i = (lis -> at(x)).begin(); i != (lis -> at(x)).end(); i++) {
                ++degree;              
            }
            return degree;  
        }

        int findDegree(int x) { // undirected graph
            return findOutDegree(x);
        }

        void bfs(int s) {
            // helper array that keeps tracks of which nodes have been visited
            bool *visited = new bool[num_nodes];
            for(int i = 0; i < num_nodes; i++) {
                visited[i] = false;
            }
            
            list<int> queue;

            //add node that was passed to list of visited and then queue it up
            visited[s] = true;
            queue.push_back(s);
            vector<vector<point>>* lis = a.get_lis();
            ofstream ofs;
            ofs.open("./results/bfs.txt", ofstream::out | ofstream::trunc); // clears text file
            ofs.close();
            while(!queue.empty()) {

                //pop vertex from the front and print it
                int ss = queue.front();
                // cout << s.getNum() << ", ";
                ofstream myfile;
                myfile.open ("./results/bfs.txt", std::ios_base::app);
                myfile << ss << endl;
                myfile.close();
                queue.pop_front();

                // Get all adjacent vertices of the dequeud vertex s. If a adjacent has not been visited, then mark it visited and enqueue it
                for (auto i = lis -> at(ss).begin(); i != lis -> at(ss).end(); i++) {
                    if(!visited[i -> getNode().getNum()]) {
                        visited[i -> getNode().getNum()] = true;
                        queue.push_back(i -> getNode().getNum());
                    }
                }
            }
        }

        int minDistanceIndex(int dist[], bool sptSet[]) { 
            int min = INF, min_index; 
        
            for (int v = 0; v < a.get_number_nodes(); v++) 
                if (sptSet[v] == false && dist[v] <= min) 
                    min = dist[v], min_index = v; 
        
            return min_index; 
        }

        void dijkstra(int src) { 
            priority_queue<pi, vector <pi>, greater<pi>> minheap;

            vector<int> dist(num_nodes, INF2); 

            minheap.push(make_pair(0, src)); 
            dist[src] = 0; 

            while (!minheap.empty()) { 
                int u = minheap.top().second; 
                minheap.pop(); 

                for (auto i = a.get_lis() -> at(u).begin(); i != a.get_lis() -> at(u).end(); ++i) { 
                    int v = i -> getNode().getNum(); 
                    int weight = i -> getWeight(); 
        
                    if (dist[v] > dist[u] + weight) { 
                        dist[v] = dist[u] + weight; 
                        minheap.push(make_pair(dist[v], v)); 
                    } 
                } 
            }

            ofstream ofs;
            ofs.open("./results/dijkstra.txt", ofstream::out | ofstream::trunc); // clears text file
            ofs.close();
            ofstream myfile;
            myfile.open ("./results/dijkstra.txt");
            for (int i = 0; i < num_nodes; i++) {
                myfile << i << " " << dist[i] << endl;
            }
            myfile.close();
        }

    // STILL NEEDS TESTING AND FORMATTING AND JUST ABOUT EVERYTHING ELSE
        void primMST(graph g) {
            //okay so I don't know how to make this without pairs because i have no clue how PQs in the
            //STL treats putting objects in heaps, with pairs you can order by weight (cost) and heap
            //recognizes the first element of the pair as the thing to prioritize
            priority_queue<pi, vector <pi>, greater<pi>> minheap; 

            vector<vector<point>>* adjlis = g.get_lis();

            //this is just the starting point, pretty sure it gets more stupid if it's negative
            int src = 0; 
        
            vector<int> cost(g.get_number_nodes(), INT_MAX); 
        
            vector<int> parent(g.get_number_nodes(), -1); 
        
            vector<bool> inMST(g.get_number_nodes(), false); 

            //this is pushing the weight-node pairs, i don't want to figure out doing it with point objects because this is just easier (i think)
            minheap.push(make_pair(0, src)); 
            cost[src] = 0; 
        
            while (!minheap.empty()) 
            { 
                int u = minheap.top().second; 
                minheap.pop(); 
        
                inMST[u] = true; 
        
                for (auto i = adjlis -> at(u).begin(); i != adjlis -> at(u).end(); i++) 
                { 
                    int v = i -> getNode().getNum(); 
                    int weight = i -> getWeight(); 
        
                    if (inMST[v] == false && cost[v] > weight) 
                    { 
                        cost[v] = weight; 
                        minheap.push(make_pair(cost[v], v)); 
                        parent[v] = u; 
                    } 
                } 
            }

            // for (int i = 1; i < g.get_number_nodes(); ++i) // I think num_nodes = # vertices = # nodes, double check though
            //     printf("%d - %d\n", parent[i], i);

            // check if I did this right
            ofstream ofs;
            ofs.open("./results/prim.txt", ofstream::out | ofstream::trunc); // clears text file
            ofs.close();
            ofstream myfile;
            myfile.open ("./results/prim.txt");
            for (int i = 0; i < g.get_number_nodes(); i++) {
                myfile << parent[i] << "-" << i << "  " << cost[i] << endl;
            }
            myfile.close();
        } 

        vector<edge> kruskal(graph g, int k) {
            vector<edge> edges = g.getEdges();
            //thank you lambda functions i love you <3
            auto compare = [](edge left, edge right) { 
                return (left.getWeight()) > (right.getWeight()); 
            };
            priority_queue<edge, vector<edge>, decltype(compare)> pq(compare);

            for(auto i = edges.begin(); i != edges.end(); i++)
                pq.push(*i);

            int clustum = 1;
            vector<edge> mst;

            while(pq.size() > 0) {
                edge current = pq.top();
                pq.pop();
                if (!g.get_list() -> at(current.getSource()).at(current.getDest()).getVisited()) {
                    mst.push_back(edge(current.getSource(), current.getDest(), current.getWeight()));
                    g.get_list() -> at(current.getSource()).at(current.getDest()).makeVisited();
                    g.get_list() -> at(current.getSource()).at(current.getDest()).setClusterNum(clustum);
                    }
                else {
                clustum++;
                if(clustum > k)
                    clustum = 1;
                }
            }
            return mst;
        }

        void krispal(graph g, int k) { // kruskal but things are in text file
            vector<edge> edges = g.getEdges();
            //thank you lambda functions i love you <3
            auto compare = [](edge left, edge right) { 
                return (left.getWeight()) > (right.getWeight()); 
            };
            priority_queue<edge, vector<edge>, decltype(compare)> pq(compare);

            for(auto i = edges.begin(); i != edges.end(); i++)
                pq.push(*i);

            int clustum = 1;
            vector<edge> mst;
            // something = g.get_ls
            // vector<vector<point>>* lis = g.get_lis();
            
            while(pq.size() > 0) {
                edge current = pq.top();
                // cout << current.getSource() << "-" << current.getDest() << "\t" << current.getWeight() << endl;
                pq.pop();
                // cerr << g.get_list() -> at(0).size() << endl;
                if (!g.get_list() -> at(current.getSource()).at(current.getDest()).getVisited()) {
                    // cerr << "kmao" << endl;
                    mst.push_back(edge(current.getSource(), current.getDest(), current.getWeight()));
                    g.get_list() -> at(current.getSource()).at(current.getDest()).makeVisited();
                    g.get_list() -> at(current.getSource()).at(current.getDest()).setClusterNum(clustum);
                }
                else {
                    // cerr << "adadf" << endl;
                    clustum++;
                    if(clustum > k)
                        clustum = 1;
                }
            }
            // cerr << "lmao3" << endl;
            ofstream myfile;
            myfile.open ("./results/kruskal.txt");
            for (auto matic = mst.begin(); matic != mst.end(); matic++) {
                // cout << matic -> getSource();
                myfile << matic -> getSource() << "-" << matic -> getDest() << "\t" << matic -> getWeight() << endl;
            }
            myfile.close();
        }

        vector<edge> kruskal(int k) {
            vector<edge> edges = a.getEdges();
            //thank you lambda functions i love you <3
            auto compare = [](edge left, edge right) { 
                return (left.getWeight()) > (right.getWeight()); 
            };
            priority_queue<edge, vector<edge>, decltype(compare)> pq(compare);

            for(auto i = edges.begin(); i != edges.end(); i++)
                pq.push(*i);

            int clustum = 1;
            vector<edge> mst;

            while(pq.size() > 0) {
                edge current = pq.top();
                pq.pop();
                if (!a.get_list() -> at(current.getSource()).at(current.getDest()).getVisited()) {
                    mst.push_back(edge(current.getSource(), current.getDest(), current.getWeight()));
                    a.get_list() -> at(current.getSource()).at(current.getDest()).makeVisited();
                    a.get_list() -> at(current.getSource()).at(current.getDest()).setClusterNum(clustum);
                    }
                else {
                clustum++;
                if(clustum > k)
                    clustum = 1;
                }
            }
            return mst;
        }

        void krispal(int k) { // kruskal but things are in text file
            vector<edge> edges = a.getEdges();
            //thank you lambda functions i love you <3
            auto compare = [](edge left, edge right) { 
                return (left.getWeight()) > (right.getWeight()); 
            };
            priority_queue<edge, vector<edge>, decltype(compare)> pq(compare);

            for(auto i = edges.begin(); i != edges.end(); i++)
                pq.push(*i);

            int clustum = 1;
            vector<edge> mst;

            while(pq.size() > 0) {
                edge current = pq.top();
                pq.pop();
                if (!a.get_list() -> at(current.getSource()).at(current.getDest()).getVisited()) {
                    mst.push_back(edge(current.getSource(), current.getDest(), current.getWeight()));
                    a.get_list() -> at(current.getSource()).at(current.getDest()).makeVisited();
                    a.get_list() -> at(current.getSource()).at(current.getDest()).setClusterNum(clustum);
                    }
                else {
                clustum++;
                if(clustum > k)
                    clustum = 1;
                }
            }
            ofstream myfile;
            myfile.open ("./results/kruskal.txt");
            for (auto matic = mst.begin(); matic != mst.end(); matic++) {
                myfile << matic -> getSource() << "-" << matic -> getDest() << "\t" << matic -> getWeight() << endl;
            }
            myfile.close();
        }




        // void primMST2(graph g) {
        //     int num_nodes = g.get_number_nodes();
        //     int parent[num_nodes];  
        //     int keys[num_nodes];  
        //     bool mstSet[num_nodes];  
        
        //     // Initialize all keys as INFINITE  
        //     for (int i = 0; i < num_nodes; i++) {
        //         keys[i] = INT_MAX, mstSet[i] = false;  
        //     }
        //     // Always include first 1st vertex in MST.  
        //     // Make key 0 so that this vertex is picked as first vertex.  
        //     keys[0] = 0;  
        //     parent[0] = -1; // First node is always root of MST  
        
        //     // The MST will have V vertices  
        //     for (int count = 0; count < num_nodes - 1; count++) {  
        //         // Pick the minimum key vertex from the  
        //         // set of vertices not yet included in MST  
        //         int min = INT_MAX, min_index;  
  
        //         for (int v = 0; v < num_nodes; v++)  
        //         if (mstSet[v] == false && keys[v] < min)  
        //             min = keys[v], min_index = v;  
  
        //         int u = min_index;  
        
        //         // Add the picked vertex to the MST Set  
        //         mstSet[u] = true;  
        
        //         // Update key value and parent index of  
        //         // the adjacent vertices of the picked vertex.  
        //         // Consider only those vertices which are not  
        //         // yet included in MST  
        //         for (int eger = 0; eger < num_nodes; eger++)  
        
        //             // graph[u][v] is non zero only for adjacent vertices of m  
        //             // mstSet[v] is false for vertices not yet included in MST  
        //             // Update the key only if graph[u][v] is smaller than key[v]  
        //             if (g.get_lis() -> at(u).at(eger).getNode().getNum() != 0 && mstSet[eger] == false && g.get_lis() -> at(u).at(eger).getNode().getNum() < keys[eger])  
        //                 parent[eger] = u, keys[eger] = g.get_lis() -> at(u).at(eger).getNode().getNum();  
        //     }  
        // }

        bool isCyclicUtil(int v, bool visited[], bool *rStack) { 
            if (visited[v] == false) { 
                visited[v] = true; 
                rStack[v] = true; 
        
                for (auto i = a.get_lis() -> at(v).begin(); i != a.get_lis() -> at(v).end(); ++i) { 
                    if (!visited[i -> getNode().getNum()] && isCyclicUtil(i -> getNode().getNum(), visited, rStack)) 
                        return true; 
                    else if (rStack[i -> getNode().getNum()]) 
                        return true; 
                } 
            } 

            rStack[v] = false;  
            return false; 
        }    

        bool isCyclicUtil(int v, bool visited[], int parent) { 
            visited[v] = true; 

            for (auto i = a.get_lis() -> at(v).begin(); i != a.get_lis() -> at(v).end(); ++i) { 
                    if (!visited[i -> getNode().getNum()] && isCyclicUtil(i -> getNode().getNum(), visited, v)) 
                        return true; 
                    else if (i -> getNode().getNum() != parent) 
                        return true; 
            } 

            return false; 
        }
        
        bool isCyclic() { 
            bool *visited = new bool[num_nodes]; 
            if (a.get_directed()){
                bool *rStack = new bool[num_nodes]; 
                for (int i = 0; i < num_nodes; i++) { 
                    visited[i] = false; 
                    rStack[i] = false; 
                } 
            
                for(int i = 0; i < num_nodes; i++) {
                    if (isCyclicUtil(i, visited, rStack)) {
                        return true; 
                    }
                }
            
                return false; 
            } else {
                for (int i = 0; i < num_nodes; i++)
                    visited[i] = false;
                for(int i = 0; i < num_nodes; i++) {
                    if (!visited[i] && isCyclicUtil(i, visited, -1)) {
                        return true;
                    }
                }

                return false;
            }

        }

        void topologicalSortUtil(int v, bool visited[], stack<int> &Stack) { 
            visited[v] = true; 
        
            for (auto i = a.get_lis() -> at(v).begin(); i != a.get_lis() -> at(v).end(); ++i) 
                if (!visited[i -> getNode().getNum()]) 
                    topologicalSortUtil(i -> getNode().getNum(), visited, Stack); 

            Stack.push(v); 
        }

        void topologicalSort() {
            if (a.get_directed() && !isCyclic()) {
                stack<int> Stack; 

                bool *visited = new bool[num_nodes]; 
                for (int i = 0; i < num_nodes; i++) 
                    visited[i] = false; 
                
                ofstream ofs;
                ofs.open("./results/topological_sort.txt", ofstream::out | ofstream::trunc); // clears text file
                ofs.close();


                for (int i = 0; i < num_nodes; i++) 
                    if (!visited[i]) 
                        topologicalSortUtil(i, visited, Stack); 

                ofstream myfile;
                myfile.open ("./results/topological_sort.txt");


                while (Stack.empty() == false) { 
                    // cout << Stack.top() << " "
                    myfile << Stack.top() << " ";
                    Stack.pop(); 
                } 
                myfile.close();
            } else {
                cout << "Graph is not a directed acyclic graph" << endl;
            }
        }

        void topologicalSortUtil(vector<int>& res, bool visited[], vector<int> indegree) {
            bool flag = false; 

            for (int i = 0; i < num_nodes; i++) { 
                if (indegree[i] == 0 && !visited[i]) { 
                    for (auto j = a.get_lis() -> at(i).begin(); j != a.get_lis() -> at(i).end(); j++) 
                        indegree[j -> getNode().getNum()]--; 

                    res.push_back(i); 
                    visited[i] = true; 
                    topologicalSortUtil(res, visited, indegree); 

                    visited[i] = false; 
                    res.erase(res.end() - 1); 
                    for (auto j = a.get_lis() -> at(i).begin(); j != a.get_lis() -> at(i).end(); j++) 
                        indegree[j -> getNode().getNum()]++; 

                    flag = true; 
                } 
            } 

            if (!flag) { 
                for (int i = 0; i < res.size(); i++) 
                    cout << res[i] << " "; 
                cout << endl; 
            }
        }

        void allTopologicalSort() {
            if (a.get_directed() && !isCyclic(a)) {
                vector<int> indegree;
                for (int i = 0; i < num_nodes; i++) {
                    indegree.push_back(0);
                }
                bool *visited = new bool[num_nodes]; 
                for (int i = 0; i < num_nodes; i++) {
                    indegree[i] = findInDegree(a, i);
                }
                for (int i = 0; i < num_nodes; i++) 
                    visited[i] = false; 
            
                vector<int> res; 
                topologicalSortUtil(a, res, visited, indegree);
            } else {
                cout << "Graph is not a directed acyclic graph" << endl;
            }
        }

        void shortestPath(int s) { 
            if (a.get_directed() && !isCyclic()) {
                stack<int> Stack; 
                int dist[num_nodes]; 

                bool *visited = new bool[num_nodes]; 
                for (int i = 0; i < num_nodes; i++) 
                    visited[i] = false; 

                for (int i = 0; i < num_nodes; i++) 
                    if (visited[i] == false) 
                        topologicalSortUtil(i, visited, Stack); 

                for (int i = 0; i < num_nodes; i++) 
                    dist[i] = INF; 
                dist[s] = 0; 

                ofstream ofs;
                ofs.open("./results/shortest_path_directed_acyclic_graph.txt", ofstream::out | ofstream::trunc); // clears text file
                ofs.close();

                while (Stack.empty() == false) { 
                    int u = Stack.top(); 
                    Stack.pop(); 

                    if (dist[u] != INF) { 
                        for (auto i = a.get_lis() -> at(u).begin(); i != a.get_lis() -> at(u).end(); ++i) {
                            if (dist[i -> getNode().getNum()] > dist[u] + i -> getWeight()) 
                                dist[i -> getNode().getNum()] = dist[u] + i -> getWeight(); 
                        }
                    } 
                } 

                ofstream myfile;
                myfile.open ("./results/shortest_path_directed_acyclic_graph.txt");
            
                for (int i = 0; i < num_nodes; i++) {
                    (dist[i] == INF)? myfile << "INF ": myfile << dist[i] << " "; 
                }
                myfile.close();
            } else {
                cout << "Graph is not a directed acyclic graph" << endl;
            }
        }

        int FordFulkerson(int &source, int &sink) {
            int maxflow = 0;
            vector<vector<point>>* lis = a.get_lis();

            int res_lis[num_nodes][num_nodes];

            for(int i = 0; i < num_nodes; i++) {
                for(auto j = lis -> at(i).begin(); j != lis -> at(i).end(); j++) {
                    res_lis[i][j -> getNode().getNum()] = j -> getWeight();
                }
            } 
 
            int parent[num_nodes];

            while(FordFulkersonUtil(num_nodes, *res_lis, source, sink, parent)) {
                int pathflow = INF;

                int v = sink;
                for (v = sink; v != source; v = parent[v]) { 
                    int u = parent[v]; 
                    pathflow = min(pathflow, res_lis[u][v]); 
                } 

                for (v = sink; v != source; v = parent[v]) {
                    int u = parent[v];
                    res_lis[u][v] -= pathflow;
                    res_lis[v][u] += pathflow;
                }

                maxflow += pathflow;
            }
            
            return maxflow;
        }

        int countEdges(graph b) {
            return b.get_edges();
        }

         int findInDegree(graph b, int x) { // x = node number
            int degree = 0;          
            vector<vector<point>>* lis = b.get_lis2();
            for (auto i = (lis -> at(x)).begin(); i != (lis -> at(x)).end(); i++) {
                ++degree;              
            }
            return degree;  
        }

        int findOutDegree(graph b, int x) {
            int degree = 0;          
            vector<vector<point>>* lis = b.get_lis();
            for (auto i = (lis -> at(x)).begin(); i != (lis -> at(x)).end(); i++) {
                ++degree;              
            }
            return degree;  
        }

        int findDegree(graph b, int x) { // undirected graph
            return findOutDegree(b, x);
        }

        void bfs(graph b, int s) {
            int num_nodes = b.get_number_nodes();
            bool *visited = new bool[num_nodes];
            for(int i = 0; i < num_nodes; i++) {
                visited[i] = false;
            }
            
            list<int> queue;

            //add node that was passed to list of visited and then queue it up
            visited[s] = true;
            queue.push_back(s);
            vector<vector<point>>* lis = b.get_lis();
            ofstream ofs;
            ofs.open("./results/bfs.txt", ofstream::out | ofstream::trunc); // clears text file
            ofs.close();
            ofstream myfile;
            myfile.open ("./results/bfs.txt");

            while(!queue.empty()) {

                //pop vertex from the front and print it
                int ss = queue.front();
                // cout << ss << ", ";
                myfile << ss << endl;
                queue.pop_front();

                // Get all adjacent vertices of the dequeud vertex s. If a adjacent has not been visited, then mark it visited and enqueue it
                for (auto i = lis -> at(ss).begin(); i != lis -> at(ss).end(); i++) {
                    if(!visited[i -> getNode().getNum()]) {
                        visited[i -> getNode().getNum()] = true;
                        queue.push_back(i -> getNode().getNum());
                    }
                }
            }
            myfile.close();
        }

        int minDistanceIndex(graph b, int dist[], bool sptSet[]) { 
            int min = INT_MAX, min_index; 
        
            for (int v = 0; v < b.get_number_nodes(); v++) 
                if (sptSet[v] == false && dist[v] <= min) 
                    min = dist[v], min_index = v; 
        
            return min_index; 
        }

        void dijkstra(graph b, int src) { 
            int num_nodes = b.get_number_nodes();
            priority_queue<pi, vector <pi>, greater<pi>> minheap;

            vector<int> dist(num_nodes, INF2); 

            minheap.push(make_pair(0, src)); 
            dist[src] = 0; 

            while (!minheap.empty()) { 
                int u = minheap.top().second; 
                minheap.pop(); 

                for (auto i = b.get_lis() -> at(u).begin(); i != b.get_lis() -> at(u).end(); ++i) { 
                    int v = i -> getNode().getNum(); 
                    int weight = i -> getWeight(); 
        
                    if (dist[v] > dist[u] + weight) { 
                        dist[v] = dist[u] + weight; 
                        minheap.push(make_pair(dist[v], v)); 
                    } 
                } 
            }
            ofstream ofs;
            ofs.open("./results/dijkstra.txt", ofstream::out | ofstream::trunc); // clears text file
            ofs.close();
            ofstream myfile;
            myfile.open ("./results/dijkstra.txt");
            for (int i = 0; i < num_nodes; i++) {
                myfile << i << " " << dist[i] << endl;
            }
            myfile.close();
        }

    // STILL NEEDS TESTING AND FORMATTING AND JUST ABOUT EVERYTHING ELSE
        void primMST() {
            graph g = a;
            //okay so I don't know how to make this without pairs because i have no clue how PQs in the
            //STL treats putting objects in heaps, with pairs you can order by weight (cost) and heap
            //recognizes the first element of the pair as the thing to prioritize
            priority_queue< pi, vector <pi> , greater<pi> > minheap; 

            vector<vector<point>>* adjlis = g.get_lis();

            //this is just the starting point, pretty sure it gets more stupid if it's negative
            int src = 0; 
        
            vector<int> cost(g.get_number_nodes(), INT_MAX); 
        
            vector<int> parent(g.get_number_nodes(), -1); 
        
            vector<bool> inMST(g.get_number_nodes(), false); 

            //this is pushing the weight-node pairs, i don't want to figure out doing it with point objects because this is just easier (i think)
            minheap.push(make_pair(0, src)); 
            cost[src] = 0; 
        
            while (!minheap.empty()) 
            { 
                int u = minheap.top().second; 
                minheap.pop(); 
        
                inMST[u] = true; 
        
                for (auto i = adjlis -> at(u).begin(); i != adjlis -> at(u).end(); i++) 
                { 
                    int v = i -> getNode().getNum(); 
                    int weight = i -> getWeight(); 
        
                    if (inMST[v] == false && cost[v] > weight) 
                    { 
                        cost[v] = weight; 
                        minheap.push(make_pair(cost[v], v)); 
                        parent[v] = u; 
                    } 
                } 
            }

            // for (int i = 1; i < g.get_number_nodes(); ++i) // I think num_nodes = # vertices = # nodes, double check though
            //     printf("%d - %d\n", parent[i], i);

            // check if I did this right
            ofstream ofs;
            ofs.open("./results/prim.txt", ofstream::out | ofstream::trunc); // clears text file
            ofs.close();
            ofstream myfile;
            myfile.open ("./results/prim.txt");
            for (int i = 0; i < g.get_number_nodes(); i++) {
                myfile << parent[i] << " " << i << endl;
            }
            myfile.close();
        } 

        bool isCyclicUtil(graph b, int v, bool visited[], bool *rStack) { 
            if (visited[v] == false) { 
                visited[v] = true; 
                rStack[v] = true; 
        
                for (auto i = b.get_lis() -> at(v).begin(); i != b.get_lis() -> at(v).end(); ++i) { 
                    if (!visited[i -> getNode().getNum()] && isCyclicUtil(b, i -> getNode().getNum(), visited, rStack)) 
                        return true; 
                    else if (rStack[i -> getNode().getNum()]) 
                        return true; 
                } 
            } 

            rStack[v] = false;  
            return false; 
        }    

        bool isCyclicUtil(graph b, int v, bool visited[], int parent) { 
            visited[v] = true; 

            for (auto i = b.get_lis() -> at(v).begin(); i != b.get_lis() -> at(v).end(); ++i) { 
                    if (!visited[i -> getNode().getNum()] && isCyclicUtil(b, i -> getNode().getNum(), visited, v)) 
                        return true; 
                    else if (i -> getNode().getNum() != parent) 
                        return true; 
            } 

            return false; 
        }
        
        bool isCyclic(graph b) { 
            int num_nodes = b.get_number_nodes();
            bool *visited = new bool[num_nodes]; 
            if (b.get_directed()){
                bool *rStack = new bool[num_nodes]; 
                for (int i = 0; i < num_nodes; i++) { 
                    visited[i] = false; 
                    rStack[i] = false; 
                } 
            
                for(int i = 0; i < num_nodes; i++) {
                    if (isCyclicUtil(b, i, visited, rStack)) {
                        return true; 
                    }
                }
            
                return false; 
            } else {
                for (int i = 0; i < num_nodes; i++)
                    visited[i] = false;
                for(int i = 0; i < num_nodes; i++) {
                    if (!visited[i] && isCyclicUtil(b, i, visited, -1)) {
                            return true;
                    }
                }

                return false;
            }

        }

        void topologicalSortUtil(graph b, int v, bool visited[], stack<int> &Stack) { 
            visited[v] = true; 
        
            for (auto i = b.get_lis() -> at(v).begin(); i != b.get_lis() -> at(v).end(); ++i) 
                if (!visited[i -> getNode().getNum()]) 
                    topologicalSortUtil(b, i -> getNode().getNum(), visited, Stack); 

            Stack.push(v); 
        }

        void topologicalSort(graph b) {
            if (b.get_directed() && !isCyclic(b)) {
                int num_nodes = b.get_number_nodes();

                stack<int> Stack; 

                bool *visited = new bool[num_nodes]; 
                for (int i = 0; i < num_nodes; i++) 
                    visited[i] = false; 
                
                ofstream ofs;
                ofs.open("./results/topological_sort.txt", ofstream::out | ofstream::trunc); // clears text file
                ofs.close();

                for (int i = 0; i < num_nodes; i++) 
                    if (!visited[i]) 
                        topologicalSortUtil(b, i, visited, Stack); 

                ofstream myfile;
                myfile.open ("./results/topological_sort.txt");


                while (Stack.empty() == false) { 
                    // cout << Stack.top() << " ";
                    myfile << Stack.top() << " ";
                    Stack.pop(); 
                } 
                myfile.close();
            } else {
                cout << "Graph is not a directed acyclic graph" << endl;
            }
        }

        void topologicalSortUtil(graph b, vector<int>& res, bool visited[], vector<int> indegree) {
            int num_nodes = b.get_number_nodes();
            bool flag = false; 

            for (int i = 0; i < num_nodes; i++) { 
                if (indegree[i] == 0 && !visited[i]) { 
                    for (auto j = b.get_lis() -> at(i).begin(); j != b.get_lis() -> at(i).end(); j++) 
                        indegree[j -> getNode().getNum()]--; 

                    res.push_back(i); 
                    visited[i] = true; 
                    topologicalSortUtil(b, res, visited, indegree); 

                    visited[i] = false; 
                    res.erase(res.end() - 1); 
                    for (auto j = b.get_lis() -> at(i).begin(); j != b.get_lis() -> at(i).end(); j++) 
                        indegree[j -> getNode().getNum()]++; 

                    flag = true; 
                } 
            } 

            if (!flag) { 
                for (int i = 0; i < res.size(); i++) 
                    cout << res[i] << " "; 
                cout << endl; 
            }
        }

        void allTopologicalSort(graph b) {
            if (b.get_directed() && !isCyclic(b)) {
                int num_nodes = b.get_number_nodes();
                vector<int> indegree;
                for (int i = 0; i < num_nodes; i++) {
                    indegree.push_back(0);
                }
                bool *visited = new bool[num_nodes]; 
                for (int i = 0; i < num_nodes; i++) {
                    indegree[i] = findInDegree(b, i);
                }
                for (int i = 0; i < num_nodes; i++) 
                    visited[i] = false; 
            
                vector<int> res; 
                topologicalSortUtil(b, res, visited, indegree);
            } else {
                cout << "Graph is not a directed acyclic graph" << endl;
            }
        }

        void shortestPath(graph b, int s) {
            if (b.get_directed() && !isCyclic(b)) {
                int num_nodes = b.get_number_nodes();
                stack<int> Stack; 
                int dist[num_nodes]; 

                bool *visited = new bool[num_nodes]; 
                for (int i = 0; i < num_nodes; i++) 
                    visited[i] = false; 

                for (int i = 0; i < num_nodes; i++) 
                    if (visited[i] == false) 
                        topologicalSortUtil(b, i, visited, Stack); 

                for (int i = 0; i < num_nodes; i++) 
                    dist[i] = INF; 
                dist[s] = 0; 

                ofstream ofs;
                ofs.open("./results/shortest_path_directed_acyclic_graph.txt", ofstream::out | ofstream::trunc); // clears text file
                ofs.close();

                while (Stack.empty() == false) { 
                    int u = Stack.top(); 
                    Stack.pop(); 

                    if (dist[u] != INF) { 
                        for (auto i = b.get_lis() -> at(u).begin(); i != b.get_lis() -> at(u).end(); ++i) 
                            if (dist[i -> getNode().getNum()] > dist[u] + i -> getWeight()) 
                                dist[i -> getNode().getNum()] = dist[u] + i -> getWeight(); 
                    } 
                } 

                ofstream myfile;
                myfile.open ("./results/shortest_path_directed_acyclic_graph.txt");
            
                for (int i = 0; i < num_nodes; i++) {
                    (dist[i] == INF)? myfile << "INF ": myfile << dist[i] << " "; 
                }
                myfile.close();
            } else {
                cout << "Graph is not a directed acyclic graph" << endl;
            }
        }

        bool FordFulkersonUtil(int num_nodes, int* res_lis2, int &src, int &sink, int parent[]) {
            int (*res_lis)[num_nodes][num_nodes] = (int (*)[num_nodes][num_nodes]) res_lis2;
            bool visited[num_nodes];
            fill_n(visited, num_nodes, false);

            queue<int> q;

            q.push(src);
            visited[src] = true;
            parent[src] = -1;

            while (!q.empty()) {
                int u = q.front();
                q.pop();

                for (int v = 0; v < num_nodes; v++) {
                    int capacity = (*res_lis)[u][v];
                    if (!visited[v] && capacity > 0) {
                        q.push(v);
                        parent[v] = u;
                        visited[v] = true;
                    }
                }
            }

            return (visited[sink] == true); 
        }  

        int FordFulkerson(graph b, int source, int sink) {
            int maxflow = 0;
            vector<vector<point>>* lis = b.get_lis();

            int num_nodes = b.get_number_nodes();

            int res_lis[num_nodes][num_nodes];

            for(int i = 0; i < num_nodes; i++) {
                for(auto j = lis -> at(i).begin(); j != lis -> at(i).end(); j++) {
                    res_lis[i][j -> getNode().getNum()] = j -> getWeight();
                }
            } 
 
            int parent[num_nodes];

            while(FordFulkersonUtil(num_nodes, *res_lis, source, sink, parent)) {
                int pathflow = INF;

                int v = sink;
                for (v = sink; v != source; v = parent[v]) { 
                    int u = parent[v]; 
                    pathflow = min(pathflow, res_lis[u][v]); 
                } 

                for (v = sink; v != source; v = parent[v]) {
                    int u = parent[v];
                    res_lis[u][v] -= pathflow;
                    res_lis[v][u] += pathflow;
                }

                maxflow += pathflow;
            }
            
            return maxflow;
        } 
};

#endif
