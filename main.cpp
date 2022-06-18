#include "analyser.hpp"
#include "mcl.hpp"
#include <chrono>
#include <eigen3/Eigen/Core>

using namespace std;

graph graph_ud_uw(false,false);
graph graph_ud_w(false,true);
graph graph_d_uw(true,false);
graph graph_d_w(false,true);



void read_ud_uw (graph g){
    ifstream infile("./data/facebook_ud_uw.txt");
    int a, b;
    while (infile >> a >> b) {
        // cout << a << b << endl;
        g.addEdge(a,b);
    }
    
}


void read_d_uw(graph g){
    ifstream infile(""); //change file path
    int a, b;
    while (infile >> a >> b) {
        g.addEdge(a,b);
    }
    
}

void read_ud_w(graph g){
    ifstream infile(""); //change file path
    int a, b, c;
    while (infile >> a >> b >> c) {
        g.addEdge(a, b, c);
    }
    
}


void read_d_w(graph g){
    ifstream infile(""); //change file path
    int a, b, c;
    while (infile >> a >> b >> c) {
        g.addEdge(a, b, c);
    }
    
}

void testing() {
    graph g = graph(false, false, 3);
    g.addEdge(0, 1);
    g.addEdge(0, 2);
    analyser a = analyser();
    cout << a.findDegree(g, 0) << endl;
    cout << a.findDegree(g, 1) << endl;

    graph f = graph(true, false, 4);
    cerr << "s1" << endl;
    f.addEdge(0, 3);
    f.addEdge(0, 1);
    f.addEdge(0, 2);
    f.addEdge(1, 3);
    f.addEdge(1, 2);
    f.addEdge(3, 0);
    f.addEdge(2, 3);
    cerr << "s2" << endl;
    cout << a.findInDegree(f, 0) << a.findOutDegree(f, 3) << a.findOutDegree(f, 1) << a.findOutDegree(f, 2) << endl;
}

void test_clustering() {
    graph g = graph(false, false, 12);
    g.addEdge(0,5);
    g.addEdge(0,6);
    g.addEdge(0,1);
    g.addEdge(0,9);
    g.addEdge(1,4);
    g.addEdge(1,2);
    g.addEdge(2,4);
    g.addEdge(2,3);
    g.addEdge(3,7);
    g.addEdge(3,8);
    g.addEdge(3,10);
    g.addEdge(4,6);
    g.addEdge(4,7);
    g.addEdge(5,9);
    g.addEdge(6,9);
    g.addEdge(7,8);
    g.addEdge(7,10);
    g.addEdge(8,10);
    g.addEdge(8,11);
    g.addEdge(10,11);

    // g.print_graph();
    mcl a = mcl();
    Eigen::SparseMatrix<double> s = a.mcl_unweighted(g, 2, 2, 15, true);
    cout << Eigen::MatrixXd(s) << endl;
}

void test_bfs() {
    graph g = graph(false, false, 4);
    g.addEdge(0, 1); 
    g.addEdge(0, 2); 
    g.addEdge(1, 2); 
    g.addEdge(2, 0); 
    g.addEdge(2, 3); 
    g.addEdge(3, 3); 
    analyser ana = analyser();
    ana.bfs(g, 2);
}

void test_cyclic() {
    // graph g = graph(true, false, 4);
    // g.addEdge(0, 1); 
    // g.addEdge(0, 2); 
    // g.addEdge(1, 2); 
    // g.addEdge(2, 0); 
    // g.addEdge(2, 3); 
    // g.addEdge(3, 3); 
    // analyser ana = analyser();
    // if(ana.isCyclic(g)) 
    //     cout << "Graph contains cycle"; 
    // else
    //     cout << "Graph doesn't contain cycle";
    graph g1 = graph(false, false, 5);
    g1.addEdge(1, 0); 
    g1.addEdge(0, 2); 
    g1.addEdge(2, 1); 
    g1.addEdge(0, 3); 
    g1.addEdge(3, 4); 
    analyser ana = analyser();
    if(ana.isCyclic(g1)) 
        cout << "Graph contains cycle" << endl; 
    else
        cout << "Graph doesn't contain cycle" << endl;
  
    graph g2 = graph(false, false, 3);
    g2.addEdge(0, 1); 
    g2.addEdge(1, 2); 
    if(ana.isCyclic(g2)) 
        cout << "Graph contains cycle" << endl; 
    else
        cout << "Graph doesn't contain cycle" << endl;
}

void test_prims() {
    graph g = graph(true, true, 5);
    g.addEdge(0, 1, 2);
    g.addEdge(0, 3, 6);
    g.addEdge(1, 0, 2);
    g.addEdge(1, 2, 3);
    g.addEdge(1, 3, 8);
    g.addEdge(1, 4, 5);
    g.addEdge(2, 1, 3);
    g.addEdge(2, 4, 7);
    g.addEdge(3, 0, 6);
    g.addEdge(3, 1, 8);
    g.addEdge(3, 4, 9);
    g.addEdge(4, 1, 5);
    g.addEdge(4, 2, 7);
    g.addEdge(4, 3, 9);

    analyser ana = analyser();
    ana.primMST(g);
}

void test_shortestPathDAG() {
    graph g = graph(true, true, 6);
    g.addEdge(0, 1, 5); 
    g.addEdge(0, 2, 3); 
    g.addEdge(1, 3, 6); 
    g.addEdge(1, 2, 2); 
    g.addEdge(2, 4, 4); 
    g.addEdge(2, 5, 2); 
    g.addEdge(2, 3, 7); 
    g.addEdge(3, 4, -1); 
    g.addEdge(4, 5, -2);
    analyser ana = analyser(g);
    int i = 0;
    ana.shortestPath(i);
}

void test_topologicalSort() {
    graph g = graph(true, false, 6);
    g.addEdge(5, 2);
    g.addEdge(5, 0); 
    g.addEdge(4, 0); 
    g.addEdge(4, 1); 
    g.addEdge(2, 3); 
    g.addEdge(3, 1);
    analyser ana = analyser();
    ana.topologicalSort(g);
}

void test_allTopologicalSort() {
    graph g = graph(true, false, 6);
    g.addEdge(5, 2);
    g.addEdge(5, 0); 
    g.addEdge(4, 0); 
    g.addEdge(4, 1); 
    g.addEdge(2, 3); 
    g.addEdge(3, 1);
    analyser ana = analyser();
    ana.allTopologicalSort(g);
}

void test_fordFulkerson() {
    graph g = graph(true, true, 6);
    g.addEdge(0, 1, 16); 
    g.addEdge(0, 2, 13); 
    g.addEdge(1, 2, 10); 
    g.addEdge(1, 3, 12); 
    g.addEdge(2, 1, 4); 
    g.addEdge(2, 0, 0);
    g.addEdge(2, 4, 14); 
    g.addEdge(3, 2, 9); 
    g.addEdge(3, 5, 20); 
    g.addEdge(4, 3, 7);
    g.addEdge(4, 5, 4);
    g.addEdge(5, 5, 0);
    g.addEdge(5, 4, 0);
    analyser ana = analyser();
    int source = 0;
    int sink = 4;

    cout << "The max flow from " << source << " to " << sink << " is: " << endl;
    cout << ana.FordFulkerson(g, source, sink) << endl;
}

void test_kruskal() {
    graph g = graph (false, true, 9);
    g.addEdge(0, 1, 4); 
    g.addEdge(0, 7, 8); 
    g.addEdge(1, 2, 8); 
    g.addEdge(1, 7, 11); 
    g.addEdge(2, 3, 7); 
    g.addEdge(2, 8, 2); 
    g.addEdge(2, 5, 4); 
    g.addEdge(3, 4, 9); 
    g.addEdge(3, 5, 14); 
    g.addEdge(4, 5, 10); 
    g.addEdge(5, 6, 2); 
    g.addEdge(6, 7, 1); 
    g.addEdge(6, 8, 6); 
    g.addEdge(7, 8, 7); 
    analyser ana = analyser();

    ana.krispal(g, 2);
}

void runtime_testing() {
    // GRAPH LOADING
    // change as needed
    cout << "Loading" << endl;
    graph g(false, false, 4039);
    read_ud_uw(g);
    cout << "Finished Loading" << endl;
    // g.print_graph();

    // CHOOSE WHICH TESTS TO RUN
    bool degree = false;
    bool dijkstra = false;
    bool bfs = false; // bfs is broken
    bool cyclic = false; // is this broken
    bool prims = false; // change this once implemented
    bool mclb = false; //mcl
    bool ford_fulkerson = false;
    bool topological_sort = false;
    bool kruskal = true;

    analyser ana = analyser();

    if(degree) {
        // Parameters for find degree
        int node = 4;

        if(g.get_directed()) {
            auto a_StartTime = std::chrono::system_clock::now();
            int inDegree = ana.findInDegree(g, node);
            auto a_EndTime = std::chrono::system_clock::now();
            auto a_dur = std::chrono::duration_cast<std::chrono::milliseconds>(a_EndTime - a_StartTime).count();
            cout << "Find in degree took " << a_dur << " ms" << endl;

            auto b_StartTime = std::chrono::system_clock::now();
            int outDegree = ana.findOutDegree(g, node);
            auto b_EndTime = std::chrono::system_clock::now();
            auto b_dur = std::chrono::duration_cast<std::chrono::milliseconds>(b_EndTime - b_StartTime).count();
            cout << "Find out degree took " << a_dur << " ms" << endl;
        } else {
            auto c_StartTime = std::chrono::system_clock::now();
            int degree = ana.findDegree(g, node);
            auto c_EndTime = std::chrono::system_clock::now();
            auto c_dur = std::chrono::duration_cast<std::chrono::milliseconds>(c_EndTime - c_StartTime).count();
            cout << "Find degree took " << c_dur << " ms" << endl;
        }
    }

    if(dijkstra) {
        // Parameters for Dijkstra
        int node = 8;

        auto d_StartTime = std::chrono::system_clock::now();
        ana.dijkstra(g, node);
        auto d_EndTime = std::chrono::system_clock::now();
        auto d_dur = std::chrono::duration_cast<std::chrono::milliseconds>(d_EndTime - d_StartTime).count();
        cout << "Dijkstra (with file saving) took " << d_dur << " ms" << endl;
    }

    if(bfs) {
        // Parameters for bfs
        int node = 3980;

        auto e_StartTime = std::chrono::system_clock::now();
        ana.bfs(g, node);
        auto e_EndTime = std::chrono::system_clock::now();
        auto e_dur = std::chrono::duration_cast<std::chrono::milliseconds>(e_EndTime - e_StartTime).count();
        cout << "BFS took " << e_dur << " ms" << endl;
    }

    if(cyclic) {
        auto f_StartTime = std::chrono::system_clock::now();
        bool isCyclic = ana.isCyclic(g);
        auto f_EndTime = std::chrono::system_clock::now();
        auto f_dur = std::chrono::duration_cast<std::chrono::milliseconds>(f_EndTime - f_StartTime).count();
        cout << "Find cyclic took " << f_dur << " ms" << endl;
        if(isCyclic) {
            cout << "------ Graph is cyclic" << endl;
        } else {
            cout << "------ Graph is not cyclic" << endl;
        }
    }

    if(prims) {
        auto g_StartTime = std::chrono::system_clock::now();
        ana.primMST(g);
        auto g_EndTime = std::chrono::system_clock::now();
        auto g_dur = std::chrono::duration_cast<std::chrono::milliseconds>(g_EndTime - g_StartTime).count();
        cout << "Prims (with file saving) took " << g_dur << " ms" << endl;
    }

    if(mclb) {
        mcl a = mcl();

        // Parameters for MCL
        int trial_to_convergence = 1;
        int expand_power = 2;
        int inflation_power = 2;
        bool add_self_loops = true;


        auto h_StartTime = std::chrono::system_clock::now();
        Eigen::SparseMatrix<double> s = a.mcl_unweighted(g, expand_power, inflation_power, trial_to_convergence, add_self_loops);
        auto h_EndTime = std::chrono::system_clock::now();
        auto h_dur = std::chrono::duration_cast<std::chrono::milliseconds>(h_EndTime - h_StartTime).count();
        cout << "MCL took " << h_dur << " ms" << endl;

        // still need to work on interpreting final matrix and getting clusters
        // only works for small matrices, was going to do spectral clustering but that would take even more time

        // cout << Eigen::MatrixXd(s) << endl;
    }

    if(ford_fulkerson) {
        // Parameters for Ford Fulkerson
        int source = 0;
        int sink = 10;

        auto i_StartTime = std::chrono::system_clock::now();
        int max_flow = ana.FordFulkerson(g, source, sink);
        auto i_EndTime = std::chrono::system_clock::now();
        auto i_dur = std::chrono::duration_cast<std::chrono::milliseconds>(i_EndTime - i_StartTime).count();
        cout << "Ford Fulkerson took " << i_dur << " ms" << endl;
        cout << "-----Max Flow was " << max_flow << endl;
    }

    if(topological_sort) {
        auto j_StartTime = std::chrono::system_clock::now();
        ana.topologicalSort(g);
        auto j_EndTime = std::chrono::system_clock::now();
        auto j_dur = std::chrono::duration_cast<std::chrono::milliseconds>(j_EndTime - j_StartTime).count();
        cout << "Topological Sorting took " << j_dur << " ms" << endl;
    }

    if(kruskal) {
        // Parameters for Kruskal
        int clusters = 3;

        auto k_StartTime = std::chrono::system_clock::now();
        ana.krispal(g, clusters);
        auto k_EndTime = std::chrono::system_clock::now();
        auto k_dur = std::chrono::duration_cast<std::chrono::milliseconds>(k_EndTime - k_StartTime).count();
        cout << "Kruskal took " << k_dur << " ms" << endl;
    }
}

int main() {
    runtime_testing();
    // test_cyclic();
    // test_bfs();
    // test_topologicalSort();
    // test_allTopologicalSort();
    // test_shortestPathDAG();
    // test_cyclic();
    // test_fordFulkerson();
    // test_prims();
    // test_kruskal();
}
