#ifndef POINT_HPP_
#define POINT_HPP_

#include "node.hpp"

using namespace std;

class point {
    public:
        point() {
            nod = node();
            weight = 1;
            visited = false;
        }
        point(int num, int w) {
            nod = node(num);
            weight = w;
            visited = false;
        }
        point(int num) {
            nod = node(num);
            weight = 1;
            visited = false;
        }
        void setNode(node n) {
            nod = n;
        }
        void setWeight(int w) {
            weight = w;
        }
        node getNode() {
            return nod;
        }
        int getWeight() {
            return weight;
        }
        bool getVisited() {
            return visited;
        }
        void makeVisited() {
            visited = true;
        }
        int getClusterNum() {
            return cluster_number;
        }
        void setClusterNum(int groupnum) {
            cluster_number = groupnum;
        }

    private:
        node nod;
        int weight;
        bool visited;
        int cluster_number;
};

#endif
