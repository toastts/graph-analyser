#ifndef EDGE_HPP_
#define EDGE_HPP_

using namespace std;

class edge {
    public:
        edge(){
            v1=0;
            v2=0;
            weight = 0;
        }

        edge(int source, int dest, int weight) {
            v1 = source;
            v2 = dest;
            this -> weight = weight;
        }

        int getWeight() {
            return this -> weight;
        }

        int getSource() {
            return this -> v1;
        }

        int getDest(){
            return this -> v2;
        }

        bool operator()(edge e1, edge e2) {
            return (e1.getWeight() > e2.getWeight());
        }


    private:
        int v1;
        int v2;
        int weight;

};




#endif
