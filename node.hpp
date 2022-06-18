#ifndef NODE_HPP_
#define NODE_HPP_

#include <cstdlib>
#include <iostream>
#include <time.h>
#include <random>

using namespace std;

class node {
    public:
        node() {
            num = 0;
        }
        node(int x) {
            num = x;
        }
        void setNum(int x) {
            num = x;
        }
        int getNum() {
            return num;
        }
    private:
        int num;

};

#endif