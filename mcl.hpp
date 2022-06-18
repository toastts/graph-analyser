#ifndef MCL_HPP_
#define MCL_HPP_

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Core>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <iostream>
#include "graph.hpp"

class mcl {
    public:
        mcl() {

        }
        
        void normize(Eigen::SparseMatrix<double>& mat) {
            for(int i = 0; i < mat.cols(); i++) {
                double sum = 0;
                for(int j = 0; j < mat.rows(); j++) {
                    sum += mat.coeff(j,i); 
                }
                for(int j = 0; j < mat.rows(); j++) {
                    // std::cout << mat.coeff(j,i) << std::endl;
                    mat.coeffRef(j,i) = mat.coeff(j,i)/sum;
                    // std::cout << mat.coeff(j,i) << std::endl;
                }
            }
        }

        void expand(Eigen::SparseMatrix<double>& mat, double expand_factor = 2) {
            Eigen::SparseMatrix<double> a(mat);
            for(int i = 1; i < expand_factor; i++) {
                a = a * mat;
            }
            mat = a;
        }

        void inflate(Eigen::SparseMatrix<double>& mat, double inflate_factor = 2) {
            Eigen::SparseMatrix<double> a(mat);
            // std::cout << "INFLATING START" << std::endl;
            // std::cout << a << std::endl;
            for(int i = 1; i < inflate_factor; i++) {
                a = a.cwiseProduct(mat);
            }
            // std::cout << a << std::endl;
            // std::cout << "INFLATING END" << std::endl;
            mat = a;
        }

        Eigen::SparseMatrix<double> mcl_unweighted(graph g, double expand_factor = 2, double inflate_factor = 2, double max_loop = 10, bool self_loops = false) {
            if(g.get_weighted()) {
                std::cout << "Given a weighted graph" << std::endl;
            } else {
                std::vector<std::vector<point>>* lis;
                if(g.get_directed()) {
                    lis = g.get_lis2();
                } else {
                    lis = g.get_lis();
                }
                Eigen::SparseMatrix<double> trans_mat;
                trans_mat.resize(g.get_number_nodes(), g.get_number_nodes());
                trans_mat.reserve(g.get_edges());
                int row = 0;
                for(auto i = lis -> begin(); i != lis -> end(); i++) {
                    // std::cout << row << std::endl;
                    double sum;
                    if(self_loops) {
                        sum = 1;
                    } else {
                        sum = 0;
                    }
                    for(auto j = i -> begin(); j != i -> end(); j++) {
                        sum += j -> getWeight();
                    }
                    // std::cout << sum << std::endl;
                    for(auto j = i -> begin(); j != i -> end(); j++) {
                        trans_mat.insert((j -> getNode()).getNum(), row) = (j -> getWeight())/sum;
                        // std::cout << (j -> getWeight())/sum << std::endl;
                    }
                    if(self_loops) {
                        trans_mat.insert(row, row) = 1.0/sum;
                    }
                    ++row;
                }
                std::cout << Eigen::MatrixXd(trans_mat) << std::endl;
                for(int t = 0; t < max_loop; t++) { // run til convergence, dont know how to yet
                    // std::cout << trans_mat << std::endl;
                    std::cout << t + 1 << " iteration" << std::endl;
                    expand(trans_mat);
                    inflate(trans_mat);
                    normize(trans_mat);
                    // std::cout << trans_mat << std::endl;
                    
                    // std::cout << trans_mat << std::endl;
                }

                // expand(trans_mat);
                // inflate(trans_mat);
                // std::cout << "Inflating" << std::endl;
                // std::cout << trans_mat << std::endl;
                // normize(trans_mat);
                // std::cout << "Normizing" << std::endl;
                // std::cout << trans_mat << std::endl;
                // expand(trans_mat);
                // std::cout << "Expanding 2" << std::endl;
                // std::cout << trans_mat << std::endl;
                // inflate(trans_mat);
                // std::cout << "Inflating 2" << std::endl;
                // std::cout << trans_mat << std::endl;
                // normize(trans_mat);
                // std::cout << "Normizing 2" << std::endl;
                // std::cout << trans_mat << std::endl;
                // for(int i = 0; i < g.get_number_nodes(); i++) {
                //     for(int j = 0; j < g.get_number_nodes(); j++) {
                //         std::cout << "ajdflasiodfja" << std::endl;
                //         std::cout << trans_mat.coeff(i,j);
                //         if(trans_mat.coeff(i,j) < 0.001) {
                //             trans_mat.prune(i,j);
                //         }
                //     }
                //     std::cout << std::endl;
                // }
                return trans_mat; // interpret this to find clusters
            }
            return Eigen::SparseMatrix<double>();
        }
};

#endif
