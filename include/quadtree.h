
#pragma once
#include <iostream>
#include <Eigen/Dense>
#include "IntervalVec.hpp"
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include "Node.hpp"
#include <queue>
#include <vector>
#include <cmath>
typedef kv::interval<double> itv;
class quadtree
{
public:
    bool status;
private:
    Eigen::MatrixXd m_mat_f1;
    Eigen::MatrixXd m_mat_f2;
    Eigen::MatrixXd m_mat_f1x;
    Eigen::MatrixXd m_mat_f1y;
    Eigen::MatrixXd m_mat_f2x;
    Eigen::MatrixXd m_mat_f2y;

    std::vector<IntervalVec<itv>> Plist;
    std::queue<IntervalVec<itv> *> Tlist;
    double m_min_size;
private:
    /* data */

    


    void Bisection(IntervalVec<itv> *p, const IntervalMat<itv> & Fprime);
    itv estimate_range(const Eigen::MatrixXd &src, const itv &xrange, const itv &yrange);
    IntervalVec<itv> F_value(const itv &xrange, const itv &yrange);
    IntervalMat<itv> F_prime_value(const itv &xrange, const itv &yrange);
    Eigen::MatrixXd diff_about_x(const Eigen::MatrixXd &src);
    Eigen::MatrixXd diff_about_y(const Eigen::MatrixXd &src);
    int count;
public:
    quadtree(const Eigen::MatrixXd &f1, const Eigen::MatrixXd &f2, const IntervalVec<itv> &range);
    ~quadtree();
    void search_procedure(IntervalVec<itv> range);
};
