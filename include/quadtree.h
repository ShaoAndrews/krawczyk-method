
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
#include "Global.h"
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

    


    bool Bisection(IntervalVec<itv> *&p, const IntervalMat<itv> & Fprime);//如果是真，说明还需要细分；如果是假，说明已经达到了指定的精度。
    itv estimate_range(const Eigen::MatrixXd &src, const itv &xrange, const itv &yrange);//给定系数矩阵，估计范围
    IntervalVec<itv> F_value(const itv &xrange, const itv &yrange);//估计F(X)的范围
    IntervalMat<itv> F_prime_value(const itv &xrange, const itv &yrange);//估计F'(X)的范围
    Eigen::MatrixXd diff_about_x(const Eigen::MatrixXd &src);//根据系数矩阵对x进行求导
    Eigen::MatrixXd diff_about_y(const Eigen::MatrixXd &src);//根据系数矩阵对y进行求导
    int count;
    void clear_Tlist();
public:
    quadtree(const Eigen::MatrixXd &f1, const Eigen::MatrixXd &f2, const IntervalVec<itv> &range);
    ~quadtree();
    void search_procedure(IntervalVec<itv> range);
};
