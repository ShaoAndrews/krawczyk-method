//存储多项式的系数矩阵

#pragma once
#include<Eigen/Dense>
#include<kv/interval.hpp>
#include<kv/rdouble.hpp>
#include"IntervalVec.hpp"
using namespace std;
using namespace Eigen;
typedef kv::interval<double> itv;
class ExCoefficients
{
private:
    MatrixXd m_mat_f1;
    MatrixXd m_mat_f2;
    MatrixXd m_mat_f1x;
    MatrixXd m_mat_f1y;
    MatrixXd m_mat_f2x;
    MatrixXd m_mat_f2y;
private:
    MatrixXd diff_about_x(const MatrixXd & src);
    MatrixXd diff_about_y(const MatrixXd & src);
public:
    ExCoefficients(const Eigen::MatrixXd & f1_mat, const Eigen::MatrixXd & f2_mat);
    ~ExCoefficients();
    itv estimate_range(const Eigen::MatrixXd &src, const itv &xrange, const itv &yrange);
    IntervalVec<itv> F_value(const itv &xrange, const itv &yrange);
    IntervalMat<itv> F_prime_value(const itv &xrange, const itv &yrange);
    IntervalVec<itv> calculate_K(const IntervalVec<itv> &initial);
};

