#include "ExCoefficients.h"

ExCoefficients::ExCoefficients(const Eigen::MatrixXd &f1_mat, const Eigen::MatrixXd &f2_mat) : m_mat_f1(f1_mat), m_mat_f2(f2_mat)
{
    m_mat_f1x = diff_about_x(m_mat_f1);
    m_mat_f1y = diff_about_y(m_mat_f1);
    m_mat_f2x = diff_about_x(m_mat_f2);
    m_mat_f2y = diff_about_y(m_mat_f2);
}

ExCoefficients::~ExCoefficients()
{
}
itv ExCoefficients::estimate_range(const Eigen::MatrixXd &src, const itv &xrange, const itv &yrange)
{
    itv result(0, 0);
    for (int i = 0; i < src.rows(); ++i)
    {
        for (int j = 0; j < src.cols(); ++j)
        {
            result += pow(xrange, i) * pow(yrange, j) * src(i, j);
        }
    }
    return itv(result.lower(), result.upper());
}
IntervalVec<itv> ExCoefficients::F_value(const itv &xrange, const itv &yrange)
{
    return IntervalVec<itv>(estimate_range(m_mat_f1, xrange, yrange), estimate_range(m_mat_f2, xrange, yrange));
}

MatrixXd ExCoefficients::diff_about_x(const MatrixXd &src)
{
    if (src.rows() == 1)
    {
        MatrixXd result(1, 1);
        result << 0;
        return result;
    }
    MatrixXd result(src.rows() - 1, src.cols());
    for (int i = 1; i < src.rows(); i++)
    {
        for (int j = 0; j < src.cols(); j++)
        {
            result(i - 1, j) = src(i, j) * i;
        }
    }
    return result;
}

MatrixXd ExCoefficients::diff_about_y(const MatrixXd &src)
{
    if (src.cols() == 1)
    {
        MatrixXd result(1, 1);
        result << 0;
        return result;
    }
    MatrixXd result(src.rows(), src.cols() - 1);
    for (int i = 0; i < src.rows(); ++i)
    {
        for (int j = 1; j < src.cols(); ++j)
        {
            result(i, j - 1) = src(i, j) * j;
        }
    }
    return result;
}

IntervalMat<itv> ExCoefficients::F_prime_value(const itv &xrange, const itv &yrange)
{
    return IntervalMat<itv>(estimate_range(m_mat_f1x, xrange, yrange),
                            estimate_range(m_mat_f1y, xrange, yrange),
                            estimate_range(m_mat_f2x, xrange, yrange),
                            estimate_range(m_mat_f2y, xrange, yrange));
}

IntervalVec<itv> ExCoefficients::calculate_K(const IntervalVec<itv> &initial)
{
    IntervalMat<itv> F_prime_value_result = F_prime_value(initial.x, initial.y);
    Region Y = F_prime_value_result.mid_value().inverse_mat();
    vec2d mX = vec2d(mid(initial.x), mid(initial.y));
    vec2d FmX = F_value(itv(mX.x, mX.x), itv(mX.y, mX.y)).mid_vec();
    IntervalMat<itv> R = Region(1.0, 0.0, 0.0, 1.0) - Y*F_prime_value_result;
    IntervalVec<itv> delta = initial - mX;
    IntervalVec<itv> K = (mX - Y*FmX) + R*delta;
    return K;
}