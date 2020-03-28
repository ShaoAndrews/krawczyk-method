
#include "quadtree.h"
quadtree::quadtree(const Eigen::MatrixXd &f1, const Eigen::MatrixXd &f2, const IntervalVec<itv> &range) : m_mat_f1(f1), m_mat_f2(f2), m_min_size(1e-3), count(0)
{
    m_mat_f1x = diff_about_x(m_mat_f1);
    m_mat_f1y = diff_about_y(m_mat_f1);
    m_mat_f2x = diff_about_x(m_mat_f2);
    m_mat_f2y = diff_about_y(m_mat_f2);
}

quadtree::~quadtree()
{
    IntervalVec<itv> *current = nullptr;
    //清空Tlist
    while (!Tlist.empty())
    {
        current = Tlist.front();
        delete current;
        current = nullptr;
        Tlist.pop();
    }
}
void quadtree::search_procedure(IntervalVec<itv> range)
{
    IntervalVec<itv> *current = nullptr;
    //清空Tlist
    while (!Tlist.empty())
    {
        current = Tlist.front();
        delete current;
        current = nullptr;
        Tlist.pop();
    }
    Plist.clear();
    current = new IntervalVec<itv>(range);
    while (1)
    {
        IntervalVec<itv> F_X = F_value(current->x, current->y);
        if (!F_X.contain_zero())
        {
            delete current;
            current = nullptr;
            if (Tlist.empty())
            {
                break;
            }
            current = Tlist.front();
            Tlist.pop();
        }
        else
        {
            IntervalMat<itv> F_prime_value_result = F_prime_value(current->x, current->y);
            Region m_F_prime = F_prime_value_result.mid_value();
            if (m_F_prime.is_reversible())
            {
                //可逆
                Region Y = m_F_prime.inverse_mat();
                IntervalMat<itv> R = Region(1.0, 0.0, 0.0, 1.0) - Y * F_prime_value_result;
                vec2d mX = vec2d(mid(current->x), mid(current->y));
                IntervalVec<itv> delta = *current - mX;
                vec2d FmX = F_value(itv(mX.x, mX.x), itv(mX.y, mX.y)).mid_vec();
                IntervalVec<itv> K = (mX - Y * FmX) + R * delta;
                if (K.is_empty_for_intersection(*current))
                {
                    //交集为空
                    delete current;
                    current = nullptr;
                    if (Tlist.empty())
                    {
                        break;
                    }
                    current = Tlist.front();
                    Tlist.pop();
                }
                else if (K.is_contained_in(*current))
                {
                    if (R.get_norm() < 1)
                    {
                        std::cout << *current << std::endl;
                        delete current;
                        current = nullptr;
                        while (!Tlist.empty())
                        {
                            current = Tlist.front();
                            delete current;
                            current = nullptr;
                            Tlist.pop();
                        }

                        break;
                    }
                    else
                    {
                        delete current;
                        current = nullptr;
                        search_procedure(K);

                        while (!Tlist.empty())
                        {
                            current = Tlist.front();
                            delete current;
                            current = nullptr;
                            Tlist.pop();
                        }

                        break;
                    }
                }
                else
                {
                    Bisection(current, F_prime_value_result);
                }
            }
            else
            {
                Bisection(current, F_prime_value_result);
            }
        }
    }

    // while (!Tlist.empty())
    // {
    //     current = Tlist.front();
    //     Tlist.pop();
    //     IntervalVec<itv> F_X = F_value(current->x, current->y);
    //     if (!F_X.contain_zero())
    //     {
    //         delete current;
    //         current = nullptr;
    //         continue;
    //     }
    //     IntervalMat<itv> F_prime_value_result = F_prime_value(current->x, current->y);
    //     Region m_F_prime = F_prime_value_result.mid_value();
    //     if (m_F_prime.is_reversible())
    //     {
    //         //可逆
    //         Region Y = m_F_prime.inverse_mat();
    //         IntervalMat<itv> R = Region(1.0, 0.0, 0.0, 1.0) - Y * F_prime_value_result;
    //         vec2d mX = vec2d(mid(current->x), mid(current->y));
    //         IntervalVec<itv> delta = *current - mX;
    //         vec2d FmX = F_value(itv(mX.x, mX.x), itv(mX.y, mX.y)).mid_vec();
    //         IntervalVec<itv> K = (mX - Y * FmX) + R * delta;
    //         if (K.is_empty_for_intersection(*current))
    //         {
    //             //交集为空
    //             delete current;
    //             current = nullptr;
    //             continue;
    //         }
    //         else if (K.is_contained_in(*current))
    //         {
    //             std::cout << R.get_norm() << std::endl;
    //             if (R.get_norm() < 1)
    //             {
    //                 std::cout << *current << std::endl;
    //                 delete current;
    //                 current = nullptr;
    //                 //清空Tlist
    //                 while (!Tlist.empty())
    //                 {
    //                     current = Tlist.front();
    //                     delete current;
    //                     current = nullptr;
    //                     Tlist.pop();
    //                 }
    //                 break;
    //             }
    //             else
    //             {
    //                 delete current;
    //                 current = nullptr;
    //                 search_procedure(K);
    //             }
    //         }
    //         else
    //         {

    //             Bisection(current, F_prime_value_result);
    //             delete current;
    //             current = nullptr;
    //         }
    //     }
    //     else
    //     {
    //         //不可逆
    //         Bisection(current, F_prime_value_result);
    //         delete current;
    //         current = nullptr;
    //     }
    // }
}

itv quadtree::estimate_range(const Eigen::MatrixXd &src, const itv &xrange, const itv &yrange)
{
    itv result(0, 0);
    for (int i = 0; i < src.rows(); ++i)
    {
        for (int j = 0; j < src.cols(); ++j)
        {
            itv temp;
            temp = pow(xrange, i) * pow(yrange, j) * src(i, j);
            double temt = src(i, j);
            result += pow(xrange, i) * pow(yrange, j) * src(i, j);
        }
    }
    return itv(result.lower(), result.upper());
}

IntervalVec<itv> quadtree::F_value(const itv &xrange, const itv &yrange)
{
    return IntervalVec<itv>(estimate_range(m_mat_f1, xrange, yrange), estimate_range(m_mat_f2, xrange, yrange));
}

IntervalMat<itv> quadtree::F_prime_value(const itv &xrange, const itv &yrange)
{
    return IntervalMat<itv>(estimate_range(m_mat_f1x, xrange, yrange),
                            estimate_range(m_mat_f1y, xrange, yrange),
                            estimate_range(m_mat_f2x, xrange, yrange),
                            estimate_range(m_mat_f2y, xrange, yrange));
}

Eigen::MatrixXd quadtree::diff_about_x(const Eigen::MatrixXd &src)
{
    if (src.rows() == 1)
    {
        Eigen::MatrixXd result(1, 1);
        result << 0;
        return result;
    }
    Eigen::MatrixXd result(src.rows() - 1, src.cols());
    for (int i = 1; i < src.rows(); i++)
    {
        for (int j = 0; j < src.cols(); j++)
        {
            result(i - 1, j) = src(i, j) * i;
        }
    }
    return result;
}

Eigen::MatrixXd quadtree::diff_about_y(const Eigen::MatrixXd &src)
{
    if (src.cols() == 1)
    {
        Eigen::MatrixXd result(1, 1);
        result << 0;
        return result;
    }
    Eigen::MatrixXd result(src.rows(), src.cols() - 1);
    for (int i = 0; i < src.rows(); ++i)
    {
        for (int j = 1; j < src.cols(); ++j)
        {
            result(i, j - 1) = src(i, j) * j;
        }
    }
    return result;
}

void quadtree::Bisection(IntervalVec<itv> *p, const IntervalMat<itv> &Fprime)
{
    if (width(p->x) < m_min_size && width(p->y) < m_min_size)
    {
        Plist.push_back(*p);
        delete p;
        p = nullptr;
        return;
    }
    std::pair<int, int> index = Fprime.find_mini_ij();

    Eigen::MatrixXd Fprime_ij_mat;
    if (index.first == 1 && index.second == 1)
    {
        Fprime_ij_mat = m_mat_f1x;
    }
    if (index.first == 1 && index.second == 2)
    {
        Fprime_ij_mat = m_mat_f1y;
    }
    if (index.first == 2 && index.second == 1)
    {
        Fprime_ij_mat = m_mat_f2x;
    }
    if (index.first == 2 && index.second == 2)
    {
        Fprime_ij_mat = m_mat_f2y;
    }
    IntervalVec<itv> *first_try1 = nullptr;
    IntervalVec<itv> *first_try2 = nullptr;
    IntervalVec<itv> *second_try1 = nullptr;
    IntervalVec<itv> *second_try2 = nullptr;
    if (count % 2 == 0)
    {
        //先对x进行分割
        first_try1 = new IntervalVec<itv>(itv(p->x.lower(), mid(p->x)), p->y);
        first_try2 = new IntervalVec<itv>(itv(mid(p->x), p->x.upper()), p->y);

        second_try1 = new IntervalVec<itv>(p->x, itv(p->y.lower(), mid(p->y)));
        second_try2 = new IntervalVec<itv>(p->x, itv(mid(p->y), p->y.upper()));
    }
    else
    {
        //先对y进行分割

        first_try1 = new IntervalVec<itv>(p->x, itv(p->y.lower(), mid(p->y)));
        first_try2 = new IntervalVec<itv>(p->x, itv(mid(p->y), p->y.upper()));

        second_try1 = new IntervalVec<itv>(itv(p->x.lower(), mid(p->x)), p->y);
        second_try2 = new IntervalVec<itv>(itv(mid(p->x), p->x.upper()), p->y);
    }

    itv first_try1_range = estimate_range(Fprime_ij_mat, first_try1->x, first_try1->y);
    itv first_try2_range = estimate_range(Fprime_ij_mat, first_try2->x, first_try2->y);
    itv second_try1_range = estimate_range(Fprime_ij_mat, second_try1->x, second_try1->y);
    itv second_try2_range = estimate_range(Fprime_ij_mat, second_try2->x, second_try2->y);

    double try1_width = width(itv::hull(first_try1_range, first_try2_range));
    double try2_width = width(itv::hull(second_try1_range, second_try2_range));
    if (try1_width <= try2_width)
    {
        IntervalVec<itv> F_range1 = F_value(first_try1->x, first_try1->y);
        IntervalVec<itv> F_range2 = F_value(first_try2->x, first_try2->y);
        double temp1 = fabs(mid(F_range1.x)) + fabs(mid(F_range1.y));
        double temp2 = fabs(mid(F_range2.x)) + fabs(mid(F_range2.y));
        if (temp1 > temp2)
        {
            delete p;
            p = nullptr;
            p = first_try2;
            Tlist.push(first_try1);
        }
        else
        {
            delete p;
            p = nullptr;
            p = first_try1;
            Tlist.push(first_try2);
        }
        delete second_try1;
        delete second_try2;
        second_try1 = nullptr;
        second_try2 = nullptr;
    }
    else
    {

        IntervalVec<itv> F_range1 = F_value(second_try1->x, second_try1->y);
        IntervalVec<itv> F_range2 = F_value(second_try2->x, second_try2->y);
        double temp1 = fabs(mid(F_range1.x)) + fabs(mid(F_range1.y));
        double temp2 = fabs(mid(F_range2.x)) + fabs(mid(F_range2.y));
        if (temp1 > temp2)
        {
            delete p;
            p = nullptr;
            p = second_try2;
            Tlist.push(second_try1);
        }
        else
        {
            delete p;
            p = nullptr;
            p = second_try1;
            Tlist.push(second_try2);
        }
        delete first_try1;
        delete first_try2;
        first_try1 = nullptr;
        first_try2 = nullptr;
    }
}