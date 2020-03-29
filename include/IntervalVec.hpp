#include "Region.h"
#include <ostream>
class Interval
{
    //输出到控制台
    friend std::ostream &operator<<(std::ostream &out, Interval &c)
    {
        out << "[ " << c.lower << ", " << c.upper << " ]" << std::endl;
    }
    //数乘
    friend Interval operator*(double lambda, const Interval &iv)
    {
        return Interval(lambda * iv.lower, lambda * iv.upper);
    }

    //与实数相减
    friend Interval operator-(double value, const Interval &iv)
    {
        return Interval(value - iv.lower, value - iv.upper);
    }

    //与实数相加
    friend Interval operator+(double value, const Interval &iv)
    {
        return Interval(value + iv.lower, value + iv.upper);
    }

    //两个区间相加
    friend Interval operator+(const Interval &iv1, const Interval &iv2)
    {
        return Interval(iv1.lower + iv2.lower, iv1.upper + iv2.upper);
    }

    //两个区间相减
    friend Interval operator-(const Interval &iv1, const Interval &iv2)
    {
        return Interval(iv1.lower - iv2.upper, iv1.upper - iv2.lower);
    }

public:
    Interval(double a, double b) : lower(a), upper(b)
    {
    }
    Interval(const Interval &I)
    {
        lower = I.lower;
        upper = I.upper;
    }
    Interval() : lower(0), upper(1) {}
    ~Interval()
    {
    }
    double length() const
    {
        return (this->upper - this->lower);
    }
    bool is_subset_of(const Interval &iv) const
    {
        if (lower >= iv.lower && upper <= iv.upper)
        {
            return true;
        }
        return false;
    }
    bool is_empty_by_intersection(const Interval &iv) const
    {
        if (lower > iv.upper || upper < iv.lower)
        {
            return true;
        }
        return false;
    }
    double mid_value() const
    {
        return (lower + upper) / 2.0;
    }

public:
    double lower;
    double upper;
};

template <typename T>
class IntervalVec
{

    //输出到控制台
    friend std::ostream &operator<<(std::ostream &out, IntervalVec &c)
    {
        out << " "<<c.x << "*"<< c.y;
    }
    //实数矩阵与该区间向量相乘
    friend IntervalVec operator*(const Region &region, const IntervalVec &c)
    {
        return IntervalVec(region.a11 * c.x + region.a12 * c.y, region.a21 * c.x + region.a22 * c.y);
    }
    //该区间向量相减
    friend IntervalVec operator-(const IntervalVec &v1, const IntervalVec &v2)
    {

        return IntervalVec(v1.x - v2.x, v1.y - v2.y);
    }
    //该区间向量相加
    friend IntervalVec operator+(const IntervalVec &v1, const IntervalVec &v2)
    {

        return IntervalVec(v1.x + v2.x, v1.y + v2.y);
    }
    friend IntervalVec operator-(const IntervalVec &v1, const vec2d &v2)
    {
        return IntervalVec(v1.x - v2.x, v1.y - v2.y);
    }

    friend IntervalVec operator+(const IntervalVec &v1, const vec2d &v2)
    {
        return IntervalVec(v1.x + v2.x, v1.y + v2.y);
    }

    friend IntervalVec operator-(const vec2d &v1, const IntervalVec &v2)
    {
        return IntervalVec(v1.x - v2.x, v1.y - v2.y);
    }

    friend IntervalVec operator+(const vec2d &v1, const IntervalVec &v2)
    {
        return IntervalVec(v1.x + v2.x, v1.y + v2.y);
    }

public:
    IntervalVec(const T &x_, const T &y_) : x(x_), y(y_)
    {
    }
    ~IntervalVec()
    {
    }
    vec2d mid_vec()
    {
        return vec2d((x.lower() + x.upper()) / 2.0, (y.lower() + y.upper()) / 2.0);
    }
    bool is_contained_in(IntervalVec<T> &target)
    {
        if (subset(x, target.x) && subset(y, target.y))
        {
            return true;
        }
        return false;
    }

    bool is_empty_for_intersection(IntervalVec<T> &target)
    {
        if ((!overlap(this->x, target.x) || (!overlap(this->y, target.y))))
        {
            return true;
        }
        return false;
    }

    double get_width() const
    {
        double x_width = width(this->x);
        double y_width = width(this->y);
        return x_width > y_width ? x_width : y_width;
    }
    double get_norm() const
    {
        double xl = fabs(x.lower());
        double xr = fabs(x.upper());
        double yl = fabs(y.lower());
        double yr = fabs(y.upper());
        double xmax = xl > xr ? xl : xr;
        double ymax = yl > yr ? yl : yr;
        return xmax > ymax ? xmax : ymax;
    }
    bool contain_zero() const
    {
        if (in(0, this->x) && in(0, this->y))
        {
            return true;
        }
        return false;
    }

public:
    T x;
    T y;
};

template <typename T>
class IntervalMat
{
    friend IntervalMat operator-(const Region &region, const IntervalMat &invmat)
    {
        return IntervalMat(region.a11 - invmat.a11, region.a12 - invmat.a12, region.a21 - invmat.a21, region.a22 - invmat.a22);
    }
    friend IntervalMat operator*(const Region &region, const IntervalMat &invmat)
    {
        return IntervalMat(region.a11 * invmat.a11 + region.a12 * invmat.a21,
                           region.a11 * invmat.a12 + region.a12 * invmat.a22,
                           region.a21 * invmat.a11 + region.a22 * invmat.a21,
                           region.a21 * invmat.a12 + region.a22 * invmat.a22);
    }
    friend IntervalVec<T> operator*(const IntervalMat &mat, const IntervalVec<T> &vec)
    {
        return IntervalVec<T>(mat.a11 * vec.x + mat.a12 * vec.y, mat.a21 * vec.x + mat.a22 * vec.y);
    }
    friend std::ostream &operator<<(std::ostream &out, IntervalMat &c)
    {
        out << "{"
            << "\n"
            << "\t" << c.a11 << "\t" << c.a12 << "\t"
            << "\t" << c.a21 << "\t" << c.a22 << "\n"
            << "}" << std::endl;
    }

public:
    IntervalMat(T a11_, T a12_, T a21_, T a22_) : a11(a11_), a12(a12_), a21(a21_), a22(a22_)
    {
    }
    Region mid_value()
    {
        return Region((a11.lower() + a11.upper()) / 2.0, (a12.lower() + a12.upper()) / 2.0, (a21.lower() + a21.upper()) / 2.0, (a22.lower() + a22.upper()) / 2.0);
    }
    double get_norm()
    {
        double line1 = norm(a11) + norm(a12);
        double line2 = norm(a21) + norm(a22);
        return line1 > line2 ? line1 : line2;
    }
    std::pair<int, int> find_mini_ij() const
    {
        double a11_numeric = width(a11);
        double a12_numeric = width(a12);
        double a21_numeric = width(a21);
        double a22_numeric = width(a22);
        if(a11_numeric >= a12_numeric && a11_numeric >= a21_numeric && a11_numeric >= a22_numeric){
            return std::make_pair(1, 1);
        }
        if(a12_numeric >= a11_numeric && a12_numeric >= a21_numeric && a12_numeric >= a22_numeric){
            return std::make_pair(1, 2);
        }
        if(a21_numeric >= a11_numeric && a21_numeric >= a12_numeric && a21_numeric >= a22_numeric){
            return std::make_pair(2, 1);
        }
        if(a22_numeric >= a11_numeric && a22_numeric >= a12_numeric && a22_numeric >= a21_numeric){
            return std::make_pair(2, 2);
        }
    }
    T &operator()(int i, int j) const
    {
        assert(i >= 1 && i <= 2 && j >= 1 && j <= 2);
        if(i==1 && j == 1){
            return a11;
        }
        if(i==1 && j == 2){
            return a12;
        }
        if(i==2 && j == 1){
            return a21;
        }
        if(i==2 && j == 2){
            return a22;
        }
    }
public:
    T a11;
    T a12;
    T a21;
    T a22;
};