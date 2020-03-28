#pragma once
#include <iostream>
#include <assert.h>
#include <cmath>
class vec2d
{
    friend std::ostream &operator<<(std::ostream &out, vec2d &vec){
        out << "[ "<<vec.x<<", "<<vec.y<<" ]"<<std::endl;
    }
    friend vec2d operator-(const vec2d &v1, const vec2d &v2)
    {
        return vec2d(v1.x - v2.x, v1.y - v2.y);
    }

    friend vec2d operator+(const vec2d &v1, const vec2d &v2)
    {
        return vec2d(v1.x + v2.x, v1.y + v2.y);
    }

public:
    double x;
    double y;
    vec2d(double x_, double y_) : x(x_), y(y_){};
    vec2d(const vec2d &vec)
    {
        x = vec.x;
        y = vec.y;
    }
    vec2d(){};
};

class Region
{
public:
    friend std::ostream &operator<<(std::ostream &out, Region &A);
    friend vec2d operator*(const Region &region, const vec2d &vec)
    {
        return vec2d(region.a11 * vec.x + region.a12 * vec.y, region.a21 * vec.x + region.a22 * vec.y);
    }

public:
    /* data */
    double a11;
    double a12;
    double a21;
    double a22;

public:
    Region(double xl, double xr, double yl, double yr);
    Region();
    const Region operator*(const Region &region) const;
    const Region operator-(const vec2d &vec) const;
    ~Region();
    double maximum_norm();                                //无穷范数
    vec2d measure_m();                                    //每个分量的宽度
    bool is_contained_in(const Region &region);           //是否包含于region中
    bool is_empty_for_intersection(const Region &region); //两个区域相交是否为空
    Region inverse_mat();                                 //逆矩阵
    bool is_reversible();                                 //是否可逆
};
