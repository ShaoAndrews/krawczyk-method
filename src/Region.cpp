#include "Region.h"

Region::Region(double xl, double xr, double yl, double yr) : a11(xl),
                                                             a12(xr),
                                                             a21(yl),
                                                             a22(yr)
{
}
Region::Region() : a11(0),
                   a12(1),
                   a21(0),
                   a22(1)
{
}
Region::~Region()
{
}

double Region::maximum_norm()
{

    double delta_x = fabs(this->a12) + fabs(this->a11);
    double delta_y = fabs(this->a22) + fabs(this->a21);
    return delta_x > delta_y ? delta_x : delta_y;
}
const Region Region::operator*(const Region &region) const
{
    return Region(this->a11 * region.a11 + this->a12 * region.a21,
                  this->a11 * region.a12 + this->a12 * region.a22,
                  this->a21 * region.a11 + this->a22 * region.a21,
                  this->a21 * region.a12 + this->a22 * region.a22);
}

vec2d Region::measure_m()
{
    vec2d result;
    result.x = this->a12 - this->a11;
    result.y = this->a22 - this->a21;
    return result;
}

const Region Region::operator-(const vec2d &vec) const
{
    return Region(this->a11 - vec.x, this->a12 - vec.x, this->a21 - vec.y, this->a22 - vec.y);
}

bool Region::is_contained_in(const Region &region)
{
    if (this->a11 >= region.a11 && this->a12 <= region.a12 && this->a21 >= region.a21 && this->a22 <= region.a22)
    {
        return true;
    }
    return false;
}
bool Region::is_empty_for_intersection(const Region &region)
{
    if ((this->a11 > region.a12 || this->a12 < region.a11) || (this->a21 > region.a22 || this->a22 < region.a21))
    {
        return true;
    }
    return false;
}
Region Region::inverse_mat()
{
    double delta = a11* a22 - a12*a21;
    return Region(a22/delta, -a12/delta, -a21/delta, a11/delta);
}
bool Region::is_reversible(){
    double delta = a11* a22 - a12*a21;
    if(fabs(delta) <1e-10){
        return false;
    }
    return true;
}