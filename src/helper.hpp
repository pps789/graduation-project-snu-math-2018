#ifndef __HELPER__H__
#define __HELPER__H__

#include<algorithm>
const double DOUBLE_EPS = 1e-6;

double nonzero(double x){
    return std::abs(x)>DOUBLE_EPS;
}

#endif
