#include "point.hpp"

template<typename T>
struct knode{
    point<T> _point;
    knode<T> * _left;
    knode<T> * _right;
};