#include "point.hpp"

template<typename coordinate>
struct knode{

    typedef point<coordinate, DIM> point_type;

    knode(const point_type& pt):
        _point{pt}, _left{nullptr}, _right{nullptr}, _axis{-1} {}

    coordinate get(std::size_t index) const{
        return _point.get(index);
    }

    double distance(const point_type& pt) const{
        return _point.distance(pt);
    }

    point_type _point;
    int _axis;
    knode<coordinate> * _left;
    knode<coordinate> * _right;
};

template<typename T>
struct knode_cmp{
    knode_cmp(std::size_t index): 
        _index{index} {}
    bool operator()(const knode<T>& n1, const knode<T>& n2) const{
        return n1._point.get(_index) < n2._point.get(_index);
    }

    std::size_t _index;
};