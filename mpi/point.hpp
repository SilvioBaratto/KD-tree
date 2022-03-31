#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#define DIM 2

/**
*   @brief this class define a point of the knode data structure
*   @param array array of n dimension x: position[0], y: position[1] etc.
*/
template<typename coordinate, std::size_t dimension>
class point{
    public:
        point(std::array<coordinate, dimension> points):
            _points{points} {}
        point(std::initializer_list<coordinate> list){
            std::size_t n = std::min(dimension, list.size());
            std::copy_n(list.begin(), n, _points.begin());
        }

        point():
            _points{0, 0} {}

        coordinate get(std::size_t index) const{
            return _points[index];
        }

        // return the x value
        coordinate x() const {return _points[0];}

        // return the y value
        coordinate y() const {return _points[1];}

        std::string print_point(int axis);

        double distance(const point& pt) const {
            double dist = 0;
            for (size_t i = 0; i < dimension; ++i) {
                double d = get(i) - pt.get(i);
                dist += d * d;
            }
            return dist;
        }

    private:
        std::array<coordinate, dimension> _points;
       
};

template<typename coordinate, std::size_t dimension>
std::ostream& operator<<(std::ostream& os, const point<coordinate, dimension>& pt) {
    os << '(';
    for (size_t i = 0; i < dimension; ++i) {
        if (i > 0)
            os << ", ";
        os << pt.get(i);
    }
    os << ')';
    return os;
}

template<typename coordinate, std::size_t dimension>
std::string point<coordinate, dimension>::print_point(int axis){
    std::string point = "";
    point += "[";

    for(auto i = 0; i < dimension; i++){
        point += std::to_string(_points[i]);
        point += ",";
    }
    point += ";";
    point += std::to_string(axis);
    point += "]";

    return point;
}