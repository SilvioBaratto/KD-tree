#include "utility.hpp"

template<typename coordinate, std::size_t dimension>
class kdtree{

    public:

        kdtree(std::string filename):
        _knodes{getValues(filename)} {
            auto time = omp_get_wtime();
            _root = make_tree(0, _knodes.size(), 0);
            time = omp_get_wtime() - time; 
            std::cout << omp_get_max_threads() << "," << time << std::endl;
        }

        knode<coordinate> * make_tree(std::size_t begin, std::size_t end, std::size_t index);
        knode<coordinate> * get_root() {return _root;}
        knode<coordinate> * make_tree_parallel(std::size_t begin, std::size_t end, std::size_t index);

        std::vector<knode<coordinate>> getValues(std::string filename);
        void nearest(knode<coordinate> * root, const point<coordinate, dimension>& point, size_t index);
        const point<coordinate, dimension>& nearest(const point<coordinate, dimension>& pt);

        bool empty() const {return _knodes.empty();}
        std::size_t visited() const {return _visited;}
        double distance() const {return std::sqrt(_best_dist);}

    private:
        knode<coordinate> * _root = nullptr;
        std::vector<knode<coordinate>> _knodes;
        knode<coordinate> * _best = nullptr;
        std::size_t _visited;
        double _best_dist = 0;
};

template<typename coordinate, std::size_t dimension>
knode<coordinate> * kdtree<coordinate, dimension>::make_tree(std::size_t begin, std::size_t end, std::size_t index){
    if(end <= begin) return nullptr;

    std::size_t med = begin + (end - begin) / 2;
    auto i = _knodes.begin();
    std::nth_element(i + begin, i + med, i + end, knode_cmp<coordinate>(index));
    _knodes[med]._axis = index;
    index = (index + 1) % DIM;

    _knodes[med]._left = make_tree(begin, med, index);

    _knodes[med]._right = make_tree(med + 1, end, index);

    return &_knodes[med];
}

template<typename coordinate, std::size_t dimension>
void kdtree<coordinate, dimension>::nearest(knode<coordinate> * root, const point<coordinate, dimension>& point, 
                                    size_t index) {
    if (root == nullptr)
        return;
    ++_visited;
    double d = root->distance(point);
    if (_best == nullptr || d < _best_dist) {
        _best_dist = d;
        _best = root;
    }
    if (_best_dist == 0)
        return;
    double dx = root->get(index) - point.get(index);
    index = (index + 1) % dimension;
    nearest(dx > 0 ? root-> _left : root-> _right, point, index);
    if (dx * dx >= _best_dist)
        return;
    nearest(dx > 0 ? root->_right : root->_left, point, index);
}

template<typename coordinate, std::size_t dimension>
const point<coordinate, dimension>& kdtree<coordinate, dimension>::nearest(const point<coordinate, dimension>& pt) {
    if (_root == nullptr)
        throw std::logic_error("tree is empty");
    _best = nullptr;
    _visited = 0;
    _best_dist = 0;
    nearest(_root, pt, 0);
    return _best->_point;
}

template<typename coordinate, std::size_t dimension>
std::vector<knode<coordinate>> kdtree<coordinate, dimension>::getValues(std::string filename){
    std::string X, Y;

    int size = getDataRows(filename);
    std::vector<knode<coordinate>> data;

    std::ifstream coeff(filename);
    int i = 0;
    if (coeff.is_open()) {
        while(!coeff.eof()) {
            getline(coeff, X, ',');
            getline(coeff, Y, '\n');
            if(!X.empty() && !Y.empty()){
                #ifdef int_data
                    data.push_back(knode<coordinate>(point<coordinate, dimension>{{stoi(X), stoi(Y)}}));
                #endif
                #ifdef double_data
                    data.push_back(knode<coordinate>(point<coordinate, dimension>{{stod(X), stod(Y)}}));
                #endif
            }
            i++;
        }
        coeff.close();
    } else std::cout << "Unable to open file" << std::endl;

    return data;
}


template<typename T>
void print_tree(knode<T> * node, const std::string &prefix, bool isLeft) {
                    
    if(node != nullptr){

        std::cout << prefix;
         std::cout << (isLeft ? "├──" : "└──" );

        // print the value of the node
        std::cout << "(" << node->_point.x() << ", " << node->_point.y() << ")" << std::endl;

        // enter the next tree level - left and right branch
        if (node->_left) print_tree(node->_left, prefix + (isLeft ? "│   " : "    "), true);

        if (node->_right) print_tree(node->_right, prefix + (isLeft ? "│   " : "    "), false);
    }
}

template<typename T>
void print_tree(knode<T> * node) {
    print_tree(node, "", false);
}

