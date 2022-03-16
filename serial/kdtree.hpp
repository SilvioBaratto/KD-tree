#include "utility.hpp"

/**
*   @brief this class perform a serial KD-TREE 
*/
template<typename coordinate, std::size_t dimension>
class kdtree{

    public:

        /**
        *@brief the constructor takes in input a filename and save in a vector of type
        *       knode the values of the filename to be processed
        *@param filename
        */
        explicit kdtree(std::string filename):
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

/**
*   @brief main function to produce the tree in serial. 
*   @param begin the beginnin of the dataset
*   @param end the end of the dataset
*   @param index the index axis of the node
*   @return the kdtree
*/

template<typename coordinate, std::size_t dimension>
knode<coordinate> * kdtree<coordinate, dimension>::make_tree(std::size_t begin, std::size_t end, std::size_t index){
    //stop criterion
    if(end <= begin) return nullptr;
    //take the median point in the dataset
    std::size_t med = begin + (end - begin) / 2;
    //save the first position in memory of the dataset
    auto i = _knodes.begin();
    // - nth_element is a partial sorting algorithm that rearranges elementes in [first, last) such that:
    //      a. The element at the nth position is the one which should be at the position if we sort the list.
    //      b. It does not sort the list, just all the elements, which precede the nth element are not greater than it,
    // - nth_element algorithm is implemented using introselect.
    //      a. introselect is a hybrid of quickselect and median of medians algorithm
    //          1. quickselect is used to find the kth smalles number in an unsorted array
    //          2. median of medians is a median selection algorithm for better pivot selection maily used in quickselect
    // 
    // In our case we use this method to sort the array up to the median and take the value at the median point
    std::nth_element(i + begin, i + med, i + end, knode_cmp<coordinate>(index));
    _knodes[med]._axis = index;
    index = (index + 1) % DIM;
    // build left part
    _knodes[med]._left = make_tree(begin, med, index);
    // build right part
    _knodes[med]._right = make_tree(med + 1, end, index);
    // return the pointer 
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

/**
*   @brief this function get the values from the filename and stores it to a vector of knode.
*          the function is templated so can store both integer and floating point numbers
*   @param filename
*   @return knodes vector of integers or floating point numbers
*/

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
                    data.push_back(knode<coordinate>(point<coordinate, dimension>{{stof(X), stof(Y)}}));
                #endif
            }
            i++;
        }
        coeff.close();
    } else std::cout << "Unable to open file" << std::endl;

    return data;
}

/**
*   @brief debug function to output the kdtree complete in a simple way
*/

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

/**
*   @brief debug function to output the kdtree complete in a simple way
*   @param node tree created using make_tree function
*/

template<typename T>
void print_tree(knode<T> * node) {
    print_tree(node, "", false);
}

