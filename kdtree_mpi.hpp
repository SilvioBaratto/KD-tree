#include "utility.hpp"

#define DIM 2;

template<typename T>
class kdtree{
    public:      
        kdtree(std::string filename): 
            _values{getValues(filename)}{
                _root = make_tree_parallel(0, _values.size(), 0);
            }
    
        std::vector<knode<T>> get_data() {return _values;}
        void printData();
        void printRootMemory();
        void sortData();
        int get_size() {return _values.size();}
        knode<T> * get_vector_data() {return _values.data();}
        knode<T> * get_root() {return _root;}
        knode<T> * getMedian(std::size_t begin, std::size_t med, std::size_t end, std::size_t index);
        knode<T> * make_tree_parallel(int begin, int end, int index);

        //print root
        void print_tree(knode<T> * node, const std::string &prefix, bool isLeft);
        void print_tree(knode<T> * node);
        void printRoot() {print_tree(_root);}

        std::vector<knode<T>> getValues(std::string filename);
        knode<T> * make_tree(int begin, int end, int index);


    private:
        knode<T> * _root;
        std::vector<knode<T>> _values;
};

template<typename T>
struct node_cmp{
    node_cmp(size_t index) : _index(index) {}
    bool operator()(knode<T>& n1, knode<T>& n2) const {
        return n1._point.index(_index) < n2._point.index(_index);
    }
    size_t _index;
};


template<typename T>
knode<T> * kdtree<T>::make_tree(int begin, int end, int index){
    if(end <= begin) return 0;
    size_t med = begin + (end - begin) / 2;
    knode<T> * tree;
    tree = getMedian(begin, med, end, index);
    index = (index + 1) % DIM
    #pragma omp task
    tree->_left = make_tree(begin, med, index);
    #pragma omp task
    tree->_right = make_tree(med + 1, end, index);  
    return tree;
}

template<typename T>
knode<T> * kdtree<T>::make_tree_parallel(int begin, int end, int index){
    knode<T> * root;
    #pragma omp parallel
    { 
        #pragma omp single 
        {
            #pragma omp task shared(root)
            root = make_tree(begin, end, index);
        }
    }
    #pragma omp barrier
    return root;
}

template<typename T>
void kdtree<T>::printData(){
    for(int i = 0; i < _values.size(); i++){
        std::cout << "(" << _values[i]._point.x() << ", " << _values[i]._point.y() << ")";
    }
    std::cout << std::endl;
}

template<typename T>
void kdtree<T>::printRootMemory(){
    std::cout << std::endl;
    for(int i = 0; i <_values.size(); i++){
        std::cout << "(" << _root[i]._left << ", " << _root[i]._right << ")" << std::endl;
    }
}

template<typename T>
knode<T> * kdtree<T>::getMedian(std::size_t begin, std::size_t med, std::size_t end, std::size_t index){
    auto i = _values.begin();
    std::sort(i + begin, i + end, node_cmp<T>(index));
    return &_values[med];
}


template<typename T>
std::vector<knode<T>> kdtree<T>::getValues(std::string filename){
    std::string X, Y;

    int size = getDataRows(filename);
    std::vector<knode<T>> data;

    std::ifstream coeff(filename);
    int i = 0;
    if (coeff.is_open()) {
        while(!coeff.eof()) {
            getline(coeff, X, ',');
            getline(coeff, Y, '\n');
            if(!X.empty() && !Y.empty()){
                data.push_back({{stoi(X), stoi(Y)}, nullptr, nullptr});
            }
            i++;
        }
        coeff.close();
    } else std::cout << "Unable to open file";
    return data;
}

template<typename T>
void kdtree<T>::print_tree(knode<T> * node, const std::string &prefix, bool isLeft) {
                    
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

/*
    Start the recursion from the root node.
*/
template<typename T>
void kdtree<T>::print_tree(knode<T> * node) {
    print_tree(node, "", false);
}