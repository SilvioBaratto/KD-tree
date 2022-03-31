#include "utility.hpp"

/**
*   @brief this class perform a MPI KD-tree construction 
*/

template<typename coordinate, std::size_t dimension>
class kdtree{

    public:

        /**
        *@brief the constructor takes in input a filename  and save in a vector of type
        *       knode the values of the filename to be processed. From the functions
        *       MPI_Comm_size and MPI_Comm_rank the constructor in the MPI version takes also
        *       the number of process and the rank
        *@param filename
        */

        kdtree(std::string filename, int nprocs, int rank):
        _knodes{getValues(filename)} {
            double mpi_time;
            mpi_time = MPI_Wtime(); 
            _root = make_tree_parallel_pointer(0, _knodes.size(), 0, nprocs, 0, MPI_COMM_WORLD, 1);
            // _root = make_tree(0, _knodes.size(), 0);
            mpi_time = MPI_Wtime() - mpi_time;
            if(rank == 0){
                std::cout << nprocs << ", " << mpi_time << std::endl;              
            }
        }

        knode<coordinate> * make_tree(std::size_t begin, std::size_t end, std::size_t index);
        knode<coordinate> * get_root() {return _root;}
        knode<coordinate> * make_tree_parallel_string(std::size_t begin, std::size_t end, 
                                                std::size_t index, int size, int depth, MPI_Comm comm, int which);
        knode<coordinate> * make_tree_parallel_pointer(std::size_t begin, std::size_t end, 
                                                std::size_t index, int nprocs, int depth, MPI_Comm comm, int next);
        knode<coordinate> * buffer(int begin, int end);

        std::vector<knode<coordinate>> getValues(std::string filename);
        void nearest(knode<coordinate> * root, const point<coordinate, dimension>& point, size_t index);
        const point<coordinate, dimension>& nearest(const point<coordinate, dimension>& pt);

        std::size_t get_size() {return _knodes.size();}

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
*   @brief This function take a pointer to a knode and return a serialize string of the knode
*   @param knode
@   @return knode to string
*/

// ******************** guglielmo try
//this function serializes the tree in an std::vector
template<typename coordinate, std::size_t dimension>
coordinate * serialize_pointer(knode<coordinate> * tree, std::vector<coordinate>& data){
	for (int i = 0; i < dimension; i++){
		data.push_back(tree->_point.get(i));
	}
	if(tree->_left != NULL) {
		serialize_pointer<coordinate, dimension>(tree->_left, data);
	}
    if(tree->_right != NULL) {
        serialize_pointer<coordinate, dimension>(tree->_right, data);
    }
    return data.data();
}

template<typename coordinate, std::size_t dimension>
point<coordinate, dimension> get_point(std::size_t begin, coordinate * data){
    point<coordinate, dimension> point;
    if(begin == 0){
        point = {data[begin], data[begin + 1]};
    }
    else{
        point = {data[begin], data[begin + 2]};
    }

    return point;
}

//this function deserializes an std::vector in  a tree

/*
template<typename coordinate, std::size_t dimension>
knode<coordinate> * deserialize_pointer(std::size_t begin, std::size_t end, 
                                coordinate * data){

	int size = end - begin;
	int rightsize = size / 2;
	int leftsize = size - rightsize - 1;

    //point<coordinate, dimension> point = get_point<coordinate, dimension>(begin, data);
    point<coordinate, dimension> point{{data[0], data[0]}};
	knode<coordinate> * tree= new knode<coordinate>{point};

	std::swap(tree->_point, point);

	if(leftsize != 0){
        tree->_left=deserialize_pointer<coordinate, dimension>(begin + rightsize + 1, end, data);
	}

	if(rightsize!=0){
		tree->_right=deserialize_pointer<coordinate, dimension>(begin + 1, begin + rightsize + 1, data);
	}

	return tree;
}
*/

template<typename coordinate, std::size_t dimension>
knode<coordinate> * deserialize_pointer(std::size_t begin, std::size_t end, 
                                        std::size_t x, std::size_t y, coordinate * data){ 

    std::size_t med = begin + (end - begin) / 2;

    int size = end - begin;

	int rightsize = size / 2;

	int leftsize = size - rightsize - 1;

    std::size_t half = end - begin;

    // point<coordinate, dimension> point = get_point<coordinate, dimension>(begin, data); 

    point<coordinate, dimension> point{{data[x], data[y]}};

    if(end <= begin){
        return nullptr;
    }

    knode<coordinate> * tree= new knode<coordinate>{point};

    // std::cout << "LEFT: (" << x << ", " << y << ")" << std::endl; 
    tree->_left = deserialize_pointer<coordinate, dimension>(begin + 1, med, x + 2, y + 2, data);
    // std::cout << "RIGHT: (" << x << ", " << y << ")" << std::endl;
    tree->_right = deserialize_pointer<coordinate, dimension>(med + 1, end, med + 1, med + 2, data);

	return tree;
}

// ******************** end guglielmo

template<typename coordinate>
std::string serialize_node(knode<coordinate> * node){

    std::string node_str = "";
        if( node -> _axis != -1) {
            node_str = node -> _point.print_point( node -> _axis);
            if ( node -> _left != NULL)
                node_str += "(" + serialize_node(node -> _left) + ")";
            if ( node -> _right != NULL)
                node_str += "(" + serialize_node(node -> _right) + ")";
        }

    return node_str;
}

/**
*   @brief this function take a serialize string of a knode and deserialize returning a pointer to a knode
*   @param string
*   @return knode pointer
*/

template<typename coordinate, std::size_t dimension>
knode<coordinate> * deserialize_node(const std::string node){

    const std::size_t size = node.size();
    int j = 0, i;
    coordinate node_point[dimension];
    std::string temp_str;
    #ifdef int_data
        int temp_val;
    #endif

    #ifdef double_data
        float temp_val;
    #endif

    if ( size == 0 ) return nullptr; 
    
    if ( node[0] == ')') return nullptr; 
    
    while ( j < size && node[j] != '(' ) 
        j ++;
    
    if ( node[1] == '[' ) temp_str = node.substr(2, j - 5);
    
    else temp_str = node.substr(1, j - 5);

    std::istringstream iss(temp_str);
    for(i = 0; i < dimension; ++i) {
        iss >> temp_val;
        node_point[i] = temp_val;
        if (iss.peek() == ',') iss.ignore();
    }
    
    point<coordinate, dimension> point{{node_point[0], node_point[1]}};

    knode<coordinate> * root = new knode<coordinate>{point};

    temp_val = (int)node[j - 2] - 48; //the numbers starts from 48 in the ASCII code
    if(temp_val > 1 || temp_val < 0) temp_val = 0;

    root -> _axis = temp_val; //where axis is saved

    int left = 0;
    i = j;

    // find separation between left and right definition
    while ( i < size ) {
        if ( node[i] == '(' ) left ++;

        if ( node[i] == ')') left --;

        if ( left == 0 ) break;

        i ++;
    }

    if ( j < size - 1) {       
        #pragma omp task  
        root-> _left = deserialize_node<coordinate, dimension>(node.substr(j + 1, i - 1 - j));
    }
    if ( i + 1 < size - 1) {      
        #pragma omp task
        root-> _right = deserialize_node<coordinate, dimension>(node.substr(i + 2, size - i - 2));   
    }  

    return root;
}

/**
*   @brief function to deserialize the string in parallel using OpenMP
*   @param string
*   @return knode pointer
*/

template<typename coordinate, std::size_t dimension>
knode<coordinate> * deserialize_node_parallel(std::string data){
    knode<coordinate> * root;
    #pragma omp parallel shared(root)
    {
        #pragma omp single
        root = deserialize_node<coordinate, dimension>(data);
    }
    #pragma omp barrier
    return root;
}

template<typename coordinate, std::size_t dimension>
knode<coordinate> * deserialize_pointer_parallel(std::size_t begin, std::size_t end, 
                                        std::size_t x, std::size_t y, coordinate * data){
    knode<coordinate> * root;
    #pragma omp parallel shared(root)
    {
        #pragma omp single
        root = deserialize_pointer<coordinate, dimension>(begin, end, x, y, data);
    }
    #pragma omp barrier
    return root;
}

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
    #pragma omp task
    _knodes[med]._left = make_tree(begin, med, index);
    // build right part
    #pragma omp task
    _knodes[med]._right = make_tree(med + 1, end, index);
    // return the pointer
    return &_knodes[med];
}

/*
knode<coordinate> * kdtree<coordinate, dimension>::make_tree_parallel(std::size_t begin, std::size_t end, 
                                                std::size_t index, int nprocs, int depth, MPI_Comm comm, int next){

    MPI_Request request;
    MPI_Status status;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    depth = initdepth(rank);

    if(end <= begin) return nullptr;

    std::size_t med = begin + (end - begin) / 2;
    auto i = _knodes.begin();
    std::nth_element(i + begin, i + med, i + end, knode_cmp<coordinate>(index));
    _knodes[med]._axis = index;
    index = (index + 1) % DIM;

    if(rank != 0){
        MPI_Recv( , dimension*sizerank(rank, n, depth), MPI_FLOAT, father(rank), 0, MPI_COMM_WORLD, &status);
    }

    if(child(rank, depth) < nprocs && sizerank(child(rank, depth), end, depth + 1) > 2){
        MPI_Isend()
    }
}
*/




template<typename coordinate, std::size_t dimension>
knode<coordinate> * kdtree<coordinate, dimension>::make_tree_parallel_string(std::size_t begin, std::size_t end, 
                                                std::size_t index, int nprocs, int depth, MPI_Comm comm, int next){
    
    int rank;
    MPI_Comm_rank(comm, &rank);
    MPI_Status status;

    if(end <= begin) return nullptr;

    std::size_t med = begin + (end - begin) / 2;
    auto i = _knodes.begin();
    std::nth_element(i + begin, i + med, i + end, knode_cmp<coordinate>(index));
    _knodes[med]._axis = index;
    index = (index + 1) % DIM;

    // vector<float>  

    if(rank != 0){
        if(nprocs / 2 != pow(2, depth)){
            depth = depth + 1;
            _knodes[med]._left = make_tree_parallel_string(begin, med, index, nprocs, depth, comm, next);
            next = next + 2;
            _knodes[med]._right = make_tree_parallel_string(med + 1, end, index, nprocs, depth, comm, next);

        }else{
            if(rank == next){
                _knodes[med]._left = make_tree(begin, med, index);
                std::string kdtree_str = serialize_node(_knodes[med]._left);            
                MPI_Send(kdtree_str.c_str(), kdtree_str.length(), MPI_CHAR, 0, 10, comm);
                // deleteTree<coordinate, dimension>(_knodes[med]._left);
            }

            next = next < nprocs - 1 ? next + 1 : 1;

            if(rank == next){
                _knodes[med]._right = make_tree(med + 1, end, index);
                std::string kdtree_str = serialize_node(_knodes[med]._right);
                MPI_Send(kdtree_str.c_str(), kdtree_str.length(), MPI_CHAR, 0, 20, comm);
                // deleteTree<coordinate, dimension>(_knodes[med]._right);
                
            }
        }
    }else if(rank == 0){
        if(nprocs / 2 != pow(2, depth)){
            depth = depth + 1;
            _knodes[med]._left = make_tree_parallel_string(begin, med, index, nprocs, depth, comm, next);
            next = next + 2;
            _knodes[med]._right = make_tree_parallel_string(med + 1, end, index, nprocs, depth, comm, next);

        }else{
            int count;         
            MPI_Probe(next, 10, comm, &status);
            // CHAR
            
            MPI_Get_count(&status, MPI_CHAR, &count);
            char * buf1 = new char[count];
            MPI_Recv(buf1, count, MPI_CHAR, next, 10, comm, &status);
            std::string bla1(buf1, count);
            delete[] buf1;

            _knodes[med]._left = deserialize_node_parallel<coordinate, dimension>(bla1);

            next = next < nprocs - 1 ? next + 1 : 1;

            MPI_Probe(next, 20, comm, &status);
            // CHAR
            
            MPI_Get_count(&status, MPI_CHAR, &count);
            char * buf2 = new char[count];
            MPI_Recv(buf2, count, MPI_CHAR, next, 20, comm, &status);
            std::string bla2(buf2, count);
            delete[] buf2;
            _knodes[med]._right = deserialize_node_parallel<coordinate, dimension>(bla2);

        }
    }

    return &_knodes[med];
}

template<typename coordinate, std::size_t dimension>
knode<coordinate> * kdtree<coordinate, dimension>::make_tree_parallel_pointer(std::size_t begin, std::size_t end, 
                                                std::size_t index, int nprocs, int depth, MPI_Comm comm, int next){
    
    int rank;
    MPI_Comm_rank(comm, &rank);
    MPI_Status status;
    MPI_Request request;

    if(end <= begin) return nullptr;

    std::size_t med = begin + (end - begin) / 2;
    auto i = _knodes.begin();
    std::nth_element(i + begin, i + med, i + end, knode_cmp<coordinate>(index));
    _knodes[med]._axis = index;
    index = (index + 1) % DIM;

    // vector<float>  

    if(rank != 0){
        if(nprocs / 2 != pow(2, depth)){
            depth = depth + 1;
            _knodes[med]._left = make_tree_parallel_pointer(begin, med, index, nprocs, depth, comm, next);
            next = next + 2;
            _knodes[med]._right = make_tree_parallel_pointer(med + 1, end, index, nprocs, depth, comm, next);

        }else{
            if(rank == next){
                _knodes[med]._left = make_tree(begin, med, index); 
                std::vector<coordinate> data;
                float * ser = serialize_pointer<coordinate, dimension>(_knodes[med]._left, data);
                MPI_Send(ser, data.size(), MPI_FLOAT, 0, 10, comm);
            }

            next = next < nprocs - 1 ? next + 1 : 1;

            if(rank == next){
                _knodes[med]._right = make_tree(med + 1, end, index);
                std::vector<coordinate> data;
                float * ser = serialize_pointer<coordinate, dimension>(_knodes[med]._right, data);
                MPI_Send(ser, data.size(), MPI_FLOAT, 0, 20, comm);
                
            }
        }
    }else if(rank == 0){
        if(nprocs / 2 != pow(2, depth)){
            depth = depth + 1;
            _knodes[med]._left = make_tree_parallel_pointer(begin, med, index, nprocs, depth, comm, next);
            next = next + 2;
            _knodes[med]._right = make_tree_parallel_pointer(med + 1, end, index, nprocs, depth, comm, next);

        }else{
            int count;         
            MPI_Probe(next, 10, comm, &status);
            
            MPI_Get_count(&status, MPI_FLOAT, &count);
            coordinate * buf1 = new coordinate[count];
            MPI_Recv(buf1, count, MPI_FLOAT, next, 10, comm, &status);

            _knodes[med]._left = deserialize_pointer<coordinate, dimension>(0, count, 0, 1, buf1);
            delete[] buf1;

            next = next < nprocs - 1 ? next + 1 : 1;

            MPI_Probe(next, 20, comm, &status);
            
            MPI_Get_count(&status, MPI_FLOAT, &count);
            coordinate * buf2 = new coordinate[count];
            MPI_Recv(buf2, count, MPI_FLOAT, next, 20, comm, &status);

            _knodes[med]._right = deserialize_pointer<coordinate, dimension>(0, count, 0, 1, buf2);
            delete[] buf2;

        }
    }

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
                    data.push_back(knode<coordinate>(point<coordinate, dimension>{{stof(X), 
                                                                                    stof(Y)}}));
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

template<typename T>
void print_points(std::vector<knode<T>> node){
    for(int i = 0; i < node.size(); i++){
        std::cout << "(" << node[i]._point.x() << ", " << node[i]._point.y() << ")";
    }
    std::cout << std::endl;
}
