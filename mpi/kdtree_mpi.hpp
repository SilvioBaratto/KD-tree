#include "utility.hpp"

template<typename coordinate, std::size_t dimension>
class kdtree{

    public:

        kdtree(std::string filename, int size, int rank):
        _knodes{getValues(filename)} {
            double mpi_time;
            mpi_time = MPI_Wtime(); 
            _root = make_tree_parallel(0, _knodes.size(), 0, size, 0, MPI_COMM_WORLD, 1);
            //_root = make_tree(0, _knodes.size(), 0);
            mpi_time = MPI_Wtime() - mpi_time;
            if(rank == 0){
                std::cout << size << ", " << mpi_time << std::endl;              
                #ifdef DEBUG
                    knode<coordinate> * root = get_root();
                    #ifdef double_data
                        point<coordinate, dimension> n = nearest({{0.361, 0.674}});
                    #endif
                    #ifdef int_data
                        point<coordinate, dimension> n = nearest({{9, 2}});
                    #endif
                    std::cout << filename << "\n";
                    std::cout << "nearest point: " << n << '\n';
                    std::cout << "distance: " << distance() << '\n';
                    std::cout << "nodes visited: " << visited() << '\n';
                    #ifdef int_data
                    if(filename == ("../datasets/integer/benchmark.csv")){
                        print_tree(root);
                    }
                    #endif
                    #ifdef double_data
                        if(filename == ("../datasets/float/benchmark.csv")){
                            print_tree(root);
                        }
                    #endif
                #endif
            }
        }

        knode<coordinate> * make_tree(std::size_t begin, std::size_t end, std::size_t index);
        knode<coordinate> * get_root() {return _root;}
        knode<coordinate> * make_tree_parallel(std::size_t begin, std::size_t end, 
                                                std::size_t index, int size, int depth, MPI_Comm comm, int which);
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

template<typename coordinate>
std::string serialize_node(knode<coordinate> * node){

    std::string s ("");

        if( node -> _axis != -1) {

            s = node -> _point.save_kpoints( node -> _axis);

            if ( node -> _left != NULL){
                s += "(" + serialize_node(node -> _left) + ")";
            }

            if ( node -> _right != NULL){  
                s += "(" + serialize_node(node -> _right) + ")";
            }
        }

    return s;
}

template<typename coordinate, std::size_t dimension>
knode<coordinate> * deserialize_node(const std::string node){

    const std::size_t size = node.size();
    int j = 0;
    int i;
    coordinate arr[dimension];
    std::string temp_str;
    #ifdef int_data
        int temp_val;
    #endif

    #ifdef double_data
        float temp_val;
    #endif

    if ( size == 0 ){
        return nullptr; 
    }
    
    if ( node[0] == ')'){ 
        return nullptr; 
    }
    
    while ( j < size && node[j] != '(' ) 
        j ++;
    
    if ( node[1] == '[' ){     
        temp_str = node.substr(2, j - 5);
    }else{
        temp_str = node.substr(1, j - 5);
    }

    std::istringstream iss(temp_str);
    for(i = 0; i < dimension; ++i) {
        iss >> temp_val;
        arr[i] = temp_val;
        if (iss.peek() == ',')
            iss.ignore();
    }
    
    point<coordinate, dimension> point{{arr[0], arr[1]}};

    knode<coordinate> * root = new knode<coordinate>{point};

    temp_val = (int)node[j - 2] - 48; //the numbers starts from 48 in the ASCII code
    if(temp_val > 1 || temp_val < 0) temp_val = 0;

    root -> _axis = temp_val; //where axis is saved

    int left = 0;
    i = j;

    while ( i < size ) {
        if ( node[i] == '(' ){ 
            left ++;
        }

        if ( node[i] == ')'){
            left --;
        }

        if ( left == 0 ) {
            break;
        }
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
knode<coordinate> * kdtree<coordinate, dimension>::make_tree(std::size_t begin, std::size_t end, std::size_t index){
    if(end <= begin) return nullptr;

    std::size_t med = begin + (end - begin) / 2;
    auto i = _knodes.begin();
    std::nth_element(i + begin, i + med, i + end, knode_cmp<coordinate>(index));
    _knodes[med]._axis = index;
    index = (index + 1) % DIM;

    #pragma omp task
    _knodes[med]._left = make_tree(begin, med, index);
    #pragma omp task
    _knodes[med]._right = make_tree(med + 1, end, index);

    return &_knodes[med];
}


template<typename coordinate, std::size_t dimension>
knode<coordinate> * kdtree<coordinate, dimension>::make_tree_parallel(std::size_t begin, std::size_t end, 
                                                std::size_t index, int size, int depth, MPI_Comm comm, int which){
    
    int rank;
    MPI_Comm_rank(comm, &rank);
    MPI_Status status;

    if(end <= begin) return nullptr;

    std::size_t med = begin + (end - begin) / 2;
    auto i = _knodes.begin();
    std::nth_element(i + begin, i + med, i + end, knode_cmp<coordinate>(index));
    _knodes[med]._axis = index;
    index = (index + 1) % DIM;

    if(rank != 0){
        if(size / 2 != pow(2, depth)){
            depth = depth + 1;
            _knodes[med]._left = make_tree_parallel(begin, med, index, size, depth, comm, which);
            which = which + 2;
            _knodes[med]._right = make_tree_parallel(med + 1, end, index, size, depth, comm, which);

        }else{
            if(rank == which){
                #pragma omp parallel
                #pragma omp single 
                _knodes[med]._left = make_tree(begin, med, index);
                std::string kdtree_str = serialize_node(_knodes[med]._left);            
                //MPI_Send(kdtree_str.c_str(), kdtree_str.length(), MPI_CHAR, 0, 10, comm);
                MPI_Send(kdtree_str.c_str(), kdtree_str.length(), MPI_CHAR, 0, 10, comm);
            }
            if(which < size - 1) which = which + 1;

            else which = 1;

            if(rank == which){
                #pragma omp parallel
                #pragma omp single
                _knodes[med]._right = make_tree(med + 1, end, index);
                std::string kdtree_str = serialize_node(_knodes[med]._right);
                // MPI_Send( kdtree_str.c_str() , kdtree_str.length() , MPI_CHAR , 0 , 20 , comm);
                MPI_Send(kdtree_str.c_str(), kdtree_str.length(), MPI_CHAR, 0, 20, comm);
            }
        }
    }else if(rank == 0){
        if(size / 2 != pow(2, depth)){
            depth = depth + 1;
            _knodes[med]._left = make_tree_parallel(begin, med, index, size, depth, comm, which);
            which = which + 2;
            _knodes[med]._right = make_tree_parallel(med + 1, end, index, size, depth, comm, which);

        }else{
            int flag = 0, count;
            
            MPI_Probe(which, 10, comm, &status);
            MPI_Get_count(&status, MPI_CHAR, &count);
            char * buf1 = new char[count];
            MPI_Recv(buf1, count, MPI_CHAR, which, 10, comm, &status);
            std::string bla1(buf1, count);
            delete[] buf1;

            _knodes[med]._left = deserialize_node_parallel<coordinate, dimension>(bla1);

            if(which < size - 1)
                which = which + 1;
            else
                which = 1;

            flag = 0;

            MPI_Probe(which, 20, comm, &status);
            MPI_Get_count(&status, MPI_CHAR, &count);

            char * buf2 = new char[count];
            MPI_Recv(buf2, count, MPI_CHAR, which, 20, comm, &status);
            std::string bla2(buf2, count);
            delete[] buf2;

            _knodes[med]._right = deserialize_node_parallel<coordinate, dimension>(bla2);

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
