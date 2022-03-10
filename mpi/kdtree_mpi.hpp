#include "utility.hpp"

template<typename coordinate, std::size_t dimension>
class kdtree{

    public:

        kdtree(std::string filename, int size):
        _knodes{getValues(filename)} {
            // auto start = std::chrono::high_resolution_clock::now();
            //double mpi_time;
            //mpi_time = MPI_Wtime();
            _root = make_tree_parallel(0, _knodes.size(), 0, size, 0, MPI_COMM_WORLD, 1);
            //_root = make_tree(0, _knodes.size(), 0);
            //mpi_time = MPI_Wtime() - mpi_time;
            //std::cout << mpi_time << std::endl;
            // auto stop = std::chrono::high_resolution_clock::now();
            // _root = make_tree(0, _knodes.size(), 0);
            // auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
            // std::cout << "Time to making the tree: " << duration.count() << std::endl;
        }

        knode<coordinate> * make_tree(std::size_t begin, std::size_t end, std::size_t index);
        knode<coordinate> * get_root() {return _root;}
        knode<coordinate> * make_tree_parallel(std::size_t begin, std::size_t end, 
                                                std::size_t index, int size, int level, MPI_Comm comm, int which);

        std::vector<knode<coordinate>> getValues(std::string filename);
        void printData();
        void printRootMemory();
        void nearest(knode<coordinate> * root, const point<coordinate, dimension>& point, size_t index);
        const point<coordinate, dimension>& nearest(const point<coordinate, dimension>& pt);

        std::size_t get_size() {return _knodes.size();}

        std::vector<knode<coordinate>> get_left_knode(std::size_t begin, std::size_t end);

        bool empty() const {return _knodes.empty();}
        std::size_t visited() const {return _visited;}
        double distance() const {return std::sqrt(_best_dist);}

        void print_knodes(){
            for(int i = 0; i < _knodes.size(); i++){
                std::cout << _knodes[i]._point;
            }
            std::cout << std::endl;
        }

    private:
        knode<coordinate> * _root = nullptr;
        std::vector<knode<coordinate>> _knodes;
        knode<coordinate> * _best = nullptr;
        std::size_t _visited;
        double _best_dist = 0;
};

template<typename T>
std::string serialize_node(knode<T> * node){

    std::string s ("");

    try{

        if( node -> _axis != -1) {
                
                s = node -> _point.save_kpoints( node -> _axis);

                if ( node -> _left != NULL)
                s += "(" + serialize_node(node -> _left) + ")";
                
                if ( node -> _right != NULL)
                s += "(" + serialize_node(node -> _right) + ")";

            }
    }catch(const std::exception& e){
        std::cout<<e.what();
    }

    return s;
}

template<typename T, std::size_t dimension>
knode<T> * deserialize_node(std::string data){

    std::string s = data;

    

    if ( s.size() == 0 ){
        // std::cout << "size 0" << std::endl;
        return nullptr; 
    }
    
    if ( s[0] == ')'){ 
        // std::cout << "return nullptr" << std::endl;
        return nullptr; 
    }

    int j = 0;

    while ( j < s.size() && s[j] != '(' ) 
        j ++;

    T arr[dimension];
    std::string temp_str;
    if ( s[1] == '[' ){     
        std::cout << "s[1] == '['" << std::endl; 
        temp_str = s.substr(2, j - 5);

    }else{
        // T arr1[N_DIM]{std::stoi(s.substr(1, j-5))}; //only the number of the array
        // for(int i =0; i<N_DIM; ++i)
        //     arr[i] = arr1[i];
        temp_str = s.substr(1, j - 5);

    }

    // std::cout << "temp str: " << temp_str << std::endl;

    int temp_val;
    std::istringstream iss(temp_str);
    for(int i = 0; i < dimension; ++i) {
        iss >> temp_val;
        arr[i] = temp_val;
        if (iss.peek() == ',')
            iss.ignore();
    }
    
    

    point<T, dimension> point{{arr[0], arr[1]}};

    // std::cout << point.x() << ", " << point.y() << ", ";

    // struct kpoint<T> point(arr);
    knode<T> * root = new knode<T>{point};
    // struct kdnode<T> * root = new kdnode<T>;

    temp_val = (int)s[j - 2] - 48;//the numbers starts from 48 in the ASCII code
    if(temp_val > 1) temp_val = 0;
    // std::cout << temp_val << std::endl;
    root -> _axis = temp_val; //where axis is saved

    int left = 0, i = j;

    // find separation between left and right definition
    while ( i < s.size() ) {
        if ( s[i] == '(' ) 
            left ++;
        else if ( s[i] == ')' ) 
            left --;

        if ( left == 0 ) {
            break;
        }
        i ++;
    }

    if ( j < s.size() - 1) {
        root-> _left = deserialize_node<T, dimension>(s.substr(j + 1, i - 1 - j));
    }
    if ( i + 1 < s.size() - 1) {
        root-> _right = deserialize_node<T, dimension>(s.substr(i + 2, s.size() - i - 2));   
    }
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

    // #pragma omp task
    _knodes[med]._left = make_tree(begin, med, index);
    // #pragma omp task
    _knodes[med]._right = make_tree(med + 1, end, index);

    return &_knodes[med];
}


template<typename coordinate, std::size_t dimension>
knode<coordinate> * kdtree<coordinate, dimension>::make_tree_parallel(std::size_t begin, std::size_t end, 
                                                std::size_t index, int size, int level, MPI_Comm comm, int which){
    
    int rank;
    double mpi_time;
    MPI_Comm_rank(comm, &rank);
    MPI_Status status;
    MPI_Request request;

    #ifdef DEBUG
        std::cout << "size: " << size << " rank: " << rank << std::endl;
    #endif

    if(end <= begin) return nullptr;

    std::size_t med = begin + (end - begin) / 2;
    auto i = _knodes.begin();
    std::nth_element(i + begin, i + med, i + end, knode_cmp<coordinate>(index));
    _knodes[med]._axis = index;
    index = (index + 1) % DIM;

    if(rank != 0){
        std::cout << "rank[" << rank << "]" << " which: " << which << std::endl;
        #ifdef DEBUG
            std::cout << "rank != 0" << std::endl;
        #endif
        if(size / 2 != pow(2, level)){
            std::cout << "Make tree parallel" << std::endl;
            level = level + 1;
            _knodes[med]._left = make_tree_parallel(begin, med, index, size, level, comm, which);
            which = which + 2;
            _knodes[med]._right = make_tree_parallel(med + 1, end, index, size, level, comm, which);
        }else{
            if(rank == which){
                std::cout << "left part go serial" << std::endl;
            #ifdef DEBUG
                std::cout << "rank[" << rank << "]" << " == which left" << std::endl;
            #endif
                _knodes[med]._left = make_tree(begin, med, index);
                // mpi_time = MPI_Wtime();
                std::string kdtree_str = serialize_node(_knodes[med]._left);
                // mpi_time = MPI_Wtime() - mpi_time;
                // std::cout << "rank: " << rank << " Time to serialize: " << mpi_time << std::endl;              
                MPI_Send(kdtree_str.c_str(), kdtree_str.length(), MPI_CHAR, 0, 10, comm);
            }
            if(which < size - 1) which = which + 1;

            else which = 1;

            if(rank == which){
                std::cout << "right part go serial" << std::endl;
            #ifdef DEBUG
                std::cout << "rank[" << rank << "]" << " == which right" << std::endl;
            #endif
                _knodes[med]._right = make_tree(med + 1, end, index);
                // mpi_time = MPI_Wtime();
                std::string kdtree_str = serialize_node(_knodes[med]._right);
                // mpi_time = MPI_Wtime() - mpi_time;
                // std::cout << "rank: " << rank << " Time to serialize: " << mpi_time << std::endl;
                MPI_Send( kdtree_str.c_str() , kdtree_str.length() , MPI_CHAR , 0 , 20 , comm );
            }
            #ifdef DEBUG
                std::cout << "\n Processor n: " << rank << " AFTER SEND \n";
            #endif
        }
    }else if(rank == 0){
        std::cout << "merge the knode in the master process" << std::endl;
        #ifdef DEBUG
            std::cout << "rank == 0" << std::endl;
        #endif
        if(size / 2 != pow(2, level)){
            std::cout << "make tree parallel 2" << std::endl;
            level = level + 1;
            _knodes[med]._left = make_tree_parallel(begin, med, index, size, level, comm, which);
            which = which + 2;
            _knodes[med]._right = make_tree_parallel(med + 1, end, index, size, level, comm, which);
        }else{
            int flag = 0, count;
            MPI_Probe(which, 10, comm, &status);
            MPI_Get_count(&status, MPI_CHAR, &count);
            char * buf1 = new char[count];
            MPI_Recv(buf1, count, MPI_CHAR, which, 10, comm, &status);
            std::string bla1(buf1, count);
            delete[] buf1;
            
            
            // mpi_time = MPI_Wtime();
            _knodes[med]._left = deserialize_node<coordinate, dimension>(bla1);
            // mpi_time = MPI_Wtime() - mpi_time;
            // std::cout << "rank: " << rank << " Time to deserialize left: " << mpi_time << std::endl;
            if(which < size - 1)
                which = which + 1;
            else
                which = 1;

            flag = 0;


            MPI_Probe(which, 20, comm, &status);
            MPI_Get_count(&status, MPI_CHAR, &count);
            std::cout << count << std::endl;
            char * buf2 = new char[count];
            MPI_Recv(buf2, count, MPI_CHAR, which, 20, comm, &status);
            std::string bla2(buf2, count);
            delete[] buf2;

            // mpi_time = MPI_Wtime();
            _knodes[med]._right = deserialize_node<coordinate, dimension>(bla2);
            // mpi_time = MPI_Wtime() - mpi_time;
            // std::cout << "rank: " << rank << " Time to deserialize right: " << mpi_time << std::endl;
        }
    }
    return &_knodes[med];
}

template<typename coordinate, std::size_t dimension>
std::vector<knode<coordinate>> kdtree<coordinate, dimension>:: get_left_knode(std::size_t begin, std::size_t end){
    std::size_t med = begin + (end - begin) / 2;
    std::vector<knode<coordinate>> left_points(_knodes.begin(), _knodes.begin() + med);
    return left_points;
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
void kdtree<coordinate, dimension>::printData(){
    int i = 0;
    for(auto x: _knodes){
        std::cout << x._left << ", ";
    }
    std::cout << std::endl;
}

template<typename coordinate, std::size_t dimension>
void kdtree<coordinate, dimension>::printRootMemory(){
    std::cout << std::endl;
    for(int i = 0; i <_knodes.size(); i++){
        std::cout << "(" << &_root[i]._left << ", " << &_root[i]._right << ")" << std::endl;
    }
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
                data.push_back(knode<coordinate>(point<coordinate, dimension>{{stoi(X), stoi(Y)}}));
            }
            i++;
        }
        coeff.close();
    } else std::cout << "Unable to open file";

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
