#include "kdtree_mpi.hpp"

int main(int argc, char * argv[]) {
    std::string filename = "dataset_bigger.csv";

    int size, rank;
    MPI_Init(&argc, &argv);

    MPI_Comm_size( MPI_COMM_WORLD, &size );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    kdtree<int, 2> tree(filename, size, rank);

    #ifdef DEBUG
        knode<int> * root = tree.get_root();
        std::string root_str = serialize_node(root);
        knode<int> * root_des = deserialize_node<int, 2>(root_str);
        print_tree(root_des);
    #endif

    MPI_Finalize();

    return 0;

}