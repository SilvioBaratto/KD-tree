#include "kdtree_mpi.hpp"

int main(int argc, char * argv[]) {

    #ifdef int_data
        const std::string filename =
        argc > 1 ? argv[1] : "../datasets/integer/benchmark.csv";
    #endif

    #ifdef double_data
        const std::string filename =
        argc > 1 ? argv[1] : "../datasets/float/benchmark.csv";
    #endif

    int size, rank;
    MPI_Init(&argc, &argv);

    MPI_Comm_size( MPI_COMM_WORLD, &size );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    #ifdef int_data
        kdtree<int, 2> tree(filename, size, rank);
    #endif

    #ifdef double_data
        kdtree<float, 2> tree(filename, size, rank);
    #endif  

    MPI_Finalize();

    return 0;

}