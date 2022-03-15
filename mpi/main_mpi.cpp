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

    int nprocs, rank;
    MPI_Init(&argc, &argv);

    MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    #ifdef int_data
        kdtree<int, 2> tree(filename, nprocs, rank);
        if(rank == 0){
            #ifdef int_data_DEBUG
            knode<int> * root = tree.get_root();
            point<int, 2> n = tree.nearest({{9, 2}});
            std::cout << filename << "\n";
            // this are utility and debug functions
            std::cout << "nearest point: " << n << '\n';
            std::cout << "distance: " << tree.distance() << '\n';
            std::cout << "nodes visited: " << tree.visited() << '\n';
            // just for debug porpuse is not reccomend to print if the size is bigger than 1000
            if(filename == ("../datasets/integer/benchmark.csv")){
                print_tree(root);
            }
            #endif
        }
    #endif

    #ifdef double_data
        kdtree<float, 2> tree(filename, nprocs, rank);
        if(rank == 0){
            #ifdef double_data_DEBUG
                knode<float> * root = tree.get_root();
                point<float, 2> n = tree.nearest({{0.361, 0.674}});
                std::cout << filename << "\n";
                // this are utility and debug functions
                std::cout << "nearest point: " << n << '\n';
                std::cout << "distance: " << tree.distance() << '\n';
                std::cout << "nodes visited: " << tree.visited() << '\n';
                // just for debug porpuse is not reccomend to print if the size is bigger than 1000
                if(filename == ("../datasets/float/benchmark.csv")){
                    print_tree(root);
                }
            #endif
        }
    #endif  
        
    MPI_Finalize();

    return 0;

}