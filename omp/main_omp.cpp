#include "kdtree_omp.hpp"

int main(int argc, char * argv[]) {

    // the main could take two tipe of knode data structure: float and integer
    // both take as input a file from a folder choose by the user in the terminal
    // if the user does not pass any file from default the functions takes the default one

    #ifdef int_data
        const std::string filename =
        argc > 1 ? argv[1] : "../datasets/integer/benchmark.csv";
        kdtree<int, 2> tree(filename);
        #ifdef int_data_DEBUG
        knode<int> * root = tree.get_root();
        point<int, 2> n = tree.nearest({{9, 2}});
        std::cout << filename << "\n";
        std::cout << "nearest point: " << n << '\n';
        std::cout << "distance: " << tree.distance() << '\n';
        std::cout << "nodes visited: " << tree.visited() << '\n';
        // just for debug porpuse is not reccomend to print if the size is bigger than 1000
        if(filename == ("../datasets/integer/benchmark.csv")){
            print_tree(root);
        }
        #endif
    #endif

    #ifdef double_data
        const std::string filename =
        argc > 1 ? argv[1] : "../datasets/float/benchmark.csv";
        kdtree<double, 2> tree(filename);
        #ifdef double_data_DEBUG
            knode<double> * root = tree.get_root();
            point<double, 2> n = tree.nearest({{0.361, 0.674}});
            std::cout << filename << "\n";
            std::cout << "nearest point: " << n << '\n';
            std::cout << "distance: " << tree.distance() << '\n';
            std::cout << "nodes visited: " << tree.visited() << '\n';
            // just for debug porpuse is not reccomend to print if the size is bigger than 1000
            if(filename == ("../datasets/float/benchmark.csv")){
                print_tree(root);
            }
        #endif
    #endif

    return 0;

}