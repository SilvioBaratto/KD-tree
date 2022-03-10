#include "kdtree_omp.hpp"

int main(int argc, char * argv[]) {
    std::string filename = "dataset_bigger.csv";

    kdtree<int, 2> tree(filename);

    return 0;

}