#include "kdtree_mpi.hpp"

int main() {
    std::string filename = "dataset.csv";
    auto start = std::chrono::high_resolution_clock::now();
    kdtree<int> tree(filename);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Time to making the tree: " << duration.count() << std::endl;
    // tree.printRoot();
}