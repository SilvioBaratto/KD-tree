// write_csv.cpp
// updated 8/12/16

// the purpose of this code is to demonstrate how to write data
// to a csv file in C++

// inlcude iostream and string libraries
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <iomanip>


int main(int argc, char** argv)
{
    std::ofstream myfile;
    int size;
    std::cout << "write the size of the dataset: \n";
    std::cin >> size;
    myfile.open ("dataset.csv");
    for(int i = 0; i < size; i++){
        std::random_device rd;
        std::default_random_engine eng(rd());
        std::uniform_real_distribution<> distr(0, 1);

        myfile << std::setprecision(3) << distr(eng) << "," << distr(eng) << "\n";   
    }
    myfile.close();
    return 0;  
}
