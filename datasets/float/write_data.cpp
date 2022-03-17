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

#define SIZE 100000000

int main()
{
    std::ofstream myfile;
    myfile.open ("dataset.csv");
    int i;
    for(i = 0; i < SIZE; i++){
        std::random_device rd;
        std::default_random_engine eng(rd());
        std::uniform_real_distribution<> distr(0, 1);

        myfile << std::setprecision(3) << distr(eng) << "," << distr(eng) << "\n";   
    }
    myfile.close();
    return 0;  
}
