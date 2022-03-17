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

int main(int argc, char** argv)
{
    std::ofstream myfile;
    myfile.open ("dataset.csv");
    for(int i = 0; i < SIZE; i++){
        int x = (rand() % 1000) + 1;
        int y = (rand() % 1000) + 1;

        myfile << x << "," << y << "\n";   
    }
    myfile.close();
    return 0;  
}
