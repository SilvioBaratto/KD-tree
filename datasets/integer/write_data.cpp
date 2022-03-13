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
    std::cout << "Write the name of the dataset: \n";
    std::string data;
    std::cin >> data;
    data = data + ".csv";
    myfile.open (data);
    for(int i = 0; i < size; i++){
        int x = (rand() % 1000) + 1;
        int y = (rand() % 1000) + 1;

        myfile << x << "," << y << "\n";   
    }
    myfile.close();
    return 0;  
}
