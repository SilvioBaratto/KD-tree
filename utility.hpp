#pragma once
#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED
#include <iostream>
#include <vector>
#include <chrono>
#include <cstring>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <filesystem>
#include <omp.h>
#include "knode.hpp"

int getDataRows(std::string filename){
    int rows=0;
    std::ifstream file(filename);
    std::string line;

    while (getline(file, line)){ 
        rows++;
    }

    return rows;
}

#endif