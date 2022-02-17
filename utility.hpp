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

int getDataRows(std::string filename);

std::vector<point<int>> getData(std::string filename);

template<typename T> 
void printKD(const std::string& prefix, const knode<T> * node, bool isLeft);

template<typename T>  
void printKD(const knode<T> * node);


#endif