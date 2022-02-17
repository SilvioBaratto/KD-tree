#include "utility.hpp"

// ******************** GET DATA SECTION ***********************
int getDataRows(std::string filename){
    int rows=0;
    std::ifstream file(filename);
    std::string line;

    while (getline(file, line)){ 
        rows++;
    }

    return rows;
}


// knode<T> * getData(std::string filename)
std::vector<point<int>> getData(std::string filename){
    std::string X, Y;

    int size = getDataRows(filename);
    std::vector<point<int>> data; 

    std::ifstream coeff(filename);
    int i = 0;
    if (coeff.is_open()) {
        while(!coeff.eof()) {
            getline(coeff, X, ',');
            getline(coeff, Y, '\n');
            if(!X.empty() && !Y.empty()){
                data.push_back(point<int>{stoi(X), stoi(Y)});
            }
            i++;
        }
        coeff.close();
    } else std::cout << "Unable to open file";
    return data;
}


// ******************** END DATA SECTION ***********************

// ******************** VARIOUS PRINT SECTION ***********************

template<typename T>
void printKD(const std::string& prefix, const knode<T> * node, bool isLeft){
    if( node != nullptr ){
        std::cout << prefix;

        std::cout << (isLeft ? "├──" : "└──" );

    // print the value of the node
    std::cout << "(" << node->x[0] << "," << node->x[1] << ")" << std::endl;

    if (node->left) printKD(prefix + (isLeft ? "│   " : "    "), node->left, true);
    if (node->right) printKD(prefix + (isLeft ? "│   " : "    "), node->right, false);

    }
}

template<typename T> 
void printKD(const knode<T> * node){
    printKD("", node, false);    
}
