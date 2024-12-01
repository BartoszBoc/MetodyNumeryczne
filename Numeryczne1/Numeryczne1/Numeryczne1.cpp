#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <fstream>

using namespace std;

template<typename T>
T func(T x) {
    return sin(x * x * x);
}

template<typename T>
T fderiv(T x) {
    return 3 * x * x * cos(x * x * x);
}

template<typename T>
T fderivA(T x, T h) {
    return (func(x + h) - func(x)) / h;
}

template<typename T>
T fderivB(T x, T h) {
    return (func(x + h) - func(x - h)) / (2 * h);
}

template<typename T>
T error(T x, T h, function<T(T, T)> compderiv) {
    return abs(compderiv(x, h) - fderiv(x));
}

int main()
{
    vector<vector<double>> Table(5);
    int size = 0;

    for (double h = 10E-18; h < 1.01; h *= pow(10, 0.2)) {
        Table[0].push_back(error<double>(0.2, h, fderivA<double>));
        Table[1].push_back(error<double>(0.2, h, fderivB<double>));
        Table[2].push_back(error<float>(0.2, h, fderivA<float>));
        Table[3].push_back(error<float>(0.2, h, fderivB<float>));
        Table[4].push_back(h);
        size++;
    }

    ofstream FileMaker("FileTable");
    
    for (int i = 0; i < size; i++) {
        FileMaker << scientific << Table[0][i];
        FileMaker << " ";
        FileMaker << scientific << Table[1][i];
        FileMaker << " ";
        FileMaker << scientific << Table[2][i];
        FileMaker << " ";
        FileMaker << scientific << Table[3][i];
        FileMaker << " ";
        FileMaker << scientific << Table[4][i];
        FileMaker << "\n";
    }
}