#include <iostream>
#include <cmath>
#include <vector>
#include <numbers>
#include "eigen-3.4.0/eigen-3.4.0/Eigen/Dense"

#define M 6

using namespace std;
using namespace Eigen;

double phi(int i, double x) {
    switch (i) {
    case 1:
        return exp(4.0 * x);
    case 2:
        return sin(4.0 * x * x);
    case 3:
        return log(abs(log(abs(2.0 * x) + 0.01)) + 1.0);
    case 4:
        return 1 / (abs(8.0 * x * x * x) + 1.0);
    case 5:
        return atan(x);
    case 6:
        return exp(-abs(sin(2.0 * x))) - 0.5;
    default:
        cout << "Error: i out of range.\n";
        return 0;
    }
}

double exercisefunction(double x, VectorXd acoefficient) {
    double result = 0.0;
    for (int i = 0; i < M; i++) result += acoefficient[i] * phi(i + 1, x);
    return result;
}

double normaldistribution(double y, double sigma) { return (1.0 / sqrt(2.0 * numbers::pi * sigma * sigma)) * exp(-(y * y) / (2.0 * sigma * sigma)); }

void vectorWriteout(VectorXd v) {
    cout << "(";
    for (int i = 0; i < v.size(); i++) {
        cout << v(i);
        if (i < v.size() - 1) cout << ", ";
        else cout << ")\n";
    }
}

void program(int n, VectorXd acoefficient, double sigma) {
    if (acoefficient.size() != M) {
        cout << "Error: Amount of coefficients doesn't match with the amount of functions.\n";
        return;
    }

    VectorXd xset(n);
    VectorXd yset(n);

    for (int i = 0; i < n; i++) {
        xset(i) = -1.0 + (2.0 * i) / (n - 1);
        yset(i) = exercisefunction(xset(i), acoefficient);
        yset(i) += normaldistribution(yset(i), sigma);
    }

    MatrixXd Amatrix(n, M);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < M; j++) {
            Amatrix(i, j) = phi(j + 1, xset(i));
        }
    }

    VectorXd pcoefficient(M);
    pcoefficient = Amatrix.jacobiSvd(ComputeThinU | ComputeThinV).solve(yset);

    cout << "Coefficients inserted: \n";
    vectorWriteout(acoefficient);
    cout << "\nCoefficients reached: \n";
    vectorWriteout(pcoefficient);
}

int main(){
    VectorXd acoefficient(M); 
    acoefficient << -0.05, 1.5, 3.0, -2.0, 1.0, -5.0;
    program(15, acoefficient, 0.5);
}