#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <numbers>
#include <random>
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

double normalDistribution(double sigma) { 
    random_device randomdevice;
    mt19937 generator(randomdevice());
    normal_distribution<> normaldistribution(0.0, sigma);
    return normaldistribution(generator);
}

void vectorWriteout(VectorXd v) {
    cout << "(";
    for (int i = 0; i < v.size(); i++) {
        cout << v(i);
        if (i < v.size() - 1) cout << ", ";
        else cout << ")\n";
    }
}

void drawFunction(string filename, VectorXd coefficients) {
    ofstream funct(filename);

    for (double j = -1.0; j <= 1.0; j += 0.001) {
        double result = exercisefunction(j, coefficients);
        funct << scientific << result << " " << j << "\n";
    }

    funct.close();
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
        yset(i) += normalDistribution(sigma);
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

    drawFunction("precise", acoefficient);
    drawFunction("approx", pcoefficient);

    ofstream pointPlotter("points");
    for (int i = 0; i < n; i++) pointPlotter << scientific << yset[i] << " " << xset[i] << "\n";
    pointPlotter.close();
}

int main(){
    VectorXd acoefficient(M); 
    acoefficient << -0.05, 1.5, 3.0, -2.0, 1.0, -5.0;
    program(400, acoefficient, 3.0);
}