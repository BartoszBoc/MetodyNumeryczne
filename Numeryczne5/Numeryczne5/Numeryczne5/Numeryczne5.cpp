#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include "eigen-3.4.0/eigen-3.4.0/Eigen/Dense"

using namespace std;

#define SIZE 200

struct Matrix {
    vector<double> lowerdiag;
    vector<double> lowdiag;
    vector<double> diag;
    vector<double> highdiag;
    vector<double> higherdiag;

    Matrix(double d) {
        lowerdiag.resize(SIZE - 2, 0.1);
        lowdiag.resize(SIZE - 1, 0.5);
        diag.resize(SIZE, d);
        highdiag.resize(SIZE - 1, 0.5);
        higherdiag.resize(SIZE - 2, 0.1);
    }
};

double normdiff(vector<double> x, vector<double> x_prev) {
    double sum = 0.0;
    for (int i = 0; i < SIZE; i++) {
        sum += (x[i] - x_prev[i]) * (x[i] - x_prev[i]);
    }
    return sqrt(sum);
}

double norm(vector<double> x) {
    double sum = 0.0;
    for (int i = 0; i < SIZE; i++) {
        sum += x[i] * x[i];
    }
    return sqrt(sum);
}

double normE(Eigen::VectorXd x) {
    double sum = 0.0;
    for (int i = 0; i < SIZE; i++) {
        sum += x(i) * x(i);
    }
    return sqrt(sum);
}

vector<double> GaussSeidel(Matrix matrix, vector<double> b, vector<double> x0, bool makeCompGraph, bool makeIteGraph) {
    vector<double> x(SIZE, 0.0);
    vector<double> x_prev = x0;

    int k = 0;
    vector<double> normdifftable, normtable;
    normdifftable.resize(40000, 0.0);
    normtable.resize(40000, 0.0);
 
    while (1) {
        if (matrix.diag[0] < (matrix.lowerdiag[0] + matrix.lowdiag[0] + matrix.highdiag[0] + matrix.higherdiag[0])) {
            if(k == 0) cout << "Algorytm Gaussa-Seidela dla tej macierzy (o d = " << matrix.diag[0] << ") jest rozbiezny. Wywolujemy 20 iteracji.\n";
            if(k == 20) break;
        }

        for (int i = 0; i < SIZE; i++) {
            double sum = b[i];

            if (i > 0) sum -= matrix.lowdiag[i - 1] * x[i - 1];
            if (i < SIZE - 1) sum -= matrix.highdiag[i] * x_prev[i + 1];
            if (i > 1) sum -= matrix.lowerdiag[i - 2] * x[i - 2];
            if (i < SIZE - 2) sum -= matrix.higherdiag[i] * x_prev[i + 2];

            x[i] = sum / matrix.diag[i];
        }

        if (normdiff(x, x_prev) < 1e-7) {
            cout << "Zbieznosc po iteracji dla Gaussa-Seidela dla d = (" << matrix.diag[0] << "): " << k << "\n";
            break;
        }
        
        if (makeCompGraph) normdifftable[k] = normdiff(x, x_prev);
        if (makeIteGraph) normtable[k] = norm(x);
        x_prev = x;
        k++;
    }

    if (makeCompGraph) {
        int i = 0;
        ofstream FileMaker("GSFileTable");
        while (normdifftable[i] != 0.0) {
            FileMaker << scientific << normdifftable[i] << " " << i + 1 << "\n";
            i++;
        }
        FileMaker.close();
    }

    if (makeIteGraph) {
        int i = 0;
        ofstream FileMaker("GSCompTable");
        while (normtable[i] != 0.0) {
            FileMaker << scientific << normtable[i] << " " << i + 1 << "\n";
            i++;
        }
        FileMaker.close();
    }

    return x;
}

vector<double> Jacobi(Matrix matrix, vector<double> b, vector<double> x0, bool makeCompGraph, bool makeIteGraph) {
    vector<double> x(SIZE, 0.0);
    vector<double> x_prev = x0;

    int k = 0;
    vector<double> normdifftable, normtable;
    normdifftable.resize(40000, 0.0);
    normtable.resize(40000, 0.0);

    while (1) {
        if (matrix.diag[0] < (matrix.lowerdiag[0] + matrix.lowdiag[0] + matrix.highdiag[0] + matrix.higherdiag[0])) {
            if(k == 0) cout << "Algorytm Jacobiego dla tej macierzy (o d = " << matrix.diag[0] << ") jest rozbiezny. Wywolujemy 20 iteracji.\n";
            if(k == 20) break;
        }

        for (int i = 0; i < SIZE; i++) {
            double sum = b[i];

            if (i > 0) sum -= matrix.lowdiag[i - 1] * x_prev[i - 1];
            if (i < SIZE - 1) sum -= matrix.highdiag[i] * x_prev[i + 1];
            if (i > 1) sum -= matrix.lowerdiag[i - 2] * x_prev[i - 2];
            if (i < SIZE - 2) sum -= matrix.higherdiag[i] * x_prev[i + 2];

            x[i] = sum / matrix.diag[i];
        }

        double vectornorm = normdiff(x, x_prev);

        if (vectornorm < 1e-7) {
            cout << "Zbieznosc po iteracji dla Jacobiego dla d = (" << matrix.diag[0] << "): " << k << "\n";
            break;
        }
        
        if (makeCompGraph) normdifftable[k] = normdiff(x, x_prev);
        if (makeIteGraph) normtable[k] = norm(x);
        x_prev = x;
        k++;
    }

    if (makeCompGraph) {
        int i = 0;
        ofstream FileMaker("JFileTable");
        while (normdifftable[i] != 0.0) {
            FileMaker << scientific << normdifftable[i] << " " << i + 1 << "\n";
            i++;
        }
        FileMaker.close();
    }

    if (makeIteGraph) {
        int i = 0;
        ofstream FileMaker("JCompTable");
        while (normtable[i] != 0.0) {
            FileMaker << scientific << normtable[i] << " " << i + 1 << "\n";
            i++;
        }
        FileMaker.close();
    }

    return x;
}

void printVector(vector<double> a) {
    cout << "Wynik: \n\n(";
    for (int i = 0; i < SIZE; i++) {
        cout << a[i];
        if (i != 0 && i % 10 == 9 && i != SIZE - 1) cout << ",\n";
        else if (i == SIZE - 1) continue;
        else cout << ", ";
    }
    cout << ")\n\n\n";
}

void printVectorEigen(Eigen::VectorXd a) {
    cout << "(";
    for (int i = 0; i < SIZE; i++) {
        cout << a(i);
        if (i != 0 && i % 10 == 9 && i != SIZE - 1) cout << ",\n";
        else if (i == SIZE - 1) continue;
        else cout << ", ";
    }
    cout << ")\n\n\n";
}

void program(Eigen::MatrixXd matrixe, Eigen::VectorXd bE, Matrix matrix, vector<double> b, bool makeCompGraph) {
    vector<double> x01;
    x01.resize(SIZE, 0);

    vector<double> x02;
    x02.resize(SIZE, 1);

    vector<double> x03 = b;

    vector<double> x04;
    x04.resize(SIZE);
    srand(time(NULL));
    for (int i = 0; i < SIZE; i++) {
        x04[i] = (double)(rand() % 10000) / 1000.0;
    }

    vector<double> Jx1 = Jacobi(matrix, b, x01, makeCompGraph, false);
    cout << "Punkt startowy dla powyzszego algorytmu: wektor z samymi 0.\n";
    printVector(Jx1);

    vector<double> GSx1 = GaussSeidel(matrix, b, x01, makeCompGraph, false);
    cout << "Punkt startowy dla powyzszego algorytmu: wektor z samymi 0.\n";
    printVector(GSx1);

    vector<double> Jx2 = Jacobi(matrix, b, x02, false, false);
    cout << "Punkt startowy dla powyzszego algorytmu: wektor z samymi 1.\n";
    printVector(Jx2);

    vector<double> GSx2 = GaussSeidel(matrix, b, x02, false, false);
    cout << "Punkt startowy dla powyzszego algorytmu: wektor z samymi 1.\n";
    printVector(GSx2);

    vector<double> Jx3 = Jacobi(matrix, b, x03, false, false);
    cout << "Punkt startowy dla powyzszego algorytmu: wektor wynikowy.\n";
    printVector(Jx3);

    vector<double> GSx3 = GaussSeidel(matrix, b, x03, false, false);
    cout << "Punkt startowy dla powyzszego algorytmu: wektor wynikowy.\n";
    printVector(GSx3);

    vector<double> Jx4 = Jacobi(matrix, b, x04, false, false);
    cout << "Punkt startowy dla powyzszego algorytmu: wektor z losowymi wartosciami.\n";
    printVector(Jx4);

    vector<double> GSx4 = GaussSeidel(matrix, b, x04, false, false);
    cout << "Punkt startowy dla powyzszego algorytmu: wektor z losowymi wartosciami.\n";
    printVector(GSx4);

    Eigen::VectorXd result(SIZE);
    result = matrixe.colPivHouseholderQr().solve(bE);
    cout << "Wynik za pomoca biblioteki Eigen:\n\n";
    printVectorEigen(result);

    ofstream FileMaker("Result");
    double resultnorm = normE(result);
    for (int i = 0; i < 10; i++) {
        FileMaker << scientific << resultnorm << " " << i + 1 << "\n";
    }
    FileMaker.close();
}

int main(){
    vector<double> b;
    Eigen::VectorXd bE(SIZE);
    b.resize(SIZE);
    for (int i = 0; i < SIZE; i++) {
        b[i] = i + 1;
        bE(i) = i + 1;
    }

    Matrix macierz1(0.5);
    Matrix macierz2(1.201);
    Matrix macierz3(20.0);

    Eigen::MatrixXd macierz1e(SIZE, SIZE);
    Eigen::MatrixXd macierz2e(SIZE, SIZE);
    Eigen::MatrixXd macierz3e(SIZE, SIZE);

    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            if (i == j) {
                macierz1e(i, j) = 0.5;
                macierz2e(i, j) = 1.201;
                macierz3e(i, j) = 20.0;
            }
            else if (i + 1 == j || i - 1 == j) {
                macierz1e(i, j) = 0.5;
                macierz2e(i, j) = 0.5;
                macierz3e(i, j) = 0.5;
            }
            else if (i + 2 == j || i - 2 == j) {
                macierz1e(i, j) = 0.1;
                macierz2e(i, j) = 0.1;
                macierz3e(i, j) = 0.1;
            }
            else {
                macierz1e(i, j) = 0.0;
                macierz2e(i, j) = 0.0;
                macierz3e(i, j) = 0.0;
            }
        }
    }

    program(macierz1e, bE, macierz1, b, false);
    program(macierz2e, bE, macierz2, b, false);
    program(macierz3e, bE, macierz3, b, true);
}