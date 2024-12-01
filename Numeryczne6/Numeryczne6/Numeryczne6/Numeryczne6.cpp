#include <iostream>
#include "eigen-3.4.0/eigen-3.4.0/Eigen/Dense"
#include <cmath>
#include <vector>
#include <fstream>
#include <string>

using namespace std;
using namespace Eigen;

void PrintVector(VectorXd v) {
    cout << "(";
    for (int i = 0; i < 4; i++) {
        cout << v(i);
        if (i < 3) cout << ", ";
    }
    cout << ")\n";
}

void PrintMatrix(MatrixXd a) {
    cout << "(";
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            cout << a(i, j);
            if (j < 3) cout << ", ";
        }
        if(i < 3) cout << "\n";
    }
    cout << ")\n";
}

void GraphResults(string filename, vector<double> table) {
    ofstream FileMaker(filename);
    int i = 0;
    while (table[i] != 0.0) {
        FileMaker << scientific << table[i] << " " << i + 1 << "\n";
        i++;
    }
    FileMaker.close();
}

void SendMatrices(string filename, vector<MatrixXd> Ai) {
    ofstream FileMaker(filename);
    for (MatrixXd A : Ai) {
        for (int j = 0; j < 4; j++){
            for (int k = 0; k < 4; k++) {
                FileMaker << A(j, k) << " ";
            }
            FileMaker << "\n";
        }
        FileMaker << "\n\n";
    }
    FileMaker.close();
}

pair<double, VectorXd> PowerMethod(MatrixXd M, double tolerance) {
    double eigenvalue = 0.0, prev_eigenvalue = 0.0;
    int i = 0;
    vector<double> eigenvaluetable;
    eigenvaluetable.resize(20, 0.0);
    VectorXd v(4);
    v << 1, 1, 1, 1;
    v.normalize();

    EigenSolver<MatrixXd> solver(M);
    VectorXd eigenvalues = solver.eigenvalues().real();
    double exacteigenvalue = eigenvalues[0];

    while (1) {
        VectorXd v_next = M * v;
        v_next.normalize();
        prev_eigenvalue = eigenvalue;

        eigenvalue = v_next.dot(M * v_next);
        v = v_next;
        eigenvaluetable[i] = abs(exacteigenvalue - eigenvalue);

        if (abs(eigenvalue - prev_eigenvalue) < tolerance) break;
        i++;
    }

    cout << "Zbieznosc algorytmu metody potegowej osiagnieta po " << i << " iteracjach.\n\n";

    GraphResults("PowerMethodEigenvalues", eigenvaluetable);

    return { eigenvalue, v };
}

MatrixXd QRAlgorithm(MatrixXd M, double tolerance) {
    MatrixXd lambda(100, 4);

    EigenSolver<MatrixXd> solver(M);
    VectorXd eigenvalues = solver.eigenvalues().real();

    vector<MatrixXd> Ai;
    Ai.resize(100);

    int i = 0;

    while (1) {
        Ai[i] = M;

        for (int j = 0; j < 4; j++) {
            lambda(i, j) = abs(M(j, j) - eigenvalues[j]);
            if (lambda(i, j) < tolerance) lambda(i, j) = 0;
        }

        HouseholderQR<MatrixXd> qr(M);
        MatrixXd Q = qr.householderQ();
        MatrixXd R = qr.matrixQR().triangularView<Upper>();
        M = R * Q;

        MatrixXd Mdiag = M.diagonal().asDiagonal();
        MatrixXd Mnodiag = M - Mdiag;

        if (lambda(i, 0) < tolerance && lambda(i, 1) < tolerance && lambda(i, 2) < tolerance && lambda(i, 3) < tolerance) break;
        i++;
    }

    vector<double> lambda1(100, 0), lambda2(100, 0), lambda3(100, 0), lambda4(100, 0);

    for (int j = 0; j <= i; j++) {
        lambda1[j] = lambda(j, 0);
        lambda2[j] = lambda(j, 1);
        lambda3[j] = lambda(j, 2);
        lambda4[j] = lambda(j, 3);
    }

    GraphResults("QREigenvalues1", lambda1);
    GraphResults("QREigenvalues2", lambda2);
    GraphResults("QREigenvalues3", lambda3);
    GraphResults("QREigenvalues4", lambda4);

    Ai.resize(i + 1);
    SendMatrices("Matrices Ai.txt", Ai);

    cout << "Zbieznosc algorytmu QR osiagnieta po " << i + 1 << " iteracjach.\n\n";
    return M;
}

int main()
{
    MatrixXd M(4, 4);
    M << 9, 2, 0, 0,
         2, 4, 1, 0,
         0, 1, 3, 1,
         0, 0, 1, 2;

    double tolerance = 1e-12;

    pair<double, VectorXd> eigen = PowerMethod(M, tolerance);
    double eigenvalue = eigen.first;
    VectorXd eigenvector = eigen.second;

    cout << "Najwieksza wartosc wlasna (metoda potegowa): " << eigenvalue << "\n\n";
    cout << "Odpowiadajacy jej wektor wlasny (metoda potegowa): \n";
    PrintVector(eigenvector);
    cout << "\n\n";

    MatrixXd A_max(4, 4);
    A_max = QRAlgorithm(M, tolerance);
    cout << "Finalna macierz Ai:\n";
    PrintMatrix(A_max);

    EigenSolver<MatrixXd> solver(M);
    VectorXd eigenvalues = solver.eigenvalues().real();
    MatrixXd eigenvectors = solver.eigenvectors().real();

    cout << "\nWartosci wlasne obliczone za pomoca biblioteki Eigen:\n";
    PrintVector(eigenvalues);
    cout << "\nZnormalizowane wektory wlasne obliczone za pomoca biblioteki Eigen zlozone w macierz:\n";
    PrintMatrix(eigenvectors);
}