#include <iostream>
#include "eigen-3.4.0/eigen-3.4.0/Eigen/Dense"
#include <vector>
#include <chrono>
#include <fstream>

using namespace std;

#define SIZE 120

vector<double> ShermanMorrison(vector<double> diag, vector<double> highdiag, vector<double> b, vector<double> u, size_t N) {
    vector<double> result(N);
    vector<double> bSolution(N);
    vector<double> uSolution(N);
    vector<double> newhighdiag(N - 1);

    vector<double> diag_copy = diag;

    for (int i = 0; i < N - 1; i++) {
        newhighdiag[i] = highdiag[i] / diag_copy[i];
        diag_copy[i + 1] -= newhighdiag[i] * highdiag[i];
    }

    bSolution[0] = b[0] / diag[0];
    uSolution[0] = u[0] / diag[0];

    for (int i = 1; i < N; i++) {
        bSolution[i] = (b[i] - highdiag[i - 1] * bSolution[i - 1]) / diag[i];
        uSolution[i] = (u[i] - highdiag[i - 1] * uSolution[i - 1]) / diag[i];
    }

    double sum_bs = 0.0, sum_us = 0.0;
    for (int i = 0; i < N; i++) {
        sum_bs += bSolution[i];
        sum_us += uSolution[i];
    }

    double factor = sum_bs / (1.0 + sum_us);

    for (int i = 0; i < N; i++) {
        result[N - 1 - i] = (bSolution[i] - uSolution[i] * factor);
    }

    return result;
}

void program(size_t N, bool shouldPrint) {
    vector<double> diag(N, 4.0);
    vector<double> highdiag(N - 1, 2.0);
    vector<double> u(N, 1.0);
    vector<double> b1(N, 2.0);

    vector<double> result2(N);
    result2 = ShermanMorrison(diag, highdiag, b1, u, N);

    if (shouldPrint) {
        cout << "\n\nWynik za pomoca wzoru Shermanna-Morrisona: \n\n";

        cout << "(";
        for (int i = 0; i < N; i++) {
            cout << result2[i];
            if (i != 0 && i % 5 == 4 && i != N - 1) cout << ",\n";
            else if (i == N - 1) continue;
            else cout << ", ";
        }
        cout << ")\n";
    }
}

void eigenprogram(size_t N) {
    Eigen::MatrixXd matrixA(N, N);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) matrixA(i, j) = 5.0;
            else if (i + 1 == j) matrixA(i, j) = 3.0;
            else matrixA(i, j) = 1.0;
        }
    }

    Eigen::VectorXd b(N);

    for (int i = 0; i < N; i++) {
        b(i) = 2.0;
    }

    Eigen::VectorXd result1(N);
    result1 = matrixA.colPivHouseholderQr().solve(b);

    cout << "Wynik za pomoca biblioteki Eigen: \n\n";

    cout << "(";
    for (int i = 0; i < N; i++) {
        cout << result1[i];
        if (i != 0 && i % 5 == 4 && i != N - 1) cout << ",\n";
        else if (i == N - 1) continue;
        else cout << ", ";
    }
    cout << ")";
}

int main()
{
    eigenprogram(SIZE);
    program(SIZE, true);

    int i = 100;
    vector<pair<int, double>> time;
    while (1) {
        using namespace chrono;

        auto start = high_resolution_clock::now();
        program(i, false);
        auto finish = high_resolution_clock::now();

        auto duration = duration_cast<microseconds>(finish - start);
        time.push_back({ i, (double)(duration.count() / 1000.0) });

        if (duration > microseconds(50000) || i > 300000) break;

        i += 100;
    }
    ofstream FileMaker("FileTable");

    for (auto value : time) {
        FileMaker << value.first << " " << value.second << "\n";
    }
}