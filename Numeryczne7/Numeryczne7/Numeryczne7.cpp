#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

double y(double x) {
    return 1 / (1 + 10 * x * x);
}

double interpolacjaLagrange(double x, vector<double> x_i, vector<double> y_i, int n) {
    double f = 0.0;

    for (int i = 0; i < n; i++) {
        double L = 1.0;
        for (int j = 0; j < n; j++) {
            if (i != j) L *= (x - x_i[j]) / (x_i[i] - x_i[j]);
        }

        f += y_i[i] * L;
    }

    return f;
}

vector<vector<double>> szesciennySpline(vector<double> x_i, vector<double> y_i, int n) {
    double h = x_i[1] - x_i[0];

    vector<double> righthandside(n, 0.0);
    for (int i = 1; i < n - 1; i++) {
        righthandside[i] = 6 * ((y_i[i + 1] - y_i[i]) / h - (y_i[i] - y_i[i - 1]) / h) / h;
    }

    vector<double> diagonal(n, 4.0), subdiagonals(n - 1, 1.0);
    for (int i = 1; i < n - 1; i++) {
        double gaussElimination = subdiagonals[i - 1] / diagonal[i - 1];
        diagonal[i] -= gaussElimination * subdiagonals[i - 1];
        righthandside[i] -= gaussElimination * righthandside[i - 1];
    }

    vector<double> xi_i(n, 0.0);
    for (int i = n - 2; i > 0; i--) {
        xi_i[i] = (righthandside[i] - subdiagonals[i] * xi_i[i + 1]) / diagonal[i];
    }

    vector<vector<double>> abcd(4, vector<double>(n - 1));
    for (int i = 0; i < n - 1; i++) {
        abcd[0][i] = (xi_i[i + 1] - xi_i[i]) / (6 * h);
        abcd[1][i] = xi_i[i] / 2.0;
        abcd[2][i] = (y_i[i + 1] - y_i[i]) / h - h * (2 * xi_i[i] + xi_i[i + 1]) / 6;
        abcd[3][i] = y_i[i];
    }

    return abcd;
}

double wartoscSpline(vector<vector<double>> abcd, vector<double> x_i, int n, double x) {
    int i = 0;
    while (i < n - 1 && x > x_i[i + 1]) i++;

    if (x < x_i[0] || x > x_i[n - 1]) {
        cerr << "Error: x out of domain: " << x << "\n";
        return 0.0;
    }

    double xprime = x - x_i[i];
    return abcd[0][i] * xprime * xprime * xprime + abcd[1][i] * xprime * xprime + abcd[2][i] * xprime + abcd[3][i];
}

void program(double (*func)(double), string name, int n) {
    vector<double> x_i(n), y_i(n);
    for (int i = 0; i < n; i++) {
        x_i[i] = -1.0 + (2.0 * i) / (n - 1);
        y_i[i] = func(x_i[i]);
    }

    ofstream PointGraph("LagrangePunkty - " + name);
    cout << "Lagrange\n\n";
    for (double x : x_i) {
        double f = interpolacjaLagrange(x, x_i, y_i, n);
        cout << "f(" << x << ") = " << f << "\n";
        PointGraph << scientific << f << " " << x << "\n";
    }
    PointGraph.close();

    vector<double> differenceLagrange(2001);
    int i = 0;
    ofstream LagrangeGraph("LagrangeWielomian - " + name);
    for (double x = -1.0; x <= 1.0; x += 0.001) {
        double f = interpolacjaLagrange(x, x_i, y_i, n);
        LagrangeGraph << scientific << f << " " << x << "\n";
        differenceLagrange[i] = abs(f - func(x));
        i++;
    }
    LagrangeGraph.close();

    vector<vector<double>> abcd = szesciennySpline(x_i, y_i, n);

    vector<double> differenceSpline(2001);
    i = 0;
    ofstream SplineGraph("SplineWykres - " + name);
    for (double x = -1.0; x <= 1.0; x += 0.001) {
        double f = wartoscSpline(abcd, x_i, n, x);
        SplineGraph << scientific << f << " " << x << "\n";
        differenceSpline[i] = abs(f - func(x));
        i++;
    }
    SplineGraph.close();
    i = 0;

    ofstream Delta("Difference - " + name);
    for (double x = -1.0; x <= 1.0; x += 0.001) {
        Delta << scientific << differenceLagrange[i] << " " << differenceSpline[i] << " " << x << "\n";
        i++;
    }
    Delta.close();
}

int main(){
    program(y, "y", 10);
}