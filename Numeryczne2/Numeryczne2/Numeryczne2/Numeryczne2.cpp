#include <iostream>
#include "eigen-3.4.0/eigen-3.4.0/Eigen/Dense"
#include <Random>

using namespace std;
using namespace Eigen;

int main()
{
	MatrixXd A_1(5, 5);
	MatrixXd A_2(5, 5);
    VectorXd B(5);

    A_1 << 5.8267103432, 1.0419816676, 0.4517861296, -0.2246976350, 0.7150286064,
        1.0419816676, 5.8150823499, -0.8642832971, 0.6610711416, -0.3874139415,
        0.4517861296, -0.8642832971, 1.5136472691, -0.8512078774, 0.6771688230,
        -0.2246976350, 0.6610711416, -0.8512078774, 5.3014166511, 0.5228116055,
        0.7150286064, -0.3874139415, 0.6771688230, 0.5228116055, 3.5431433879;

    A_2 << 5.4763986379, 1.6846933459, 0.3136661779, -1.0597154562, 0.0083249547,
        1.6846933459, 4.6359087874, -0.6108766748, 2.1930659258, 0.9091647433,
        0.3136661779, -0.6108766748, 1.4591897081, -1.1804364456, 0.3985316185,
        -1.0597154562, 2.1930659258, -1.1804364456, 3.3110327980, -1.1617171573,
        0.0083249547, 0.9091647433, 0.3985316185, -1.1617171573, 2.1174700695;

    B << -2.8634904630, -4.8216733374, -4.2958468309, -0.0877703331, -2.0223464006;

    VectorXd x_1 = A_1.colPivHouseholderQr().solve(B);
    VectorXd x_2 = A_2.colPivHouseholderQr().solve(B);

    cout << "A_1 * x_1 = B\n\nx_1:\n\n" << x_1 << "\n";
    cout << "\n\nA_2 * x_2 = B\n\nx_2:\n\n" << x_2 << "\n";

    for (int i = 1; i <= 3; i++) {

        VectorXd deltaB = VectorXd::Random(5);
        deltaB.normalize();
        deltaB *= 1e-6;

        VectorXd totB = B + deltaB;

        VectorXd x_1 = A_1.colPivHouseholderQr().solve(totB);
        VectorXd x_2 = A_2.colPivHouseholderQr().solve(totB);
        
        cout << "\nITERATION " << i << "\n\n";
        cout << "A_1 * x_1 = B + db\n\nx_1:\n\n" << x_1 << "\n";
        cout << "\n\nA_2 * x_2 = B + db\n\nx_2:\n\n" << x_2 << "\n";
        cout << "\nThe error: \n\n" << deltaB << "\n";
    }
}