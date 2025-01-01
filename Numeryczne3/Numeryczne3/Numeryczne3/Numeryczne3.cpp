#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>

#define SIZE 300

using namespace std;

struct FourDiagonalMatrix {			// Oznaczenia L i U pomocne dla macierzy LU która zostanie zwrócona przez funkcje LUDecomposition
	vector<double> U_diag;			// (S¹ to z³¹czone dwie macierze L i U, gdzie diagonala i elementy nad diagonal¹ s¹ macierz¹ U...
	vector<double> U_highdiag;		// ... a elementy pod diagonal¹ s¹ macierz¹ L, gdzie diagonala zawsze bêdzie wype³niona jedynkami)
	vector<double> U_higherdiag;
	vector<double> L_lowdiag;
};

FourDiagonalMatrix LUDecomposition(FourDiagonalMatrix Matrix, int size) {

	FourDiagonalMatrix LU;

	LU.L_lowdiag.resize(size - 1);
	LU.U_diag.resize(size);
	LU.U_highdiag.resize(size - 1);
	LU.U_higherdiag.resize(size - 2);

	LU.U_diag = Matrix.U_diag;
	LU.U_highdiag = Matrix.U_highdiag;
	LU.U_higherdiag = Matrix.U_higherdiag;

	for (int i = 1; i < size; ++i) {
		
		LU.L_lowdiag[i - 1] = Matrix.L_lowdiag[i - 1] / LU.U_diag[i - 1];	// Obliczenie diagonali pod g³ówn¹ diagonal¹ dla macierzy L (g³ówna diagonala bêdzie zawieraæ tylko 1-ki)

		LU.U_diag[i] -= LU.L_lowdiag[i - 1] * LU.U_highdiag[i - 1];	// Obliczanie g³ównej przek¹tnej w macierzy U

		if (i < (size - 1)) LU.U_highdiag[i] -= LU.L_lowdiag[i - 1] * Matrix.U_higherdiag[i - 1];	// Obliczanie przek¹tnej nad g³ówn¹ w macierzy U (a¿ do jej koñca)

		if (i < (size - 2)) LU.U_higherdiag[i] = Matrix.U_higherdiag[i]; // Obliczanie przek¹tnej nad przek¹tn¹ nad g³ówn¹ w macierzy U
	}

	return LU;
}

vector<double> Forward_Substitution(FourDiagonalMatrix LU, vector<double> x) {

	int vectorsize = x.size();
	vector<double> z(vectorsize);
	z[0] = x[0];							// Pierwszy element jest równy dla obu wektorów, gdy¿ przek¹tna L ma same jedynki
	
	for (int i = 1; i < vectorsize; ++i) z[i] = x[i] - LU.L_lowdiag[i - 1] * z[i - 1];

	return z;
}

vector<double> Backward_Substitution(FourDiagonalMatrix LU, vector<double> z) {

	int vectorsize = z.size();
	vector<double> y(vectorsize);
	y[vectorsize - 1] = z[vectorsize - 1] / LU.U_diag[vectorsize - 1];		// Analogicznie, lecz od koñca (backward_substitution)

	for (int i = vectorsize - 2; i >= 0; i--) y[i] = (z[i] - LU.U_highdiag[i] * y[i + 1] - ((i < vectorsize - 2) ? (LU.U_higherdiag[i]) * y[i + 2] : 0)) / LU.U_diag[i];

	return y;
}

double determinant(FourDiagonalMatrix LU) {
	double determinant = 1.0;
	for (double U : LU.U_diag) determinant *= U;	// Wyznacznik macierzy A jest równy wyznacznikom macierzy L i U (det(L) = 1, wiêc tylko wyznacznikowi macierzy U)
	return determinant;								// Wyznacznik macierzy U jest równy iloczynowi elementów na diagonali
}	

void program(int size, bool shouldPrint) {
	FourDiagonalMatrix Macierz;

	Macierz.U_diag.resize(size);
	Macierz.U_highdiag.resize(size - 1);
	Macierz.U_higherdiag.resize(size - 2);
	Macierz.L_lowdiag.resize(size - 1);

	for (int i = 0; i < size; ++i) {
		Macierz.U_diag[i] = 1.01;
		if (i < size - 1) Macierz.U_highdiag[i] = 0.2 / (i + 1);
		if (i < size - 2) Macierz.U_higherdiag[i] = 0.15 / ((i + 1) * (i + 1) * (i + 1));
		if (i > 0) Macierz.L_lowdiag[i - 1] = 0.3;
	}

	vector<double> x(size);

	for (int i = 0; i < size; ++i) {
		x[i] = i + 1;
	}

	FourDiagonalMatrix LU = LUDecomposition(Macierz, size);

	vector<double> z(size);
	z = Forward_Substitution(LU, x);
	vector<double> y(size);
	y = Backward_Substitution(LU, z);
	double det = determinant(LU);
	
	if (shouldPrint) {
		cout << "Rozwiazanie to:\n(";
		for (int i = 0; i < size; ++i) {
			cout << y[i];
			if (i != size - 1) cout << "; ";
			if (i % 10 == 9 && i != 0 && i != size - 1) cout << "\n";
		}
		cout << ")\n\n" << "Wyznacznik macierzy A to: " << det << "\n";
	}
}

int main()
{
	program(SIZE, true);

	int i = 100;
	vector<pair<int, double>> time;

	while(1) {
		using namespace chrono;

		auto start = high_resolution_clock::now();
		program(i, false);
		auto finish = high_resolution_clock::now();

		auto duration = duration_cast<microseconds>(finish - start);
		time.push_back({i, (double)(duration.count()/1000.0)});

		if (duration > microseconds(100000) || i > 1000000) break;

		i += 100;
	}

	ofstream FileMaker("FileTable");
	for (auto value : time) {
		FileMaker << value.first << " " << value.second << "\n";
	}
	FileMaker.close();

	system("gnuplot DisplayGraph.txt");
}