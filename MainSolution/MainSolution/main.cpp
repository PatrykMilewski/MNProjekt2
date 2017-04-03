#include <iostream>
#include <math.h>
#include <time.h>

using namespace std;


void createMatrix(double **tab, int tabSize, int third, int fourth) {
	double a1 = 5 + fourth;
	double a2 = -1, a3 = -1;
	double tabPom[5] = { a3, a2, a1, a2, a3 };

	int delta = -2;
	for (int i = 0; i < tabSize; i++) {
		for (int j = 0; j < 5; j++) {
			if (j + delta < tabSize && j + delta >= 0)
				tab[i][j + delta] = tabPom[j];
		}

		delta++;
		tab[i][tabSize] = sin((double)i * ((double)third + 1) / 50);
	}
}

void printMatrix(double **tab, int tabSize) {
	for (int i = 0; i < tabSize; i++) {
		for (int j = 0; j < tabSize + 1; j++)
			cout << tab[i][j] << " ";

		cout << endl;
	}

}

void gaussSeidl(double **matrix, double *results, int size) {
	double temp;

	for (int i = 1; i < size; i++) {
		for (int k = i; k < size; k++) {
			if (matrix[i - 1][i - 1] != 0)
				temp = matrix[k][i - 1] / matrix[i - 1][i - 1];
			else
				continue;

			for (int j = i - 1; j < size + 1; j++)
					matrix[k][j] = matrix[k][j] - (temp * matrix[i - 1][j]);

		}
	}

	for (int i = size - 1; i >= 0; i--) {
		temp = matrix[i][size];
		for (int j = size - 1; j >= i + 1; j--)
			temp -= matrix[i][j] * results[j];

		results[i] = temp / matrix[i][i];
	}

}

void jacobi(double **matrix, double **inverseMatrix, double *vector, double **results, int size, double epsilon) {
	double temp;

	// calculate diagonal^(-1)
	for (int i = 0; i < size; i++) {
		if (matrix[i][i] != 0)
			vector[i] = 1 / matrix[i][i];
		else
			throw new exception;
	}


	// calculate M = -1/D * (L + U)
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (i == j)
				inverseMatrix[i][j] = 0;
			else
				inverseMatrix[i][j] = (-1) * (matrix[i][j] * vector[i]);
		}
	}

	// calculate results
	bool resultInEpsilon;
	while (true) {
		resultInEpsilon = true;
		for (int i = 0; i < size; i++) {
			results[i][2] = results[i][1];					// x3 = x2
			results[i][1] = vector[i] * matrix[i][size];	// x2 = 1/D * Bi

			for (int j = 0; j < size; j++) {
				results[i][1] += inverseMatrix[i][j] * results[i][0];	// x2 = 1/M * x1
			}
		}
		for (int i = 0; i < size; i++) {
			results[i][0] = results[i][1];		// x1 = x2

			if (fabs(results[i][0] - results[i][2]) > epsilon)	// calculate if is not inside: -epsilon < delta < epsilon
				resultInEpsilon = false;
		}
		if (resultInEpsilon)
			break;
	}

}

void copyExtendedMatrix(double **matrix, double **destiny, int size) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size + 1; j++)
			destiny[i][j] = matrix[i][j];
	}
}

int main() {
	int size = 5;

	double **matrixGaussSeidl = new double*[size];
	double **matrixJacobi = new double*[size];
	double **matrixGauss = new double*[size];
	double **matrix = new double*[size];
	double *resultsGaussSeidl = new double[size];
	double **resultsJacobi = new double*[size];
	double *resultsGauss = new double[size];
	double *vector = new double[size];
	for (int i = 0; i < size; i++) {
		matrixGaussSeidl[i] = new double[size + 1];
		matrixJacobi[i] = new double[size + 1];
		matrixGauss[i] = new double[size + 1];
		matrix[i] = new double[size];
		resultsJacobi[i] = new double[2];
	}

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size + 1; j++) {
			matrixGaussSeidl[i][j] = 0;
			matrix[i][j] = 0;
		}
		vector[i] = 0;
		resultsJacobi[i][0] = 0;
		resultsJacobi[i][1] = 0;
		resultsJacobi[i][2] = 0;
	}


	createMatrix(matrixGaussSeidl, size, 0, 5);
	copyExtendedMatrix(matrixGaussSeidl, matrixJacobi, size);
	copyExtendedMatrix(matrixGaussSeidl, matrixGauss, size);


	//printMatrix(matrixJacobi, size);
	clock_t Start;

	Start = clock();
	jacobi(matrixJacobi, matrix, vector, resultsJacobi, size, 1e-12);
	cout << "Time Difference: " << clock() - Start << endl;
	cout << endl;
	printMatrix(matrix, size);

	cout << endl;

	for (int i = 0; i < size; i++)
		cout << resultsJacobi[i][0] << endl;

	//printMatrix(matrix, size);
	//cout << endl;
	Start = clock();
	gaussSeidl(matrixGaussSeidl, resultsGaussSeidl, size);
	//printMatrix(matrix, size);
	cout << "Time Difference: " << clock() - Start << endl;
	cout << endl;

	for (int i = 0; i < size; i++)
		cout << resultsGaussSeidl[i] << endl;

	for (int i = 0; i < size; i++) {
		delete[] matrixGaussSeidl[i];
		delete[] matrixJacobi[i];
		delete[] matrixGauss[i];
	}
	delete[] matrixGaussSeidl;
	delete[] matrixJacobi;
	delete[] matrixGauss;
	delete[] resultsGaussSeidl;
	delete[] resultsJacobi;
	delete[] resultsGauss;
	return 0;
}