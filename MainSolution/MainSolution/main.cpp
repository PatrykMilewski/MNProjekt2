#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>

using namespace std;


void createMatrix(double **tab, int tabSize, int third, int fourth) {
	double a1 = 3;
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

void printVector(double *tab, int size) {
	for (int i = 0; i < size; i++)
		cout << tab[i] << endl;
}

void gaussJordan(double **matrix, double *results, int size) {
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

double calculateResiduum(double **matrix, double *vector, int size) {
	double *tab = new double[size];
	double sum = 0;

	for (int i = 0; i < size; i++)
		tab[i] = 0;

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++)
			tab[i] += matrix[i][j] * vector[j];
		tab[i] -= matrix[i][size];
		tab[i] *= tab[i];
	}

	for (int i = 0; i < size; i++)
		sum += tab[i];

	return sqrt(sum);
}


int jacobi(double **matrix, double **inverseMatrix, double *diagonal, double **results, int size, double epsilon) {
	// calculate 1/D and M = -1/D * (L + U)
	for (int i = 0; i < size; i++) {
		diagonal[i] = 1 / matrix[i][i];		// 1/D

		for (int j = 0; j < size; j++) {
			if (i == j)
				inverseMatrix[i][j] = 0;
			else
				inverseMatrix[i][j] = (-1) * (matrix[i][j] * diagonal[i]);
		}
	}

	// calculate results
	bool resultInEpsilon;
	int iterationsAmount = 0;
	while (true) {
		resultInEpsilon = true;
		for (int i = 0; i < size; i++) {
			results[i][2] = results[i][1];					// x3 = x2
			results[i][1] = diagonal[i] * matrix[i][size];	// x2 = 1/D * Bi

			for (int j = 0; j < size; j++) {
				results[i][1] += inverseMatrix[i][j] * results[j][0];	// x2 = 1/M * x1
			}
		}
		for (int i = 0; i < size; i++) {
			results[i][0] = results[i][1];					// x1 = x2

		if (fabs(results[i][0] - results[i][2]) > epsilon)	// check if results outside of epsilon
				resultInEpsilon = false;
		}
		iterationsAmount++;
		if (resultInEpsilon)
			break;
		if (iterationsAmount > 1000)						// max iterations break
			break;
	}
	return iterationsAmount;
}

int gaussSeidel(double **matrix, double **inverseMatrix, double *diagonal, double **results, int size, double epsilon) {
	// calculate 1/D ; 1/D * b ; 1/D * L ; 1/D * U
	for (int i = 0; i < size; i++) {
		diagonal[i] = 1 / matrix[i][i];			// 1/D
		matrix[i][size] *= diagonal[i];			// 1/D * b

		for (int j = 0; j < size; j++) {
			if (i != j)
				inverseMatrix[i][j] = matrix[i][j] * diagonal[i];	// 1/D * L && 1/D * U
			else
				inverseMatrix[i][j] = matrix[i][j];					// A
		}
	}

	bool resultInEpsilon;
	int iterationsAmount = 0;

	while (true) {
		resultInEpsilon = true;
		for (int i = 0; i < size; i++) {
			results[i][0] = matrix[i][size];	// x1 = 1/D * b

			for (int j = 0; j < i; j++)
				results[i][0] -= inverseMatrix[i][j] * results[j][0];	// 1/D * L * x
			for (int j = i + 1; j < size; j++)
				results[i][0] -= inverseMatrix[i][j] * results[j][0];	// 1/D * U * x

			
		}
		for (int i = 0; i < size; i++) {
			if (fabs(results[i][0] - results[i][1]) > epsilon)	// check if results outside of epsilon
				resultInEpsilon = false;

			results[i][1] = results[i][0];						// x2 = x1
		}
		iterationsAmount++;
		if (resultInEpsilon)
			break;
		if (iterationsAmount > 1000)							// max iterations break
			break;
	}
	return iterationsAmount;
}

void copyExtendedMatrix(double **matrix, double **destiny, int size) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size + 1; j++)
			destiny[i][j] = matrix[i][j];
	}
}

// function made for checking if all results are in given epsilon
int checkResults(double *ideal, double **check, int size, double epsilon) {
	int counter = 0;
	for (int i = 0; i < size; i++) {
		if (fabs(ideal[i] - check[i][0]) > epsilon)
			counter++;
	}
	return counter;
}

int main() {
	// variables and tables initializations
	const int testsAmount = 11;
	int tabTime[testsAmount * 3];
	int tabIterations[testsAmount * 3];
	double residuum;
	int size = 903;
	clock_t Start, Stop;
	double epsilon = 1e-9;

	for (int i = 0; i < testsAmount; i++) {

		double **matrixGaussSeidel = new double*[size];
		double **matrixJacobi = new double*[size];
		double **matrixGaussJordan = new double*[size];
		double **matrix = new double*[size];
		double **resultsGaussSeidel = new double*[size];
		double **resultsJacobi = new double*[size];
		double *resultsGaussJordan = new double[size];
		double *vector = new double[size];
		for (int i = 0; i < size; i++) {
			matrixGaussSeidel[i] = new double[size + 1];
			matrixJacobi[i] = new double[size + 1];
			matrixGaussJordan[i] = new double[size + 1];
			matrix[i] = new double[size];
			resultsJacobi[i] = new double[3];
			resultsGaussSeidel[i] = new double[2];
		}

		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size + 1; j++) {
				matrixGaussSeidel[i][j] = 0;
				matrix[i][j] = 0;
			}
			vector[i] = 0;
			resultsJacobi[i][0] = 0;
			resultsJacobi[i][1] = 0;
			resultsJacobi[i][2] = 0;

			resultsGaussSeidel[i][0] = 0;
			resultsGaussSeidel[i][1] = 0;
		}


		createMatrix(matrixGaussSeidel, size, 0, 5);
		copyExtendedMatrix(matrixGaussSeidel, matrixJacobi, size);
		copyExtendedMatrix(matrixGaussSeidel, matrixGaussJordan, size);

		Start = clock();
		tabIterations[i] = jacobi(matrixJacobi, matrix, vector, resultsJacobi, size, epsilon);
		Stop = clock() - Start;
		tabTime[i] = (int)Stop;


		Start = clock();
		tabIterations[i + testsAmount] = gaussSeidel(matrixGaussSeidel, matrix, vector, resultsGaussSeidel, size, epsilon);
		Stop = clock() - Start;
		tabTime[i + testsAmount] = (int)Stop;

		Start = clock();
		gaussJordan(matrixGaussJordan, resultsGaussJordan, size);
		tabIterations[i + 2 * testsAmount] = 1;
		Stop = clock() - Start;
		residuum = calculateResiduum(matrixGaussJordan, resultsGaussJordan, size);
		tabTime[i + 2 * testsAmount] = (int)Stop;

		for (int i = 0; i < size; i++) {
			delete[] matrixGaussSeidel[i];
			delete[] matrixJacobi[i];
			delete[] matrixGaussJordan[i];
		}
		delete[] matrix;
		delete[] matrixGaussSeidel;
		delete[] matrixJacobi;
		delete[] matrixGaussJordan;
		delete[] resultsGaussSeidel;
		delete[] resultsJacobi;
		delete[] resultsGaussJordan;
		delete[] vector;
	}
	int count = 0;
	fstream output("results.txt", fstream::out);
	for (int i = 0; i < testsAmount * 3; i++) {
		if (i % testsAmount != 0) {										// delete first result from every test
			output << tabTime[i] << "\t" << tabIterations[i] << endl;
			count++;
		}
	}
	output << "Residuum: " << residuum << endl;
	output.flush();
	output.close();

	return 0;
}