/*
 ============================================================================
 Name        : SVD.c
 Author      : Jiahua Liang, Kevin Chen, Xiangyun Wang, Yinuo Wang
 Version     : 1.0
 Copyright   : Your copyright notice
 Description : Matrix Solver in C, Ansi-style
 ============================================================================


 ============================================================================
 This is the original version of the SVD solver.
 Note that since the algorithm of finding eigenvalue in c is not applicable
 due to segmentation faults, the computed eigen matrix and values are declared
 directly in the code.
 ============================================================================
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

int size = 5;

static double eigenvector[5][5] = { { -0.37850527, 0.22754263, -0.18262938,
		-0.87840887, 0.0015737 }, { 0.23120478, 0.49272987, -0.78356889,
		0.19133505, 0.23059125 }, { -0.03987157, -0.24071506, -0.41177417,
		0.03886613, -0.87715334 }, { -0.39238419, 0.74061481, 0.32566007,
		0.29263578, -0.32532163 }, { -0.80481241, -0.31462293, -0.27758576,
		0.3234849, 0.2675688 }, };

static double eigenvalues[5] = { 216.37561647, 100.08410035, 49.38342505,
		34.38342505, 0.35340114 };

static double b[5] = { 1, 2, 5, 3, 1 };

double** dotTranspose(double** A, int size) {
	// create a transpose matrix
	int i = 0;
	double** A_transpose = malloc(size * sizeof(double*));
	for (i = 0; i < size; i++) {
		A_transpose[i] = malloc(size * sizeof(double));
	}
	int j = 0;
	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			A_transpose[i][j] = A[j][i];
		}
	}
	// allocate a new matrix for returning
	double** product = malloc(size * sizeof(double*));
	for (i = 0; i < size; i++) {
		product[i] = malloc(size * sizeof(double));
	}
	// do multiplication
	int k = 0;
	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {

			double m = 0;
			for (k = 0; k < size; k++) {
				m = m + A_transpose[i][k] * A[k][j];
			}

			product[i][j] = m;
		}
	}
	return product;
}

void printMatrix(double** a, int size) {
	int i = 0;
	int j = 0;

	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			printf("%f ", a[i][j]);
			if (j == size - 1) {
				printf("\n");
			}
		}
	}
	printf("\n");
}

double** readMatrix(int size) {
	int i = 0;
	int j = 0;

	double** mat = malloc(size * sizeof(double*));
	for (i = 0; i < size; i++) {
		mat[i] = malloc(size * sizeof(double));
	}

	FILE *file;
	file = fopen("data.txt", "r");

	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			if (!fscanf(file, "%lf", &mat[i][j]))
				break;
		}

	}
	fclose(file);
	return mat;

}

double** transpose(double** a, int size) {
	int i = 0;
	double** mat = malloc(size * sizeof(double*));
	for (i = 0; i < size; i++) {
		mat[i] = malloc(size * sizeof(double));
	}

	int j = 0;
	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			mat[i][j] = a[j][i];
		}
	}

	return mat;

}

double** MdotV(double** a, double** v, int size) {
	int i, j;
	double** product = malloc(size * sizeof(double*));
	for (i = 0; i < size; i++) {
		product[i] = malloc(size * sizeof(double));
	}
	// do multiplication
	int k = 0;
	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {

			double m = 0;
			for (k = 0; k < size; k++) {
				m = m + a[i][k] * v[k][j];
			}

			product[i][j] = m;
		}
	}
	return product;

}

double** normalizeVectors(double** v, int size) {
	int i, j, k;
	double** result = malloc(size * sizeof(double*));
	for (i = 0; i < size; i++) {
		result[i] = malloc(size * sizeof(double));
	}

	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {

			double pythag = 0;
			for (k = 0; k < size; k++) {
				pythag += v[k][j] * v[k][j];
			}

			result[i][j] = v[i][j] / sqrt(pythag);
		}
	}

	return result;
}

double* sqrtVector(double* v, int size) {
	int i;
	double* result = malloc(size * sizeof(double));

	for (i = 0; i < size; i++) {
		result[i] = sqrt(v[i]);
	}
	return result;
}

double** matrixDotDiag(double** a, double* v, int size) {
	int i, j;
	double** result = malloc(size * sizeof(double*));
	for (i = 0; i < size; i++) {
		result[i] = malloc(size * sizeof(double));
	}

	for (j = 0; j < size; j++) {
		double m = v[j];
		for (i = 0; i < size; i++) {

			result[i][j] = a[i][j] * m;
		}
	}

	return result;
}

double* SVDSolver(int size) {

	// V matrix is simply the normalized eigen vectors
	int i, j;
	double** V = malloc(size * sizeof(double*));
	for (i = 0; i < size; i++) {
		V[i] = malloc(size * sizeof(double));
	}

	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			V[i][j] = eigenvector[i][j];
		}
	}

	//sigma (diagonal matrix in expressed in vector form is the square root of the original eigenvalues)
	double* sigma = sqrtVector(eigenvalues, size);

	// U is original matrix * v vector, and then normalized
	double** A = readMatrix(size);
	double** U = MdotV(A, V, size);
	U = normalizeVectors(U, size);
	//printMatrix(U, size);

	// now using svd components to solve Ax = b
	// x = v*[diag(1/w)]*Ut*b

	// first, get the reciprocal of the diagonal matrix (vector)
	for (i = 0; i < size; i++) {
		sigma[i] = 1.0 / sigma[i];

	}

	// allocate some intermediate variable and a vector X to hold final solution
	double** param;
	double* x = malloc(size * sizeof(double));
	for (i = 0; i < size; i++) {
		x[i] = 0;
	}

	// this part calculates  v*[diag(1/w)]
	param = matrixDotDiag(V, sigma, size);

	// this part calculates  v*[diag(1/w)]*Ut
	param = MdotV(param, transpose(U, size), size);

	// this part calculates  v*[diag(1/w)]*Ut*b
	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			x[i] += param[i][j] * b[j];
		}

	}

	return x;

}

int main() {
	double** A = readMatrix(size);
	printMatrix(A, size);

	A = dotTranspose(A, size);

	printMatrix(A, size);

	printf("--------------------------after eigen----------------------\n");

	int i = 0;
	double* a;
	clock_t begin = clock();
	for (int i = 0; i < 100000; i++) {
		//a = QRSolver(input, b, n);
		a = SVDSolver(size);
		//a = choleskySolver(input, b, n);
	}
	clock_t end = clock();
	double time = (double) (end - begin) / ((double) CLOCKS_PER_SEC) * 1000;
	printf("SVD uses %f ms\n", time);
	// see the solution vector
	for (i = 0; i < size; i++) {
		printf("%.11f \n", a[i]);
	}

}
