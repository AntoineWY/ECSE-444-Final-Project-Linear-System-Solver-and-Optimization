/*
 ============================================================================
 Name        : float_mem_SVD.c
 Author      : Jiahua Liang, Kevin Chen, Xiangyun Wang, Yinuo Wang
 Version     : 1.2
 Copyright   : Your copyright notice
 Description : Matrix Solver in C, Ansi-style
 ============================================================================


 ============================================================================
 This is the original version of the SVD solver.
 Note that since the algorithm of finding eigenvalue in c is not applicable
 due to segmentation faults, the computed eigen matrix and values are declared
 directly in the code.

 Feature notice:
 1. All double type are replace with float type to boost the performance.
 2. At the end of the iteration, intermediate variables with dynamic allocated
 memory are freed. So that the 100000 iteration will not deplete system's memory.
 ============================================================================
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

int size = 5;

static float eigenvector[5][5] = { { (float) -0.37850527, (float) 0.22754263,
		(float) -0.18262938, (float) -0.87840887, (float) 0.0015737 }, {
		(float) 0.23120478, (float) 0.49272987, (float) -0.78356889,
		(float) 0.19133505, (float) 0.23059125 }, { (float) -0.03987157,
		(float) -0.24071506, (float) -0.41177417, (float) 0.03886613,
		(float) -0.87715334 }, { (float) -0.39238419, (float) 0.74061481,
		(float) 0.32566007, (float) 0.29263578, (float) -0.32532163 }, {
		(float) -0.80481241, (float) -0.31462293, (float) -0.27758576,
		(float) 0.3234849, (float) 0.2675688 }, };

static float eigenvalues[5] = { (float) 216.37561647, (float) 100.08410035,
		(float) 49.38342505, (float) 34.38342505, (float) 0.35340114 };


static float b[5] = { 1, 2, 5, 3, 1 };

float** dotTranspose(float** A, int size) {
	// create a transpose matrix
	int i = 0;
	float** A_transpose = malloc(size * sizeof(float*));
	for (i = 0; i < size; i++) {
		A_transpose[i] = malloc(size * sizeof(float));
	}
	int j = 0;
	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			A_transpose[i][j] = A[j][i];
		}
	}
	// allocate a new matrix for returning
	float** product = malloc(size * sizeof(float*));
	for (i = 0; i < size; i++) {
		product[i] = malloc(size * sizeof(float));
	}
	// do multiplication
	int k = 0;
	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {

			float m = 0;
			for (k = 0; k < size; k++) {
				m = m + (float)A_transpose[i][k] * A[k][j];
			}

			product[i][j] = m;
		}
	}
	return product;
}

void printMatrix(float** a, int size) {
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

float** readMatrix(int size) {
	int i = 0;
	int j = 0;

	float** mat = malloc(size * sizeof(float*));
	for (i = 0; i < size; i++) {
		mat[i] = malloc(size * sizeof(float));
	}

	FILE *file;
	file = fopen("data.txt", "r");

	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			if (!fscanf(file, "%f", &mat[i][j]))
				break;
		}

	}
	fclose(file);
	return mat;

}

float** transpose(float** a, int size) {
	int i = 0;
	float** mat = malloc(size * sizeof(float*));
	for (i = 0; i < size; i++) {
		mat[i] = malloc(size * sizeof(float));
	}

	int j = 0;
	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			mat[i][j] = a[j][i];
		}
	}

	return mat;

}

float** MdotV(float** a, float** v, int size) {
	int i, j;
	float** product = malloc(size * sizeof(float*));
	for (i = 0; i < size; i++) {
		product[i] = malloc(size * sizeof(float));
	}
	// do multiplication
	int k = 0;
	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {

			float m = 0;
			for (k = 0; k < size; k++) {
				m = m + (float)(a[i][k] * v[k][j]);
			}

			product[i][j] = m;
		}
	}
	return product;

}

float** normalizeVectors(float** v, int size) {
	int i, j, k;
	float** result = malloc(size * sizeof(float*));
	for (i = 0; i < size; i++) {
		result[i] = malloc(size * sizeof(float));
	}

	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {

			float pythag = 0;
			for (k = 0; k < size; k++) {
				pythag += (float)v[k][j] * v[k][j];
			}

			result[i][j] = v[i][j] / (float)sqrt(pythag);
		}
	}

	return result;
}

float* sqrtVector(float* v, int size) {
	int i;
	float* result = malloc(size * sizeof(float));

	for (i = 0; i < size; i++) {
		result[i] = sqrt(v[i]);
	}
	return result;
}

float** matrixDotDiag(float** a, float* v, int size) {
	int i, j;
	float** result = malloc(size * sizeof(float*));
	for (i = 0; i < size; i++) {
		result[i] = malloc(size * sizeof(float));
	}

	for (j = 0; j < size; j++) {
		float m = v[j];
		for (i = 0; i < size; i++) {

			result[i][j] = (float)a[i][j] * m;
		}
	}

	return result;
}

float* SVDSolver(int size) {

	// V matrix is simply the normalized eigen vectors
	int i, j;
	float** V = malloc(size * sizeof(float*));
	for (i = 0; i < size; i++) {
		V[i] = malloc(size * sizeof(float));
	}

	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			V[i][j] = eigenvector[i][j];
		}
	}

	//sigma (diagonal matrix in expressed in vector form is the square root of the original eigenvalues)
	float* sigma = sqrtVector(eigenvalues, size);

	// U is original matrix * v vector, and then normalized
	float** A = readMatrix(size);
	float** U = MdotV(A, V, size);
	U = normalizeVectors(U, size);
	//printMatrix(U, size);

	// now using svd components to solve Ax = b
	// x = v*[diag(1/w)]*Ut*b

	// first, get the reciprocal of the diagonal matrix (vector)
	for (i = 0; i < size; i++) {
		sigma[i] = (float)(1.0 / sigma[i]);

	}

	// allocate some intermediate variable and a vector X to hold final solution
	float** param;
	float* x = malloc(size * sizeof(float));
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
			x[i] += (float)(param[i][j] * b[j]);
		}

	}

	free(A);
	free(U);
	free(V);
	free(sigma);

	return x;

}

int main() {
	float** A = readMatrix(size);
	printMatrix(A, size);

	A = dotTranspose(A, size);

	printMatrix(A, size);

	printf("--------------------------after eigen----------------------\n");

	int i = 0;
	float* a;
	clock_t begin = clock();
	for (int i = 0; i < 100000; i++) {
		//a = QRSolver(input, b, n);
		a = SVDSolver(size);
		//a = choleskySolver(input, b, n);
	}
	clock_t end = clock();
	float time = (float) (end - begin) / ((float) CLOCKS_PER_SEC) * 1000;
	printf("SVD uses %f ms\n", time);
	// see the solution vector
	for (i = 0; i < size; i++) {
		printf("%.11f \n", a[i]);
	}

}
