/*
 ============================================================================
 Name        : ECSE444_MatrixSolver.c
 Author      : Jiahua Liang, Kevin Chen, Xiangyun Wang, Yinuo Wang
 Version     : 1.2
 Copyright   : Your copyright notice
 Description : Matrix Solver in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//----------General function-------------
/*
 * Forward substitution function
 */
double *forward_optimized(double *L, double *b, int n){
	double *x = (double*)calloc(n,sizeof(double));
	for(int i = 0; i<n;i++){
		x[i] = b[i];
		for (int j = 0; j<i; j++){
			x[i] -= L[i*n+j]*x[j];
		}
		if(L[i*n+i] == 0) return NULL;
		x[i] /= L[i*n+i];
	}
	return x;
}

/*
 * Backward substitution function
 */
double *backward_optimized(double *U, double *b, int n){
	double *x = (double*)calloc(n,sizeof(double));
	for(int i = n-1; i>(-1);i--){
		x[i] = b[i];
		for (int j = n-1; j>i; j--){
			x[i] -= U[i*n+j]*x[j];
		}
		x[i] /= U[i*n+i];
		if(U[i*n+i] == 0) return NULL;
	}
	return x;
}

/*
 * Compute the dot product of a matrix and a vector
 */
double* dotproduct(int m, double* mat, double* vec){
	double* out = (double*)calloc(m, sizeof(double));
	for (int i = 0; i<m; i++){
		for (int j = 0; j<m;j++){
			out[i] += mat[i*m+j]*vec[j];
		}
	}
	return out;
}

/*
 * Print matrix in a desired format
 */
void printM (double *A, int n){
	for (int i = 0; i < n; i++) {
		printf("[");
		for (int j = 0; j < n-1; j++){
			printf("%.2f, ",A[i*n+j]);
		}
		printf("%.2f] \n",A[i*n+n-1]);
	}
}

/*
 * Print vector in a desired format
 */
void printV (double *B, int n){
	printf("[");
	for (int i = 0; i<n-1; i++){
		printf("%.4f; ",B[i]);
	}
	printf("%.4f",B[n-1]);
	printf("] \n");
}

/*
 * Copy a matrix to a block of newly allocated memory
 */
double* copymatrix(double* A, int n){
	double* out = (double*)calloc(n * n, sizeof(double));
	for (int i = 0; i<n; i++){
		for (int j = 0; j<n; j++){
			out[i*n+j] = A[i*n+j];
		}
	}
	return out;
}

/*
 * Compute the forward of the matrix problem
 * (Problem has form Ax=b, compute Ax-b after getting x)
 */
double computeError(double* A, double* b, double* x, int n){
	double error = 0;
	double* bb = dotproduct(n, A, x);
	for (int i = 0; i<n; i++){
		error += fabs(b[i]-bb[i]);
	}
	return error;
}

//-----------------General function end----------------------

//------------------Cholesky Decomposition--------------------
/*
 * Given matrix A, want A = L*transpose(L)
 * This function returns L
 */
double *cholesky_optimized(double *A, int n) {
    double *L = (double*)calloc(n * n, sizeof(double));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < (i+1); j++) {
        	double term = 0;
            for (int k = 0; k < j; k++) term += L[i * n + k] * L[j * n + k];
            if(i==j){
            	L[i * n + j] = sqrt(A[i * n + i] - term);
            }
            else{
            	L[i * n + j] = 1.0 / L[j * n + j] * (A[i * n + j] - term);
            }
            if(L[i * n + j] == 0)	return NULL;
        }
    }
    return L;
}

/*
 * An optimized backward substitution only for Cholesky Decomposition Solver.
 * Backward substitution is usually for upper triangular matrix.
 * Cholesky requires to backward substitute the transpose of a lower triangular matrix.
 * Instead of transpose a matrix and do backward substitution, we can directly do backward
 * substitution for a lower triangular matrix, with some changes on the for loop order.
 */
double *backward_cholesky(double *U, double *b, int n){
	double *x = (double*)calloc(n,sizeof(double));
	for(int j = n-1; j>(-1);j--){
		x[j] = b[j];
		for (int i = n-1; i>j; i--){
			x[j] -= U[i*n+j]*x[i];
		}
		x[j] /= U[j*n+j];
		if(U[j*n+j] == 0) return NULL;
	}
	return x;
}

/*
 * Solve matrix using Cholesky decomposition
 */
double *choleskySolver(double *A, double *b, int n){
	double *L = cholesky_optimized(A, n);
	if(L==NULL) return NULL;
	double *y = forward_optimized(L, b, n);
	if(y==NULL) return NULL;
	double *x = backward_cholesky(L, y, n);
	free(L);
	free(y);
	return x;
}

/*
 * Check if the matrix is symmetric or not
 */
int checkSymetric(double *input, int n){
	for (int i = 0; i<n; i++){
		for (int j = 0; j<i; j++){
			if(input[i*n+j]!=input[j*n+i]){
				return 0;
			}
		}
	}
	return 1;
}
//-------------Cholesky Decomposition End--------------
//-------------QR deconposition -----------------------
/*
 * Compute the dot product of two columns of a matrix
 */
double dot(int m, double* mat, int col1, int col2){
	double output = 0;
	for(int i = 0 ; i < m; i++){
		output += mat[i*m+col1]*mat[i*m+col2];
	}
	return output;
}

/*
 * Compute the norm of a column of a matrix
 */
double norm(int m, double* mat,int col){
	double output = 0;
	for(int i = 0 ; i < m; i++){
		output += mat[i*m+col]*mat[i*m+col];
	}
	if(output<=0) return 0;
	return sqrt(output);
}

/*
 * transpose a matrix
 */
void transpose(int m, double* mat){
	double temp;
	for(int i = 0 ; i < m; i++){
		for(int j = i+1 ; j < m; j++){
			temp = mat[i*m+j];
			mat[i*m+j] = mat[j*m+i];
			mat[j*m+i] = temp;
		}
	}
}

/*
 * compute the norm vector of a column of a matrix
 */
int unit(int m, double* mat, int col){
	double matNorm = norm(m,mat,col);
	if (matNorm == 0) return 0;
	for(int j = 0 ; j < m; j++){
		mat[j*m+col] = mat[j*m+col]/matNorm;
	}
	return 1;
}

/*
 * Compute Q in the QR decomposition
 */
double* decompQ(int m, double* mat){
	double* output = copymatrix(mat, m);
	int check = unit(m,output,0);
	if(check == 0) return NULL;
	double product = 0;
	for(int i = 1 ; i < m; i++){
		for(int j = 0; j < i; j++){
			product = dot(m,output,i,j);
			for(int k = 0; k < m; k++){
				output[k*m+i] -= output[k*m+j]*product;
			}
			unit(m,output,i);
		}
	}
	return output;
}

/*
 * Compute R in the QR decomposition
 */
double* decompR(int m, double* mat,double* Q){
	double* R = (double*)calloc(m * m, sizeof(double));
	double* tempQ = copymatrix(Q, m);
	transpose(m,tempQ);
	double val;
	for(int i = 0 ; i < m; i++){
		for(int j = 0;j < m; j++){
			val = 0;
			for(int k = 0 ; k < m ; k++){
				val += tempQ[i*m+k]*mat[k*m+j];
			}
			R[i*m+j] = val;
		}
	}
	return R;
}

/*
 * Solve matrix using QR decomposition
 */
double* QRSolver(double* mat, double* b, int m){
	double* Q = decompQ(m,mat);
	if(Q==NULL) return NULL;
	double* R = decompR(m,mat,Q);
	transpose(m, Q);
	double* y = dotproduct(m,Q,b);
	double* x = backward_optimized(R, y, m);
	free(Q);
	free(R);
	free(y);
	return x;
}
//----------QR Decomposition end---------------
//---------------LU Decomposition Start----------------
/*
 * Compute L and U for LU decomposition
 */
int doolittleLU(int n, double* A, double* L, double* U){ // n is the size of the input matrix
    for (int i = 0; i < n; i++) {  		// for row i in U, i.e i th iteration
        for (int k = i; k < n; k++) { 	// for k th entry in row i of U
            double sum = 0;
            for (int m = 0; m < i; m++){
                sum = sum + L[i*n+m]*U[m*n+k];
            }
            U[i*n+k] = A[i*n+k] - sum;
        }
        if (i<n-1){ 					// fill in the non diagonal entries of L
            for (int j = i+1; j < n; j++){ 	// for row j in L
                if (U[i*n+i] == 0 ){ 		// check if divisor is 0
                    return 0;
                }
                double sum = 0;
                for (int m = 0; m < i; m++){
                	sum = sum + L[j*n+m]*U[m*n+i];
                }
                L[j*n+i] = (A[j*n+i] - sum)/U[i*n+i];
            }
        }else{
            for (int c = 0; c < n; c++){
                L[c*n+c] = 1; 				//diagonal entries of L are all 1s
            }
        }
    }
    return 1;
}

/*
 * Use L and U with forward and backward substitution to solve the matrix
 */
double* LUSolver(double* A, double*b, int n){
	double* L = (double*)calloc(n * n, sizeof(double));
	double* U = (double*)calloc(n * n, sizeof(double));
	int check = doolittleLU(n, A, L, U);
	if (check == 0) return NULL;
	double* y = forward_optimized(L, b, n);
	if(y==NULL) return NULL;
	double* x = backward_optimized(U, y, n);
	free(L);
	free(U);
	free(y);
	return x;
}
//----------------LU decomposition End--------------------
//----------Matrix Solver----------------------
int main(void) {
	setbuf(stdout, NULL);
	printf("The matrix problem is in the form of: Ax=b \n1. 'A' is a square matrix. \n2. 'b' is the result wanted. \n3. 'x' is the solution to the linear system. \n\n");
	int n = 0;
	printf("Size of matrix: ");
	scanf("%d", &n);	//user input for matrix size
	double *input = (double*)calloc(n * n, sizeof(double));
	double *b = (double*)calloc(n, sizeof(double));
	printf("A: \n");
	for (int i = 0; i<n; i++){				//user input for matrix to be solved
		for (int j = 0; j < n; j++){
			printf("[%d][%d]= ", i, j);
			scanf("%lf",&input[i*n+j]);
		}
	}
	printf("b: \n");						//user input for the right hand side of the problem
	for (int i = 0; i<n; i++){
		printf("[%d]= ", i);
		scanf("%lf", &b[i]);
	}
	int condition = 0;			//condition number is predefined to )
	printf("Is the matrix well-conditioned? (0/1) (0 for 'no', 1 for 'yes'; if not sure, type '0'): \n");
	scanf("%d", &condition);	//ask user for condition number
	printf("A = \n");			//print problem given
	printM(input, n);
	printf("\nb = \n");
	printV(b,n);
	printf("\n");

	int symetric = checkSymetric(input, n);		//check if symmetric
	double* answer = NULL;

	if(condition){		//if problem well conditioned, use LU decomposition
		answer = LUSolver(input, b, n);
		if(answer != NULL){
			printf("LU Decomposition Method is used. \n\nAnswer x = \n");
			printV(answer,n);
			return 0;
		}
	}else{			//if not well conditioned, use cholesky or QR
		if(symetric) answer = choleskySolver(input, b, n);		//if symmetrix, try Cholesky
		if(answer != NULL){
			printf("Cholesky Decomposition Method is used. \n\nAnswer x = \n");
			printV(answer,n);
			return 0;
		}		//if not symmetric, or Cholesky failed, use QR
		answer = QRSolver(input, b, n);
		if(answer != NULL && computeError(input, b, answer, n) > 0.1) answer = NULL;
		if(answer != NULL){
			printf("QR Decomposition Method is used. \n\nAnswer x= \n");
			printV(answer,n);
			return 0;
		}
	}
	printf("Matrix cannot be solved...\n"); //if matrix is not solvable
	return 0;
}
