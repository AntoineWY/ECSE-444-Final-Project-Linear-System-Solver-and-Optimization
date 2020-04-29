/*
 ============================================================================
 Name        : ECSE444_MatrixSolver.c
 Author      :
 Version     : 1.0
 Copyright   : Your copyright notice
 Description : Matrix Solver Tester in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

//----------General method-------------
double *forward_optimized(double *L, double *b, int n){
	double *x = (double*)calloc(n,sizeof(double));
	for(int i = 0; i<n;i++){
		x[i] = b[i];
		for (int j = 0; j<i; j++){
			x[i] -= L[i*n+j]*x[j];
		}
		x[i] /= L[i*n+i];
	}
	return x;
}

double *backward_optimized(double *U, double *b, int n){
	double *x = (double*)calloc(n,sizeof(double));
	for(int i = n-1; i>(-1);i--){
		x[i] = b[i];
		for (int j = n-1; j>i; j--){
			x[i] -= U[i*n+j]*x[j];
		}
		x[i] /= U[i*n+i];
	}
	return x;
}

double* copymatrix(double* A, int n){
	double* out = (double*)calloc(n * n, sizeof(double));
	for (int i = 0; i<n; i++){
		for (int j = 0; j<n; j++){
			out[i*n+j] = A[i*n+j];
		}
	}
	return out;
}

void printM (double *A, int n){
	for (int i = 0; i < n; i++) {
		printf("[");
		for (int j = 0; j < n-1; j++){
			printf("%.2f, ",A[i*n+j]);
		}
		printf("%.2f] \n",A[i*n+n-1]);
	}
}

void printV (double *B, int n){
	printf("[");
	for (int i = 0; i<n-1; i++){
		printf("%.4f; ",B[i]);
	}
	printf("%.4f",B[n-1]);
	printf("] \n");
}
//-----------------General function end----------------------

//------------------Cholesky Decomposition--------------------
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

double *backward_cholesky(double *U, double *b, int n){
	double *x = (double*)calloc(n,sizeof(double));
	for(int j = n-1; j>(-1);j--){
		x[j] = b[j];
		for (int i = n-1; i>j; i--){
			x[j] -= U[i*n+j]*x[i];
		}
		x[j] /= U[j*n+j];
	}
	return x;
}

double *choleskySolver(double *A, double *b, int n){
	double *L = cholesky_optimized(A, n);
	if(L==NULL) return NULL;
	double *y = forward_optimized(L, b, n);
	double *x = backward_cholesky(L, y, n);
	free(L);
	free(y);
	return x;
}

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
double dot(int m, double* mat, int col1, int col2){
	double output = 0;
	for(int i = 0 ; i < m; i++){
		output += mat[i*m+col1]*mat[i*m+col2];
	}
	return output;
}

double norm(int m, double* mat,int col){
	double output = 0;
	for(int i = 0 ; i < m; i++){
		output += mat[i*m+col]*mat[i*m+col];
	}
	if(output<=0) return 0;
	return sqrt(output);
}

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

int unit(int m, double* mat, int col){
	double matNorm = norm(m,mat,col);
	if (matNorm == 0) return 0;
	for(int j = 0 ; j < m; j++){
		mat[j*m+col] = mat[j*m+col]/matNorm;
	}
	return 1;
}

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


double* dotproduct(int m, double* mat, double* vec){
	double* out = (double*)calloc(m, sizeof(double));
	for (int i = 0; i<m; i++){
		for (int j = 0; j<m;j++){
			out[i] += mat[i*m+j]*vec[j];
		}
	}
	return out;
}

double *QRSolver(double* mat, double* b, int m){
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
//----------- LU solver --------------
int doolittleLU(int n, double* A, double* L, double* U){ // n is the size of the input matrix
    // fill up the entries of L and U
    for (int i = 0; i < n; i++) {  // for row i in U, i.e i th iteration
        for (int k = i; k < n; k++) { // for k th entry in row i of U
            double sum = 0;
            for (int m = 0; m < i; m++){
                sum = sum + L[i*n+m]*U[m*n+k];
            }
            U[i*n+k] = A[i*n+k] - sum;
        }
        if (i<n-1){ // fill in the non diagonal entries of L
            for (int j = i+1; j < n; j++){ // for row j in L
                if (U[i*n+i] == 0 ){ // check if divisor is 0
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
                L[c*n+c] = 1; //diagonal entries of L are all 1s
            }
        }
    }
    return 1;
}

double* LUSolver(double* A, double*b, int n){
	double* L = (double*)calloc(n * n, sizeof(double));
	double* U = (double*)calloc(n * n, sizeof(double));
	int check = doolittleLU(n, A, L, U);
	if (check == 0) return NULL;
	double* y = forward_optimized(L, b, n);
	double* x = backward_optimized(U, y, n);
	free(L);
	free(U);
	free(y);
	return x;
}
//---------- LU solver end---------------
double computeError(double* A, double* b, double* x, int n){
	double error = 0;
	double* bb = dotproduct(n, A, x);
	for (int i = 0; i<n; i++){
		error += fabs(b[i]-bb[i]);
	}
	return error;
}

//--------- time comparison-------------------

int main(void){
	int n = 4;
	//------------------------spd--------------
	double input[] = {18, 22, 54, 42, 22, 70, 86, 62, 54, 86, 174, 134, 42, 62, 134, 106};
	//double input[] = {38, 18, 35, 20, 42, 18, 41, 28, 18, 40, 35, 28, 42, 24, 43, 20, 18, 24, 15, 24, 42, 40, 43, 24, 66};
	//double input[] = {44, 32, 22, 24, 14, 42, 32, 43, 31, 19, 10, 49, 22, 31, 45,  9, 13, 40, 24, 19,  9, 15,  7, 24, 14, 10, 13,  7, 15, 20, 42, 49, 40, 24, 20, 65};
	//double input[] = {54, 18, 45, 24, 21, 45, 21, 18, 38, 34,  9, 21, 12, 18, 45, 34, 62, 20, 27, 36, 32, 24,  9, 20, 22, 23, 28, 11, 21, 21, 27, 23, 35, 26, 18, 45, 12, 36, 28, 26, 51, 17, 21, 18, 32, 11, 18, 17, 20};
	//double input[] = {61, 55, 43, 34, 33, 17, 60, 42, 55, 64, 44, 43, 33, 37, 65, 47, 43, 44, 53, 29, 26, 23, 51, 39, 34, 43, 29, 44, 28, 35, 35, 39, 33, 33, 26, 28, 30, 19, 30, 24, 17, 37, 23, 35, 19, 44, 29, 28, 60, 65, 51, 35, 30, 29, 78, 48, 42, 47, 39, 39, 24, 28, 48, 56};
	//double input[] = {27, 22, 28, 26, 17, 40, 26, 25, 22, 22, 30, 18, 32, 12, 40, 23, 20, 23, 28, 18, 59, 33, 30, 55, 36, 47, 53, 26, 32, 33, 67, 30, 58, 31, 29, 51, 17, 12, 30, 30, 34, 36, 17, 22, 40, 40, 40, 55, 58, 36, 87, 47, 48, 66, 26, 23, 36, 31, 17, 47, 39, 38, 39, 25, 20, 47, 29, 22, 48, 38, 44, 47, 22, 23, 53, 51, 40, 66, 39, 47, 82};
	//double input[] = {32, 32, 30, 32, 30, 28, 36, 17, 35, 36, 32, 62, 50, 32, 47, 38, 49, 33, 41, 62, 30, 50, 55, 29, 50, 38, 43, 35, 48, 64, 32, 32, 29, 38, 22, 24, 36, 19, 35, 37, 30, 47, 50, 22, 65, 45, 46, 34, 46, 62, 28, 38, 38, 24, 45, 50, 50, 34, 41, 51, 36, 49, 43, 36, 46, 50, 73, 28, 39, 62, 17, 33, 35, 19, 34, 34, 28, 47, 34, 48, 35, 41, 48, 35, 46, 41, 39, 34, 60, 57, 36, 62, 64, 37, 62, 51, 62, 48, 57, 83};
	//-----------------------well conditioned----------------------
	//double input[] = {36,13,8,40,13,16,2,1,9,39,11,19,13,11,34,35};
	//double input[] = {25,31,28,2,35,1,34,26,23,0,15,12,26,1,47,23,19,41,14,2,33,26,0,32,8};
	//double input[] = {19,7,35,36,2,31,12,23,28,3,46,43,39,42,1,47,6,22,16,44,31,29,21,18,14,26,4,30,46,15,2,12,18,40,3,40};
	//double input[] = {13,0,19,10,35,2,1,9,17,19,41,29,32,16,18,10,40,24,30,40,5,15,32,45,36,30,45,42,28,32,20,15,46,42,23,28,4,39,27,14,20,43,16,4,50,20,28,35,50};
	//double input[] = {60,15,18,23,4,37,33,20,25,47,10,29,13,33,47,20,32,27,38,18,1,50,25,11,22,47,26,31,10,4,39,41,41,6,20,3,10,16,39,25,20,45,36,9,35,46,26,28,41,13,46,10,12,28,12,4,27,21,24,60,41,33,24,31};
	//-----------------------ill conditioned---------------------
	//double input[] = {98,0,9,500,95,85,12,505,111,95,10,499,101,88,11,495};
	//double input[] = {89,3,19,8131,169,87,18120,1,88,65,94,36,50,89,714,20,284,295,27,98,56,3,578,19,320};
	//double input[] = {89,3,19,56,8131,169,87,18120,327,88,65,1,94,36,50,6,2852,714,20,9284,295,36,27,98,56,3,578,3,19,320,121,0,2,287,84950,99};
	//double input[] = {1000,15,18,23,4,37,33,25,100000,10,29,13,33,47,100,27,38,18,98989,50,25,200,47,26,31,100,9999,39,41,6,20,3,10,16,39,20,45,36,9,99999,46,26,41,13,46,10,12,28,12,8888,21,202,60,41,33,24};
	//double input[] = {100,15,18,23,4,37,33,453,25,100,10,29,13,33,47,67,100,27,38,18,98989,50,25,4,200,47,9869,31,100,9999,39,1,41,8888,20,3,10,16,39,95634,20,45,36,9,999,46,26,34576,41,233,0,10,12,28,12,9666,8888,21,202,60,41,33,24,23,909,453,324,876,34,76,9887,9999};
	double b[] = {1,2,3,4,5,6,7,8,9,0};
	printM(input,n);
	//double* a = QRSolver(input, b, n);
	double* a = choleskySolver(input, b, n);
	//double* a = LUSolver(input, b, n);
	double er = computeError(input, b, a, n);
	printf("The error is: %.10f\n", er);
	printV(a,n);
	clock_t begin = clock();
	for(int i = 0; i< 100000; i++){
		//a = QRSolver(input, b, n);
		a = LUSolver(input, b, n);
		//a = choleskySolver(input, b, n);
	}

	clock_t end = clock();
	double time = (double)(end-begin) / ((double)CLOCKS_PER_SEC)*1000;
	printf("%f ms", time);
	return 0;
}


//----------Matrix Solver----------------------
/*
int main(void) {
	setbuf(stdout, NULL);
	int n = 0;
	printf("Size of matrix: ");
	scanf("%d", &n);
	double *input = (double*)calloc(n * n, sizeof(double));
	double *b = (double*)calloc(n, sizeof(double));
	printf("A: \n");
	for (int i = 0; i<n; i++){
		for (int j = 0; j < n; j++){
			printf("[%d][%d]= ", i, j);
			scanf("%f",&input[i*n+j]);
		}
	}
	printf("b: \n");
	for (int i = 0; i<n; i++){
		printf("[%d]= ", i);
		scanf("%f", &b[i]);
	}
	int condition = 0;
	printf("Is the matrix well-conditioned? (0/1) (0 for 'no', 1 for 'yes'; if not sure, type '0')\n");
	scanf("%d", &condition);
	//int n = 4;
	//double input[] = {10,9,8,15,13,12,8,9,6,1,3,7,15,23,16,4};
	//double b[] = {1,2,3,4};
	printf("A = \n");
	printM(input, n);
	printf("b = \n");
	printV(b,n);

	int symetric = checkSymetric(input, n);
	double* answer = NULL;

	if(condition){
		if(symetric) answer = choleskySolver(input, b, n);
		if(answer != NULL){
			printf("Cholesky Decomposition Method is used. \nAnswer x= \n");
			printV(answer,n);
			return 0;
		}

		//answer = LUSolver(input, b, n);
		//if(answer != NULL){
			//printf("LU Decomposition Method is used. \nAnswer x= \n");
			//printV(answer,n);
			//return 0;
		//}

	//}
	answer = QRSolver(input, b, n);
	if(answer == NULL){
		printf("Cannot solve the matrix...");
	}else{
		printf("QR Decomposition Method is used. \nAnswer x= \n");
		printV(answer,n);
	}
	return 0;
}*/
