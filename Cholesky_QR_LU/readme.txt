Welcome to the Square Matrix Solver!!!

The attached C file is designed to be exceted using gcc in the Trottier Linux Environment. 

Instructions on How to Run (* file name is 'ECSE444_MatrixSolver')

1.  In the Linux terminal, navigate to the folder where you put the C file

2. Compile the C file using gcc, make sure you use the c99 mode of gcc. 
/*--------------------------------------------------------------------------------------
** If you are too lazy, just copy and paste the following command to the terminal:
** Default level optimization command: gcc ECSE444_MatrixSolver.c -o main -std=c99 -lm
** Level 0 optimization command (if you need): gcc -O0 ECSE444_MatrixSolver.c -o main -std=c99 -lm
*/--------------------------------------------------------------------------------------

3. Execute the executable file you just created. 
(The executable file is called "main" if you copied and pasted the command above...
just type in "./main" in your terminal and press enter..)

4. You will see a short description of the problem format and term explanations

5. Type in the size of the square matrix that you want to solve (the size is at least 1).
(Example: If you want to solve a 4 by 4 matrix, just enter '4')

6. Type in each term of the matrix that you want to solve in the row-major order. 

7. Type in the left-hand side vector of the matrix problem. 

8. The program will ask if the matrix entered is well-conditioned or not.
   If it is well-conditioned, enter '1'; if not or you do not know, type '0'.

9. There will be two outcomes:
	a. The answer will be displayed in a line, numbers are seperated by ';'.
	b. 'Matrix cannot be solved' will be displayed, meaning the matrix has no solution.

//----------------------------------------------------------------------------------------
Example:
The problem has the format: Ax = b, where 'A' is 3 by 3, 
and you know that the matrix is well-conditioned (according to other programs?).

    [1 0 0]      [1]
A = [0 1 0], b = [2]   
    [0 0 1]      [3]

You wanna solve for x...

The entire inputs and results in the terminal should look like this: 
________________________________________________________________________________________
|The matrix problem is in the form of: Ax=b						|
|1. 'A' is a square matrix.								|
|2. 'b' is the result wanted.								|
|3. 'x' is the solution to the linear system						|
|											|
|Size of matrix: 3									|
|A:											|
|[0][0]= 1										|
|[0][1]= 0										|
|[0][2]= 0										|
|[1][0]= 0										|
|[1][1]= 1										|
|[1][2]= 0										|
|[2][0]= 0										|
|[2][1]= 0										|
|[2][2]= 1										|
|b:											|
|[0]= 1											|
|[1]= 2											|
|[2]= 3											|
|Is the matrix well-conditioned? (0/1) (0 for 'no', 1 for 'yes'; if not sure, type'0'):	|
|1											|
|A =											|
|[1.00, 0.00, 0.00]									|
|[0.00, 1.00, 0.00]									|
|[0.00, 0.00, 1.00]									|
|											|
|b = 											|
|[1.0000; 2.0000; 3.0000]								|
|											|
|LU Decomposition is used.								|
|											|
|Answer x = 										|
|[1.0000; 2.0000; 3.0000]								|
|_______________________________________________________________________________________|

