# ECSE-444-Final-Project-Linear-System-Solver-and-Optomization
Instructions (for SVD):
SVD files: SVD.c. floatSVD.c, float_mem_SVD.c and float_mem_vectorized_SVD.c
To see demos, simply compile the corresponding file with command "gcc -o" and run the .exe. Each file is a stand-alone program with its own set of implementation and independent drivers. Example demo:  
Two matrix printed are original matrix and the covariance matrix (its transpose dot itself). You can use printMatrix(size) to see more n-by-n matrices.

Eigen value files:
Those two files included in folder “eigen utils” are deprecated. They are canned routines to perform eigen decomposition but not successfully integrated in the environment of this C project.

Python:
There is a python notebook included in SVD folder to compare the result in the demo to the correct result. You can visit https://github.com/AntoineWY/ECSE-444-Final-Project-Linear-System-Solver-and-Optimization/blob/master/SVD/Eigen%20Calculation%20and%20answer%20validation.ipynb to see example output corresponding to the demo.

Data:
For testing data of the entire group please refer to spread sheet “ecse-444-test-data” in red section for SVD.
Exploring new matrix?
As explained in the report, the implementation of the algorithm is interrupted by the python injection of eigen decomposition result. Thus, every time a new matrix wanted to be tried out, you should perform eigen decomposition in python notebook first, get the result copied in C code and run again.
