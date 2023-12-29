#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <png.h>
#include <string.h>
#include "functions.h"

int main(int argc, char *argv[])
{

    // Input handling

    if (argc < 2)
    {
        printf("ERROR: filename was not provided. Please try again with the following format: ./main filename.mtx plot \nThe last input 'plot' is optional if user wants to plot the sparsity pattern.\n");
        return 1;
    }
    
    const char *filename = argv[1];

    CSRMatrix A;

    // Reading the .mtx file, if it is read as triangular assume it is meant to be symmetric and fill in the missing data.

    ReadMMtoCSR(filename, &A);
    printf("\nFinished reading %s file.\n", argv[1]);

    const int file_read_non_zeros = A.num_non_zeros;

    if (check_symmetry(&A) == 0)
    {
        fill_symmetry(&A);
    } 

    // Print some matrix data

    printf("\nMatrix data from %s file. Note that if this matrix was given as triangular it was assumed symmetric and the relevant data was copied to the other triangle.\n", argv[1]);
    printf("\nMatrix Dimensions: %d by %d", A.num_rows, A.num_cols);
    printf("\nNumber of File Read Non-Zeros: %d", file_read_non_zeros);
    printf("\nNumber of Total Non-Zeros: %d", A.num_non_zeros);

    // For small matrices output the data read from the file in CSR format.

    if (A.num_rows <= 15)
    {
        printf("\nRow Pointer: ");
        for (int i = 0; i <= A.num_rows; i++)
        {
            printf("%d ", A.row_ptr[i]);
        }

        printf("\nColumn Index: ");
        for (int i = 0; i < A.num_non_zeros; i++)
        {
            printf("%d ", A.col_ind[i]);
        }

        printf("\nValue: ");
        for (int i = 0; i < A.num_non_zeros; i++)
        {
            printf("%.4f ", A.csr_data[i]);
        }
    }
    
    printf("\n");

    // Calls for creating a sparsity pattern, two different functions depending on the matrix's size. These functions are both from ChatGPT.

    if (A.num_rows < 2000 && argc == 3) 
    {
        printf("\nVisualizing matrix sparsity pattern for a small matrix, file will be created called 'sparsity_pattern.png'.\n");
        visualize_matrix_small(&A, "sparsity_pattern.png");
    }
    else if (argc ==3) 
    {
        printf("\nVisualizing matrix sparsity pattern for a large matrix, file will be created called 'sparsity_pattern.png'.\n");
        visualize_matrix_large(&A, "sparsity_pattern.png");
    }

    // Initializing all the needed column vectors.

    double *b = (double *)malloc(A.num_rows * sizeof(double));
    double *x = (double *)malloc(A.num_rows * sizeof(double));
    double *r = (double *)malloc(A.num_rows * sizeof(double));
    double norm = 0;

    if (b == NULL || x == NULL || r == NULL)
    {
        printf("ERROR: memory allocation failed in main. Exiting the program.\n");
        return(1);
    }

    // Set all elements of b to 1.

    for (int i = 0; i < A.num_rows; ++i)
    {
        b[i] = 1.0;
    }

    // Defining parameters for the solver, these could easily be implemented as user inputs.

    int max_iteration = 25000;
    double stop_criteria = 1e-16;

    // Seed rand() and start clock to calculate CPU Time.

    srand(time(NULL));
    clock_t start_time, end_time;
    double cpu_time_used;
    start_time = clock();

    // Solve Ax=b

    solver(&A, b, x, max_iteration, stop_criteria);

    // Print the CPU time taken for the code.

    end_time = clock();
    cpu_time_used = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    printf("CPU time taken to solve Ax=b: %f seconds\n", cpu_time_used);

    // Compute the residual and norm to determine the accuracy of the solver.

    compute_residual(&A, x, b, r);
    compute_norm(&A, r, &norm);

    printf("\nResidual Norm: %.4e\n\n", norm);

    // Deallocate memory

    free(b);
    free(x);
    free(r);
    free(A.csr_data);
    free(A.col_ind);
    free(A.row_ptr);
}