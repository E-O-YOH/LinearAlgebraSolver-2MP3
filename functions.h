#ifndef FUNCTIONS_H
#define FUNCTIONS_H


typedef struct {
    double *csr_data;   // Array of non-zero values
    int *col_ind;       // Array of column indices
    int *row_ptr;       // Array of row pointers
    int num_non_zeros;  // Number of non-zero elements
    int num_rows;       // Number of rows in matrix
    int num_cols;       // Number of columns in matrix
} CSRMatrix;

// Structure used in for the qsort algorithm.

typedef struct {
    int row_data;
    int col_data;
    double value_data;
} SortData;

// The purpose of each function is explained in functions.c

int compare(const void *a, const void *b);
void ReadMMtoCSR(const char *filename, CSRMatrix *matrix);
void spmv_csr(const CSRMatrix *A, const double *x, double *y);
int check_symmetry(CSRMatrix *A);
void fill_symmetry(CSRMatrix *A);
void solver(const CSRMatrix *A, const double *b, double *x, int max_iteration, double stop_criteria);
void compute_residual(const CSRMatrix *A, const double *x, const double *b, double *r);
void compute_norm(const CSRMatrix *A, const double *r, double *norm);
void visualize_matrix_small(CSRMatrix *matrix, const char *output_filename);
void visualize_matrix_large(CSRMatrix *matrix, const char *output_filename);

#endif
