#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <png.h>
#include "functions.h"

// Comparison function required for qsort algorithm.

int compare(const void *a, const void *b)
{
    return ((SortData*)a)->row_data - ((SortData*)b)->row_data;
}

// Reads a file in matrix market format and stores it in memory in CSR format.

void ReadMMtoCSR(const char *filename, CSRMatrix *matrix)
{
    FILE *file = fopen(filename,"r");
    if (file == NULL)
    {
        printf("ERROR: could not open the file. Exiting the program.\n");
        exit(1);
    }

    // The following few lines for skipping the comments in the .mtx file are from ChatGPT.

    char line[256]; 
    while (fgets(line, sizeof(line), file) != NULL)
    {
        if (line[0] != '%')
        {
            break;
        }
    }

    sscanf(line, "%d %d %d", &matrix->num_rows, &matrix->num_cols, &matrix->num_non_zeros);

    // Allocate memory for the arrays in matrix and a few additional arrays which will be used for sorting later in the function.

    matrix->csr_data = malloc(matrix->num_non_zeros * sizeof(double));
    matrix->col_ind = malloc(matrix->num_non_zeros * sizeof(int));
    matrix->row_ptr = malloc((matrix->num_rows + 1) * sizeof(int));
    SortData *temp = malloc(matrix->num_non_zeros * sizeof(SortData));

    if (matrix->csr_data == NULL || matrix->col_ind == NULL || matrix->row_ptr == NULL || temp == NULL)
    {
        printf("ERROR: memory allocation failed in ReadMMtoCSR. Exiting the program.\n");
        exit(1);
    }

    int i,j;

    // Set all values of row_ptr to zero.

    for (i = 0; i <= matrix->num_rows; i++)
    {
        matrix->row_ptr[i] = 0;
    }

    // Read all other values from the file.

    for (i = 0; i < matrix->num_non_zeros; i++)
    {
        int row, col;
        double value;
        if (fscanf(file, "%d %d %lf", &row, &col, &value) != 3)
        {
            printf("ERROR: expected value was missing when reading the file. Exiting the program.\n");
            exit(1);
        }

        col--;
        temp[i].row_data = row;
        temp[i].col_data = col;
        temp[i].value_data = value;

        // Populate row_ptr with the correct values.

        for (j = row; j <= matrix->num_rows; j++)
        {
            matrix->row_ptr[j]++;
        }
    }

    fclose(file);

    // Sort col_ind and csr_data in order of increasing rows to align with the row_ptr using the qsort algorithm.

    qsort(temp, matrix->num_non_zeros, sizeof(SortData), compare);

    // Writes this sorted data to the col_ind and csr_data arrays in the CSRMatrix structure.

    for (i = 0; i < matrix->num_non_zeros; i++)
    {
        matrix->col_ind[i] = temp[i].col_data;
        matrix->csr_data[i] = temp[i].value_data;
    }

    // Frees this temporary structure used for qsort sorting.

    free(temp);
}

// Function performs a multiplication between a sparse matrix in CSR format and a column vector.

void spmv_csr(const CSRMatrix *A, const double *x, double *y)
{
    int i,j;

    // Ensure all values of solution vector y are initalized to zero.

    for (i = 0; i < A->num_cols; i++)
    {
        y[i] = 0;
    }

    // Multiple each row by the column vector.

    for (i = 0; i < A->num_rows; i++)
    {
        for (j = A->row_ptr[i]; j < A->row_ptr[i+1]; j++)
        {
            y[i] += A->csr_data[j] * x[A->col_ind[j]];
        }
    }
}

// Function checks whether the given matrix is upper or lower triangular and if it is we assume it is symmetric as indicated by the MM document given.

int check_symmetry(CSRMatrix *A)
{
    int i,j;
    int code = 0; 
    
    // Checks the upper triangle by seeing if there are any non-zero values stored above the diagonal.

    for (i = 0; i < A->num_rows; i++)
    {
        for (j = A->row_ptr[i]; j < A->row_ptr[i+1]; j++)
        {
            if (A->col_ind[j] > i)
            {
                code = 1;
            }
        }
    }

    // If there were no values in the upper triangle we can assume it is lower triangular and therefore symmetric.

    if (code == 0)
    {
        printf("This matrix is lower triangular assuming this means it is symmetric, program will continue.\n");
        return 0;
    }

    // Checks to see if there are any values below the diagonal.

    for (i = 0; i < A->num_rows; i++)
    {
        for (j = A->row_ptr[i]; j < A->row_ptr[i+1]; j++)
        {
            if (A->col_ind[j] < i && code == 1)
            {
                printf("This matrix is not triangular assuming nonsymmetrical, program will continue.\n");
                return 1;
            }
        }
    }

    printf("This matrix is upper triangular assuming this means it is symmetric, program will continue.\n");
    return 0;
}

// This function will recieve a triangular matrix that it assumes is symmetric and will fill in the missing data in the other triangle.

void fill_symmetry(CSRMatrix *A)
{
    // Allocate memory for temporary structure/arrays for sorting purposes, same as ReadMMtoCSR.

    SortData *temp = malloc(A->num_non_zeros * sizeof(SortData));
    int num_diagonal = 0;
    int i,j;

    if (temp == NULL)
    {
        printf("ERROR: memory allocation failed in fill_symmetry. Exiting the program.\n");
        exit(1);
    }

    // Find number of diagonal entries and get row data of each entry.

    for (i = 0; i < A->num_rows; i++)
    {
        for (j = A->row_ptr[i]; j < A->row_ptr[i+1]; j++)
        {
            temp[j].row_data = i;
            if (i == A->col_ind[j])
            {
                num_diagonal++;
            }
        }
    }

    // Reallocate memory to make room for the values which will be duplicated on the other side of the diagonal.

    int old_num_non_zeros = A->num_non_zeros;
    A->num_non_zeros = A->num_non_zeros * 2 - num_diagonal;

    A->col_ind = (int *)realloc(A->col_ind, A->num_non_zeros * sizeof(int));
    A->csr_data = (double *)realloc(A->csr_data, A->num_non_zeros * sizeof(double));
    temp = realloc(temp, A->num_non_zeros * sizeof(SortData));

    if (A->col_ind == NULL || A->csr_data == NULL || temp == NULL)
    {
        printf("ERROR: memory reallocation failed in fill_symmetry. Exiting the program.\n");
        exit(1);
    }

    // Copy the current data into the temporary arrays so that it can participate in the sorting later.

    for (i = 0; i < old_num_non_zeros; i++)
    {
        temp[i].col_data = A->col_ind[i];
        temp[i].value_data = A->csr_data[i];
    }

    // Add all the duplicated data to temp_col_ind, temp_csr_data, and row_data.

    int num_new = 0;

    for (i = 0; i < A->num_rows; i++)
    {
        for (j = A->row_ptr[i]; j < A->row_ptr[i+1]; j++)
        {
            if (i != A->col_ind[j])
            {
                temp[old_num_non_zeros + num_new].col_data = i;
                temp[old_num_non_zeros + num_new].value_data = A->csr_data[j];
                temp[old_num_non_zeros + num_new].row_data = A->col_ind[j];
                num_new++;
            }
        }
    }

    // Update the row_ptr to account for all of the new duplicated data.

    for (i = old_num_non_zeros; i < A->num_non_zeros; i++)
    {
        for (j = temp[i].row_data; j < A->num_rows; j++)
        {
            A->row_ptr[j+1]++;
        }
    }

    // Sort all of the col_ind and csr_data values based on row data to align with the row_ptr as was done in the ReadMMtoCSR function above.

    qsort(temp, A->num_non_zeros, sizeof(SortData), compare);

    // Writes this sorted data to the col_ind and csr_data arrays in the CSRMatrix structure.

    for (i = 0; i < A->num_non_zeros; i++)
    {
        A->col_ind[i] = temp[i].col_data;
        A->csr_data[i] = temp[i].value_data;
    }

    // Deallocate all of the temporary arrays.

    free(temp);
}

// Function which will use the Gauss-Seidel method to solve Ax=b for x.

void solver(const CSRMatrix *A, const double *b, double *x, int max_iteration, double stop_criteria)
{
    // Create another x vector to determine change between generations to see how close we are to convergence and the solution.
    
    double *tempx = (double *)malloc(A->num_rows * sizeof(double));
    int i,j,diagonal;
    double sum;

    if (tempx == NULL)
    {
        printf(" ERROR: memory allocation failed in solver. Exiting the program.\n");
        exit(1);
    }

    for (i = 0; i < A->num_rows; i++)
    {
        x[i] = 0.0;
        tempx[i] = 0.0;
    }

    int iteration = 0;
    double error = 1.0;

    // The iterative part of the method.

    while (iteration < max_iteration && error > stop_criteria)
    {

        // Going through the matrix row by row and calculating the next guesses for x based on the previous guess.

        for (i = 0; i < A->num_rows; i++)
        {
            sum = b[i];
            diagonal = -1;
            for (j = A->row_ptr[i]; j < A->row_ptr[i+1]; j++)
            {
                if (i != A->col_ind[j])
                {
                    // This is the difference between the Jacobi methond and the Gauss-Seidel method, as we use the most up to date value for x[i] rather than the value from the last iteration. This will allow the values to converge faster.
                    
                    sum -= A->csr_data[j] * tempx[A->col_ind[j]];
                }
                else {diagonal = j;}   
            }

            // Checks to see if a diagonal entry is missing as if so the method will not work.

            if (diagonal == -1)
            {
                printf("ERROR: there was a value missing from the diagonal. Solver failed. Exiting program.\n");
                exit(1);
            }

            // Determine the new x[i] value for that row.

            tempx[i] = sum / A->csr_data[diagonal];
        }

        // Checks the difference between this iteration and the last and moves to the next iteration.

        error = 0.0;

        for (i = 0; i < A->num_rows; i++)
        {
            error += fabs(x[i] - tempx[i]);
            x[i] = tempx[i];
        }

        iteration++;
    }

    // Deallocate the array as it is no longer needed.

    free(tempx);

    // Checks how we left the while loop to give the user more feedback.

    if (iteration == max_iteration)
    {
        printf("\nGauss-Seidel solver terminated by reaching the maximum iterations of %d.\n", iteration);
    }
    else 
    {
        printf("\nGauss-Seidel solver terminated by reaching the stop criteria of %.1e in %d iterations.\n", stop_criteria, iteration);
    }
}

// Function computes the value of the residual vector which is r = Ax - b.

void compute_residual(const CSRMatrix *A, const double *x, const double *b, double *r)
{

    spmv_csr(A, x, r);

    for (int i = 0; i < A->num_cols; i++)
    {
        r[i] -= b[i];
    }
}

// Function computes the norm/magnitude of a vector.

void compute_norm(const CSRMatrix *A, const double *r, double *norm)
{
    *norm = 0;
    for (int i = 0; i < A->num_cols; i++)
    {
        *norm += r[i]*r[i];
    }
    *norm = sqrt(*norm);
}

// This function was entirely made using ChatGPT and is used to visualize the sparsity of a matrix, specifically small matrices as it is more precise.

void visualize_matrix_small(CSRMatrix *matrix, const char *output_filename) {
    FILE *fp = fopen(output_filename, "wb");
    if (!fp) {
        fprintf(stderr, "Error: Could not open file for writing\n");
        return;
    }

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) {
        fprintf(stderr, "Error: Could not create PNG write structure\n");
        fclose(fp);
        return;
    }

    png_infop info = png_create_info_struct(png);
    if (!info) {
        fprintf(stderr, "Error: Could not create PNG info structure\n");
        png_destroy_write_struct(&png, NULL);
        fclose(fp);
        return;
    }

    if (setjmp(png_jmpbuf(png))) {
        fprintf(stderr, "Error during PNG creation\n");
        png_destroy_write_struct(&png, &info);
        fclose(fp);
        return;
    }

    png_init_io(png, fp);

    png_set_IHDR(
        png,
        info,
        matrix->num_cols, matrix->num_rows,
        8,  // Bit depth
        PNG_COLOR_TYPE_RGBA,  // Color type
        PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_DEFAULT,
        PNG_FILTER_TYPE_DEFAULT
    );

    png_bytep *row_pointers = (png_bytep *)malloc(matrix->num_rows * sizeof(png_bytep));
    for (int i = 0; i < matrix->num_rows; i++) {
        row_pointers[i] = (png_bytep)malloc(4 * matrix->num_cols * sizeof(png_byte));
    }

    for (int i = 0; i < matrix->num_rows; i++) {
        for (int j = 0; j < matrix->num_cols; j++) {
            int index = matrix->row_ptr[i];
            int end = matrix->row_ptr[i + 1];

            // Find the matrix entry in the current row
            int color = 0;  // Assume zero
            for (; index < end; index++) {
                if (matrix->col_ind[index] == j) {
                    color = 255;  // Non-zero entry
                    break;
                }
            }

            // Fill RGBA values
            row_pointers[i][4 * j] = color;        // Red
            row_pointers[i][4 * j + 1] = color;    // Green
            row_pointers[i][4 * j + 2] = color;    // Blue
            row_pointers[i][4 * j + 3] = 255;      // Alpha (fully opaque)
        }
    }

    png_set_rows(png, info, row_pointers);
    png_write_png(png, info, PNG_TRANSFORM_IDENTITY, NULL);

    for (int i = 0; i < matrix->num_rows; i++) {
        free(row_pointers[i]);
    }
    free(row_pointers);

    png_destroy_write_struct(&png, &info);
    fclose(fp);
}

// This function was entirely made using ChatGPT and is used to visualize the sparsity of a matrix, specifically large matrices as it is less precise.

void visualize_matrix_large(CSRMatrix *matrix, const char *output_filename) {
    const int downsample_factor = 10;  // Adjust this factor based on your needs

    FILE *fp = fopen(output_filename, "wb");
    if (!fp) {
        fprintf(stderr, "Error: Could not open file for writing\n");
        return;
    }

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) {
        fprintf(stderr, "Error: Could not create PNG write structure\n");
        fclose(fp);
        return;
    }

    png_infop info = png_create_info_struct(png);
    if (!info) {
        fprintf(stderr, "Error: Could not create PNG info structure\n");
        png_destroy_write_struct(&png, NULL);
        fclose(fp);
        return;
    }

    if (setjmp(png_jmpbuf(png))) {
        fprintf(stderr, "Error during PNG creation\n");
        png_destroy_write_struct(&png, &info);
        fclose(fp);
        return;
    }

    png_init_io(png, fp);

    // Downsample the image size
    int downsampled_width = matrix->num_cols / downsample_factor;
    int downsampled_height = matrix->num_rows / downsample_factor;

    png_set_IHDR(
        png,
        info,
        downsampled_width, downsampled_height,
        8,  // Bit depth
        PNG_COLOR_TYPE_RGBA,  // Color type
        PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_DEFAULT,
        PNG_FILTER_TYPE_DEFAULT
    );

    png_bytep *row_pointers = (png_bytep *)malloc(downsampled_height * sizeof(png_bytep));
    for (int i = 0; i < downsampled_height; i++) {
        row_pointers[i] = (png_bytep)malloc(4 * downsampled_width * sizeof(png_byte));
    }

    for (int i = 0; i < downsampled_height; i++) {
        for (int j = 0; j < downsampled_width; j++) {
            int original_row = i * downsample_factor;
            int original_col = j * downsample_factor;

            int index = matrix->row_ptr[original_row];
            int end = matrix->row_ptr[original_row + 1];

            int color = 0;
            for (; index < end; index++) {
                if (matrix->col_ind[index] == original_col) {
                    color = 255;  // Non-zero entry
                    break;
                }
            }

            row_pointers[i][4 * j] = color;        // Red
            row_pointers[i][4 * j + 1] = color;    // Green
            row_pointers[i][4 * j + 2] = color;    // Blue
            row_pointers[i][4 * j + 3] = 255;      // Alpha (fully opaque)
        }
    }

    png_set_rows(png, info, row_pointers);
    png_write_png(png, info, PNG_TRANSFORM_IDENTITY, NULL);

    for (int i = 0; i < downsampled_height; i++) {
        free(row_pointers[i]);
    }
    free(row_pointers);

    png_destroy_write_struct(&png, &info);
    fclose(fp);
}