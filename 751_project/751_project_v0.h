#pragma once

//// select the solver
#define SCHUR_COMPLEMENT
//#define ADDITIVE_SCHWARZ

#define MAX_ITRN_COUNT 1000

//// defines to set the print verbosity
//#define FLOW_PRINT   
//#define DEBUG_PRINT

#ifdef FLOW_PRINT
//#define myprint_f(x) printf("FLOW : "); printf x
#define myprint_f(x) printf x
#else 
#define myprint_f(x)  
#endif

#ifdef DEBUG_PRINT
//#define myprint_d(x) printf("DEBUG : "); printf x
#define myprint_d(x) printf x
#else 
#define myprint_d(x)  
#endif

//// file name defines 

#define INPUT_CONFIG_FILE "input_config.txt"
#define INPUT_DATA_A "input_data_a.txt"
#define INPUT_DATA_J "input_data_j.txt"
#define INPUT_DATA_E "input_data_e.txt"
#define INPUT_DATA_F "input_data_f.txt"
#define INPUT_DATA_VZ "input_data_vz.txt"

#define SCHWARZ_INPUT_DATA_A "schwarz_input_data_a.txt"
#define SCHWARZ_INPUT_DATA_V "schwarz_input_data_v.txt"


//// function forward declarations

// generic functions
int parse_size_configs(int *K);
int parse_data_configs(int *size_of_subdomain, int K);
int read_matrix(FILE *fp, double **mat, double m, double n); 

int myprint_vect(const char *vect_name, double *vect, int m);
int myprint_matrix(const char *mat_name, double **mat, int m, int n);

int fprint_vect(const char *vect_name, FILE* fp, double *vect, int m);
int fprint_matrix(const char *mat_name, FILE* fp, double **mat, int m, int n);

void addMatrix(double** matrix1, double** matrix2, double** matrix3, int m);
void subMatrix(double** matrix1, double** matrix2, double** matrix3, int m);
void multiplyMatrix(double** matrix1, double** matrix2, double** matrix3, int m, int q, int n);
void LUdecomposition_nonpivot(double** A, double** L, double** U, int m);
void LUdecomposition(double** A, double** L, double** U, int *P, int m);
void LU_solve(double** L, double** U, double* b, double* x, int *P, int m);
void LU_solve_nonpivot(double** L, double** U, double* b, double* x, int m);
void vector_multiply(double** A, double* x, double* b, int m, int q);
void vector_add(double* a, double* b, double*c, int m);
void vector_sub(double*a, double* b, double* c, int m);
void lin_solve(double** A, double*x, double* b, int m);
void Matrix_LUsolve_nonpivot(double** L, double** U, double**C, double** B, int m, int n);
void Matrix_LUsolve(double** L, double** U, double**C, double** B, int *P, int m, int n);
double l2_norm(double* u, int m);

// schur functions
int read_input_data(double ***A, double **J, double ***E, double ***F, double **v, double *z, int K, int *size_of_subdomain);
int myprint(double ***A, double **J, double ***E, double ***F, double **v, double *z, int K, int *size_of_subdomain); 

// schwarz functions
int schwarz_parse_data_configs(int **range_of_subdomain, int K);
int schwarz_read_input_data(double **A, double *v, int K, int *size_of_subdomain, int size_of_system);
