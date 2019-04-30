#include "pch.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "751_project_v0.h"

using namespace std; 

/////////////////////////////////////////////////////////////////////////////////////////

int parse_size_configs(int *K)
{

	myprint_f(("	Parsing the size configs\n"));

	FILE *fp_input_config;
	fp_input_config = fopen(INPUT_CONFIG_FILE, "r");

	char name_of_variable[100];

	while (fscanf(fp_input_config, "%s", name_of_variable)==1)
	{
		if (strcmp(name_of_variable, "No_of_sub-domains_K") == 0)
		{
			fscanf(fp_input_config, " = %d", K);
		}
	}

	fclose(fp_input_config); 

	return 0; 
}


/////////////////////////////////////////////////////////////////////////////////////////

int parse_data_configs(int *size_of_subdomain, int K)
{
	myprint_f(("	Parsing the data configs\n"));

	FILE *fp_input_config;
	fp_input_config = fopen(INPUT_CONFIG_FILE, "r");

	char name_of_variable[100];

	while (fscanf(fp_input_config, "%s", name_of_variable)==1)
	{
		if (strcmp(name_of_variable, "Size_of_sub-domains") == 0)
		{
			fscanf(fp_input_config, " =");
			for (int i = 0; i < K + 1; i++)
			{
				fscanf(fp_input_config, " %d", &size_of_subdomain[i]);
			}
		}
	}

	fclose(fp_input_config);

	return 0; 
}

int schwarz_parse_data_configs(int **range_of_subdomain, int K)
{
	myprint_f(("	Parsing the data configs\n"));

	FILE *fp_input_config;
	fp_input_config = fopen(INPUT_CONFIG_FILE, "r");

	char name_of_variable[100];

	while (fscanf(fp_input_config, "%s", name_of_variable) == 1)
	{
		if (strcmp(name_of_variable, "Schwarz_Range_of_sub-domains") == 0)
		{
			fscanf(fp_input_config, " =");
			for (int i = 0; i < K; i++)
			{
				fscanf(fp_input_config, " %d %d", &range_of_subdomain[i][0], &range_of_subdomain[i][1]);
			}
		}
	}

	fclose(fp_input_config);

	return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////

int read_matrix(FILE *fp, double **mat, double m, double n)
{
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			fscanf(fp, "%lf", &mat[i][j]);
		}
	}
	return 0; 
}


/////////////////////////////////////////////////////////////////////////////////////////

int read_input_data(double ***A, double **J, double ***E, double ***F, double **v, double *z, int K, int *size_of_subdomain)
{
	myprint_f(("	Reading the data files for A's and J\n"));

	FILE* fp_a;
	fp_a = fopen(INPUT_DATA_A, "r");
	for (int i = 0; i < K; i++)
	{
		read_matrix(fp_a, A[i], size_of_subdomain[i], size_of_subdomain[i]);
	}
	fclose(fp_a);

	FILE* fp_j;
	fp_j = fopen(INPUT_DATA_J, "r");
	read_matrix(fp_j, J, size_of_subdomain[K], size_of_subdomain[K]);
	fclose(fp_j);

	myprint_f(("	Reading the data files for E's\n"));
	FILE* fp_e;
	fp_e = fopen(INPUT_DATA_E, "r");
	for (int i = 0; i < K; i++)
	{
		read_matrix(fp_e, E[i], size_of_subdomain[K], size_of_subdomain[i]);
	}
	fclose(fp_e);

	myprint_f(("	Reading the data files for F's\n"));
	FILE* fp_f;
	fp_f = fopen(INPUT_DATA_F, "r");
	for (int i = 0; i < K; i++)
	{
		read_matrix(fp_e, F[i], size_of_subdomain[i], size_of_subdomain[K]);
	}
	fclose(fp_f);

	myprint_f(("	Reading the data files for v's and z\n"));

	FILE* fp_vz;
	fp_vz = fopen(INPUT_DATA_VZ, "r");

	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < size_of_subdomain[i]; j++)
		{
			fscanf(fp_vz, "%lf", &v[i][j]);
		}
	}

	for (int i = 0; i < size_of_subdomain[K]; i++)
	{
		fscanf(fp_vz, "%lf", &z[i]);
	}

	fclose(fp_vz);

	return 0;
}

int schwarz_read_input_data(double **A, double *v, int K, int *size_of_subdomain, int size_of_system)
{
	myprint_f(("	Reading the data files for A's\n"));

	FILE* fp_a;
	fp_a = fopen(SCHWARZ_INPUT_DATA_A, "r");
		read_matrix(fp_a, A, size_of_system, size_of_system);
	fclose(fp_a);

	FILE* fp_v;
	fp_v = fopen(SCHWARZ_INPUT_DATA_V, "r");
	for (int i = 0; i < size_of_system; i++)
	{
		fscanf(fp_v, "%lf", &v[i]);
	}
	fclose(fp_v);

	return 0; 
}


/////////////////////////////////////////////////////////////////////////////////////////

int myprint_matrix(const char *mat_name, double **mat, int m, int n)
{
	if (strcmp(mat_name, "") != 0)
		myprint_d(("\n%s : \n", mat_name));
	for (int i = 0; i < m; i++)
	{
		myprint_d(("\nrow %d : \n", i));
		for (int j = 0; j < n; j++)
		{
			myprint_d(("%lf ", mat[i][j]));
		}
	}
	myprint_d(("\n"));

	return 0; 
}

int myprint_vect(const char *vect_name, double *vect, int m)
{
	if (strcmp(vect_name, "") !=0)
		myprint_d(("\n%s : \n", vect_name));
	for (int i = 0; i < m; i++)
	{
		myprint_d(("%lf ", vect[i]));
	}
	myprint_d(("\n")); 

	return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////

int myprint(double ***A, double **J, double ***E, double ***F, double **v, double *z, int K, int *size_of_subdomain)
{
	myprint_d(("Printing data: \n"));

	for (int i = 0; i < K; i++)
	{
		myprint_d(("\n\nA[%d] : \n", i));
		for (int j = 0; j < size_of_subdomain[i]; j++)
		{
			myprint_d(("\nrow = %d: ",j));
			for (int k = 0; k < size_of_subdomain[i]; k++)
			{
				myprint_d(("%lf ", A[i][j][k]));
			}
		}
	}

	myprint_d(("\n\nJ : \n"));
	for (int j = 0; j < size_of_subdomain[K]; j++)
	{
		myprint_d(("\nrow = %d: ", j));
		for (int k = 0; k < size_of_subdomain[K]; k++)
		{
			myprint_d(("%lf ", J[j][k]));
		}
	}

	for (int i = 0; i < K; i++)
	{
		myprint_d(("\n\nE[%d] : \n", i));
		for (int j = 0; j < size_of_subdomain[K]; j++)
		{
			myprint_d(("\nrow = %d: ", j));
			for (int k = 0; k < size_of_subdomain[i]; k++)
			{
				myprint_d(("%lf ", E[i][j][k]));
			}
		}
	}

	for (int i = 0; i < K; i++)
	{
		myprint_d(("\n\nF[%d] : \n", i));
		for (int j = 0; j < size_of_subdomain[i]; j++)
		{
			myprint_d(("\nrow = %d: ", j));
			for (int k = 0; k < size_of_subdomain[K]; k++)
			{
				myprint_d(("%lf ", F[i][j][k]));
			}
		}
	}

	for (int j = 0; j < K; j++)
	{
		myprint_d(("\n\nv[%d] : \n", j));
		for (int k = 0; k < size_of_subdomain[j]; k++)
		{
			myprint_d(("%lf ", v[j][k]));
		}
	}

	myprint_d(("\n\nz: \n"));
	for (int k = 0; k < size_of_subdomain[K]; k++)
	{
		myprint_d(("%lf ", z[k]));
	}
	myprint_d(("\n"));

	return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////

int fprint_matrix(const char *mat_name, FILE* fp, double **mat, int m, int n)
{
	if (strcmp(mat_name, "") != 0)
		fprintf(fp, "\n%s : \n", mat_name);
	for (int i = 0; i < m; i++)
	{
		fprintf(fp, "\nrow %d : \n", i);
		for (int j = 0; j < n; j++)
		{
			fprintf(fp, "%lf ", mat[i][j]);
		}
	}
	fprintf(fp, "\n");

	return 0;
}

int fprint_vect(const char *vect_name, FILE* fp, double *vect, int m)
{
	if (strcmp(vect_name, "") != 0)
		fprintf(fp, "\n%s : \n", vect_name);
	for (int i = 0; i < m; i++)
	{
		fprintf(fp, "%lf ", vect[i]);
	}
	fprintf(fp, "\n");

	return 0;
}

