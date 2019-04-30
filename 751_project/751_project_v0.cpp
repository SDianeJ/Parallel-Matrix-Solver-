// 751_project_v0.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


#include "751_project_v0.h"

using namespace std;

#ifdef SCHUR_COMPLEMENT 
int main()
{
	printf("Starting the serial Schurs Complement solver... \n") ; 


	/////////////////////////////////////////////////////////////////////////////////////////
	
	myprint_f(("Parsing the input config file\n")) ;

	int K; 
	int *size_of_subdomain;

	parse_size_configs(&K); 
	
	myprint_d(("K = %d\n", K)); 

	size_of_subdomain = new int[K + 1];

	parse_data_configs(size_of_subdomain, K);

	myprint_d(("size_of_subdomain = "));

	printf("\n");
	printf("\n===========================================================\n");
	printf("no of subdomains K = %d\n", K);
	printf("sizes of the subdomains : \n");
	for (int i = 0; i < K ; i++)
	{
		myprint_d((" %d", size_of_subdomain[i]));
		printf(" %d", size_of_subdomain[i]);
	}
	myprint_d(("\n"));
	printf("\n"); 
	printf("Size of the interconnect = %d\n", size_of_subdomain[K]); 


	/////////////////////////////////////////////////////////////////////////////////////////

	// declare all the variables : 

	double ***A,  **J ;
	double ***E, ***F ; 
	double  **x,   *y ; 
	double  **v,   *z ; 

	double ***LA, ***UA ; // L and U for A's
	int **P;
	double **Aiv, ***AiF ; // v/A  and f/A
	double **EAiv, ***EAiF; // E(v/A) and E(F/A)
	double *sEAiv, **sEAiF; // summations
	double  *b_y, **A_y ;   // b_y = z - sEAiv, A_y = J - sEAiF 

	double **Fiy; 
	double **bi_x;  //  = vi - Fiy

	clock_t start_t, end_t, total_t ;

	/////////////////////////////////////////////////////////////////////////////////////////

	// allocate memory for all the variables

	A = new double **[K];
	E = new double **[K];
	F = new double **[K];
	x = new double *[K];
	v = new double *[K];

	for (int i = 0; i < K; i++)
	{
		A[i] = new double*[size_of_subdomain[i]];
		F[i] = new double*[size_of_subdomain[i]];
		for (int j = 0; j < size_of_subdomain[i]; j++)
		{
			A[i][j] = new double[size_of_subdomain[i]];
			F[i][j] = new double[size_of_subdomain[K]];
		}

		E[i] = new double*[size_of_subdomain[K]];
		for (int j = 0; j < size_of_subdomain[K]; j++)
		{
			E[i][j] = new double[size_of_subdomain[i]];
		}

		x[i] = new double[size_of_subdomain[i]];
		v[i] = new double[size_of_subdomain[i]];
	}	

	J = new double *[size_of_subdomain[K]];
	y = new double[size_of_subdomain[K]]; 
	z = new double[size_of_subdomain[K]]; 
	
	for (int i = 0; i < size_of_subdomain[K]; i++)
	{
		J[i] = new double[size_of_subdomain[K]];
	}

	LA = new double **[K];
	UA = new double **[K];
	P = new int *[K];
	Aiv = new double *[K];
	AiF = new double **[K];
	EAiv = new double *[K];
	EAiF = new double **[K];

	for (int i = 0; i < K; i++)
	{
		LA[i] = new double*[size_of_subdomain[i]];
		UA[i] = new double*[size_of_subdomain[i]];
		P[i] = new int[size_of_subdomain[i]];
		AiF[i] = new double*[size_of_subdomain[i]];
		EAiF[i] = new double*[size_of_subdomain[K]];
		for (int j = 0; j < size_of_subdomain[i]; j++)
		{
			LA[i][j] = new double[size_of_subdomain[i]];
			UA[i][j] = new double[size_of_subdomain[i]];
			AiF[i][j] = new double[size_of_subdomain[K]];
		}

		for (int j = 0; j < size_of_subdomain[K]; j++)
		{
			EAiF[i][j] = new double[size_of_subdomain[K]];
		}

		Aiv[i] = new double[size_of_subdomain[i]];
		EAiv[i] = new double[size_of_subdomain[K]];
	}

	sEAiv = new double[size_of_subdomain[K]];
	sEAiF = new double *[size_of_subdomain[K]];
	b_y = new double[size_of_subdomain[K]];
	A_y = new double *[size_of_subdomain[K]];

	Fiy = new double *[K];
	bi_x = new double *[K];

	for (int i = 0; i < size_of_subdomain[K]; i++)
	{
		sEAiF[i] = new double[size_of_subdomain[K]];
		A_y[i] = new double[size_of_subdomain[K]];
	}

	for (int i = 0; i < K; i++)
	{
		Fiy[i] = new double[size_of_subdomain[i]];
		bi_x[i] = new double[size_of_subdomain[i]];
	}


	/////////////////////////////////////////////////////////////////////////////////////////

	myprint_f(("Loading the matrix data from the input data files\n"));

	read_input_data(A, J, E, F, v, z, K, size_of_subdomain); 
	myprint(A, J, E, F, v, z, K, size_of_subdomain); 


	/////////////////////////////////////////////////////////////////////////////////////////


	myprint_f(("Computing y\n")); 

	start_t = clock(); 

	for (int i = 0; i < K; i++)
	{
		LUdecomposition(A[i], LA[i], UA[i], P[i], size_of_subdomain[i]); 
	}

	//myprint_matrix("LA[1]", LA[1], size_of_subdomain[1], size_of_subdomain[1]);
	//myprint_matrix("UA[1]", UA[1], size_of_subdomain[1], size_of_subdomain[1]);

	for (int i = 0; i < K; i++)
	{
		LU_solve(LA[i], UA[i], v[i], Aiv[i], P[i], size_of_subdomain[i]); 
		myprint_vect("Aiv[i]", Aiv[i], size_of_subdomain[i]);

	}

	myprint_vect("Aiv[1]", Aiv[1], size_of_subdomain[1]);


	for (int i = 0; i < K; i++)
	{
		Matrix_LUsolve(LA[i], UA[i], F[i], AiF[i], P[i], size_of_subdomain[i], size_of_subdomain[K]);
		myprint_matrix("AiF[i]", AiF[i], size_of_subdomain[i], size_of_subdomain[K]);
	}

	for (int i = 0; i < K; i++)
	{
		vector_multiply(E[i], Aiv[i], EAiv[i], size_of_subdomain[K], size_of_subdomain[i]);
		myprint_vect("EAiv[i]", EAiv[i], size_of_subdomain[K]);
	}

	for (int i = 0; i < K; i++)
	{
		multiplyMatrix(E[i], AiF[i], EAiF[i], size_of_subdomain[K], size_of_subdomain[i], size_of_subdomain[K]);
		myprint_matrix("EAIFi", EAiF[i], size_of_subdomain[K], size_of_subdomain[K]);
	}

	for (int i = 0; i < size_of_subdomain[K]; i++)
	{
		sEAiv[i] = EAiv[0][i];
		for (int j = 0; j < size_of_subdomain[K]; j++)
		{
			sEAiF[i][j] = EAiF[0][i][j]; 
		}
	}
	for (int i = 1; i < K; i++)
	{
		vector_add(EAiv[i], sEAiv, sEAiv, size_of_subdomain[K]);
	}

	myprint_vect("sEAIv1", sEAiv, size_of_subdomain[K]);
	for (int i = 1; i < K; i++)
	{
		addMatrix(EAiF[i], sEAiF, sEAiF, size_of_subdomain[K]);
	}
	myprint_matrix("SEAIF1", sEAiF, size_of_subdomain[K], size_of_subdomain[K]);
	vector_sub(z, sEAiv, b_y, size_of_subdomain[K]);
	subMatrix(J, sEAiF, A_y, size_of_subdomain[K]);

	myprint_vect("b_y", b_y, size_of_subdomain[K]);
	myprint_matrix("A_y", A_y, size_of_subdomain[K], size_of_subdomain[K]);

	lin_solve(A_y, b_y, y, size_of_subdomain[K]); 


	myprint_vect("y", y, size_of_subdomain[K]); 


	/////////////////////////////////////////////////////////////////////////////////////////

	myprint_f(("solving for the xi's")); 

	for (int i = 0; i < K; i++)
	{
		vector_multiply(F[i], y, Fiy[i], size_of_subdomain[i], size_of_subdomain[K]);
	}

	for (int i = 0; i < K; i++)
	{
		vector_sub(v[i], Fiy[i], bi_x[i], size_of_subdomain[i]);
	}

	for (int i = 0; i < K; i++)
	{
		LU_solve(LA[i], UA[i], bi_x[i], x[i], P[i], size_of_subdomain[i]);
		myprint_d(("\nx[%d] = \n", i)); 
		myprint_vect("", x[i], size_of_subdomain[i]);
	}


	end_t = clock(); 

	printf("\n===========================================================\n");
	printf("Solution : \n");
	printf("x = \n");
	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < size_of_subdomain[i]; j++)
		{
			printf("%f  ", x[i][j]);
		}
		printf("\n");
	}
	for (int j = 0; j < size_of_subdomain[K]; j++)
	{
		printf("%f  ", y[j]);
	}
	printf("\n");
	printf("\n===========================================================\n");

	total_t = (double) (end_t - start_t) / CLOCKS_PER_SEC ; 
	printf("Total cycles taken by CPU: %f\n",end_t - start_t  );
	printf("Total time taken by CPU: %f\n", total_t  );
	printf("\n");
	printf("\n===========================================================\n");
   	printf("Exiting of the program...\n");

	return 0; 
}
#endif
