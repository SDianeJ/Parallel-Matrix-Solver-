#include "pch.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


#include "751_project_v0.h"

using namespace std;

#ifdef ADDITIVE_SCHWARZ
int main()
{
	myprint_f(("Starting the serial Schwarz solver... \n"));


	/////////////////////////////////////////////////////////////////////////////////////////

	myprint_f(("Parsing the input config file\n"));

	int K;
	int **range_of_subdomain;
	int  *size_of_subdomain;
	int	 size_of_system;

	parse_size_configs(&K);

	myprint_d(("K = %d\n", K));

	range_of_subdomain = new int*[K];
	for (int i = 0; i < K; i++)
		range_of_subdomain[i] = new int[2];
	size_of_subdomain = new int[K];

	schwarz_parse_data_configs(range_of_subdomain, K);

	for (int i = 0; i < K; i++)
	{
		size_of_subdomain[i] = range_of_subdomain[i][1] - range_of_subdomain[i][0] + 1;
	}

	size_of_system = range_of_subdomain[K - 1][1] + 1;

	myprint_d(("range_of_subdomain = "));
	for (int i = 0; i < K; i++)
	{
		myprint_d((" %d %d; ", range_of_subdomain[i][0], range_of_subdomain[i][1]));
		printf(" %d", size_of_subdomain[i]);
	}
	myprint_d(("\n"));
	printf("\n");


	/////////////////////////////////////////////////////////////////////////////////////////

	// declare all the variables : 

	double **A, **RAi, **Ai;
	double  *x;
	double  *v;

	double ***R, ***RT;
	double ***LA, ***UA; // L and U for A's
	int **P; 
	double *x_curr, *x_prev, *x_diff;
	double *Ax_curr, *e_curr;  // ej = v - Ax_curr
	double **ei_curr, **Aiei_curr; // ei_curr = Ri * e_curr
	double **xei, *sxei;  // xei = RTi * Aiei_curr

	double convergence_condn;
	int itrn_count;
	clock_t start_t, end_t, total_t ;


	/////////////////////////////////////////////////////////////////////////////////////////

	// allocate memory for all the variables


	A = new double *[size_of_system];
	x = new double[size_of_system];
	v = new double[size_of_system];
	R = new double **[K];
	RT = new double **[K];
	LA = new double **[K];
	UA = new double **[K];
	P = new int *[K]; 
	x_curr = new double[size_of_system];
	x_prev = new double[size_of_system];
	x_diff = new double[size_of_system];
	Ax_curr = new double[size_of_system];
	e_curr = new double[size_of_system];
	ei_curr = new double *[K];
	Aiei_curr = new double *[K];
	xei = new double *[K];
	sxei = new double[size_of_system];

	for (int i = 0; i < K; i++)
	{
		R[i] = new double*[size_of_subdomain[i]];
		RT[i] = new double*[size_of_system];
		LA[i] = new double*[size_of_subdomain[i]];
		UA[i] = new double*[size_of_subdomain[i]];

		for (int j = 0; j < size_of_subdomain[i]; j++)
		{
			R[i][j] = new double[size_of_system];
			LA[i][j] = new double[size_of_subdomain[i]];
			UA[i][j] = new double[size_of_subdomain[i]];
		}

		for (int j = 0; j < size_of_system; j++)
		{
			RT[i][j] = new double[size_of_subdomain[i]];
		}

		P[i] = new int[size_of_subdomain[i]]; 

		ei_curr[i] = new double[size_of_subdomain[i]];
		Aiei_curr[i] = new double[size_of_subdomain[i]];
		xei[i] = new double[size_of_system];
	}

	for (int i = 0; i < size_of_system; i++)
	{
		A[i] = new double[size_of_system];
	}


	/////////////////////////////////////////////////////////////////////////////////////////

	myprint_f(("Loading the matrix data from the input data files\n"));

	schwarz_read_input_data(A, v, K, size_of_subdomain, size_of_system);
	myprint_matrix("A", A, size_of_system, size_of_system); 
	myprint_vect("v", v, size_of_system); 


	///////////////////////////////////////////////////////////////////////////////////////////

	start_t = clock(); 

	myprint_f(("configuring R and RT matrices\n"));

	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < size_of_subdomain[i]; j++)
		{
			for (int k = 0; k < size_of_system; k++)
			{
				R[i][j][k] = (k == (range_of_subdomain[i][0] + j));
			}
		}	myprint_matrix("R", R[i], size_of_subdomain[i], size_of_system);
	}


	// R - R[o], R[1], each R[i] = size_of_subdomain[i] x size_of_system  4x5, 3x5

	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < size_of_system; j++)
		{
			for (int k = 0; k < size_of_subdomain[i]; k++)
			{
				RT[i][j][k] = R[i][k][j];
			}
		}
	}
	// RT - RT[o], RT[1], each RT[i] = size_of_system x size_of_subdomain[i]  

	myprint_f(("LU factorize Ai matrices\n"));
	for (int i = 0; i < K; i++)
	{
		RAi = new double *[size_of_subdomain[i]];   // R * A  , 4x5, 3x5
		Ai = new double *[size_of_subdomain[i]];    // Ai[0], Ai[1] 4x4,3x3
		for (int j = 0; j < size_of_subdomain[i]; j++)
		{
			RAi[j] = new double[size_of_system];
			Ai[j] = new double[size_of_subdomain[i]];
		}

		multiplyMatrix(R[i], A, RAi, size_of_subdomain[i], size_of_system, size_of_system);
		multiplyMatrix(RAi, RT[i], Ai, size_of_subdomain[i], size_of_system, size_of_subdomain[i]);
		//LUdecomposition(Ai, LA[i], UA[i], size_of_subdomain[i]);   /// LA[0/1],  UA[0/1]  - 4x4,3x3 
		LUdecomposition(Ai, LA[i], UA[i], P[i], size_of_subdomain[i]); 
		//myprint_matrix("LA1", LA[1], size_of_subdomain[1], size_of_subdomain[1]);
		//myprint_matrix("UA1", UA[1], size_of_subdomain[1], size_of_subdomain[1]);
		//multiplyMatrix(LA[1], UA[1], RAi, size_of_subdomain[0], size_of_subdomain[0], size_of_subdomain[0]);
		//myprint_matrix("RAi", RAi, size_of_subdomain[0], size_of_subdomain[0]);
		//myprint_d(("\nP = ")); 
		//for (int j = 0; j < size_of_subdomain[i]; j++)
		//	myprint_d(("%d ", P[1][j])); 
		//myprint_d(("\n")); 

		for (int j = 0; j < size_of_subdomain[i]; j++)
		{
			delete[]RAi[j];
			delete[]Ai[j];
		}
		delete[]RAi;
		delete[]Ai;
	}

	myprint_f(("Guess an initial value for x\n"));


	//x_curr - 5
	for (int i = 0; i < size_of_system; i++)
	{
		x_curr[i] = 0;
		x_prev[i] = 0;
	}



	/////////////////////////////////////////////////////////////////////////////////////////////

	myprint_f(("Starting the iterations... \n"));

	FILE *fp_out; 
	fp_out = fopen("out.txt", "w"); 

	itrn_count = 0; convergence_condn = 0;
	while ((!convergence_condn) && (itrn_count < MAX_ITRN_COUNT))
	{
		myprint_f(("Iteration no : %d\n", itrn_count));

		vector_multiply(A, x_curr, Ax_curr, size_of_system, size_of_system);    // Ax_curr - 5

		vector_sub(v, Ax_curr, e_curr, size_of_system);   // e_curr - 5

		for (int i = 0; i < K; i++)
		{
			vector_multiply(R[i], e_curr, ei_curr[i], size_of_subdomain[i], size_of_system);   // ei_curr[0/1]  4, 3
			//myprint_vect("ei curri", ei_curr[i], size_of_subdomain[i]);
		}
		for (int i = 0; i < K; i++)
		{
			//LU_solve(LA[i], UA[i], ei_curr[i], Aiei_curr[i], size_of_subdomain[i]);    // Aiei_curr[0/1] - 4,3
			LU_solve(LA[i], UA[i], ei_curr[i], Aiei_curr[i], P[i], size_of_subdomain[i]);
			//myprint_vect("Aiei_curr[i]", Aiei_curr[i], size_of_subdomain[i]);
		}

		for (int i = 0; i < K; i++)
		{
			//myprint_matrix("RT", RT[i], size_of_system, size_of_subdomain[i]);
			vector_multiply(RT[i], Aiei_curr[i], xei[i], size_of_system, size_of_subdomain[i]);  //xei[0/1] 5, 5 
			//myprint_vect("xei[i]", xei[i], size_of_system);
		}

		for (int i = 0; i < K; i++)
		{
			for (int j = 0; j < size_of_system; j++)
			{
				if (i < (K - 1))
				{
					xei[i][j] = ((j <= range_of_subdomain[i][1]) && (j >= range_of_subdomain[i + 1][0])) ? xei[i][j] / 2 : xei[i][j]; 
				}
				if (i > 0)
				{
					xei[i][j] = ((j >= range_of_subdomain[i][0]) && (j <= range_of_subdomain[i - 1][1])) ? xei[i][j] / 2 : xei[i][j];
				}
			}
		}
		
		for (int i = 0; i < size_of_system; i++)
		{
			sxei[i] = xei[0][i];  
		}  //sxei   -   5

		for (int i = 1; i < K; i++)
		{
			vector_add(xei[i], sxei, sxei, size_of_system);
		}    //sxei   -   5
		//myprint_vect("sxei", sxei, size_of_system);

		//// damping factor
		//for (int i = 0; i < size_of_system; i++)
		//{
		//	sxei[i] = 2 * sxei[i];
		//}

		vector_add(x_prev, sxei, x_curr, size_of_system);   	//x_curr - 5


		myprint_vect("Xcurr_next", x_curr, size_of_system);
		fprint_vect("Xcurr_next", fp_out, x_curr, size_of_system);
		//convergence_condn = (l2_norm(sxei, size_of_system) / l2_norm(x_prev, size_of_system))*100;   
		convergence_condn = (l2_norm(e_curr, size_of_system) / l2_norm(v, size_of_system)) * 100;
		cout << "\nThis is......................................" << convergence_condn;
		if (convergence_condn < 0.1)
			convergence_condn = 1;
		else
			convergence_condn = 0; 
		cout << "\nThis is......................................"<<convergence_condn; 

		//cout << "\nThis is......................................" << itrn_count;

		itrn_count++;
		for (int i = 0; i < size_of_system; i++)
		{
			x_prev[i] = x_curr[i];
		}

	}

	
	for (int i=0; i<size_of_system; i++)
	{
			printf("%f  ", x_curr[i]); 
	}
	total_t = (double) (end_t - start_t) / CLOCKS_PER_SEC ; 
	printf("Total cycles taken by CPU: %lf\n",end_t - start_t  );
	printf("Total time taken by CPU: %lf\n", total_t  );
   	printf("Exiting of the program...\n");


	return 0;
}
#endif
