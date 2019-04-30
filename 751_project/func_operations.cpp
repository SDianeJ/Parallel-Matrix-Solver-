#include "pch.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include<vector>
#include "751_project_v0.h"

using namespace std;


void addMatrix(double** matrix1, double** matrix2, double** matrix3, int m)
{
	for (int row = 0; row < m; row++)
		for (int col = 0; col < m; col++)
			matrix3[row][col] = matrix1[row][col] + matrix2[row][col];

}


void subMatrix(double** matrix1, double** matrix2, double** matrix3, int m)
{
	for (int row = 0; row < m; row++)
		for (int col = 0; col < m; col++)
			matrix3[row][col] = matrix1[row][col] - matrix2[row][col];

}


void multiplyMatrix(double** matrix1, double** matrix2, double** matrix3, int m, int q, int n)
{
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
		{
			matrix3[i][j] = 0;
			for (int k = 0; k < q; ++k)
				matrix3[i][j] = matrix3[i][j] + (matrix1[i][k] * matrix2[k][j]);
		}
}



void LUdecomposition_nonpivot( double** A, double** L, double** U, int m)
{
	int i = 0, j = 0, k = 0;
	for (i = 0; i < m; i++) 
	{
		for (j = 0; j < m; j++) 
		{
			if (j < i)
			L[j][i] = 0;
			else
			{
				L[j][i] = A[j][i];
				for (k = 0; k < i; k++)
				{
					L[j][i] = L[j][i] - L[j][k] * U[k][i];
				}
			}
		}

		for (j = 0; j < m; j++) 
		{
			if (j < i)
				U[i][j] = 0;
			else if (j == i)
				U[i][j] = 1;
			else 
			{
				U[i][j] = A[i][j] / L[i][i];
				for (k = 0; k < i; k++)
				{
					U[i][j] = U[i][j] - ((L[i][k] * U[k][j]) / L[i][i]);
				}
			}
		}
	}
}



void LUdecomposition(double** A, double** L, double** U, int *P, int m)
{
	int i = 0, j = 0, k = 0;
	int max_j = 0; 
	double temp; 
	int temp_p; 

	double ** A_temp; 
	A_temp = new double *[m]; 
	for (i = 0; i < m; i++)
	{
		A_temp[i] = new double[m]; 
	}

	// row pivoting
	for ( i = 0; i < m; i++)
	{
		P[i] = i;
		for (j = 0; j < m; j++)
			A_temp[i][j] = A[i][j]; 
	}
	for ( i = 0; i < m; i++)
	{
		max_j = i;
		for (j = i; j < m; j++)
		{
			if (A_temp[j][i] > A_temp[max_j][i])
				max_j = j; 
		}
		if (max_j != i)
		{
			temp_p = P[i]; 
			P[i] = P[max_j]; 
			P[max_j] = temp_p; 
			for (k = 0; k < m; k++)
			{
				temp = A_temp[i][k];
				A_temp[i][k] = A_temp[max_j][k];
				A_temp[max_j][k] = temp;
			}
		}
	}

	for (i = 0; i < m; i++)
	{
		for (j = 0; j < m; j++)
		{
			if (j < i)
				L[j][i] = 0;
			else
			{
				L[j][i] = A_temp[j][i];
				for (k = 0; k < i; k++)
				{
					L[j][i] = L[j][i] - L[j][k] * U[k][i];
				}
			}
		}

		for (j = 0; j < m; j++)
		{
			if (j < i)
				U[i][j] = 0;
			else if (j == i)
				U[i][j] = 1;
			else
			{
				U[i][j] = A_temp[i][j] / L[i][i];
				for (k = 0; k < i; k++)
				{
					U[i][j] = U[i][j] - ((L[i][k] * U[k][j]) / L[i][i]);
				}
			}
		}
	}	
	for (i = 0; i < m; i++)
	{
		delete[]A_temp[i];
	}
	delete[]A_temp;

}


void LU_solve_nonpivot(double** L, double** U, double* b, double* x, int m) //Ax=b
{
	int i, j, k; 
	float temp;
	//double** v;
	//v = new double*[m];

	vector<double> v(m);
	for (i = 0; i < m; ++i)
		v[i] = 0;
	for (i = 0; i < m; i++)
	{
		temp = 0;
	for (j = 0; j < i; j++)
	{
		temp = temp + L[i][j] * v[j];
	}
	 v[i] = (b[i] - temp) / L[i][i];
	 /*cout << "This is L" << L[i][i]<<"\n";
	 cout << "This is b" << b[i]<<"\n";
	 cout << v[i];*/

	}

	for ( i = m - 1; i >= 0; i--)
	{
		temp = 0; x[i] = 0;
		for ( k = i + 1; k < m; k++)
		{
			temp = temp + U[i][k] * x[k];
		}
			x[i] = (v[i] - temp) / U[i][i];
		
			/*cout << "This is v" << v[i] << "\n";
			cout << "This is U" << U[i][i] << "\n";
			cout << "X " << x[i] << "\n";*/

	}

}


void LU_solve(double** L, double** U, double* b, double* x, int *P, int m) //Ax=b
{
	int i, j, k;
	float temp;
	double* b_temp;
	b_temp = new double[m];

	for (i = 0; i < m; i++)
	{
		b_temp[i] = b[P[i]]; 
	}

	vector<double> v(m);
	for (i = 0; i < m; ++i)
		v[i] = 0;
	for (i = 0; i < m; i++)
	{
		temp = 0;
		for (j = 0; j < i; j++)
		{
			temp = temp + L[i][j] * v[j];
		}
		v[i] = (b_temp[i] - temp) / L[i][i];
		/*cout << "This is L" << L[i][i]<<"\n";
		cout << "This is b" << b[i]<<"\n";
		cout << v[i];*/

	}

	for (i = m - 1; i >= 0; i--)
	{
		temp = 0; x[i] = 0;
		for (k = i + 1; k < m; k++)
		{
			temp = temp + U[i][k] * x[k];
		}
		x[i] = (v[i] - temp) / U[i][i];

		/*cout << "This is v" << v[i] << "\n";
		cout << "This is U" << U[i][i] << "\n";
		cout << "X " << x[i] << "\n";*/

	}

}


void lin_solve(double** A, double*b, double* x, int m)
{
	double **L; double** U;
	int *P; 
		L = new double* [m];
		U = new double* [m];
		P = new int[m]; 
		for (int i = 0; i < m; i++)
		{
			L[i] = new double [m];
			U[i] = new double [m];
		}

		LUdecomposition(A, L, U, P,  m);
		LU_solve(L, U, b, x, P, m);

}


void vector_add(double* a, double* b, double*c, int m)
{
	for (int i = 0; i < m; i++)
		c[i] = a[i] + b[i];
}

void vector_sub(double* a, double* b, double*c, int m)
{
	for (int i = 0; i < m; i++)
		c[i] = a[i] - b[i];
}


void vector_multiply(double** A, double* b, double*c, int m, int q)
{
	for (int i = 0; i < m; ++i)
	{
		c[i] = 0;
		for (int j = 0; j < q; ++j)
		{
				c[i] = c[i] + (A[i][j] * b[j]);
		}
	}
}

void Matrix_LUsolve(double**L, double** U, double** B, double** C, int *P, int m, int n) //AB=C
{
	int i, j, k;
	double temp;

	/*double** V;
	V = new double*[m];
	for (i = 0; i < m; i++)
		V[i] = new double[n];
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			V[i][j] = 0;
*/
	double* b; double* c;
	b = new double[m];
	c = new double[m];
	//myprint_matrix("B", B, m, n);
	for (j = 0; j < n; j++)
	{
			for (i = 0; i< m; i++)
			{
				b[i] = B[i][j];
			}

			//myprint_vect("B", b, m);
			LU_solve(L, U, b, c, P, m);
			
			//myprint_vect("C", c, m);
			for (i = 0; i < m; i++)
				{
					C[i][j] = c[i];
				}


	}


/*
	for (j = 0; j < n; ++j)
		for (i = 0; i < m; i++)
		{
			V[i][j] = 0;
			temp = 0;
			for (k = 0; k < i; k++)
			{
				temp = temp + L[i][k] * V[k][i];
				cout << "\nVki  " << V[k][i];
			} cout << "\n Temp" << temp;
			V[i][j] = (C[i][j] - temp) / L[i][i];
			cout << "\n" << C[i][j] << "    " << L[i][i];
			cout << "\n Vij " << V[i][j];
		}
*/

	/*or (j=0; j<n;j++)
		for (i = m - 1; i >= 0; i--)
		{
			temp = 0; B[i][j] = 0;
			for (k = i + 1; k < m; k++)
			{
				temp = temp + U[i][k] * B[k][i];
			}
			B[i][j] = (V[i][j] - temp) / U[i][i];
		}
*/
		/*cout << "This is v" << v[i] << "\n";
		cout << "This is U" << U[i][i] << "\n";
		cout << "X " << x[i] << "\n";*/

		


}



double l2_norm(double* u, int m) 
{
	double accum = 0 ;
	for (int i = 0; i <m ; ++i) 
	{
		accum += u[i] * u[i];
	}
	return sqrt(accum);
}


