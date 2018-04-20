// Paulina Szreder 147418
// Uk³ady róznañ liniowych 
// NUM - 918
// A - macierz o rozmiarze NUMxNUM
// b - wektor d³ugoœci NUMxNUM
// x - wektor rozwi¹zania 

#include <iostream>
#include <cmath>
#include <time.h>
#include "Windows.h"

#define NUM 918	
#define ITER 1000
using namespace std;

#pragma region CRT

void crtMatrix(double **A, int e)
{
	for (int i = 0; i < NUM; i++)
	for (int j = 0; j < NUM; j++)
	{
		A[i][j] = 0;

		if (j == i)
			A[i][j] = e;		
		else if ((i == (j - 1) || j == (i - 1)) && j >= 0)
			A[i][j] = -1;
		else if ((i == (j - 2) || j == (i - 2)) && i >= 0)
			A[i][j] = -1;
	}



}

void crtVector_b(double *b)
{
	for (int n = 0; n < NUM; n++)
	{
		double tmp = n * 9;
		tmp = tmp / 50;
		b[n] = (double)sin(tmp); // sin(n * (f + 1)/50)
	}
		
}

void crtVector_N(double *N, double **A)
{
	for (int i = 0; i < NUM; i++)
		N[i] = (double)pow(A[i][i], -1);
}

void crtMatrix_M(double **M, double **A, double *N)
{
	for (int i = 0; i < NUM; i++)
	for (int j = 0; j < NUM; j++)
	{
		if (i == j)
			M[i][j] = 0;
		else
			M[i][j] = (-A[i][j] * N[i]);
	}
}

void  AintoLDU(double **A, int **L, int **D, int **U)
{
	for (int i = 0; i < NUM; i++)
	for (int j = 0; j < NUM; j++)
	{
		if (i < j)
			U[i][j] = A[i][j];
		else if (i > j)
			L[i][j] = A[i][j];
		else
			D[i][j] = A[i][j];
	}

}

#pragma endregion

void mnozenieMacierzy(double **A, double *b, double *x, double *R)
{
	for (int i = 0; i < NUM; i++)
	{
		for (int j = 0; j < NUM; j++)
			R[i] += A[j][i] * x[i];
			R[i] -= b[i];

	}
}

int Jacobi(double *x1, double *x2, double **A, double *b, double *N, double **M)
{
	double res = pow(10, -20);	
	double *R = (double*)malloc(sizeof(double)*NUM);
	int k;	

	for (k = 0; k < ITER; k++)
	{
		double suma = 0;
		for (int i = 0; i<NUM; i++) 
		{
			for (int j = 0; j < NUM; j++)
				R[j] = 0;

			x2[i] = N[i] * b[i];
			for (int j = 0; j < NUM; j++)
				x2[i] += M[i][j] * x1[j];

			//mnozenieMacierzy(A, b, x2, R);

			R[i] = x1[i] - x2[i];

			
		}

		for (int i = 0; i < NUM; i++)
			suma += abs(R[i]);
		

		
		for (int i = 0; i < NUM; i++)
			x1[i] = x2[i];
		//for (int i = 0; i < NUM; i++)
		//	suma += R[i];

		if (suma < res)
			return k;

	}
	return k;
	
}

int GaussaSeidla(double *x, double **A, double *b, double *N, int **L, int **U)
{
	double res = pow(10, -20);
	double *R = (double*)malloc(sizeof(double)*NUM);
	double *x1 = (double*)malloc(sizeof(double)*NUM);
	int k;

	for (int i = 0; i < NUM; i++)
		x1[i] = 0;

	double suma;

	for (k = 0; k < ITER; k++)
	{
		suma = 0;

		for (int i = 0; i < NUM; i++)
		{
			x[i] = b[i] * N[i];
			for (int j = 0; j < i; j++)
				x[i] -= L[i][j] * N[j] * x[j];
			for (int j = i + 1; j < NUM; j++)
				x[i] -= U[i][j] * N[j] * x[j];

			R[i] = x[i] - x1[i];
		}

		for (int i = 0; i < NUM; i++)
			x1[i] = x[i];

		for (int i = 0; i < NUM; i++)
			suma += abs(R[i]);

		if (suma < res)
			return k;
			
	}
	return k;

}	

double Gaussa(double **A, double *b, double *x)
{
	double **tmpA;
	double res = pow(10, -20);
	double suma = 0;


	double *R = (double*)malloc(sizeof(double)*NUM);

	for (int i = 0; i < NUM; i++)
		R[i] = 0;

	tmpA = (double**)malloc(sizeof(double)*NUM);
	for (int i = 0; i < NUM; i++)
		tmpA[i] = (double*)malloc(sizeof(double)*NUM + 1);

	for (int i = 0; i < NUM; i++)
	{
		for (int j = 0; j < NUM; j++)
		{
			tmpA[i][j] = A[i][j];
		}
		tmpA[i][NUM] = b[i];
	}

	double tmp = 0;

	for (int k = 0; k < NUM - 1; k++)
	{
		for (int i = k + 1; i < NUM; i++)
		{
			tmp = tmpA[i][k] / tmpA[k][k];
			for (int j = k; j < NUM + 1; j++)
			{
				tmpA[i][j] -= tmp * tmpA[k][j];
			}
		}
	}

	for (int k = NUM - 1; k >= 0; k--)
	{
		tmp = 0;
		for (int j = k + 1; j < NUM; j++)
		{
			tmp += tmpA[k][j] * x[j];
		}
		x[k] = (tmpA[k][NUM] - tmp) / tmpA[k][k];


	}
	
	mnozenieMacierzy(A, b, x, R);

	for (int i = 0; i < NUM; i++)
		suma += R[i];

	return suma;
}


int main()
{
	cout.precision(20);

#pragma region ZMIENNE


	double **A;
	double **C;
	double *b;
	double *N; // N = D^-1
	double **M;

	double *x1;
	double *x2;
	double *x3;
	double *x;

	int **L;
	int **D;
	int **U;

	int liczbaIteracji;

	long long start = GetTickCount();


	A = (double**)malloc(sizeof(double)*NUM);
	for (int i = 0; i < NUM; i++)
		A[i] = (double*)malloc(sizeof(double)*NUM);

	C = (double**)malloc(sizeof(double)*NUM);
	for (int i = 0; i < NUM; i++)
		C[i] = (double*)malloc(sizeof(double)*NUM);

	M = (double**)malloc(sizeof(double)*NUM);
	for (int i = 0; i < NUM; i++)
		M[i] = (double*)malloc(sizeof(double)*NUM);

	b = (double*)malloc(sizeof(double)*NUM);

	N = (double*)malloc(sizeof(double)*NUM);

	x1 = (double*)malloc(sizeof(double)*NUM);
	x2 = (double*)malloc(sizeof(double)*NUM);
	x3 = (double*)malloc(sizeof(double)*NUM);
	x = (double*)malloc(sizeof(double)*NUM);

	L = (int**)malloc(sizeof(int)*NUM);
	for (int i = 0; i < NUM; i++)
		L[i] = (int*)malloc(sizeof(int)*NUM);

	D = (int**)malloc(sizeof(int)*NUM);
	for (int i = 0; i < NUM; i++)
		D[i] = (int*)malloc(sizeof(int)*NUM);

	U = (int**)malloc(sizeof(int)*NUM);
	for (int i = 0; i < NUM; i++)
		U[i] = (int*)malloc(sizeof(int)*NUM);

	for (int i = 0; i < NUM; i++)
	{
		x1[i] = 0;
		x2[i] = 0;
		x3[i] = 0;
		x[i] = 0;
	}

#pragma endregion

	crtMatrix(A, 9);
	crtMatrix(C, 3);
	crtVector_b(b);
	crtVector_N(N, A);
	crtMatrix_M(M, A, N);

#pragma region MATRIX_A
	
	cout << "Matrix A " << endl;

	cout << "Start Jacobi " << GetTickCount() - start << endl;
	cout << endl;
	cout << Jacobi(x1, x, A, b, N, M) << endl << endl;

	cout << "Wyniki dla metody Jacobiego" << endl;
	cout << endl;
	cout << "Stop Jacobi " << GetTickCount() - start << endl;
	cout << endl;
	for (int i = 0; i < NUM; i++)
		cout <</* "x[" << i << "] = " << */x[i] << endl;

	cout << endl;
	cout << GetTickCount() - start << endl;
	cout << endl;

	AintoLDU(A, L, D, U);
	cout << "Start GaussaSeidla " << GetTickCount() - start << endl;
	cout << endl;
	cout << GaussaSeidla(x2, A, b, N, L, U) << endl << endl;

	cout << endl << endl << endl;


	cout << "Wyniki dla metody Gaussa-Seidla" << endl;
	cout << endl;
	cout << "Stop GaussaSeidla " << GetTickCount() - start << endl;
	cout << endl;
	for (int i = 0; i < NUM; i++)
		cout <</* "x[" << i << "] = " <<*/ x2[i] << endl;

	cout << "Start Gaussa " << GetTickCount() - start << endl;
	cout << endl;

	cout << Gaussa(A, b, x3) << endl;

	cout << endl << endl << endl;

	cout << "Wyniki dla metody Gaussa" << endl;
	cout << endl;
	cout << "Stop Gaussa " << GetTickCount() - start << endl;
	cout << endl;
	
	for (int i = 0; i < NUM; i++)
		cout <</* "x[" << i << "] = " <<*/ x3[i] << endl;
	

#pragma endregion


	for (int i = 0; i < NUM; i++)
	{
		x1[i] = 0;
		x2[i] = 0;
		x3[i] = 0;
		x[i] = 0;
	}

	crtVector_N(N, C);
	crtMatrix_M(M, C, N);


#pragma region MATRIX_C

	cout << "Matrix C "<< endl;

	cout << "Start Jacobi " << time(NULL) - start << endl;
	cout << endl;
	Jacobi(x1, x, C, b, N, M);

	cout << "Wyniki dla metody Jacobiego" << endl;
	cout << endl;
	cout << "Stop Jacobi " << time(NULL) - start << endl;
	cout << endl;
	for (int i = 0; i < NUM; i++)
		cout <</* "x[" << i << "] = " <<*/ x[i] << endl;

	cout << endl;
	cout << time(NULL) - start << endl;
	cout << endl;

	AintoLDU(C, L, D, U);

	cout << "Start GaussaSeidla " << time(NULL) - start << endl;
	cout << endl;
	GaussaSeidla(x2, C, b, N, L, U);

	cout << endl << endl << endl;


	cout << "Wyniki dla metody Gaussa-Seidla" << endl;
	cout << endl;
	cout << "Stop GaussaSeidla " << time(NULL) - start << endl;
	cout << endl;
	for (int i = 0; i < NUM; i++)
		cout <</* "x[" << i << "] = " <<*/ x2[i] << endl;

	cout << "Start Gaussa " << time(NULL) - start << endl;
	cout << endl;

	cout << Gaussa(C, b, x3) << endl;

	cout << endl << endl << endl;

	cout << "Wyniki dla metody Gaussa" << endl;
	cout << endl;
	cout << "Stop Gaussa " << time(NULL) - start << endl;
	cout << endl;

	for (int i = 0; i < NUM; i++)
		cout <</* "x[" << i << "] = " <<*/ x3[i] << endl;


#pragma endregion

	return 0;
}