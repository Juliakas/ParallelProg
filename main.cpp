#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include <iomanip>

using namespace std;

int numDP = 10000; // Vietoviu skaicius (demand points, max 10000)
int numPF = 5;	   // Esanciu objektu skaicius (preexisting facilities)
int numCL = 50;	   // Kandidatu naujiems objektams skaicius (candidate locations)
int numX = 3;	   // Nauju objektu skaicius
int iters = 10120; // Iteraciju skaicius
int threadCount = 4;

double **adjacencyMatrix; // Gretimumo matrica turinti atstumus tarp tasku (vietoviu)
double **demandPoints;	  // Geografiniai duomenys
int *X;					  // Sprendinys

//=============================================================================

double getTime();
void loadDemandPoints();
void randomSolution(int *X);
double HaversineDistance(double *a, double *b);
double evaluateSolution(int *X);
void calculateAllDistances();

//=============================================================================

int main()
{
	omp_set_num_threads(threadCount);
	double ts = getTime(); // Algoritmo vykdymo pradzios laikas
	double seqT = 0;
	double parT = 0;

	loadDemandPoints(); // Nuskaitomi duomenys

	double tp = getTime();
	seqT += tp - ts;

	calculateAllDistances();

	ts = getTime();
	parT += ts - tp;

	X = new int[numX];			// Sprndinys
	double u;					// Sprendinio tikslo funkcijos reiksme
	int *bestX = new int[numX]; // Geriausias rastas sprendinys
	double bestU = -1;			// Geriausio sprendinio tikslo funkcijos reiksme

	//----- Pagrindinis ciklas ------------------------------------------------

	omp_set_num_threads(threadCount);

	for (int iter = 0; iter < iters; iter++)
	{
		// Generuojam atsitiktini sprendini ir tikrinam ar jis nera geresnis uz geriausia zinoma
		randomSolution(X);

		tp = getTime();
		seqT += tp - ts;

		u = evaluateSolution(X);

		ts = getTime();
		parT += ts - tp;

		if (u > bestU)
		{ // Jei geresnis, tai issaugojam kaip geriausia zinoma
			bestU = u;
			for (int i = 0; i < numX; i++)
				bestX[i] = X[i];
		}
	}
	//----- Rezultatu spausdinimas --------------------------------------------

	seqT += getTime() - ts;

	cout << "Geriausias sprendinys: ";
	for (int i = 0; i < numX; i++)
		cout << bestX[i] << " ";
	cout << "(" << bestU << ")" << endl
		 << "Nuosekliosios dalies trukme: " << fixed << setprecision(2) << seqT << endl
		 << "Lygiagrecios dalies trukme: " << fixed << setprecision(2) << parT << endl
		 << "Visa trukme: " << fixed << setprecision(2) << seqT + parT << endl;

	delete adjacencyMatrix;
	delete demandPoints;
}

//=============================================================================

void loadDemandPoints()
{

	//----- Load demand points ------------------------------------------------
	FILE *f;
	f = fopen("demandPoints.dat", "r");
	demandPoints = new double *[numDP];
	for (int i = 0; i < numDP; i++)
	{
		demandPoints[i] = new double[3];
		fscanf(f, "%lf%lf%lf", &demandPoints[i][0], &demandPoints[i][1], &demandPoints[i][2]);
	}
	fclose(f);
}

//=============================================================================

double HaversineDistance(double *a, double *b)
{
	double dlon = fabs(a[0] - b[0]);
	double dlat = fabs(a[1] - b[1]);
	double aa = pow((sin((double)dlon / (double)2 * 0.01745)), 2) + cos(a[0] * 0.01745) * cos(b[0] * 0.01745) * pow((sin((double)dlat / (double)2 * 0.01745)), 2);
	double c = 2 * atan2(sqrt(aa), sqrt(1 - aa));
	double d = 6371 * c;
	return d;
}

//=============================================================================

double getTime()
{
	struct timeval laikas;
	gettimeofday(&laikas, NULL);
	double rez = (double)laikas.tv_sec + (double)laikas.tv_usec / 1000000;
	return rez;
}

//=============================================================================

void randomSolution(int *X)
{
	int unique;
	for (int i = 0; i < numX; i++)
	{
		do
		{
			unique = 1;
			X[i] = (int)((double)rand() / RAND_MAX * numCL);
			for (int j = 0; j < i; j++)
				if (X[j] == X[i])
				{
					unique = 0;
					break;
				}
		} while (unique == 0);
	}
}

//=============================================================================

double evaluateSolution(int *X)
{
	double U = 0;
#pragma omp parallel
	{
		int bestPF;
		int bestX;
		double d;
#pragma omp for
		for (int i = 0; i < numDP; i++)
		{
			bestPF = 1e5;
			for (int j = 0; j < numPF; j++)
			{
				d = adjacencyMatrix[i][j];
				if (d < bestPF)
					bestPF = d;
			}
			bestX = 1e5;
			for (int j = 0; j < numX; j++)
			{
				d = adjacencyMatrix[i][X[j]];
				if (d < bestX)
					bestX = d;
			}
			if (bestX < bestPF)
#pragma omp atomic
				U += demandPoints[i][2];
			else if (bestX == bestPF)
#pragma omp atomic
				U += 0.3 * demandPoints[i][2];
		}
	}
	return U;
}

//=============================================================================

void calculateAllDistances()
{
	adjacencyMatrix = new double *[numDP];
#pragma omp parallel for
	for (int i = 0; i < numDP; i++)
	{
		adjacencyMatrix[i] = new double[numDP];
		for (int j = 0; j < numDP; j++)
		{
			adjacencyMatrix[i][j] = HaversineDistance(demandPoints[i], demandPoints[j]);
		}
	}
}