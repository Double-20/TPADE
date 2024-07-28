#pragma once
#ifndef _POPULATION_H
#define _POPULATION_H
#include <stdlib.h>
#include <math.h>
#include <time.h>
using namespace std;
#define PI 3.1415926

template<class T>
class Population
{
public:
	Population(int psize, int nn, T lb, T ub);
	virtual ~Population();

public:
	T** pop;
	double *pop_fit;
	T** newpop;
	double *newpop_fit;
	T* gbest;
	double gbest_fit;
    T* ibest_fit;
    T** ibest;

	int Popsize0;
    int Popsize;
	int Nvar;
	T Ubound;
	T Lbound;
	int cur_best;
	int cur_worst;

public:
	T** CreateMatrix(int nRow, int nCol);
	void DeleteMatrix(T **ppT, int nRow);
	inline double randval(double low, double high);
	void worst_and_best();
};

template<class T>
Population<T>::Population(int psize, int nn, T lb, T ub)
	: Nvar(nn), Popsize(psize), Lbound(lb), Ubound(ub),
	pop_fit(0), newpop_fit(0), gbest_fit(0), ibest_fit(0), cur_best(0), cur_worst(0)
{
    Popsize0 = Popsize;
	pop = CreateMatrix(Popsize0, Nvar);
	newpop = CreateMatrix(Popsize0, Nvar);
	gbest = new T[Nvar];
	pop_fit = new double[Popsize0];
	newpop_fit = new double[Popsize0];
    ibest = CreateMatrix(Popsize0, Nvar);
    ibest_fit = new double[Popsize0];
	for (int i = 0; i < Popsize0; i ++)
		for (int j = 0; j < Nvar; j ++)
			pop[i][j] = (T)randval(Lbound, Ubound);
			
	for (int i = 0; i < Popsize0; i ++)
		for (int j = 0; j < Nvar; j ++)
			newpop[i][j] = (T)randval(Lbound, Ubound);

    for (int i = 0; i < Popsize0; i ++)
        for (int j = 0; j < Nvar; j ++)
            ibest[i][j] = (T)randval(Lbound, Ubound);

	for (int i = 0; i < Nvar; i ++)
		gbest[i] = (T)randval(Lbound, Ubound);

}

template<class T>
Population<T>::~Population()
{
	DeleteMatrix(pop, Popsize0);
	DeleteMatrix(newpop, Popsize0);
    DeleteMatrix(ibest, Popsize0);
	delete [] gbest;
	delete [] pop_fit;
	delete [] newpop_fit;
    delete [] ibest_fit;
}

template<class T>
T** Population<T>::CreateMatrix(int nRow, int nCol)
{
	T **ppT = new T*[nRow];
	for(int iRow = 0; iRow < nRow; iRow ++)
	{
		ppT[iRow] = new T[nCol];
	}
	return ppT;
}

template<class T>
void Population<T>::DeleteMatrix(T **ppT, int nRow)
{
	for(int iRow = 0; iRow < nRow; iRow ++)
	{
		delete[] ppT[iRow];
	}
	delete[] ppT;
	ppT = NULL;
}

template<class T>
double Population<T>::randval(double low, double high)
{
	return ((double)(rand() % 1000) / 1000.0) * (high - low) + low;
}

template<class T>
void Population<T>::worst_and_best()
{
	int i;
	cur_best = 0;
	cur_worst = 0;
	for (i = 0; i < Popsize; i++)
	{
		if (pop_fit[i] < pop_fit[cur_best])
			cur_best = i;
		else if (pop_fit[i] > pop_fit[cur_worst])
			cur_worst = i;
	}
}


#endif
