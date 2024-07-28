#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <random>
#include <cassert>
#include <climits>
#include <iostream>
#include "Population.h"
#include "CEC17Problems.h"
#include "CEC14Problems.h"
#include "CEC22Problems.h"
using namespace std;

#ifndef _MULTIMETHOD_H
#define _MULTIMETHOD_H

typedef double (*FF) (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag);

class MultiMet : public Population<double>
{
public:
    MultiMet();
    MultiMet(int psize, int nn, double lb, double ub, FF evaluate);
    ~MultiMet();

public:
    FF EvaluFunc;
    int nfes;
    bool ini_check;

    double uF;
    double uCR;
    int NP0;
    int archive_size;
    vector<double> Mf;
    vector<double> Mcr;
    vector<vector<double>>archive;
    vector<vector<double>>archivePop;

    int func_num;
    int s_flag;
    int r_flag;
    double* M;
    double* OShift;
    int* SS;

public:

    void pop_update(int p_start, int p_end);
    void Initial(int funcnum, int CECnum);
    double randDouble();
    double gaussss(double mu, double sigma);
    double cauchy_g(double mu, double gamma);
    void TPADE(double iniCR, double iniF, double r_arc, double pmax, double pmin, int MAXNFEs);
};
#endif
