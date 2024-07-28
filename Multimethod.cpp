#pragma once
#include <cstring>
#include <iostream>
#include <random>
#include <fstream>
#include <ctime>
#include <cmath>
#include<math.h>
#include "Multimethod.h"

double MultiMet::randDouble() {
    return (double)rand() / (double) RAND_MAX;
}
double MultiMet::gaussss(double mu, double sigma){
    return mu + sigma * sqrt(-2.0 * log(randDouble())) * sin(2.0 * M_PI * randDouble());
}
double MultiMet::cauchy_g(double mu, double gamma) {
    return mu + gamma * tan(M_PI*(randDouble() - 0.5));
}

void MultiMet::pop_update(int p_start, int p_end)
{
    for (int i = p_start; i < p_end; i ++)
    {
        if (pop_fit[i] < ibest_fit[i])
        {
            for (int j = 0; j < Nvar; j ++)
                ibest[i][j] = pop[i][j];
            ibest_fit[i] = pop_fit[i];
        }
    }
}

MultiMet::MultiMet(int psize, int nn, double lb, double ub, FF evaluate)
        : Population(psize, nn, lb, ub)
{
    EvaluFunc = evaluate;
}

MultiMet::~MultiMet()
{
    delete [] M;
    delete [] OShift;

    delete [] SS;
}

void MultiMet::Initial(int funcnum, int CECnum)
{
    ini_check = true;
    nfes = 0;

    if (CECnum == 2014) {
        int i;
        int cf_num = 10;
        r_flag = 1;
        s_flag = 1;
        func_num = funcnum;
        if (func_num == 8 || func_num == 10)
            r_flag = 0;
        if (!(Nvar == 2 || Nvar == 10 || Nvar == 20 || Nvar == 30 || Nvar == 50 || Nvar == 100)) {
            printf("\nError: Test functions are only defined for D=2,10,20,30,50,100.\n");
        }
        if (Nvar == 2 && ((func_num >= 17 && func_num <= 22) || (func_num >= 29 && func_num <= 30))) {
            printf("\nError: hf01,hf02,hf03,hf04,hf05,hf06,cf07&cf08 are NOT defined for D=2.\n");
        }
        char FileName[256];

        /* Load Matrix M*/
        sprintf(FileName, "../CEC2014/input_data/M_%d_D%d.txt", func_num, Nvar);
        std::ifstream fpt;
        fpt.open(FileName);
        if (!fpt.is_open()) {
            printf("\n Error: Cannot open input file for reading \n");
        }
        if (func_num < 23) {
            M = new double[Nvar * Nvar];
            for (i = 0; i < Nvar * Nvar; i++) {
                fpt >> M[i];
            }
        } else {
            M = new double[cf_num * Nvar * Nvar];
            for (i = 0; i < cf_num * Nvar * Nvar; i++) {
                fpt >> M[i];
            }
        }
        fpt.close();

        /* Load shift_data */
        sprintf(FileName, "../CEC2014/input_data/shift_data_%d.txt", func_num);
        fpt.open(FileName);
        if (!fpt.is_open()) {
            printf("\n Error: Cannot open input file for reading \n");
        }

        if (func_num < 23) {
            OShift = new double[Nvar];
            for (i = 0; i < Nvar; i++) {
                fpt >> OShift[i];
            }
        } else {
            double value;
            OShift = new double[Nvar * cf_num];
            for (i = 0; i < cf_num - 1; i++) {
                for (int j = 0; j < 100; j++) {
                    fpt >> value;
                    if (j < Nvar) {
                        OShift[i * Nvar + j] = value;
                    }
                }
            }
            for (int j = 0; j < Nvar; j++) {
                fpt >> OShift[(cf_num - 1) * Nvar + j];
            }

        }
        fpt.close();

        /* Load Shuffle_data */
        SS=new int[Nvar];
        if (func_num>=17&&func_num<=22)
        {
            sprintf(FileName, "../CEC2014/input_data/shuffle_data_%d_D%d.txt", func_num, Nvar);
            fpt.open(FileName);
            if (!fpt.is_open())
            {
                printf("\n Error: Cannot open input file for reading \n");
            }
            SS=new int[Nvar];
            for(i=0;i<Nvar;i++)
            {
                fpt >> SS[i];
            }
            fpt.close();
        }
        else if (func_num==29||func_num==30)
        {
            sprintf(FileName, "../CEC2014/input_data/shuffle_data_%d_D%d.txt", func_num, Nvar);
            fpt.open(FileName);
            if (!fpt.is_open())
            {
                printf("\n Error: Cannot open input file for reading \n");
            }
            SS=new int [Nvar*cf_num];
            for(i=0;i<Nvar*cf_num;i++)
            {
                fpt >> SS[i];
            }
            fpt.close();
        }
    }

    if(CECnum == 2017) {
        func_num = funcnum;
        s_flag = 1;
        r_flag = 1;
        int cf_num = 10;
        char FileName[256];
        std::ifstream fpt;
        if (!(Nvar == 2 || Nvar == 10 || Nvar == 20 || Nvar == 30 || Nvar == 50 || Nvar == 100)) {
            std::cout << "\nError: Test functions are only defined for D=2,10,20,30,50,100.\n";
        }

        if (Nvar == 2 && ((func_num >= 17 && func_num <= 22) || (func_num >= 29 && func_num <= 30))) {
            std::cout << "\nError: hf01,hf02,hf03,hf04,hf05,hf06,cf07&cf08 are NOT defined for D=2.\n";
        }

        sprintf(FileName, "../CEC2017/input_data/M_%d_D%d.txt", func_num, Nvar);
        fpt.open(FileName);
        if (!fpt.is_open()) {
            std::cout << "\n Error: Cannot open input file for reading \n";
        }

        if (func_num < 20) {
            M = new double[Nvar * Nvar];
            for (int i = 0; i < Nvar * Nvar; i++)
                M[i] = 0.0;
            for (int i = 0; i < Nvar * Nvar; i++) {
                fpt >> M[i];
            }
        } else {
            M = new double[cf_num * Nvar * Nvar];
            for (int i = 0; i < cf_num * Nvar * Nvar; i++)
                M[i] = 0.0;
            for (int i = 0; i < cf_num * Nvar * Nvar; i++) {
                fpt >> M[i];
            }
        }
        fpt.close();

        sprintf(FileName, "../CEC2017/input_data/shift_data_%d.txt", func_num);
        fpt.open(FileName);
        if (!fpt.is_open()) {
            std::cout << "\n Error: Cannot open input file for reading \n";
        }

        if (func_num < 20) {
            OShift = new double[Nvar];
            for (int i = 0; i < Nvar; i++) {
                fpt >> OShift[i];
            }
        } else {
            OShift = new double[Nvar * cf_num];
            for (int i = 0; i < cf_num - 1; i++) {
                double value;
                for (int j=0;j<100;j++)
                {
                    fpt >> value;
                    if(j < Nvar) {
                        OShift[i * Nvar + j] = value;
                    }
                }
            }
            for (int j = 0; j < Nvar; j++) {
                fpt >> OShift[(cf_num - 1) * Nvar + j];
            }
        }
        fpt.close();

        SS = new int[Nvar];
        if (func_num >= 11 && func_num <= 20) {
            sprintf(FileName, "../CEC2017/input_data/shuffle_data_%d_D%d.txt", func_num, Nvar);
            fpt.open(FileName);
            if (!fpt.is_open()) {
                std::cout << "\n Error: Cannot open input file for reading \n";
            }
            SS = new int[Nvar];
            for (int i = 0; i < Nvar; i++) {
                SS[i] = 0;
            }
            for (int i = 0; i < Nvar; i++) {
                fpt >> SS[i];
            }
            fpt.close();
        } else if (func_num == 29 || func_num == 30) {
            sprintf(FileName, "../CEC2017/input_data/shuffle_data_%d_D%d.txt", func_num, Nvar);
            fpt.open(FileName);
            if (!fpt.is_open()) {
                std::cout << "\n Error: Cannot open input file for reading \n";
            }
            SS = new int[Nvar * cf_num];
            for (int i = 0; i < Nvar * cf_num; i++) {
                SS[i] = 0;
            }
            for (int i = 0; i < Nvar * cf_num; i++) {
                fpt >> SS[i];
            }
            fpt.close();
        }
    }

    if(CECnum == 2022) {
        int cf_num = 12;
        func_num = funcnum;
        r_flag = 1;
        s_flag = 1;
        FILE *fpt;
        char FileName[256];

        if (!((Nvar == 2) || (Nvar == 10) || (Nvar == 20)))
        {
            printf("\nError: Test functions are only defined for D=2,10,20.\n");
        }
        if ((Nvar == 2) && ((func_num == 6) || (func_num == 7) || (func_num == 8)))
        {
            printf("\nError:  NOT defined for D=2.\n");
        }

        /* Load Matrix M*/
        sprintf(FileName, "../CEC2022/input_data/M_%d_D%d.txt", func_num, Nvar);
        fpt = fopen(FileName, "r");
        if (fpt == NULL)
        {
            printf("\n Error: Cannot open M_%d_D%d.txt for reading \n", func_num, Nvar);
        }
        if (func_num < 9)
        {
            M = (double *)malloc(Nvar * Nvar * sizeof(double));
            if (M == NULL)
                printf("\nError: there is insufficient memory available!\n");
            for (int i = 0; i < Nvar * Nvar; i++)
            {
                fscanf(fpt, "%lf", &M[i]);
            }
        }
        else
        {
            M = (double *)malloc(cf_num * Nvar * Nvar * sizeof(double));
            if (M == NULL)
                printf("\nError: there is insufficient memory available!\n");
            for (int i = 0; i < cf_num * Nvar * Nvar; i++)
            {
                fscanf(fpt, "%lf", &M[i]);
            }
        }
        fclose(fpt);

        /* Load shift_data */
        sprintf(FileName, "../CEC2022/input_data/shift_data_%d.txt", func_num);
        fpt = fopen(FileName, "r");
        if (fpt == NULL)
        {
            printf("\n Error: Cannot open shift_data_%d.txt for reading \n", func_num);
        }

        if (func_num < 9)
        {
            OShift = (double *)malloc(Nvar * sizeof(double));
            if (OShift == NULL)
                printf("\nError: there is insufficient memory available!\n");
            for (int i = 0; i < Nvar; i++)
            {
                fscanf(fpt, "%lf", &OShift[i]);
            }
        }
        else
        {
            OShift = (double *)malloc(Nvar * cf_num * sizeof(double));
            if (OShift == NULL)
                printf("\nError: there is insufficient memory available!\n");
            for (int i = 0; i < cf_num - 1; i++)
            {
                for (int j = 0; j < Nvar; j++)
                {
                    fscanf(fpt, "%lf", &OShift[i * Nvar + j]);
                }
                fscanf(fpt, "%*[^\n]%*c");
            }
            for (int j = 0; j < Nvar; j++)
            {
                fscanf(fpt, "%lf", &OShift[Nvar * (cf_num - 1) + j]);
            }
        }
        fclose(fpt);

        /* Load Shuffle_data */
        SS = (int *)malloc(Nvar * sizeof(int));
        if ((func_num >= 6) && (func_num <= 8))
        {
            sprintf(FileName, "../CEC2022/input_data/shuffle_data_%d_D%d.txt", func_num, Nvar);
            fpt = fopen(FileName, "r");
            if (fpt == NULL)
            {
                printf("\n Error: Cannot open shuffle_data_%d_D%d.txt for reading \n", func_num, Nvar);
            }
            SS = (int *)malloc(Nvar * sizeof(int));
            if (SS == NULL)
                printf("\nError: there is insufficient memory available!\n");
            for (int i = 0; i < Nvar; i++)
            {
                fscanf(fpt, "%d", &SS[i]);
            }
            fclose(fpt);
        }
    }

    for (int i = 0; i < Popsize; i ++) {
        for (int j = 0; j < Nvar; j++) {
            newpop[i][j] = pop[i][j] = ibest[i][j] = randval(Lbound, Ubound);
        }
        newpop_fit[i] = pop_fit[i] = ibest_fit[i] = EvaluFunc(pop[i], Nvar, OShift, M, SS, s_flag, r_flag);
    }

    nfes += Popsize;

    worst_and_best();
    for (int j = 0; j < Nvar; j ++)
        gbest[j] = pop[cur_best][j];
    gbest_fit = pop_fit[cur_best];
}

void MultiMet::TPADE(double iniCR, double iniF, double r_arc, double pmax, double pmin, int MAXNFEs)
{
    if (ini_check) {
        archivePop.resize(0);
        archive.resize(0);
        archive_size = (int)(Popsize * r_arc);
        NP0 = Popsize;
        uF = iniF;
        uCR = iniCR;
        ini_check = false;
    }

    vector<int> index;
    for (int i = 0; i < Popsize; i++)
        index.push_back(i);
    double **copy_pop = CreateMatrix(Popsize, Nvar);
    double *copy_fit = new double[Popsize];
    for (int i = 0; i < Popsize; i++) {
        for (int j = 0; j < Nvar; j++)
            copy_pop[i][j] = pop[i][j];
        copy_fit[i] = pop_fit[i];
    }
    sort(index.begin(), index.end(),
         [&](const int &a, const int &b) {
             return (pop_fit[a] < pop_fit[b]);
         }
    );
    for (int i = 0; i < Popsize; i++) {
        for (int j = 0; j < Nvar; j++)
            pop[i][j] = copy_pop[index[i]][j];
        pop_fit[i] = copy_fit[index[i]];
    }

    for (int i = 0; i < Popsize; i++) {
        for (int j = 0; j < Nvar; j++)
            copy_pop[i][j] = ibest[i][j];
        copy_fit[i] = ibest_fit[i];
    }
    for (int i = 0; i < Popsize; i++) {
        for (int j = 0; j < Nvar; j++)
            ibest[i][j] = copy_pop[index[i]][j];
        ibest_fit[i] = copy_fit[index[i]];
    }

    DeleteMatrix(copy_pop, Popsize);
    delete[] copy_fit;

    double *ind = new double[Nvar];
    vector<double>delta_fcr;
    vector<double>delta_ff;
    vector<double> Sf;
    vector<double> Scr;
    double countt = 0;
    archivePop.resize(0);
    int EN = int(Popsize * (pmax - (pmax - pmin) * nfes / MAXNFEs));
    int FN = int(Popsize * (pmax - (pmax - pmin) * nfes / MAXNFEs));
    if(EN < 0)
        EN = 0;
    if(FN < 0)
        FN = 0;
    int MN = Popsize - EN - FN;
    archivePop.resize(0);
    vector<vector<double>>newarchive;
    vector<int>indexpop;
    for (int i = 0; i < Popsize; i++) {
        double F;
        do {
            F = cauchy_g(uF, 0.1);
        } while (F < 0);
        if (F > 1.0)
            F = 1.0;
        if (nfes < 0.6 * MAXNFEs)
            F = min(F, 0.7);
        double F2;
        if (nfes < 0.2 * MAXNFEs)
            F2 = 0.7 * F;
        else if (nfes < 0.4 * MAXNFEs)
            F2 = 0.8 * F;
        else
            F2 = 1.2 * F;

        double CR;
        CR = gaussss(uCR, 0.1);
        if (CR < 0)
            CR = 0;
        else if (CR > 1)
            CR = 1;

        if (i < EN) {
            int pos, mem;
            do {
                pos = rand() % Popsize;
                mem = rand() % (Popsize + archive.size());
            } while (pos == mem || pos == i);
            int x1 = rand() % Popsize;
            int x2;
            do{
                x2 = rand() % Popsize;
            }while(x2 == x1);
            F2 = F2 * 1.2;
            double w = 1.0 - pow((double)nfes/MAXNFEs,2);
            for (int j = 0; j < Nvar; j++) {
                if(mem < Popsize)
                    ind[j] = w * pop[i][j] + (1-w) * ibest[rand() % Popsize][j] + F2 * (pop[pos][j] - pop[mem][j]) + F2 * (pop[x1][j] - pop[x2][j]); // 后一个本来是F，改成F2（325）
                else
                    ind[j] = w * pop[i][j] + (1-w) * ibest[rand() % Popsize][j] + F2 * (pop[pos][j] - archive[mem - Popsize][j]) + F2 * (pop[x1][j] - pop[x2][j]); // 后一个本来是F
                if (ind[j] > Ubound)
                    ind[j] = (pop[i][j] + Ubound) / 2.0;
                else if (ind[j] < Lbound)
                    ind[j] = (pop[i][j] + Lbound) / 2.0;
            }
        }
        else if (i > EN + MN) {
            int xp = rand() % EN;
            int xl = rand() % (Popsize - EN + archive.size()) + EN;
            F2 = F2 * 0.8;
            for (int j = 0; j < Nvar; j++) {
                if (xl < Popsize)
                    ind[j] = pop[i][j] + F2 * (pop[xp][j] - pop[xl][j]);
                else
                    ind[j] = pop[i][j] + F2 * (pop[xp][j] - archive[xl - Popsize][j]);
                if (ind[j] > Ubound)
                    ind[j] = (pop[i][j] + Ubound) / 2.0;
                else if (ind[j] < Lbound)
                    ind[j] = (pop[i][j] + Lbound) / 2.0;
            }
        }
        else {
            int pos, mem;
            do {
                mem = rand() % (Popsize + archive.size());
                pos = rand() % Popsize;
            } while (pos == mem || pos == i);
            int xp = rand() % EN;
            for (int j = 0; j < Nvar; j++) {
                if(mem < Popsize)
                    ind[j] = pop[i][j] + F2 * (pop[xp][j] - pop[i][j]) + F * (pop[pos][j] - pop[mem][j]);
                else
                    ind[j] = pop[i][j] + F2 * (pop[xp][j] - pop[i][j]) + F * (pop[pos][j] - archive[mem - Popsize][j]);
                if (ind[j] > Ubound)
                    ind[j] = (pop[i][j] + Ubound) / 2.0;
                else if (ind[j] < Lbound)
                    ind[j] = (pop[i][j] + Lbound) / 2.0;
            }
        }

        int d = rand() % Nvar;
        double CR1 = pow((double)nfes/MAXNFEs,2);
        double CR2 = pow((double)nfes/MAXNFEs,2);
        for (int j = 0; j < Nvar; j++) {
            double rr = randval(0, 1);
            if (i > EN + MN) {
                if (rr < CR1 && j != d)
                    ind[j] = gbest[j];
            }
            else if (i < EN) {
                if (rr < CR2 && j != d)
                    ind[j] = ibest[i][j];
            }
            else {
                if (rr > CR && j != d)
                    ind[j] = pop[i][j];
            }
        }

        double fit = EvaluFunc(ind, Nvar, OShift, M, SS, s_flag, r_flag);
        if (fit <= pop_fit[i]) {
            for (int j = 0; j < Nvar; j++)
                newpop[i][j] = ind[j];
            newpop_fit[i] = fit;
        }
        else {
            if(i < EN){
                for (int j = 0; j < Nvar; j++)
                    newpop[i][j] = pop[i][j];
                newpop_fit[i] = pop_fit[i];
                vector<double>individual;
                for(int k = 0; k < Nvar; k++)
                    individual.push_back(ind[k]);
                individual.push_back(fit);
                if(fit <= pop_fit[Popsize - 1])
                    archivePop.push_back(individual);
            }
            else if(i >= EN && i <= EN + MN){
                vector<double>individual;
                for(int k = 0; k < Nvar; k++)
                    individual.push_back(ind[k]);
                individual.push_back(fit);
                if(fit <= pop_fit[Popsize - 1])
                    archivePop.push_back(individual);
                for (int j = 0; j < Nvar; j++)
                    newpop[i][j] = pop[i][j];
                newpop_fit[i] = pop_fit[i];
            }
            else {
                indexpop.push_back(i);
            }
        }
        if (fit < pop_fit[i]) {
            double dis = 0;
            for(int k = 0; k < Nvar; k++)
                dis += pow(pop[i][k] - ind[k],2);
            dis = sqrt(dis);
            Sf.push_back(F);
            if(i <= EN + MN && i >= EN){
                Scr.push_back(CR);
                delta_fcr.push_back(dis);
            }
            delta_ff.push_back(dis);
            vector<double> individual;
            for (int j = 0; j < Nvar; j++)
                individual.push_back(pop[i][j]);
            newarchive.push_back(individual);
        }
    }

    int point1 = 0;
    int point2 = 0;
    vector<int>ids;
    double p = pow((double)nfes/MAXNFEs,0.5);
    for(auto it = indexpop.begin(); it != indexpop.end(); ++it){
        double r = randval(0,1);
        if(r > p && point1 < archivePop.size()) {
            for (int j = 0; j < Nvar; j++)
                newpop[*it][j] = archivePop[point1][j];
            newpop_fit[*it] = archivePop[point1][Nvar];
            point1++;
            ids.push_back(point1);
        }
        else{
            for (int j = 0; j < Nvar; j++)
                newpop[*it][j] = pop[indexpop[point2]][j];
            newpop_fit[*it] = pop_fit[indexpop[point2]];
            point2 ++;
        }
    }

    for(int i = 0; i < newarchive.size(); i++)
        archive.push_back(newarchive[i]);
    int Asize = archive.size();
    if(Asize > archive_size){
        int delta_size = Asize - archive_size;
        for(int i = 0; i < archive_size; i++){
            for(int j = 0; j < Nvar; j++){
                archive[i][j] = archive[i + delta_size][j];
            }
        }
        archive.resize(archive_size);
    }

    if(delta_ff.size() > 0) {
        double c = 0.1;
        double sum_f = accumulate(delta_ff.begin(), delta_ff.end(), 0.);
        double *omega = new double[Sf.size()];
        for (int u = 0; u < Sf.size(); u++)
            omega[u] = delta_ff[u] / sum_f;
        double sum1 = 0.;
        double sum2 = 0.;
        for (int u = 0; u < Sf.size(); u++) {
            sum1 += omega[u] * pow(Sf[u], 2);
            sum2 += omega[u] * Sf[u];
        }
        uF = (1 - c) * uF + c * sum1 / sum2;
        delete [] omega;
    }
    if(delta_fcr.size() > 0) {
        double c = 0.1;
        double sum_f = accumulate(delta_fcr.begin(), delta_fcr.end(), 0.);
        double *omega = new double[Scr.size()];
        for (int u = 0; u < Scr.size(); u++)
            omega[u] = delta_fcr[u] / sum_f;
        double sum1 = 0.;
        double sum2 = 0.;
        for (int u = 0; u < Scr.size(); u++) {
            sum1 += omega[u] * pow(Scr[u], 2);
            sum2 += omega[u] * Scr[u];
        }
        uCR = (1 - c) * uCR + c * sum1 / sum2;
        delete[] omega;
    }


    nfes += Popsize;

    int min_NP = 10;
    int N = round(NP0 - (double)(NP0 - min_NP) / MAXNFEs * nfes);
    vector<int>index2;
    for(int u = 0; u < Popsize; u++)
        index2.push_back(u);
    sort(index2.begin(), index2.end(),
         [&](const int& a, const int& b) {
             return (newpop_fit[a] < newpop_fit[b]);
         }
    );
    for(int i = 0; i < N; i++) {
        for (int j = 0; j < Nvar; j++)
            pop[i][j] = newpop[index2[i]][j];
        pop_fit[i] = newpop_fit[index2[i]];
    }

    if (pop_fit[0] < gbest_fit)
    {
        for (int i = 0; i < Nvar; i ++)
            gbest[i] = pop[0][i];
        gbest_fit = pop_fit[0];
    }

    Popsize = N;
    archive_size = int(Popsize * r_arc);

    delete[] ind;
}

