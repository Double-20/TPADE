#include "Multimethod.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;
#define IRUN 30
FF select_func(int funcnum)
{
    switch (funcnum) {
        case 1:
            return zakharov_func;
        case 2:
            return rosenbrock_func;
        case 3:
            return escaffer6_func;
        case 4:
            return step_rastrigin_func;
        case 5:
            return levy_func;
        case 6:
            return hf02;
        case 7:
            return hf10;
        case 8:
            return hf06;
        case 9:
            return cf01;
        case 10:
            return cf02;
        case 11:
            return cf06;
        case 12:
            return cf07;
        default:
            printf("\nError: There are only 12 test functions in this test suite!\n");
            break;
    }
}
double fun_value(int funcnum)
{
    switch (funcnum) {
        case 1:
            return 300;
        case 2:
            return 400;
        case 3:
            return 600;
        case 4:
            return 800;
        case 5:
            return 900;
        case 6:
            return 1800;
        case 7:
            return 2000;
        case 8:
            return 2200;
        case 9:
            return 2300;
        case 10:
            return 2400;
        case 11:
            return 2600;
        case 12:
            return 2700;
        default:
            printf("\nError: There are only 12 test functions in this test suite!\n");
            break;
    }
}

int main()
{
    int NVAR = 10;
    int MAXNFEs;
    if(NVAR == 10)
        MAXNFEs = 2 * 1e5;
    else if(NVAR == 20)
        MAXNFEs = 1e6;
    else{
        cout << "NVAR is invalid." << endl;
        return 0;
    }
    int CECnum = 2022;
    double* seed_values = new double[1000];
    ifstream infile("..//CEC2022//input_data//Rand_Seeds.txt");
    if (!infile) {
        cerr << "Error opening file" << endl;
        return 1;
    }
    int u = 0;
    double value;
    while (infile >> value) {
        seed_values[u++] = value;
    }
    infile.close();
    for (int funcnum = 1; funcnum <= 12; funcnum++) {
        cout << "------- func_num = " << funcnum << " -------" << endl;
        clock_t t1, t2;
        ofstream outf1;
        outf1.open("..//2022results//data//TPADE//" + to_string(NVAR) + "D//f" + to_string(funcnum) + ".txt");
        ofstream outf2;
        outf2.open("..//2022results//results//TPADE//" + to_string(NVAR) + "D//f" + to_string(funcnum) + ".txt");
        double *result = new double[IRUN];
        for (int i = 0; i < IRUN; i++)
            result[i] = 0;
        srand(unsigned(time(0)));
        t1 = clock();
        for (int nrun = 0; nrun < IRUN; nrun++) {
            int generation = 0;
            double rNinit = 18;
            int Popsize = int(NVAR * rNinit);
            MultiMet solver(Popsize, NVAR, -100, 100, select_func(funcnum));
            solver.Initial(funcnum, CECnum);
            outf1 << solver.nfes << " " << solver.gbest_fit - fun_value(funcnum) << endl;
            generation += 1;
            while (solver.nfes < MAXNFEs)
            {
                solver.TPADE(0.5, 0.5, 2.6, 0.2, 0.1, MAXNFEs);
                solver.pop_update(0, solver.Popsize);
                outf1 << solver.nfes << " " << solver.gbest_fit - fun_value(funcnum) << endl;
                generation++;
            }
            cout << "irun: " << nrun << " , The best solution of F" << funcnum << " = "
                 << solver.gbest_fit - fun_value(funcnum) << endl;
            result[nrun] = solver.gbest_fit - fun_value(funcnum);
        }

        t2 = clock();
        double ave_result = 0;
        for (int irun = 0; irun < IRUN; irun++) {
            ave_result += result[irun];
        }
        ave_result /= IRUN;

        cout << "The ave solution = " << ave_result << endl;
        cout << "Time = " << (double) (t2 - t1) / CLOCKS_PER_SEC / IRUN << " s" << endl;

        double variance = 0;
        for (int i = 0; i < IRUN; i++)
            variance += pow(result[i] - ave_result, 2);
        variance /= (double) (IRUN - 1);
        variance = sqrt(variance);

        for (int i = 0; i < IRUN; i++)
            outf2 << result[i] << " ";
        outf2 << endl;

        outf2 << ave_result << endl;
        outf2 << variance << endl;
        outf2 << (double) (t2 - t1) / CLOCKS_PER_SEC / IRUN << endl;

        delete[] result;

        outf1.close();
        outf2.close();
    }
    return 0;
}