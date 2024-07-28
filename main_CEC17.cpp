#include "Multimethod.h"
#include <iostream>
#include <fstream>
#include <time.h>
using namespace std;
#define IRUN 51

FF select_func(int funcnum)
{
    switch (funcnum) {
        case 1:
            return bent_cigar_func;
        case 2:
            return bent_cigar_func;
        case 3:
            return zakharov_func;
        case 4:
            return rosenbrock_func;
        case 5:
            return rastrigin_func;
        case 6:
            return schaffer_F7_func;
        case 7:
            return bi_rastrigin_func;
        case 8:
            return step_rastrigin_func;
        case 9:
            return levy_func;
        case 10:
            return schwefel_func;
        case 11:
            return hf01;
        case 12:
            return hf02;
        case 13:
            return hf03;
        case 14:
            return hf04;
        case 15:
            return hf05;
        case 16:
            return hf06;
        case 17:
            return hf07;
        case 18:
            return hf08;
        case 19:
            return hf09;
        case 20:
            return hf10;
        case 21:
            return cf01;
        case 22:
            return cf02;
        case 23:
            return cf03;
        case 24:
            return cf04;
        case 25:
            return cf05;
        case 26:
            return cf06;
        case 27:
            return cf07;
        case 28:
            return cf08;
        case 29:
            return cf09;
        case 30:
            return cf10;
        default:
            printf("\nError: There are only 29 test functions in this test suite!\n");
            break;
    }
}

int main()
{
    int NVAR = 10;
    int MAXNFEs = NVAR * 1e4;
    int CECnum = 2017;
    for (int funcnum = 1; funcnum <= 30; funcnum++) {
        if (funcnum != 2) {
            cout << "------- func_num = " << funcnum << " -------" << endl;
            clock_t t1, t2;
            ofstream outf1;
            outf1.open("..//2017results//data//TPADE//" + to_string(NVAR) + "D//f" + to_string(funcnum) + ".txt");
            ofstream outf2;
            outf2.open("..//2017results//results//TPADE//" + to_string(NVAR) + "D//f" + to_string(funcnum) + ".txt");
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
                outf1 << solver.nfes << " " << solver.gbest_fit - funcnum * 100.0 << endl;
                generation += 1;
                while (solver.nfes < MAXNFEs)
                {
                    solver.TPADE(0.5, 0.5, 2.6, 0.2, 0.1, MAXNFEs);
                    solver.pop_update(0, solver.Popsize);
                    outf1 << solver.nfes << " " << solver.gbest_fit - funcnum * 100.0 << endl;
                    generation++;
                }
                cout << "irun: " << nrun << " , The best solution of F" << funcnum << " = "
                     << solver.gbest_fit - funcnum * 100.0 << endl;
                result[nrun] = solver.gbest_fit - funcnum * 100.0;
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
    }
    return 0;
}
