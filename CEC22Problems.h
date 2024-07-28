#include <math.h>
#include <time.h>
#include <random>
using namespace std; //pragma once

double ackley_func(double* var, int num, double* OShift, double* M, int* S, int s_flag, int r_flag); /* Ackley's */
double bent_cigar_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag);/* Discus */
double discus_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag); /* Bent_Cigar */
double ellips_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag);/* Ellipsoidal */
double escaffer6_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag);/* Expanded Scaffer¡¯s F6  */
double griewank_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag);/* Griewank's  */
double grie_rosen_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag);/* Griewank-Rosenbrock  */
double happycat_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag);/* HappyCat */
double hgbat_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag);/* HGBat  */
double rosenbrock_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag);/* Rosenbrock's */
double rastrigin_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag);/* Rastrigin's  */
double schwefel_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag);/* Schwefel's */
double schaffer_F7_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag);/* Schwefel's F7 */
double step_rastrigin_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag);/* Noncontinuous Rastrigin's  */
double levy_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag);/* Levy */
double zakharov_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag);/* ZAKHAROV */
double katsuura_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag);/* Katsuura */



double hf02(double* var, int num, double* OShift, double* M, int* S, int s_flag, int r_flag); /* Hybrid Function 2 */
double hf06(double* var, int num, double* OShift, double* M, int* S, int s_flag, int r_flag); /* Hybrid Function 6 */
double hf10(double* var, int num, double* OShift, double* M, int* S, int s_flag, int r_flag); /* Hybrid Function 10 */


double cf01(double* var, int num, double* OShift, double* M, int* S, int s_flag, int r_flag); /* Composition Function 1 */
double cf02(double* var, int num, double* OShift, double* M, int* S, int s_flag, int r_flag); /* Composition Function 2 */
double cf06(double* var, int num, double* OShift, double* M, int* S, int s_flag, int r_flag); /* Composition Function 6 */
double cf07(double* var, int num, double* OShift, double* M, int* S, int s_flag, int r_flag); /* Composition Function 7 */



void shiftfunc(double* var, double* varshift, int num, double* OShift);
void rotatefunc(double* var, double* varrot, int num, double* M);
void sr_func(double* var, double* sr_var, int num, double* OShift, double* M, double sh_rate, int s_flag, int r_flag); /* shift and rotate */
void asyfunc(double* var, double* varasy, int num, double beta);
void oszfunc(double* var, double* varosz, int num);
double cf_cal(double* var, int num, double* OShift, double* delta, double* bias, double* fit, int cf_num);

extern double* OShift, * M, * y, * z, * x_bound;
extern int ini_flag, n_flag, func_flag, * SS;
