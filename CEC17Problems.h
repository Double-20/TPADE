double sphere_func (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Sphere */
double ellips_func(double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Ellipsoidal */
double bent_cigar_func(double *var, int num, double *Os,double *Mr, int * S, int s_flag,int r_flag); /* Discus */
double discus_func(double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag);  /* Bent_Cigar */
double dif_powers_func(double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag);  /* Different Powers */
double rosenbrock_func (double* var, int num, double* Os,double* Mr,int* S,int s_flag, int r_flag); /* Rosenbrock's */
double schaffer_F7_func (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Schwefel's F7 */
double ackley_func (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Ackley's */
double rastrigin_func (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Rastrigin's  */
double weierstrass_func (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Weierstrass's  */
double griewank_func (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Griewank's  */
double schwefel_func (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Schwefel's */
double katsuura_func (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Katsuura */
double bi_rastrigin_func (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Lunacek Bi_rastrigin */
double grie_rosen_func (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Griewank-Rosenbrock  */
double escaffer6_func (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Expanded Scaffer¡¯s F6  */
double step_rastrigin_func (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Noncontinuous Rastrigin's  */
double happycat_func (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* HappyCat */
double hgbat_func (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* HGBat  */

/* New functions Noor Changes */
double sum_diff_pow_func(double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Sum of different power */
double zakharov_func(double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* ZAKHAROV */
double levy_func(double* var, int num, double* Os,double* Mr,int* S,int s_flag, int r_flag); /* Levy */
double dixon_price_func(double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Dixon and Price */

double hf01 (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Hybrid Function 1 */
double hf02 (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Hybrid Function 2 */
double hf03 (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_fla); /* Hybrid Function 3 */
double hf04 (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Hybrid Function 4 */
double hf05 (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Hybrid Function 5 */
double hf06 (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Hybrid Function 6 */
double hf07 (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Hybrid Function 7 */
double hf08 (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Hybrid Function 8 */
double hf09 (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Hybrid Function 9 */
double hf10 (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Hybrid Function 10 */

double cf01 (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Composition Function 1 */
double cf02 (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Composition Function 2 */
double cf03 (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Composition Function 3 */
double cf04 (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Composition Function 4 */
double cf05 (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Composition Function 5 */
double cf06 (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Composition Function 6 */
double cf07 (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Composition Function 7 */
double cf08 (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Composition Function 8 */
double cf09 (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Composition Function 9 */
double cf10 (double* var, int num, double* Os,double* Mr, int* S,int s_flag, int r_flag); /* Composition Function 10 */

void shiftfunc (double* var, int num, double* Os,double* xshift);
void rotatefunc (double* var, int num, double* Mr,double* xrot);
void sr_func(double * var, int num, double* Os, double* Mr, double sh_rate, int s_flag, int r_flag,double* sr_x);/* shift and rotate */
double cf_cal(double * var, int num, double * Os,double * delta,double * bias,double * fit,int cf_num);