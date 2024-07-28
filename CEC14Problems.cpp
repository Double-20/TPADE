#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <fstream>
#include <iostream>
#include "CEC14Problems.h"
#define INF 1.0e99
#define EPS 1.0e-14
#define E  2.7182818284590452353602874713526625
#define PI 3.1415926535897932384626433832795029

void shiftfunc (double* var,int num,double* OShift, double* xshift)
{
    int i;
    for (i=0; i<num; i++)
        xshift[i]=var[i]-OShift[i];
}

void rotatefunc (double* var, int num, double* M, double* xrot)
{
    int i,j;
    for (i=0; i<num; i++)
    {
        xrot[i]=0;
        for (j=0; j<num; j++)
        {
            xrot[i]=xrot[i]+var[j]*M[i*num+j];
        }
    }
}

void sr_func (double * var, int num, double* OShift, double* M, double sh_rate, int s_flag, int r_flag,double* sr_x)
{
    int i;
    double* y = new double[num];
    if (s_flag==1)
    {
        if (r_flag==1)
        {
            shiftfunc(var, num, OShift, y);
            for (i=0; i<num; i++)//shrink to the orginal search range
            {
                y[i]=y[i]*sh_rate;
            }
            rotatefunc(y, num, M,sr_x);
        }
        else
        {
            shiftfunc(var,num, OShift,sr_x);
            for (i=0; i<num; i++)//shrink to the orginal search range
            {
                sr_x[i]=sr_x[i]*sh_rate;
            }
        }
    }
    else
    {
        if (r_flag==1)
        {
            for (i=0; i<num; i++)//shrink to the orginal search range
            {
                y[i]=var[i]*sh_rate;
            }
            rotatefunc(y, num, M,sr_x);
        }
        else {
            for (i = 0; i < num; i++)//shrink to the orginal search range
            {
                sr_x[i] = var[i] * sh_rate;
            }
        }
    }
    delete[] y;
}

void asyfunc (double * var,int num, double beta, double* xasy)
{
    int i;
    for (i=0; i<num; i++)
    {
        if (var[i]>0)
            xasy[i]=pow(var[i],1.0+beta*i/(num-1)*pow(var[i],0.5));
    }
}

void oszfunc (double * var, int num, double* xOShiftz)
{
    int i,sx;
    double c1,c2,xx;
    for (i=0; i<num; i++)
    {
        if (i==0||i==num-1)
        {
            if (var[i]!=0)
                xx=log(fabs(var[i]));
            if (var[i]>0)
            {
                c1=10;
                c2=7.9;
            }
            else
            {
                c1=5.5;
                c2=3.1;
            }
            if (var[i]>0)
                sx=1;
            else if (var[i]==0)
                sx=0;
            else
                sx=-1;
            xOShiftz[i]=sx*exp(xx+0.049*(sin(c1*xx)+sin(c2*xx)));
        }
        else
            xOShiftz[i]=var[i];
    }
}


double cf_cal(double * var,int num, double * OShift,double * delta,double * bias,double * fit,int cf_num)
{
    int i,j;
    double *w = new double[cf_num];
    double w_max=0,w_sum=0;
    // w=(double *)malloc(cf_num * sizeof(double));
    for (i=0; i<cf_num; i++)
    {
        fit[i]+=bias[i];
        w[i]=0;
        for (j=0; j<num; j++)
        {
            w[i]+=pow(var[j]-OShift[i*num+j],2.0);
        }
        if (w[i]!=0)
            w[i]=pow(1.0/w[i],0.5)*exp(-w[i]/2.0/num/pow(delta[i],2.0));
        else
            w[i]=INF;
        if (w[i]>w_max)
            w_max=w[i];
    }

    for (i=0; i<cf_num; i++)
    {
        w_sum=w_sum+w[i];
    }
    if(w_max==0)
    {
        for (i=0; i<cf_num; i++)
            w[i]=1;
        w_sum=cf_num;
    }
    double f = 0.;
    for (i=0; i<cf_num; i++)
    {
        f=f+w[i]/w_sum*fit[i];
    }
    delete [] w;
    return f;
}

double sphere_func (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Sphere */
{
    int i;
    double f = 0.0;
    double* z = new double [num];
    sr_func (var, num, OShift, M, 1.0, s_flag, r_flag,z);  /* shift and rotate */
    for (i=0; i<num; i++)
    {
        f += z[i]*z[i];
    }
    delete [] z;
    return f;
}

double ellips_func (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Ellipsoidal */
{
    int i;
    double* z = new double[num];
    double f = 0.;
    sr_func (var, num, OShift, M,1.0, s_flag, r_flag,z); /* shift and rotate */
    for (i=0; i<num; i++)
    {
        f += pow(10.0,6.0*i/(num-1))*z[i]*z[i];
    }
    delete [] z;
    return f + 100;
}

double bent_cigar_func (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Bent_Cigar */
{
    int i;
    double* z = new double[num];
    sr_func (var, num, OShift, M,1.0, s_flag, r_flag,z); /* shift and rotate */
    double f = z[0]*z[0];
    for (i=1; i<num; i++)
    {
        f += pow(10.0,6.0)*z[i]*z[i];
    }
    delete [] z;
    return f + 200;
}

double discus_func (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Discus */
{
    int i;
    double* z = new double[num];
    sr_func (var, num, OShift, M,1.0, s_flag, r_flag,z); /* shift and rotate */
    double f = pow(10.0,6.0)*z[0]*z[0];
    for (i=1; i<num; i++)
    {
        f += z[i]*z[i];
    }
    delete [] z;
    return f + 300.;
}

double dif_powers_func (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Different Powers */
{
    int i;
    double* z = new double[num];
    double f = 0.;
    sr_func (var, num, OShift, M, 1.0, s_flag, r_flag,z);  /* shift and rotate */
    for (i=0; i<num; i++)
    {
        f += pow(fabs(z[i]),2+4*i/(num-1));
    }
    f=pow(f,0.5);
    delete [] z;
    return f;
}


double rosenbrock_func (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* ROShiftenbrock's */
{
    int i;
    double tmp1,tmp2;
    double f = 0.;
    double* z = new double[num];
    sr_func (var, num, OShift, M, 2.048/100.0, s_flag, r_flag,z); /* shift and rotate */
    z[0] += 1.0;//shift to orgin
    for (i=0; i<num-1; i++)
    {
        z[i+1] += 1.0;//shift to orgin
        tmp1=z[i]*z[i]-z[i+1];
        tmp2=z[i]-1.0;
        f += 100.0*tmp1*tmp1 +tmp2*tmp2;
    }
    delete [] z;
    return f + 400.;
}

double schaffer_F7_func (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Schwefel's 1.2  */
{
    int i;
    double tmp;
    double f = 0.;
    double* z = new double[num];
    sr_func (var, num, OShift, M, 1.0, s_flag, r_flag,z);  /* shift and rotate */
    for (i=0; i<num-1; i++)
    {
        z[i]=pow(z[i]*z[i]+z[i+1]*z[i+1],0.5);
        tmp=sin(50.0*pow(z[i],0.2));
        f += pow(z[i],0.5)+pow(z[i],0.5)*tmp*tmp ;
    }
    f = f*f/(num-1)/(num-1);
    delete [] z;
    return f;
}

double ackley_func (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Ackley's  */
{
    int i;
    double sum1, sum2;
    sum1 = 0.0;
    sum2 = 0.0;
    double* z = new double[num];
    sr_func (var, num, OShift, M, 1.0, s_flag, r_flag,z);  /* shift and rotate */

    for (i=0; i<num; i++)
    {
        sum1 += z[i]*z[i];
        sum2 += cos(2.0*PI*z[i]);
    }
    sum1 = -0.2*sqrt(sum1/num);
    sum2 /= num;
    delete [] z;
    return E - 20.0*exp(sum1) - exp(sum2) +20.0 + 500.0;
}


double weierstrass_func (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Weierstrass's  */
{
    int i,j,k_max;
    double sum,sum2, a, b;
    a = 0.5;
    b = 3.0;
    k_max = 20;
    double f = 0.;
    double* z = new double[num];
    sr_func (var, num, OShift, M, 0.5/100.0, s_flag, r_flag, z); /* shift and rotate */

    for (i=0; i<num; i++)
    {
        sum = 0.0;
        sum2 = 0.0;
        for (j=0; j<=k_max; j++)
        {
            sum += pow(a,j)*cos(2.0*PI*pow(b,j)*(z[i]+0.5));
            sum2 += pow(a,j)*cos(2.0*PI*pow(b,j)*0.5);
        }
        f += sum;
    }
    f -= num*sum2;
    delete [] z;
    return f + 600.;
}


double griewank_func (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Griewank's  */
{
    int i;
    double s, p;
    s = 0.0;
    p = 1.0;
    double* z = new double[num];
    sr_func (var, num, OShift, M, 600.0/100.0, s_flag, r_flag, z); /* shift and rotate */

    for (i=0; i<num; i++)
    {
        s += z[i]*z[i];
        p *= cos(z[i]/sqrt(1.0+i));
    }
    delete [] z;
    return 1.0 + s/4000.0 - p + 700.0;
}

double rastrigin_func (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Rastrigin's  */
{
//    string pathname = "../inputdata.txt";
//    ifstream infile;
//    infile.open(pathname.10D());
//    for (int j = 0; j < num; j++) {
//        infile >> var[j];
//    }
//    infile.close();
    int i;
    double f = 0.;
    double* z = new double[num];
    sr_func (var, num, OShift, M, 5.12/100.0, s_flag, r_flag,z); /* shift and rotate */

    for (i=0; i<num; i++)
    {
        f += (z[i]*z[i] - 10.0*cos(2.0*PI*z[i]) + 10.0);
    }
    delete [] z;
    if(r_flag == 0)
        return f + 800.0;
    else
        return f + 900.0;
}

double step_rastrigin_func (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Noncontinuous Rastrigin's  */
{
//    int i;
//    double f=0.0;
//    double* z = new double[num];
//    for (i=0; i<num; i++)
//    {
//        if (fabs(y[i]-OShift[i])>0.5)
//            y[i]=OShift[i]+floor(2*(y[i]-OShift[i])+0.5)/2;
//    }
//
//    sr_func (var, num, OShift, M, 5.12/100.0, s_flag, r_flag,z); /* shift and rotate */
//
//    for (i=0; i<num; i++)
//    {
//        f += (z[i]*z[i] - 10.0*cos(2.0*PI*z[i]) + 10.0);
//    }
//    delete [] z;
//    return f;
}

double schwefel_func (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Schwefel's  */
{
    int i;
    double tmp;
    double f=0.0;
    double* z = new double[num];
    sr_func (var, num, OShift, M, 1000.0/100.0, s_flag, r_flag,z); /* shift and rotate */

    for (i=0; i<num; i++)
    {
        z[i] += 4.209687462275036e+002;
        if (z[i]>500)
        {
            f-=(500.0-fmod(z[i],500))*sin(pow(500.0-fmod(z[i],500),0.5));
            tmp=(z[i]-500.0)/100;
            f+= tmp*tmp/num;
        }
        else if (z[i]<-500)
        {
            f-=(-500.0+fmod(fabs(z[i]),500))*sin(pow(500.0-fmod(fabs(z[i]),500),0.5));
            tmp=(z[i]+500.0)/100;
            f+= tmp*tmp/num;
        }
        else
            f-=z[i]*sin(pow(fabs(z[i]),0.5));
        if(f < 1e-20)
            f = f;
    }
    f +=4.189828872724338e+002*num;
    delete [] z;
    if(r_flag == 0)
        return f + 1000.0;
    else
        return f + 1100.0;
}

double katsuura_func (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Katsuura  */
{
    int i,j;
    double temp,tmp1,tmp2,tmp3;
    double f=1.0;
    tmp3=pow(1.0*num,1.2);
    double* z = new double[num];
    sr_func (var, num, OShift, M, 5.0/100.0, s_flag, r_flag,z); /* shift and rotate */

    for (i=0; i<num; i++)
    {
        temp=0.0;
        for (j=1; j<=32; j++)
        {
            tmp1=pow(2.0,j);
            tmp2=tmp1*z[i];
            temp += fabs(tmp2-floor(tmp2+0.5))/tmp1;
        }
        f *= pow(1.0+(i+1)*temp,10.0/tmp3);
    }
    tmp1=10.0/num/num;
    f=f*tmp1-tmp1;
    delete [] z;
    return f + 1200.0;
}

double bi_rastrigin_func (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Lunacek Bi_rastrigin Function */
{
    int i;
    double mu0=2.5,d=1.0,s,mu1,tmp,tmp1,tmp2;
    double *tmpx = new double[num];
    s=1.0-1.0/(2.0*pow(num+20.0,0.5)-8.2);
    mu1=-pow((mu0*mu0-d)/s,0.5);
    double* y = new double[num];
    double* z = new double[num];
    if (s_flag==1)
        shiftfunc(var, num, OShift, y);
    else
    {
        for (i=0; i<num; i++)//shrink to the orginal search range
        {
            y[i] = var[i];
        }
    }
    for (i=0; i<num; i++)//shrink to the orginal search range
    {
        y[i] *= 10.0/100.0;
    }

    for (i = 0; i < num; i++)
    {
        tmpx[i]=2*y[i];
        if (OShift[i] < 0.0)
            tmpx[i] *= -1.;
    }
    for (i=0; i<num; i++)
    {
        z[i]=tmpx[i];
        tmpx[i] += mu0;
    }
    tmp1=0.0;tmp2=0.0;
    for (i=0; i<num; i++)
    {
        tmp = tmpx[i]-mu0;
        tmp1 += tmp*tmp;
        tmp = tmpx[i]-mu1;
        tmp2 += tmp*tmp;
    }
    tmp2 *= s;
    tmp2 += d*num;
    tmp=0.0;
    double f = 0.;
    if (r_flag==1)
    {
        rotatefunc(z, num, M, y);
        for (i=0; i<num; i++)
        {
            tmp+=cos(2.0*PI*y[i]);
        }
        if(tmp1<tmp2)
            f = tmp1;
        else
            f = tmp2;
        f += 10.0*(num-tmp);
    }
    else
    {
        for (i=0; i<num; i++)
        {
            tmp+=cos(2.0*PI*z[i]);
        }
        if(tmp1<tmp2)
            f = tmp1;
        else
            f = tmp2;
        f += 10.0*(num-tmp);
    }

    delete [] tmpx;
    delete [] z;
    delete [] y;
    return f;
}

double grie_rosen_func (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Griewank-ROShiftenbrock  */
{
    int i;
    double temp,tmp1,tmp2;
    double f=0.0;
    double* z = new double[num];
    sr_func (var, num, OShift, M, 5.0/100.0, s_flag, r_flag,z); /* shift and rotate */

    z[0] += 1.0;//shift to orgin
    for (i=0; i<num-1; i++)
    {
        z[i+1] += 1.0;//shift to orgin
        tmp1 = z[i]*z[i]-z[i+1];
        tmp2 = z[i]-1.0;
        temp = 100.0*tmp1*tmp1 + tmp2*tmp2;
        f += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
    tmp1 = z[num-1]*z[num-1]-z[0];
    tmp2 = z[num-1]-1.0;
    temp = 100.0*tmp1*tmp1 + tmp2*tmp2;;
    f += (temp*temp)/4000.0 - cos(temp) + 1.0 ;
    delete [] z;
    return f + 1500.0;
}

double escaffer6_func (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Expanded Scaffer��s F6  */
{
    int i;
    double temp1, temp2;
    double* z = new double[num];
    sr_func (var, num, OShift, M, 1.0, s_flag, r_flag,z);  /* shift and rotate */

    double f = 0.;
    for (i=0; i<num-1; i++)
    {
        temp1 = sin(sqrt(z[i]*z[i]+z[i+1]*z[i+1]));
        temp1 =temp1*temp1;
        temp2 = 1.0 + 0.001*(z[i]*z[i]+z[i+1]*z[i+1]);
        f += 0.5 + (temp1-0.5)/(temp2*temp2);
    }
    temp1 = sin(sqrt(z[num-1]*z[num-1]+z[0]*z[0]));
    temp1 =temp1*temp1;
    temp2 = 1.0 + 0.001*(z[num-1]*z[num-1]+z[0]*z[0]);
    f += 0.5 + (temp1-0.5)/(temp2*temp2);
    delete [] z;
    return f + 1600.0;
}

double happycat_func (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* HappyCat, provdided by Hans-Georg Beyer (HGB) */
/* original global optimum: [-1,-1,...,-1] */
{
    int i;
    double alpha,r2,sum_z;
    alpha=1.0/8.0;
    double* z = new double[num];
    sr_func (var, num, OShift, M, 5.0/100.0, s_flag, r_flag,z); /* shift and rotate */

    r2 = 0.0;
    sum_z=0.0;
    for (i=0; i<num; i++)
    {
        z[i]=z[i]-1.0;//shift to orgin
        r2 += z[i]*z[i];
        sum_z += z[i];
    }
    delete[]z;
    return pow(fabs(r2-num),2*alpha) + (0.5*r2 + sum_z)/num + 0.5 + 1300.0;
}

double hgbat_func (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* HGBat, provdided by Hans-Georg Beyer (HGB)*/
/* original global optimum: [-1,-1,...,-1] */
{
    int i;
    double alpha,r2,sum_z;
    alpha=1.0/4.0;
    double* z = new double[num];
    sr_func (var, num, OShift, M, 5.0/100.0, s_flag, r_flag,z); /* shift and rotate */

    r2 = 0.0;
    sum_z=0.0;
    for (i=0; i<num; i++)
    {
        z[i]=z[i]-1.0;//shift to orgin
        r2 += z[i]*z[i];
        sum_z += z[i];
    }
    delete[]z;
    return pow(fabs(pow(r2,2.0)-pow(sum_z,2.0)),2*alpha) + (0.5*r2 + sum_z)/num + 0.5 + 1400.0;
}

double cf01 (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* CompOShiftition Function 1 */
{
    int i,cf_num=5;
    double fit[5];
    double delta[5] = {10, 20, 30, 40, 50};
    double bias[5] = {0, 100, 200, 300, 400};

    i=0;
    fit[i] = rosenbrock_func(var,num,&OShift[i*num],&M[i*num*num],S,1,r_flag);
    fit[i] -= 400.;
    fit[i]=10000*fit[i]/1e+4;
    i=1;
    fit[i] = ellips_func(var,num,&OShift[i*num],&M[i*num*num],S,1,r_flag);
    fit[i] -= 100;
    fit[i]=10000*fit[i]/1e+10;
    i=2;
    fit[i] = bent_cigar_func(var,num,&OShift[i*num],&M[i*num*num],S,1,r_flag);
    fit[i] -= + 200;
    fit[i]=10000*fit[i]/1e+30;
    i=3;
    fit[i] = discus_func(var,num,&OShift[i*num],&M[i*num*num],S,1,r_flag);
    fit[i] -= 300.;
    fit[i]=10000*fit[i]/1e+10;
    i=4;
    fit[i] = ellips_func(var,num,&OShift[i*num],&M[i*num*num],S,1,0);
    fit[i] -= 100;
    fit[i]=10000*fit[i]/1e+10;
    return cf_cal(var, num, OShift, delta,bias,fit,cf_num) + 2300.0;
}

double cf02 (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* CompOShiftition Function 2 */
{
    int i,cf_num=3;
    double fit[3];
    double delta[3] = {20,20,20};
    double bias[3] = {0, 100, 200};

    i=0;
    fit[i] = schwefel_func(var,num,&OShift[i*num],&M[i*num*num],S,1,0);
    fit[i] -= 1000.0;
    i=1;
    fit[i] = rastrigin_func(var,num,&OShift[i*num],&M[i*num*num],S,1,r_flag);
    fit[i] -= 900.0;
    i=2;
    fit[i] = hgbat_func(var,num,&OShift[i*num],&M[i*num*num],S,1,r_flag);
    fit[i] -= 1400.0;
    return cf_cal(var, num, OShift, delta,bias,fit,cf_num)+ 2400.0;
}

double cf03 (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* CompOShiftition Function 3 */
{
    int i,cf_num=3;
    double fit[3];
    double delta[3] = {10,30,50};
    double bias[3] = {0, 100, 200};
    i=0;
    fit[i] = schwefel_func(var,num,&OShift[i*num],&M[i*num*num],S,1,r_flag);
    fit[i] -= 1100.0;
    fit[i]=1000*fit[i]/4e+3;
    i=1;
    fit[i] = rastrigin_func(var,num,&OShift[i*num],&M[i*num*num],S,1,r_flag);
    fit[i] -= 900.0;
    fit[i]=1000*fit[i]/1e+3;
    i=2;
    fit[i] = ellips_func(var,num,&OShift[i*num],&M[i*num*num],S,1,r_flag);
    fit[i] -= 100;
    fit[i]=1000*fit[i]/1e+10;
    return cf_cal(var, num, OShift, delta,bias,fit,cf_num)+ 2500.0;
}

double cf04 (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* CompOShiftition Function 4 */
{
    int i,cf_num=5;
    double fit[5];
    double delta[5] = {10,10,10,10,10};
    double bias[5] = {0, 100, 200, 300, 400};
    i=0;
    fit[i] = schwefel_func(var,num,&OShift[i*num],&M[i*num*num],S,1,r_flag);
    fit[i] -= 1100.0;
    fit[i]=1000*fit[i]/4e+3;
    i=1;
    fit[i] = happycat_func(var,num,&OShift[i*num],&M[i*num*num],S,1,r_flag);
    fit[i] -= 1300.0;
    fit[i]=1000*fit[i]/1e+3;
    i=2;
    fit[i] = ellips_func(var,num,&OShift[i*num],&M[i*num*num],S,1,r_flag);
    fit[i] -= 100;
    fit[i]=1000*fit[i]/1e+10;
    i=3;
    fit[i] = weierstrass_func(var,num,&OShift[i*num],&M[i*num*num],S,1,r_flag);
    fit[i] -= 600.;
    fit[i]=1000*fit[i]/400;
    i=4;
    fit[i] = griewank_func(var,num,&OShift[i*num],&M[i*num*num],S,1,r_flag);
    fit[i] -= 700.0;
    fit[i]=1000*fit[i]/100;
    return cf_cal(var, num, OShift, delta,bias,fit,cf_num)+ 2600.0;
}

double cf05 (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* CompOShiftition Function 4 */
{
    int i,cf_num=5;
    double fit[5];
    double delta[5] = {10,10,10,20,20};
    double bias[5] = {0, 100, 200, 300, 400};
    i=0;
    fit[i] = hgbat_func(var,num,&OShift[i*num],&M[i*num*num],S,1,r_flag);
    fit[i] -= 1400.0;
    fit[i]=10000*fit[i]/1000;
    i=1;
    fit[i] = rastrigin_func(var,num,&OShift[i*num],&M[i*num*num],S,1,r_flag);
    fit[i] -= 900.0;
    fit[i]=10000*fit[i]/1e+3;
    i=2;
    fit[i] = schwefel_func(var,num,&OShift[i*num],&M[i*num*num],S,1,r_flag);
    fit[i] -= 1100.0;
    fit[i]=10000*fit[i]/4e+3;
    i=3;
    fit[i] = weierstrass_func(var,num,&OShift[i*num],&M[i*num*num],S,1,r_flag);
    fit[i] -= 600.;
    fit[i]=10000*fit[i]/400;
    i=4;
    fit[i] = ellips_func(var,num,&OShift[i*num],&M[i*num*num],S,1,r_flag);
    fit[i] -= 100;
    fit[i]=10000*fit[i]/1e+10;
    return cf_cal(var, num, OShift, delta,bias,fit,cf_num)+ 2700.0;
}

double cf06 (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* CompOShiftition Function 6 */
{
    int i,cf_num=5;
    double fit[5];
    double delta[5] = {10,20,30,40,50};
    double bias[5] = {0, 100, 200, 300, 400};
    i=0;
    fit[i] = grie_rosen_func(var,num,&OShift[i*num],&M[i*num*num],S,1,r_flag);
    fit[i] -= 1500.0;
    fit[i]=10000*fit[i]/4e+3;
    i=1;
    fit[i] = happycat_func(var,num,&OShift[i*num],&M[i*num*num],S,1,r_flag);
    fit[i] -= 1300.0;
    fit[i]=10000*fit[i]/1e+3;
    i=2;
    fit[i] = schwefel_func(var,num,&OShift[i*num],&M[i*num*num],S,1,r_flag);
    fit[i] -= 1100.0;
    fit[i]=10000*fit[i]/4e+3;
    i=3;
    fit[i] = escaffer6_func(var,num,&OShift[i*num],&M[i*num*num],S,1,r_flag);
    fit[i] -= 1600.0;
    fit[i]=10000*fit[i]/2e+7;
    i=4;
    fit[i] = ellips_func(var,num,&OShift[i*num],&M[i*num*num],S,1,r_flag);
    fit[i] -= 100.0;
    fit[i]=10000*fit[i]/1e+10;
    return cf_cal(var, num, OShift, delta,bias,fit,cf_num)+ 2800.0;
}

double cf07 (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* CompOShiftition Function 7 */
{
    int i,cf_num=3;
    double fit[3];
    double delta[3] = {10,30,50};
    double bias[3] = {0, 100, 200};
    i=0;
    fit[i] = hf01(var,num,&OShift[i*num],&M[i*num*num],&S[i*num],1,r_flag);
    fit[i] -= 1700.0;
    i=1;
    fit[i] = hf02(var,num,&OShift[i*num],&M[i*num*num],&S[i*num],1,r_flag);
    fit[i] -= 1800;
    i=2;
    fit[i] = hf03(var,num,&OShift[i*num],&M[i*num*num],&S[i*num],1,r_flag);
    fit[i] -= 1900.0;
    return cf_cal(var, num, OShift, delta,bias,fit,cf_num)+ 2900.0;
}

double cf08 (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* CompOShiftition Function 8 */
{
//    string pathname = "../inputdata.txt";
//    ifstream infile;
//    infile.open(pathname.10D());
//    for (int j = 0; j < num; j++) {
//        infile >> var[j];
//    }
//    infile.close();
    int i,cf_num=3;
    double fit[3];
    double delta[3] = {10,30,50};
    double bias[3] = {0, 100, 200};
    i=0;
    fit[i] = hf04(var,num,&OShift[i*num],&M[i*num*num],&S[i*num],1,r_flag);
    fit[i] -= 2000.0;
    i=1;
    fit[i] = hf05(var,num,&OShift[i*num],&M[i*num*num],&S[i*num],1,r_flag);
    fit[i] -= 2100.0;
    i=2;
    fit[i] = hf06(var,num,&OShift[i*num],&M[i*num*num],&S[i*num],1,r_flag);
    fit[i] -= 2200.0;
    return cf_cal(var, num, OShift, delta,bias,fit,cf_num)+ 3000.0;
}


double hf01 (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Hybrid Function 1 */
{
    int i,tmp,cf_num=3;
    double fit[3];
    int G[3],G_num[3];
    double Gp[3]={0.3,0.3,0.4};
    double* z = new double[num];
    double* y = new double[num];
    tmp=0;
    for (i=0; i<cf_num-1; i++)
    {
        G_num[i] = ceil(Gp[i]*num);
        tmp += G_num[i];
    }
    G_num[cf_num-1]=num-tmp;
    G[0]=0;
    for (i=1; i<cf_num; i++)
    {
        G[i] = G[i-1]+G_num[i-1];
    }

    sr_func (var, num, OShift, M, 1.0, s_flag, r_flag,z);  /* shift and rotate */

    for (i=0; i<num; i++)
    {
        y[i]=z[S[i]-1];
    }
    i=0;
    fit[i] = schwefel_func(&y[G[i]],G_num[i],OShift,M,S,0,0);
    fit[i] -= 1000.0;
    i=1;
    fit[i] = rastrigin_func(&y[G[i]],G_num[i],OShift,M,S,0,0);
    fit[i] -= 800.0;
    i=2;
    fit[i] = ellips_func(&y[G[i]],G_num[i],OShift,M,S,0,0);
    fit[i] -= 100.0;
    double f=0.0;
    for(i=0;i<cf_num;i++)
    {
        f += fit[i];
    }
    delete [] y;
    delete [] z;
    return f + 1700.0;
}

double hf02 (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Hybrid Function 2 */
{

//    string pathname = "../inputdata.txt";
//    ifstream infile;
//    infile.open(pathname.10D());
//    for (int j = 0; j < num; j++) {
//        infile >> var[j];
//    }
//    infile.close();

    int i,tmp,cf_num=3;
    double fit[3];
    int G[3],G_num[3];
    double Gp[3]={0.3,0.3,0.4};
    double* z = new double[num];
    double* y = new double[num];
    tmp=0;
    for (i=0; i<cf_num-1; i++)
    {
        G_num[i] = ceil(Gp[i]*num);
        tmp += G_num[i];
    }
    G_num[cf_num-1]=num-tmp;

    G[0]=0;
    for (i=1; i<cf_num; i++)
    {
        G[i] = G[i-1]+G_num[i-1];
    }

    sr_func (var, num, OShift, M, 1.0, s_flag, r_flag,z);  /* shift and rotate */

    for (i=0; i<num; i++)
    {
        y[i]=z[S[i]-1];
    }
    i=0;
    fit[i] = bent_cigar_func(&y[G[i]],G_num[i],OShift,M,S,0,0);
    fit[i] -= 200.0;
    i=1;
    fit[i] = hgbat_func(&y[G[i]],G_num[i],OShift,M,S,0,0);
    fit[i] -= 1400.0;
    i=2;
    fit[i] = rastrigin_func(&y[G[i]],G_num[i],OShift,M,S,0,0);
    fit[i] -= 800.0;
    double f=0.0;
    for(i=0;i<cf_num;i++)
    {
        f += fit[i];
    }
    delete [] y;
    delete [] z;
    return f + 1800;
}

double hf03 (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Hybrid Function 3 */
{
    int i,tmp,cf_num=4;
    double fit[4];
    int G[4],G_num[4];
    double Gp[4]={0.2,0.2,0.3,0.3};
    double* z = new double[num];
    double* y = new double[num];
    tmp=0;
    for (i=0; i<cf_num-1; i++)
    {
        G_num[i] = ceil(Gp[i]*num);
        tmp += G_num[i];
    }
    G_num[cf_num-1]=num-tmp;

    G[0]=0;
    for (i=1; i<cf_num; i++)
    {
        G[i] = G[i-1]+G_num[i-1];
    }

    sr_func (var, num, OShift, M, 1.0, s_flag, r_flag,z);  /* shift and rotate */

    for (i=0; i<num; i++)
    {
        y[i]=z[S[i]-1];
    }
    i=0;
    fit[i] = griewank_func(&y[G[i]],G_num[i],OShift,M,S,0,0);
    fit[i] -= 700.0;
    i=1;
    fit[i] = weierstrass_func(&y[G[i]],G_num[i],OShift,M,S,0,0);
    fit[i] -= 600.0;
    i=2;
    fit[i] = rosenbrock_func(&y[G[i]],G_num[i],OShift,M,S,0,0);
    fit[i] -= 400.0;
    i=3;
    fit[i] = escaffer6_func(&y[G[i]],G_num[i],OShift,M,S,0,0);
    fit[i] -= 1600.0;
    double f=0.0;
    for(i=0;i<cf_num;i++)
    {
        f += fit[i];
    }
    delete [] y;
    delete [] z;
    return f + 1900.0;
}

double hf04 (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Hybrid Function 4 */
{
    int i,tmp,cf_num=4;
    double fit[4];
    int G[4],G_num[4];
    double Gp[4]={0.2,0.2,0.3,0.3};
    double* z = new double[num];
    double* y = new double[num];
    tmp=0;
    for (i=0; i<cf_num-1; i++)
    {
        G_num[i] = ceil(Gp[i]*num);
        tmp += G_num[i];
    }
    G_num[cf_num-1]=num-tmp;

    G[0]=0;
    for (i=1; i<cf_num; i++)
    {
        G[i] = G[i-1]+G_num[i-1];
    }

    sr_func (var, num, OShift, M, 1.0, s_flag, r_flag,z);  /* shift and rotate */

    for (i=0; i<num; i++)
    {
        y[i]=z[S[i]-1];
    }
    i=0;
    fit[i] = hgbat_func(&y[G[i]],G_num[i],OShift,M,S,0,0);
    fit[i] -= 1400.0;
    i=1;
    fit[i] = discus_func(&y[G[i]],G_num[i],OShift,M,S,0,0);
    fit[i] -= 300.;
    i=2;
    fit[i] = grie_rosen_func(&y[G[i]],G_num[i],OShift,M,S,0,0);
    fit[i] -= 1500.0;
    i=3;
    fit[i] = rastrigin_func(&y[G[i]],G_num[i],OShift,M,S,0,0);
    fit[i] -= 800.0;

    double f=0.0;
    for(i=0;i<cf_num;i++)
    {
        f += fit[i];
    }
    delete [] y;
    delete [] z;
    return f + 2000.0;
}
double hf05 (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Hybrid Function 5 */
{
    int i,tmp,cf_num=5;
    double fit[5];
    int G[5],G_num[5];
    double Gp[5]={0.1,0.2,0.2,0.2,0.3};
    double* z = new double[num];
    double* y = new double[num];
    tmp=0;
    for (i=0; i<cf_num-1; i++)
    {
        G_num[i] = ceil(Gp[i]*num);
        tmp += G_num[i];
    }
    G_num[cf_num-1]=num-tmp;

    G[0]=0;
    for (i=1; i<cf_num; i++)
    {
        G[i] = G[i-1]+G_num[i-1];
    }

    sr_func (var, num, OShift, M, 1.0, s_flag, r_flag,z);  /* shift and rotate */

    for (i=0; i<num; i++)
    {
        y[i]=z[S[i]-1];
    }
    i=0;
    fit[i] = escaffer6_func(&y[G[i]],G_num[i],OShift,M,S,0,0);
    fit[i] -=1600.0;
    i=1;
    fit[i] = hgbat_func(&y[G[i]],G_num[i],OShift,M,S,0,0);
    fit[i] -=1400.0;
    i=2;
    fit[i] = rosenbrock_func(&y[G[i]],G_num[i],OShift,M,S,0,0);
    fit[i] -= 400.;
    i=3;
    fit[i] = schwefel_func(&y[G[i]],G_num[i],OShift,M,S,0,0);
    fit[i] -= 1000.0;
    i=4;
    fit[i] = ellips_func(&y[G[i]],G_num[i],OShift,M,S,0,0);
    fit[i] -= 100;

            double f=0.0;
    for(i=0;i<cf_num;i++)
    {
        f += fit[i];
    }
    delete [] y;
    delete [] z;
    return f + 2100.0;
}

double hf06 (double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Hybrid Function 6 */
{
    int i,tmp,cf_num=5;
    double fit[5];
    int G[5],G_num[5];
    double Gp[5]={0.1,0.2,0.2,0.2,0.3};
    double* z = new double[num];
    double* y = new double[num];
    tmp=0;
    for (i=0; i<cf_num-1; i++)
    {
        G_num[i] = ceil(Gp[i]*num);
        tmp += G_num[i];
    }
    G_num[cf_num-1]=num-tmp;

    G[0]=0;
    for (i=1; i<cf_num; i++)
    {
        G[i] = G[i-1]+G_num[i-1];
    }

    sr_func (var, num, OShift, M, 1.0, s_flag, r_flag,z);  /* shift and rotate */

    for (i=0; i<num; i++)
    {
        y[i]=z[S[i]-1];
    }
    i=0;
    fit[i] = katsuura_func(&y[G[i]],G_num[i],OShift,M,S,0,0);
    fit[i] -= 1200.0;
    i=1;
    fit[i] = happycat_func(&y[G[i]],G_num[i],OShift,M,S,0,0);
    fit[i] -= 1300.0;
    i=2;
    fit[i] =  grie_rosen_func(&y[G[i]],G_num[i],OShift,M,S,0,0);
    fit[i] -= 1500.0;
    i=3;
    fit[i] = schwefel_func(&y[G[i]],G_num[i],OShift,M,S,0,0);
    fit[i] -= 1000.0;
    i=4;
    fit[i] = ackley_func(&y[G[i]],G_num[i],OShift,M,S,0,0);
    fit[i] -= 500.0;

    double f=0.0;
    for(i=0;i<cf_num;i++)
    {
        f += fit[i];
    }
    delete [] y;
    delete [] z;
    return f + 2200.0;
}
