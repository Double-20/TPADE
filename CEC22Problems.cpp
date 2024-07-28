#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include "CEC22Problems.h"
#define INF 1.0e99
#define EPS 1.0e-14
#define E  2.7182818284590452353602874713526625
#define PI 3.1415926535897932384626433832795029
//(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag)
double ellips_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Ellipsoidal */
{

    int i;
    double f = 0.;
    double* z = new double[num];
    sr_func(var, z, num, OShift, M, 1.0, s_flag, r_flag); /* shift and rotate */
    for (i = 0; i < num; i++)
    {
        f += pow(10.0, 6.0 * i / (num - 1)) * z[i] * z[i];
    }
    delete[] z;
    return f;
}

double bent_cigar_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Bent_Cigar */
{
    int i;
    double f = 0.;
    double* z = new double[num];
    sr_func(var, z, num, OShift, M, 1.0, s_flag, r_flag); /* shift and rotate */
    f = z[0] * z[0];
    for (i = 1; i < num; i++)
    {
        f += pow(10.0, 6.0) * z[i] * z[i];
    }
    delete[] z;
    return f;
}

double discus_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Discus */
{
    int i;
    double* z = new double[num];
    sr_func(var, z, num, OShift, M, 1.0, s_flag, r_flag); /* shift and rotate */
    double f = pow(10.0, 6.0) * z[0] * z[0];
    for (i = 1; i < num; i++)
    {
        f += z[i] * z[i];
    }
    delete[] z;
    return f;
}

double rosenbrock_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Rosenbrock's */
{
//    ifstream inputFile("../testdata.txt", ios::binary);
//
//    if (!inputFile.is_open()) {
//        cerr << "无法打开文件" << endl;
//        return 1;
//    }
//    double temp;
//    int size = 0;
//
//    while (inputFile >> temp) {
//        double* newVar = new double[size+1]; // 创建一个新的更大的数组
//        for (int i = 0; i < size; ++i) {
//            newVar[i] = var[i]; // 将旧数据复制到新数组
//        }
//        newVar[size] = temp; // 添加新数据到新数组
//        delete[] var; // 释放旧数组内存
//        var = newVar; // 更新指针
//        size++;
//    }
//
//    inputFile.close();
//
//    cout << "从文件中读取的内容为：" << endl;
//    for (int i = 0; i < size; ++i) {
//        cout << var[i] << endl;
//    }

    int i;
    double tmp1, tmp2;
    double f = 0.;
    double* z = new double[num];
    sr_func(var, z, num, OShift, M, 2.048 / 100.0, s_flag, r_flag); /* shift and rotate */
    z[0] += 1.0; // shift to orgin
    for (i = 0; i < num - 1; i++)
    {
        z[i + 1] += 1.0; // shift to orgin
        tmp1 = z[i] * z[i] - z[i + 1];
        tmp2 = z[i] - 1.0;
        f += 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
    }
    delete[] z;
//    cout<<f<<endl;
    return f + 400.0;
}

double ackley_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Ackley's  */
{
    int i;
    double sum1, sum2;
    sum1 = 0.0;
    sum2 = 0.0;
    double f = 0.;
    double* z = new double[num];
    sr_func(var, z, num, OShift, M, 1.0, s_flag, r_flag); /* shift and rotate */

    for (i = 0; i < num; i++)
    {
        sum1 += z[i] * z[i];
        sum2 += cos(2.0 * PI * z[i]);
    }
    sum1 = -0.2 * sqrt(sum1 /num);
    sum2 /= num;
    f = E - 20.0 * exp(sum1) - exp(sum2) + 20.0;
    delete[] z;
    return f;
}

double griewank_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Griewank's  */
{
    int i;
    double s, p;
    s = 0.0;
    p = 1.0;
    double f = 0.;
    double* z = new double[num];
    sr_func(var, z, num, OShift, M, 600.0 / 100.0, s_flag, r_flag); /* shift and rotate */

    for (i = 0; i < num; i++)
    {
        s += z[i] * z[i];
        p *= cos(z[i] / sqrt(1.0 + i));
    }
    f = 1.0 + s / 4000.0 - p;
    delete[] z;
    return f;
}

double rastrigin_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Rastrigin's  */
{
    int i;
    double f = 0.;
    double* z = new double[num];
    sr_func(var, z, num, OShift, M, 5.12 / 100.0, s_flag, r_flag); /* shift and rotate */

    for (i = 0; i < num; i++)
    {
        f += (z[i] * z[i] - 10.0 * cos(2.0 * PI * z[i]) + 10.0);
    }
    delete[] z;
    return f;
}


double schwefel_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Schwefel's  */
{
    int i;
    double tmp;
    double f = 0.;
    double* z = new double[num];
    sr_func(var, z, num, OShift, M, 1000.0 / 100.0, s_flag, r_flag); /* shift and rotate */

    for (i = 0; i < num; i++)
    {
        z[i] += 4.209687462275036e+002;
        if (z[i] > 500)
        {
            f -= (500.0 - fmod(z[i], 500)) * sin(pow(500.0 - fmod(z[i], 500), 0.5));
            tmp = (z[i] - 500.0) / 100;
            f += tmp * tmp / num;
        }
        else if (z[i] < -500)
        {
            f -= (-500.0 + fmod(fabs(z[i]), 500)) * sin(pow(500.0 - fmod(fabs(z[i]), 500), 0.5));
            tmp = (z[i] + 500.0) / 100;
            f += tmp * tmp / num;
        }
        else
            f -= z[i] * sin(pow(fabs(z[i]), 0.5));
    }
    f += 4.189828872724338e+002 * num;
    delete[] z;
    return f;
}

double grie_rosen_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Griewank-Rosenbrock  */
{
    int i;
    double temp, tmp1, tmp2;
    double f = 0.;
    double* z = new double[num];
    sr_func(var, z, num, OShift, M, 5.0 / 100.0, s_flag, r_flag); /* shift and rotate */

    z[0] += 1.0; // shift to origin
    for (i = 0; i < num - 1; i++)
    {
        z[i + 1] += 1.0; // shift to origin
        tmp1 = z[i] * z[i] - z[i + 1];
        tmp2 = z[i] - 1.0;
        temp = 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
        f += (temp * temp) / 4000.0 - cos(temp) + 1.0;
    }
    tmp1 = z[num - 1] * z[num - 1] - z[0];
    tmp2 = z[num - 1] - 1.0;
    temp = 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
    f += (temp * temp) / 4000.0 - cos(temp) + 1.0;
    delete[] z;
    return f;
}


double escaffer6_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Evarpanded Scaffer's F6  */
{
    int i;
    double temp1, temp2;
    double* z = new double[num];
    sr_func(var, z, num, OShift, M, 1.0, s_flag, r_flag); /* shift and rotate */

    double f = 0.;
    for (i = 0; i < num - 1; i++)
    {
        temp1 = sin(sqrt(z[i] * z[i] + z[i + 1] * z[i + 1]));
        temp1 = temp1 * temp1;
        temp2 = 1.0 + 0.001 * (z[i] * z[i] + z[i + 1] * z[i + 1]);
        f += 0.5 + (temp1 - 0.5) / (temp2 * temp2);
    }
    temp1 = sin(sqrt(z[num - 1] * z[num - 1] + z[0] * z[0]));
    temp1 = temp1 * temp1;
    temp2 = 1.0 + 0.001 * (z[num - 1] * z[num - 1] + z[0] * z[0]);
    f += 0.5 + (temp1 - 0.5) / (temp2 * temp2);
    delete[] z;
    return f +600.0;
}

double happycat_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* HappyCat, provided by Hans-Georg Beyer (HGB) */
/* original global optimum: [-1,-1,...,-1] */
{
    int i;
    double alpha, r2, sum_z;
    alpha = 1.0 / 8.0;
    double f = 0.;
    double* z = new double[num];
    sr_func(var, z, num, OShift, M, 5.0 / 100.0, s_flag, r_flag); /* shift and rotate */

    r2 = 0.0;
    sum_z = 0.0;
    for (i = 0; i < num; i++)
    {
        z[i] = z[i] - 1.0; // shift to origin
        r2 += z[i] * z[i];
        sum_z += z[i];
    }
    f = pow(fabs(r2 - num), 2 * alpha) + (0.5 * r2 + sum_z) / num + 0.5;
    delete[] z;
    return f;
}

double hgbat_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* HGBat, provided by Hans-Georg Beyer (HGB) */
/* original global optimum: [-1,-1,...,-1] */
{
    int i;
    double alpha, r2, sum_z;
    alpha = 1.0 / 4.0;
    double f = 0.;
    double* z = new double[num];
    sr_func(var, z, num, OShift, M, 5.0 / 100.0, s_flag, r_flag); /* shift and rotate */

    r2 = 0.0;
    sum_z = 0.0;
    for (i = 0; i < num; i++)
    {
        z[i] = z[i] - 1.0; // shift to origin
        r2 += z[i] * z[i];
        sum_z += z[i];
    }
    f = pow(fabs(pow(r2, 2.0) - pow(sum_z, 2.0)), 2 * alpha) + (0.5 * r2 + sum_z) / num + 0.5;
    delete[] z;
    return f;
}

double schaffer_F7_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Schwefel's 1.2  */
{
//    ifstream inputFile("../testdata.txt", ios::binary);
//
//    if (!inputFile.is_open()) {
//        cerr << "无法打开文件" << endl;
//        return 1;
//    }
//    double temp;
//    int size = 0;
//
//    while (inputFile >> temp) {
//        double* newVar = new double[size+1]; // 创建一个新的更大的数组
//        for (int i = 0; i < size; ++i) {
//            newVar[i] = var[i]; // 将旧数据复制到新数组
//        }
//        newVar[size] = temp; // 添加新数据到新数组
//        delete[] var; // 释放旧数组内存
//        var = newVar; // 更新指针
//        size++;
//    }
//
//    inputFile.close();
//
//    cout << "从文件中读取的内容为：" << endl;
//    for (int i = 0; i < size; ++i) {
//        cout << var[i] << endl;
//    }


    int i;
    double tmp;
    double f = 0.;
    double* z = new double[num];
    double* y = new double[num];
    sr_func(var, y, num, OShift, M, 1.0, s_flag, r_flag); /* shift and rotate */  // z
    for (i = 0; i < num - 1; i++)
    {
        z[i] = pow(y[i] * y[i] + y[i + 1] * y[i + 1], 0.5);
        tmp = sin(50.0 * pow(z[i], 0.2));
        f += pow(z[i], 0.5) + pow(z[i], 0.5) * tmp * tmp;
    }
    f = f * f / (num - 1) / (num - 1);
    delete[] z;
    delete[] y;
    //cout << f;
    return f;
}

double step_rastrigin_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Noncontinuous Rastrigin's  */
{
//    ifstream inputFile("../testdata.txt", ios::binary);
//
//    if (!inputFile.is_open()) {
//        cerr << "无法打开文件" << endl;
//        return 1;
//    }
//    double temp;
//    int size = 0;
//
//    while (inputFile >> temp) {
//        double* newVar = new double[size+1]; // 创建一个新的更大的数组
//        for (int i = 0; i < size; ++i) {
//            newVar[i] = var[i]; // 将旧数据复制到新数组
//        }
//        newVar[size] = temp; // 添加新数据到新数组
//        delete[] var; // 释放旧数组内存
//        var = newVar; // 更新指针
//        size++;
//    }
//
//    inputFile.close();
//
//    cout << "从文件中读取的内容为：" << endl;
//    for (int i = 0; i < size; ++i) {
//        cout << var[i] << endl;
//    }

    int i;
    double f = 0.;
    double* y = new double[num];
    for (i = 0; i < num; i++)
    {
        if (fabs(y[i] - OShift[i]) > 0.5)
            y[i] = OShift[i] + floor(2 * (y[i] - OShift[i]) + 0.5) / 2;
    }
    double* z = new double[num];
    sr_func(var, z, num, OShift, M, 5.12 / 100.0, s_flag, r_flag); /* shift and rotate */

    for (i = 0; i < num; i++)
    {
        f += (z[i] * z[i] - 10.0 * cos(2.0 * PI * z[i]) + 10.0);
    }
    delete[] z;
    delete[] y;
    //cout<<f;
    return f + 800.0;
}

double levy_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Levy */
{
//    ifstream inputFile("../testdata.txt", ios::binary);
//
//    if (!inputFile.is_open()) {
//        cerr << "无法打开文件" << endl;
//        return 1;
//    }
//    double temp;
//    int size = 0;
//
//    while (inputFile >> temp) {
//        double* newVar = new double[size+1]; // 创建一个新的更大的数组
//        for (int i = 0; i < size; ++i) {
//            newVar[i] = var[i]; // 将旧数据复制到新数组
//        }
//        newVar[size] = temp; // 添加新数据到新数组
//        delete[] var; // 释放旧数组内存
//        var = newVar; // 更新指针
//        size++;
//    }
//
//    inputFile.close();
//
//    cout << "从文件中读取的内容为：" << endl;
//    for (int i = 0; i < size; ++i) {
//        cout << var[i] << endl;
//    }

    int i;
    double f = 0.;
    double* z = new double[num];
    sr_func(var, z, num, OShift, M, 1.0, s_flag, r_flag); /* shift and rotate */

    double* w = new double[num];

    double sum1 = 0.0;
    for (i = 0; i < num; i++)
    {
        w[i] = 1.0 + (z[i] - 0.0) / 4.0;
    }

    double term1 = pow((sin(PI * w[0])), 2);
    double term3 = pow((w[num - 1] - 1), 2) * (1 + pow((sin(2 * PI * w[num - 1])), 2));

    double sum = 0.0;

    for (i = 0; i < num - 1; i++)
    {
        double wi = w[i];
        double newv = pow((wi - 1), 2) * (1 + 10 * pow((sin(PI * wi + 1)), 2));
        sum = sum + newv;
    }

    f = term1 + sum + term3; // - 1.442600987052770; // - 1.442600987052770
    delete[] z;
    free(w);                     // ADD THIS LINE to free memory! Thanks for Dr. Janez
    return f + 900.0;
}

double zakharov_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* zakharov */
{
//
//    ifstream inputFile("../testdata.txt", ios::binary);
//
//    if (!inputFile.is_open()) {
//        cerr << "无法打开文件" << endl;
//        return 1;
//    }
//    double temp;
//    int size = 0;
//
//    while (inputFile >> temp) {
//        double* newVar = new double[size+1]; // 创建一个新的更大的数组
//        for (int i = 0; i < size; ++i) {
//            newVar[i] = var[i]; // 将旧数据复制到新数组
//        }
//        newVar[size] = temp; // 添加新数据到新数组
//        delete[] var; // 释放旧数组内存
//        var = newVar; // 更新指针
//        size++;
//    }
//
//    inputFile.close();
//
//    cout << "从文件中读取的内容为：" << endl;
//    for (int i = 0; i < size; ++i) {
//        cout << var[i] << endl;
//    }


    int i;
    double* z = new double[num];
    sr_func(var, z, num, OShift, M, 1.0, s_flag, r_flag); // shift and rotate
    double f = 0.;
    double sum1 = 0.0;
    double sum2 = 0.0;
    for (i = 0; i < num; i++)
    {
        double xi = z[i];
        sum1 = sum1 + pow(xi, 2);
        sum2 = sum2 + 0.5 * (i + 1) * xi;
    }

    f = sum1 + pow(sum2, 2) + pow(sum2, 4);
    delete[] z;
//    cout << "一轮输出：";
//    cout << f+300.0 << endl;
    return f + 300.0;
}

double katsuura_func(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Katsuura  */
{
    int i, j;
    double temp, tmp1, tmp2, tmp3;
    double f = 1.0;
   
    tmp3 = pow(1.0 * num, 1.2);
    double* z = new double[num];
    sr_func(var, z, num, OShift, M, 5.0 / 100.0, s_flag, r_flag); /* shift and rotate */

    for (i = 0; i < num; i++)
    {
        temp = 0.0;
        for (j = 1; j <= 32; j++)
        {
            tmp1 = pow(2.0, j);
            tmp2 = tmp1 * z[i];
            temp += fabs(tmp2 - floor(tmp2 + 0.5)) / tmp1;
        }
        f *= pow(1.0 + (i + 1) * temp, 10.0 / tmp3);
    }
    tmp1 = 10.0 / num / num;
    f = f * tmp1 - tmp1;
    delete[] z;
    return f;
}

double hf02(double* var, int num, double* OShift, double* M, int* S, int s_flag, int r_flag) /* Hybrid Function 2 */
{
//    ifstream inputFile("../testdata.txt", ios::binary);
//
//    if (!inputFile.is_open()) {
//        cerr << "无法打开文件" << endl;
//        return 1;
//    }
//    double temp;
//    int size = 0;
//
//    while (inputFile >> temp) {
//        double* newVar = new double[size+1]; // 创建一个新的更大的数组
//        for (int i = 0; i < size; ++i) {
//            newVar[i] = var[i]; // 将旧数据复制到新数组
//        }
//        newVar[size] = temp; // 添加新数据到新数组
//        delete[] var; // 释放旧数组内存
//        var = newVar; // 更新指针
//        size++;
//    }
//
//    inputFile.close();
//
//    cout << "从文件中读取的内容为：" << endl;
//    for (int i = 0; i < size; ++i) {
//        cout << var[i] << endl;
//    }

    int i, tmp, cf_num = 3;
    double fit[3];
    int G[3], G_num[3];
    double Gp[3] = { 0.4, 0.4, 0.2 };
    double* y = new double[num];

    tmp = 0;
    for (i = 0; i < cf_num - 1; i++)
    {
        G_num[i] = ceil(Gp[i] * num);
        tmp += G_num[i];
    }
    G_num[cf_num - 1] = num - tmp;

    G[0] = 0;
    for (i = 1; i < cf_num; i++)
    {
        G[i] = G[i - 1] + G_num[i - 1];
    }
    double* z = new double[num];
    sr_func(var, z, num, OShift, M, 1.0, s_flag, r_flag); /* shift and rotate */

    for (i = 0; i < num; i++)
    {
        y[i] = z[S[i] - 1];
    }
    i = 0;
    fit[i] = bent_cigar_func(&y[G[i]], G_num[i], OShift, M,S, 0, 0);
    i = 1;
    fit[i] = hgbat_func(&y[G[i]], G_num[i], OShift, M,S, 0, 0);
    i = 2;
    fit[i] = rastrigin_func(&y[G[i]], G_num[i], OShift, M,S, 0, 0);
    
    double f = 0.;
    for (i = 0; i < cf_num; i++)
    {
        f += fit[i];
    }
    delete[] z;
    delete[] y;
    return f + 1800.0;
}

double hf10(double* var, int num, double* OShift, double* M, int* S, int s_flag, int r_flag) /* Hybrid Function 6 */
{

//    ifstream inputFile("../testdata.txt", ios::binary);
//
//    if (!inputFile.is_open()) {
//        cerr << "无法打开文件" << endl;
//        return 1;
//    }
//    double temp;
//    int size = 0;
//
//    while (inputFile >> temp) {
//        double* newVar = new double[size+1]; // 创建一个新的更大的数组
//        for (int i = 0; i < size; ++i) {
//            newVar[i] = var[i]; // 将旧数据复制到新数组
//        }
//        newVar[size] = temp; // 添加新数据到新数组
//        delete[] var; // 释放旧数组内存
//        var = newVar; // 更新指针
//        size++;
//    }
//
//    inputFile.close();
//
//    cout << "从文件中读取的内容为：" << endl;
//    for (int i = 0; i < size; ++i) {
//        cout << var[i] << endl;
//    }

    int i, tmp, cf_num = 6;
    double fit[6];
    int G[6], G_num[6];
    double Gp[6] = { 0.1, 0.2, 0.2, 0.2, 0.1, 0.2 };
    double* y = new double[num];
    tmp = 0;
    for (i = 0; i < cf_num - 1; i++)
    {
        G_num[i] = ceil(Gp[i] * num);
        tmp += G_num[i];
    }
    G_num[cf_num - 1] = num - tmp;

    G[0] = 0;
    for (i = 1; i < cf_num; i++)
    {
        G[i] = G[i - 1] + G_num[i - 1];
    }
    double* z = new double[num];
    sr_func(var, z, num, OShift, M, 1.0, s_flag, r_flag); /* shift and rotate */

    for (i = 0; i < num; i++)
    {
        y[i] = z[S[i] - 1];
    }

    i = 0;
    fit[i] = hgbat_func(&y[G[i]], G_num[i], OShift, M,S, 0, 0);
    i = 1;
    fit[i] = katsuura_func(&y[G[i]], G_num[i], OShift, M,S, 0, 0);
    i = 2;
    fit[i] = ackley_func(&y[G[i]], G_num[i], OShift, M,S, 0, 0);
    i = 3;
    fit[i] = rastrigin_func(&y[G[i]], G_num[i], OShift, M,S, 0, 0);
    i = 4;
    fit[i] = schwefel_func(&y[G[i]], G_num[i], OShift, M,S, 0, 0);
    i = 5;
    fit[i] = schaffer_F7_func(&y[G[i]], G_num[i], OShift, M,S, 0, 0);

    double f = 0.;
    for (i = 0; i < cf_num; i++)
    {
        f += fit[i];
    }
    delete[] z;
    delete[] y;
    return f + 2000.0;
}

double hf06(double* var, int num, double* OShift, double* M, int* S, int s_flag, int r_flag) /* Hybrid Function 6 */
{
//    ifstream inputFile("../testdata.txt", ios::binary);
//
//    if (!inputFile.is_open()) {
//        cerr << "无法打开文件" << endl;
//        return 1;
//    }
//    double temp;
//    int size = 0;
//
//    while (inputFile >> temp) {
//        double* newVar = new double[size+1]; // 创建一个新的更大的数组
//        for (int i = 0; i < size; ++i) {
//            newVar[i] = var[i]; // 将旧数据复制到新数组
//        }
//        newVar[size] = temp; // 添加新数据到新数组
//        delete[] var; // 释放旧数组内存
//        var = newVar; // 更新指针
//        size++;
//    }
//
//    inputFile.close();
//
//    cout << "从文件中读取的内容为：" << endl;
//    for (int i = 0; i < size; ++i) {
//        cout << var[i] << endl;
//    }

    int i, tmp, cf_num = 5;
    double fit[5];
    int G[5], G_num[5];
    double Gp[5] = { 0.3, 0.2, 0.2, 0.1, 0.2 };
    double* y = new double[num];
    tmp = 0;
    for (i = 0; i < cf_num - 1; i++)
    {
        G_num[i] = ceil(Gp[i] * num);
        tmp += G_num[i];
    }
    G_num[cf_num - 1] = num - tmp;

    G[0] = 0;
    for (i = 1; i < cf_num; i++)
    {
        G[i] = G[i - 1] + G_num[i - 1];
    }
    double* z = new double[num];
    sr_func(var, z, num, OShift, M, 1.0, s_flag, r_flag); /* shift and rotate */

    for (i = 0; i < num; i++)
    {
        y[i] = z[S[i] - 1];
    }
    i = 0;
    fit[i] = katsuura_func(&y[G[i]], G_num[i], OShift, M,S, 0, 0);
    i = 1;
    fit[i] = happycat_func(&y[G[i]], G_num[i], OShift, M,S, 0, 0);
    i = 2;
    fit[i] = grie_rosen_func(&y[G[i]], G_num[i], OShift, M,S, 0, 0);
    i = 3;
    fit[i] = schwefel_func(&y[G[i]], G_num[i], OShift, M,S, 0, 0);
    i = 4;
    fit[i] = ackley_func(&y[G[i]], G_num[i], OShift, M,S, 0, 0);
    double f = 0.;
    for (i = 0; i < cf_num; i++)
    {
        f += fit[i];
    }
    delete[] z;
    delete[] y;
    return f + 2200.0;
}

double cf01(double* var, int num, double* OShift, double* M, int* S, int s_flag, int r_flag) /* Composition Function 1 */
{
//    ifstream inputFile("../testdata.txt", ios::binary);
//
//    if (!inputFile.is_open()) {
//        cerr << "无法打开文件" << endl;
//        return 1;
//    }
//    double temp;
//    int size = 0;
//
//    while (inputFile >> temp) {
//        double* newVar = new double[size+1]; // 创建一个新的更大的数组
//        for (int i = 0; i < size; ++i) {
//            newVar[i] = var[i]; // 将旧数据复制到新数组
//        }
//        newVar[size] = temp; // 添加新数据到新数组
//        delete[] var; // 释放旧数组内存
//        var = newVar; // 更新指针
//        size++;
//    }
//
//    inputFile.close();
//
//    cout << "从文件中读取的内容为：" << endl;
//    for (int i = 0; i < size; ++i) {
//        cout << var[i] << endl;
//    }

    int i, cf_num = 5;
    double fit[5];
    double delta[5] = { 10, 20, 30, 40, 50 };
    double bias[5] = { 0, 200, 300, 100, 400 };
    i = 0;
    fit[i] = rosenbrock_func(var, num, &OShift[i * num], &M[i * num * num], S,1, r_flag);
    fit[i] = 10000 * fit[i] / 1e+4;
    i = 1;
    fit[i] = ellips_func(var, num, &OShift[i * num], &M[i * num * num], S,1, r_flag);
    fit[i] = 10000 * fit[i] / 1e+10;
    i = 2;
    fit[i] = bent_cigar_func(var, num, &OShift[i * num], &M[i * num * num], S,1, r_flag);
    fit[i] = 10000 * fit[i] / 1e+30;
    i = 3;
    fit[i] = discus_func(var, num, &OShift[i * num], &M[i * num * num], S,1, r_flag);
    fit[i] = 10000 * fit[i] / 1e+10;
    i = 4;
    fit[i] = ellips_func(var,  num, &OShift[i * num], &M[i * num * num], S,1, 0);
    fit[i] = 10000 * fit[i] / 1e+10;
    double f = cf_cal(var, num, OShift, delta, bias, fit, cf_num);
    return f + 2300.0;
}

double cf02(double* var, int num, double* OShift, double* M, int* S, int s_flag, int r_flag) /* Composition Function 2 */
{
//    ifstream inputFile("../testdata.txt", ios::binary);
//
//    if (!inputFile.is_open()) {
//        cerr << "无法打开文件" << endl;
//        return 1;
//    }
//    double temp;
//    int size = 0;
//
//    while (inputFile >> temp) {
//        double* newVar = new double[size+1]; // 创建一个新的更大的数组
//        for (int i = 0; i < size; ++i) {
//            newVar[i] = var[i]; // 将旧数据复制到新数组
//        }
//        newVar[size] = temp; // 添加新数据到新数组
//        delete[] var; // 释放旧数组内存
//        var = newVar; // 更新指针
//        size++;
//    }
//
//    inputFile.close();
//
//    cout << "从文件中读取的内容为：" << endl;
//    for (int i = 0; i < size; ++i) {
//        cout << var[i] << endl;
//    }
    int i, cf_num = 3;
    double fit[3];
    double delta[3] = { 20, 10, 10 };
    double bias[3] = { 0, 200, 100 };

    i = 0;
    fit[i] = schwefel_func(var, num, &OShift[i * num], &M[i * num * num], S, 1, 0);
    i = 1;
    fit[i] = rastrigin_func(var, num, &OShift[i * num], &M[i * num * num], S,1, r_flag);
    i = 2;
    fit[i] = hgbat_func(var, num, &OShift[i * num], &M[i * num * num], S, 1, r_flag);
    double f = cf_cal(var, num, OShift, delta, bias, fit, cf_num);
    return f + 2400.0;
}

double cf06(double* var, int num, double* OShift, double* M, int* S, int s_flag, int r_flag) /* Composition Function 4 */
{
//    ifstream inputFile("../testdata.txt", ios::binary);
//
//    if (!inputFile.is_open()) {
//        cerr << "无法打开文件" << endl;
//        return 1;
//    }
//    double temp;
//    int size = 0;
//
//    while (inputFile >> temp) {
//        double* newVar = new double[size+1]; // 创建一个新的更大的数组
//        for (int i = 0; i < size; ++i) {
//            newVar[i] = var[i]; // 将旧数据复制到新数组
//        }
//        newVar[size] = temp; // 添加新数据到新数组
//        delete[] var; // 释放旧数组内存
//        var = newVar; // 更新指针
//        size++;
//    }
//
//    inputFile.close();
//
//    cout << "从文件中读取的内容为：" << endl;
//    for (int i = 0; i < size; ++i) {
//        cout << var[i] << endl;
//    }
    int i, cf_num = 5;
    double fit[5];
    double delta[5] = { 20, 20, 30, 30, 20 };
    double bias[5] = { 0, 200, 300, 400, 200 };
    i = 0;
    fit[i] = escaffer6_func(var, num, &OShift[i * num], &M[i * num * num], S,1, r_flag);
    fit[i] = 10000 * fit[i] / 2e+7;
    i = 1;
    fit[i] = schwefel_func(var, num, &OShift[i * num], &M[i * num * num], S,1, r_flag);
    i = 2;
    fit[i] = griewank_func(var, num, &OShift[i * num], &M[i * num * num], S,1, r_flag);
    fit[i] = 1000 * fit[i] / 100;
    i = 3;
    fit[i] = rosenbrock_func(var, num, &OShift[i * num], &M[i * num * num], S,1, r_flag);
    i = 4;
    fit[i] = rastrigin_func(var, num, &OShift[i * num], &M[i * num * num], S,1, r_flag);
    fit[i] = 10000 * fit[i] / 1e+3;
    double f = cf_cal(var, num, OShift, delta, bias, fit, cf_num);
    return f + 2600.0;
}

double cf07(double * var, int num, double * OShift,double * M, int * S,int s_flag, int r_flag) /* Composition Function 4 */
{
//    ifstream inputFile("../testdata.txt", ios::binary);
//
//    if (!inputFile.is_open()) {
//        cerr << "无法打开文件" << endl;
//        return 1;
//    }
//    double temp;
//    int size = 0;
//
//    while (inputFile >> temp) {
//        double* newVar = new double[size+1]; // 创建一个新的更大的数组
//        for (int i = 0; i < size; ++i) {
//            newVar[i] = var[i]; // 将旧数据复制到新数组
//        }
//        newVar[size] = temp; // 添加新数据到新数组
//        delete[] var; // 释放旧数组内存
//        var = newVar; // 更新指针
//        size++;
//    }
//
//    inputFile.close();
//
//    cout << "从文件中读取的内容为：" << endl;
//    for (int i = 0; i < size; ++i) {
//        cout << var[i] << endl;
//    }
    int i, cf_num = 6;
    double fit[6];
    double delta[6] = { 10, 20, 30, 40, 50, 60 };
    double bias[6] = { 0, 300, 500, 100, 400, 200 };
    i = 0;
    fit[i] = hgbat_func(var, num, &OShift[i * num], &M[i * num * num], S,1, r_flag);
    fit[i] = 10000 * fit[i] / 1000;
    i = 1;
    fit[i] = rastrigin_func(var, num, &OShift[i * num], &M[i * num * num], S,1, r_flag);
    fit[i] = 10000 * fit[i] / 1e+3;
    i = 2;
    fit[i] = schwefel_func(var, num, &OShift[i * num], &M[i * num * num], S,1, r_flag);
    fit[i] = 10000 * fit[i] / 4e+3;
    i = 3;
    fit[i] = bent_cigar_func(var, num, &OShift[i * num], &M[i * num * num], S,1, r_flag);
    fit[i] = 10000 * fit[i] / 1e+30;
    i = 4;
    fit[i] = ellips_func(var, num, &OShift[i * num], &M[i * num * num], S,1, r_flag);
    fit[i] = 10000 * fit[i] / 1e+10;
    i = 5;
    fit[i] = escaffer6_func(var, num, &OShift[i * num], &M[i * num * num], S,1, r_flag);
    fit[i] = 10000 * fit[i] / 2e+7;
    double f = cf_cal(var, num, OShift, delta, bias, fit, cf_num);
    return f + 2700.0;
}

void shiftfunc(double* var, double* varshift, int num, double* OShift)
{
    int i;
    for (i = 0; i < num; i++)
    {
        varshift[i] = var[i] - OShift[i];
    }
}

void rotatefunc(double* var, double* varrot, int num, double* M)
{
    int i, j;
    for (i = 0; i < num; i++)
    {
        varrot[i] = 0;
        for (j = 0; j < num; j++)
        {
            varrot[i] = varrot[i] + var[j] * M[i * num + j];
        }
    }
}

void sr_func(double* var, double* sr_var, int num, double* OShift, double* M, double sh_rate, int s_flag, int r_flag)
{
    int i;
    double* y = new double[num];
    if (s_flag == 1)
    {
        if (r_flag == 1)
        {
            shiftfunc(var, y, num, OShift);
            for (i = 0; i < num; i++) // shrink to the original search range
            {
                y[i] = y[i] * sh_rate;
            }
            rotatefunc(y, sr_var, num, M);
        }
        else
        {
            shiftfunc(var, sr_var, num, OShift);
            for (i = 0; i < num; i++) // shrink to the original search range
            {
                sr_var[i] = sr_var[i] * sh_rate;
            }
        }
    }
    else
    {
        if (r_flag == 1)
        {
            for (i = 0; i < num; i++) // shrink to the original search range
            {
                y[i] = var[i] * sh_rate;
            }
            rotatefunc(y, sr_var, num, M);
        }
        else
            for (i = 0; i < num; i++) // shrink to the original search range
            {
                sr_var[i] = var[i] * sh_rate;
            }
    }
    delete[] y;
}

void asyfunc(double* var, double* varasy, int num, double beta)
{
    int i;
    for (i = 0; i < num; i++)
    {
        if (var[i] > 0)
            varasy[i] = pow(var[i], 1.0 + beta * i / (num - 1) * pow(var[i], 0.5));
    }
}

void oszfunc(double* var, double* varosz, int num)
{
    int i, svar;
    double c1, c2, varvar;
    for (i = 0; i < num; i++)
    {
        if (i == 0 || i == num - 1)
        {
            if (var[i] != 0)
                varvar = log(fabs(var[i]));
            if (var[i] > 0)
            {
                c1 = 10;
                c2 = 7.9;
            }
            else
            {
                c1 = 5.5;
                c2 = 3.1;
            }
            if (var[i] > 0)
                svar = 1;
            else if (var[i] == 0)
                svar = 0;
            else
                svar = -1;
            varosz[i] = svar * exp(varvar + 0.049 * (sin(c1 * varvar) + sin(c2 * varvar)));
        }
        else
            varosz[i] = var[i];
    }
}

double cf_cal(double* var, int num, double* OShift, double* delta, double* bias, double* fit, int cf_num)
{
    int i, j;
    double w_mavar = 0, w_sum = 0;
    double* w = new double[num];

    for (i = 0; i < cf_num; i++)
    {
        fit[i] += bias[i];
        w[i] = 0;
        for (j = 0; j < num; j++)
        {
            w[i] += pow(var[j] - OShift[i * num + j], 2.0);
        }
        if (w[i] != 0)
            w[i] = pow(1.0 / w[i], 0.5) * exp(-w[i] / 2.0 / num / pow(delta[i], 2.0));
        else
            w[i] = INF;
        if (w[i] > w_mavar)
            w_mavar = w[i];
    }

    for (i = 0; i < cf_num; i++)
    {
        w_sum = w_sum + w[i];
    }
    if (w_mavar == 0)
    {
        for (i = 0; i < cf_num; i++)
            w[i] = 1;
        w_sum = cf_num;
    }
    double f = 0.;
    for (i = 0; i < cf_num; i++)
    {
        f = f + w[i] / w_sum * fit[i];
    }
    delete[] w;
    return f;
}
