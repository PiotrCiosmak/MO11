#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "calerf.h"

using namespace std;

constexpr double D = 1.0;
constexpr double LAMBDA = 1.0;// 0.4 dla kmb
constexpr double T_MAX = 2.0;
constexpr double H = 0.1;
constexpr double eps = 1e-16;
const double A = 6.0 * sqrt(D * T_MAX);
constexpr double DT = (LAMBDA * H * H) / D;
constexpr int N = static_cast<int>(T_MAX / DT);
const int M = static_cast<int>(A / H);

double **allocateMatrix()
{
    auto **matrix = new double *[N];
    for (int i = 0; i < N; i++)
        matrix[i] = new double[M];
    return matrix;
}

void deleteMatrix(double **mat)
{
    for (int i = 0; i < N; i++)
        delete[] mat[i];
    delete[] mat;
}

void saveToFile(const string &fileName, double **matrix)
{
    ofstream file(fileName);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
            file << matrix[i][j] << "\t";
        file << endl;
    }
}

double *maxErrFL(double **error)
{
    double maxVal;
    auto *outVec = new double[N];
    for (int i = 0; i < N; i++)
    {
        maxVal = fabs(error[i][0]);
        for (int j = 0; j < M; j++)
        {
            if (maxVal < fabs(error[i][j]))
                maxVal = fabs(error[i][j]);
        }
        outVec[i] = maxVal;
    }
    return outVec;
}

void set(double *dflag, double *hflag)
{
    dflag[0] = 0.0;
    for (int i = 1; i < N; i++)
        dflag[i] = i * DT;

    hflag[0] = 0.0;
    for (int i = 1; i < M; i++)
        hflag[i] = i * H;
}

void calculateError(double **err, double **analit, double **calculatedValues)
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
            err[i][j] = abs(calculatedValues[i][j] - analit[i][j]);
}

double **analiticalSolition(const double *dtVector, const double *hVector)
{
    double **flag = allocateMatrix();
    for (int i = 0; i < M; i++)
        flag[0][i] = 0.0;
    for (int i = 0; i < N; i++)
    {
        flag[i][M - 1] = 0.0;
        flag[i][0] = 1.0;
    }
    for (int i = 1; i < N; i++)
        for (int j = 1; j < M - 1; j++)
            flag[i][j] = calerf::ERFC_L(hVector[j] / (2.0 * sqrt(D * dtVector[i])));
    return flag;
}

void directMethod(double **U)
{
    for (int i = 0; i < M; i++)
        U[0][i] = 0;
    for (int i = 0; i < N; i++)
    {
        U[i][M - 1] = 0.0;
        U[i][0] = 1.0;
    }

    for (int k = 1; k < N; k++)
        for (int i = 1; i < M - 1; i++)
            U[k][i] = U[k - 1][i] + LAMBDA * (U[k - 1][i - 1] - (2 * U[k - 1][i]) + U[k - 1][i + 1]);
}

void thomasAlgorithm(double **matrix, const double *b, double *res)
{
    double *beta, *gamma;
    beta = new double[M];
    gamma = new double[M];

    beta[0] = -(matrix[0][1] / matrix[0][0]);
    gamma[0] = (b[0] / matrix[0][0]);

    for (int i = 1; i < M; i++)
    {
        if (i <= M - 2)
            beta[i] = -((matrix[i][2 + i - 1]) / ((matrix[i][i - 1] * beta[i - 1]) + matrix[i][1 + i - 1]));

        gamma[i] = (b[i] - (matrix[i][i - 1] * gamma[i - 1])) / (matrix[i][i - 1] * beta[i - 1] + matrix[i][1 + i - 1]);
    }
    beta[M - 1] = 0;

    res[M - 1] = gamma[M - 1];
    for (int i = M - 2; i >= 0; i--)
        res[i] = beta[i] * res[i + 1] + gamma[i];

    delete beta;
}

void LUdecomposition(double **matrix, const double *b, double *res)
{
    double x;
    for (int k = 0; k < M - 1; k++)
        for (int i = k + 1; i < M; i++)
        {
            x = matrix[i][k] / matrix[k][k];
            matrix[i][k] = x;
            for (int j = k + 1; j < M; j++)
                matrix[i][j] = matrix[i][j] - (x * matrix[k][j]);
        }

    double suma;
    auto *z = new double[M];

    for (int i = 0; i < M; i++)
    {
        suma = 0;
        for (int j = 0; j <= i - 1; j++)
            suma += matrix[i][j] * z[j];
        z[i] = b[i] - suma;
    }

    for (int i = M - 1; i >= 0; i--)
    {
        suma = 0;
        for (int j = i + 1; j < M; j++)
            suma += matrix[i][j] * res[j];
        res[i] = (z[i] - suma) / matrix[i][i];
    }
}

void laasonenSolution(double **laasonen, const string &function)
{
    auto *b = new double[M];
    auto *res = new double[M];

    auto **matrix = new double *[M];
    for (int i = 0; i < M; i++) matrix[i] = new double[M];

    for (int i = 0; i < M; i++)
        for (int j = 0; j < M; j++)
            matrix[i][j] = 0;

    for (int j = 0; j < M; j++)
        laasonen[0][j] = 0;

    for (int i = 0; i < N; i++)
        laasonen[i][0] = 1;
    for (int i = 0; i < N; i++)
        laasonen[i][M - 1] = 0;

    for (int k = 1; k < N; k++)
    {
        matrix[0][0] = 1;
        b[0] = 1;

        for (int i = 1; i < M - 1; i++)
        {
            matrix[i][i] = -(1 + (2 * LAMBDA));
            matrix[i][i + 1] = LAMBDA;
            matrix[i][i - 1] = LAMBDA;
            b[i] = -laasonen[k - 1][i];
        }

        b[M - 1] = 0;
        matrix[M - 1][M - 1] = 1;

        if (function == "Thomas")
            thomasAlgorithm(matrix, b, res);
        if (function == "LU")
            LUdecomposition(matrix, b, res);

        for (int i = 1; i < M - 1; i++)
            laasonen[k][i] = res[i];
    }

    for (int i = M - 1; i >= 0; i--) delete[]matrix[i];
    delete[]matrix;
    delete[] b;
}

void results(double *dtVector, double *h_l, double **analytical, const string &method)
{
    double **U = allocateMatrix();
    double **errorMatrix = allocateMatrix();

    if (method == "LU")
        laasonenSolution(U, "LU");
    else if (method == "Thomas")
        laasonenSolution(U, "Thomas");
    else if (method == "KMB")
        directMethod(U);

    calculateError(errorMatrix, analytical, U);
    saveToFile( method + "errors.txt", errorMatrix);
    saveToFile(method + "values.txt", U);

    ofstream file04(method + "_analit04.txt");
    ofstream file10(method + "_analit10.txt");
    ofstream file16(method + "_analit16.txt");

    if (method != "KMB")
    {
        for (int i = 0; i < N; i++)
        {
            if (abs(h_l[i + 1]) < eps || abs(U[40][i + 1]) < eps)
                break;
            file04 << h_l[i + 1] << "\t" << U[40][i + 1] << "\t" << analytical[40][i + 1] << endl; // d_t = 0.3
        }

        for (int i = 0; i < N; i++)
        {
            if (abs(h_l[i + 1]) < eps || abs(U[100][i + 1]) < eps)
                break;
            file10 << h_l[i + 1] << "\t" << U[100][i + 1] << "\t" << analytical[100][i + 1] << endl; // d_t = 1.0
        }

        for (int i = 0; i < N; i++)
        {
            if (abs(h_l[i + 1]) < eps || abs(U[160][i + 1]) < eps)
                break;
            file16 << h_l[i + 1] << "\t" << U[160][i + 1] << "\t" << analytical[160][i + 1] << endl; // d_t = 1.7
        }
    }
    else
    {
        for (int i = 0; i < N; i++)
        {
            if (abs(h_l[i + 1]) < eps || abs(U[100][i + 1]) < eps)
                break;
            file04 << h_l[i + 1] << "\t" << U[100][i + 1] << "\t" << analytical[100][i + 1] << endl; // d_t = 0.3
        }

        for (int i = 0; i < N; i++)
        {
            if (abs(h_l[i + 1]) < eps || abs(U[250][i + 1]) < eps)
                break;
            file10 << h_l[i + 1] << "\t" << U[250][i + 1] << "\t" << analytical[250][i + 1] << endl; // d_t = 1.0
        }

        for (int i = 0; i < N; i++)
        {
            if (abs(h_l[i + 1]) < eps || abs(U[400][i + 1]) < eps)
                break;
            file16 << h_l[i + 1] << "\t" << U[400][i + 1] << "\t" << analytical[400][i + 1] << endl; // d_t = 1.7
        }
    }

    double *maxerrors = maxErrFL(errorMatrix);
    ofstream filemaxerrors;
    string filename6 = method + "-time.txt";
    filemaxerrors.open(filename6);
    for (int i = 0; i < N - 1; i++)
    {
        if (abs(dtVector[i + 1]) < eps || (maxerrors[i + 1]) < eps)
            break;
        filemaxerrors << dtVector[i + 1] << "\t" << maxerrors[i + 1] << endl;
    }

    cout << method << " " << H << " " << maxerrors[M - 1] << endl;

    deleteMatrix(U);
    deleteMatrix(errorMatrix);
    delete[] maxerrors;
}

int main()
{
    double *dtVector{new double[N]};
    auto *hVector{new double[M]};
    double **analiticMatrix = allocateMatrix();

    set(dtVector, hVector);

    analiticMatrix = analiticalSolition(dtVector, hVector);

    saveToFile("analytical_wartosci.txt", analiticMatrix);

    //results(dt_l, hVector, analiticMatrix, "KMB"); //LAMBDE = 0.4

    results(dtVector, hVector, analiticMatrix, "Thomas"); //LAMBDA = 1.0

    results(dtVector, hVector, analiticMatrix, "LU"); //LAMBDA = 1.0
}