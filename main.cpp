#include <iostream>
#include <cmath>
#include <limits>
#include <vector>

using namespace std;

struct Matrix
{
    int rows;
    int cols;
    vector<double> matrix;

    Matrix(int rows, int cols, double *mtr = nullptr) : rows(rows), cols(cols), matrix(vector<double>(rows * cols, 0))
    {
        if (mtr)
            for (unsigned i = 0; i < matrix.size(); i++)
                matrix[i] = mtr[i];
    }

    double& operator()(int i, int j) { return matrix[i * cols + j]; }
};

double scalar_product(const vector<double> &vec1, const vector<double> &vec2)
{

}

int axpby(vector<double> &vec1, const vector<double> &vec2, double a, double b)
{
    const size_t ss = (vec1.size() / 4) * 4;
    // float time = omp_get_wtime();

//#pragma omp parallel for
    for (size_t i = 0; i < ss; i += 4)
    {
        vec1[i] = a * vec1[i] + b * vec2[i];
        vec1[i + 1] = a * vec1[i + 1] + b * vec2[i + 1];
        vec1[i + 2] = a * vec1[i + 2] + b * vec2[i + 2];
        vec1[i + 3] = a * vec1[i + 3] + b * vec2[i + 3];
    }

    // #pragma omp parallel for
    for (size_t i = ss; i < vec1.size(); i++)
        vec1[i] = a * vec1[i] + b * vec2[i];

    // #pragma omp parallel for
    // for (int i = 0; i < vec1.size; i++)
    //     vec1[i] = a * vec1[i] + b * vec2[i];

    // time = omp_get_wtime() - time;
    // globtime += time;
    // cout << "> Time of computation = " << (time) << endl;

    return 0;
}

vector<double> multiplication(const Matrix &matrix, const vector<double> &vec)
{
    c = new double[rowsA];
    memset(c, 0, rowsA * sizeof(double));

    if (colsA != size)
    {
        cout << "Incorrect sizes!" << endl;
        return 1;
    }

    starttime = MPI_Wtime();

    if (rowsA >= colsA)
    {
        int part = rowsA / proc_amount;
        int start = part * proc_rank;

        for (int i = start; i < start + part; i++)
            for (int j = 0; j < colsA; j++)
                vectorC[i] += matrixA[i][j] * vectorB[j];

        if (proc_rank == proc_amount - 1)
            for (int i = part * proc_amount; i < rowsA; i++)
                for (int j = 0; j < colsA; j++)
                    vectorC[i] += matrixA[i][j] * vectorB[j];
    }

    if (colsA > rowsA)
    {
        int part = colsA / proc_amount;
        int start = part * proc_rank;

        for (int i = 0; i < rowsA; i++)
            for (int j = start; j < start + part; j++)
                vectorC[i] += matrixA[i][j] * vectorB[j];

        if (proc_rank == proc_amount - 1)
            for (int i = 0; i < rowsA; i++)
                for (int j = part * proc_amount; j < colsA; j++)
                    vectorC[i] += matrixA[i][j] * vectorB[j];
    }

    endtime = MPI_Wtime();
    restime = endtime - starttime;
    maxtime = restime;

    MPI_Reduce(&maxtime, &timemax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(vectorC.data(), c, rowsA, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (!proc_rank)
    {
        char buf[100];
        sprintf(buf, "max process time for %d x %d", rowsA, colsA);
        FILE *stat = fopen(buf, "a+");
        fprintf(stat, "%d %lf\n", proc_amount, timemax);
        fclose(stat);

        FILE *fileC = fopen(argv[3], "w");
        fprintf(fileC, "%d\n", rowsA);
        for (int i = 0; i < rowsA; i++)
            fprintf(fileC, "%lf ", c[i]);
        fclose(fileC);
    }

}

struct result
{
    vector<double> res;
    double tol;
    int maxit;
    int nit;
};

void solver(const Matrix &matrix, const vector<double> &right_part, result &data)
{
    Matrix DD(matrix);
    vector<double> temp_right_part1 = right_part;
    vector<double> temp_right_part2 = right_part;
    vector<double> PP(right_part.size(), 0), VV(right_part.size(), 0);
    double initres = sqrt(scalar_product(right_part, right_part));
    double eps = max(numeric_limits<double>::min(), data.tol * initres);
    double res = initres;
    double alphai_1 = 1,
           betai_1  = 1,
           Rhoi_2 = 1,
           wi_1 = 1;
    int i = 0;
    data.res = vector<double>(right_part.size(), 0);
    for (i = 0; i < data.maxit; i++) {
        printf("It %d: res = %e; tol = %e\n", i, res, res / initres);
        if (res < eps) break;
        if (res > initres / numeric_limits<double>::min()) {
            data.nit = -1;
            return;
        }
        double Rhoi_1 = initres * initres;
        if (i) Rhoi_1 = scalar_product(temp_right_part2, temp_right_part1);
        if (abs(Rhoi_1) < numeric_limits<double>::min()) {
            data.nit = -1;
            return;
        }
        if (i) {
            betai_1 = (Rhoi_1 *  alphai_1) / (Rhoi_2 * wi_1);
            axpby(PP, temp_right_part1, betai_1, 1.0); //p = r + betai_1 * (p - w1 * v)
            axpby(PP, VV, 1.0, -wi_1 * betai_1);
        }
        vector<double> PP2 = multiplication(DD, PP);
        VV = multiplication(matrix, PP2);
        double alphai = scalar_product(temp_right_part2, VV);
        if (abs(alphai) < numeric_limits<double>::min()) {
            data.nit = -3;
            return;
        }
        alphai = Rhoi_1 / alphai;
        vector<double> SS = temp_right_part1; //s = r - alphai * v
        axpby(SS, VV, 1.0, -alphai);
        vector<double > SS2 = multiplication(DD, SS);
        vector<double> TT = multiplication(matrix, SS2);
        double wi = scalar_product(TT, TT);
        if (abs(wi) < numeric_limits<double>::min()) {
            data.nit = -4;
            return;
        }
        wi = scalar_product(TT, SS) / wi;
        if (abs(wi) < numeric_limits<double>::min()) {
            data.nit = -5;
            return;
        }
        //x = x + alphai * p2 + wi * s2
        axpby(data.res, PP2,1.0,alphai);
        axpby(data.res, SS2,1.0,wi);
        temp_right_part1 = SS; //r = s - wi * t
        axpby(temp_right_part1, TT, 1.0, -wi);
        alphai_1 = alphai;
        Rhoi_2 = Rhoi_1;
        wi_1 = wi;
        res = sqrt(scalar_product(temp_right_part1, temp_right_part1));
    }
    printf("Solver_BiCGSTAB: outres: %g\n", res);

    data.nit = i;
}

int main()
{


    return 0;
}