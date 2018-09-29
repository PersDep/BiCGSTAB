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

result solver(const Matrix &matrix, const vector<double> &right_part, result &data)
{
    vector<double> temp_right_part1 = right_part;
    vector<double> temp_right_part2 = right_part;
    double initres = sqrt(scalar_product(right_part, right_part));
    double eps = max(numeric_limits<double>::min(), data.tol * initres);
    double res = initres;
    for(int i = 0; i < data.maxit; i++) {
        printf("It %d: res = %e; tol = %e\n", i, res, res / initres);
        if(res<eps) break;
        if(res>initres/mineps) return -1;
        if(i==0) Rhoi_1 = initres*initres;
        else Rhoi_1 = dot(RR2,RR);
        if(fabs(Rhoi_1)<RhoMin) return -1
        if(i==0) PP=RR;
        else{
            betai_1=(Rhoi_1*alphai_1)/(Rhoi_2*wi_1);
            axpby(PP, RR, betai_1, 1.0); // p=r+betai_1*(p-w1*v)
            axpby(PP, VV, 1.0, -wi_1*betai_1);
        }
        SpMV(DD,PP,PP2);
        SpMV(A,PP2,VV);
        alphai = dot(RR2,VV);
        if(fabs(alphai)<RhoMin) return -3;
        alphai = Rhoi_1/alphai;
        SS=RR; // s=r-alphai*v
        axpby(SS, VV, 1.0, -alphai);
        SpMV(DD, SS, SS2);
        SpMV(A, SS2, TT);
        wi = dot(TT, TT);
        if(fabs(wi)<RhoMin) return -4;
        wi = dot(TT, SS) / wi;
        if(fabs(wi)<RhoMin) return -5;
        // x=x+alphai*p2+wi*s2
        axpby(XX, PP2,1.0,alphai);
        axpby(XX, SS2,1.0,wi);
        RR=SS; // r=s-wi*t
        9
        axpby(RR, TT, 1.0, -wi);
        alphai_1=alphai;
        Rhoi_2=Rhoi_1;
        wi_1=wi;
        res = sqrt(dot(RR,RR));
    }
    if(info) printf("Solver_BiCGSTAB: outres: %g\n",res);
    return I;
}

int main()
{


    return 0;
}