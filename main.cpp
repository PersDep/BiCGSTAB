#include <iomanip>
#include <iostream>
#include <cmath>
#include <limits>
#include <vector>

#include <omp.h>

using namespace std;

double scalar_time, sum_time, multiplication_time, sparsed_scalar_time;
int scalars_amount, sums_amount, multiplications_amount;

int GetRand(int min, int max)
{
    return rand() % (max - min + 1) + min;
}

//double GetRand(int min = 1, int max = 10000000)
//{
//    return double(rand() % (max - min + 1) + min) / (double(max) / 10);
//}

struct SparsedMatrix
{
    int rows;
    int cols;
    vector<double> matrix;
    vector<int> rowsOrigins, colsNumbers;

    void MatrixGeneration(int a, int b)
    {
        rows = a, cols = b;
        for (int i = 0; i < rows; i++) {
            rowsOrigins.emplace_back(matrix.size());
            double diagonal_element = 0;
            int diagonal_element_index = 0;
            for (int j = 0; j < cols; j++) {
                if (GetRand(0, 3) == 0 || i == j) {
                    diagonal_element_index += (i == j) * matrix.size();
                    colsNumbers.emplace_back(j);
                    matrix.emplace_back(sin(i + j + 1));
                    diagonal_element += (i != j) * matrix.back();
                }
            }
            matrix[diagonal_element_index] = 1.2 * diagonal_element;
        }
    }

    explicit SparsedMatrix(int Nx = 10, int Ny = 10, int Nz = 10)
    {
        rows = Nx * Ny * Nz;
        cols = rows;
        for (int k = 0; k < Nz; k++)
            for (int j = 0; j < Ny; j++)
                for (int i = 0; i < Nx; i++) {
                    double diagonal_element = 0;
                    int index = k * (Nx * Ny) + j * Nx + i;
                    rowsOrigins.emplace_back(matrix.size());
                    if (k > 0) {
                        colsNumbers.emplace_back(index - Nx * Ny);
                        matrix.emplace_back(sin(index + colsNumbers.back() + 1));
                        diagonal_element += abs(matrix.back());
                    }
                    if (j > 0) {
                        colsNumbers.emplace_back(index - Nx);
                        matrix.emplace_back(sin(index + colsNumbers.back() + 1));
                        diagonal_element += abs(matrix.back());
                    }
                    if (i > 0) {
                        colsNumbers.emplace_back(index - 1);
                        matrix.emplace_back(sin(index + colsNumbers.back() + 1));
                        diagonal_element += abs(matrix.back());
                    }
                    colsNumbers.emplace_back(index);
                    matrix.emplace_back(0);
                    int counter = 1;
                    if (i < Nx - 1) {
                        colsNumbers.emplace_back(index + 1);
                        matrix.emplace_back(sin(index + colsNumbers.back() + 1));
                        diagonal_element += abs(matrix.back());
                        counter++;
                    }
                    if (j < Ny - 1) {
                        colsNumbers.emplace_back(index + Nx);
                        matrix.emplace_back(sin(index + colsNumbers.back() + 1));
                        diagonal_element += abs(matrix.back());
                        counter++;
                    }
                    if (k < Nz - 1) {
                        colsNumbers.emplace_back(index + Nx * Ny);
                        matrix.emplace_back(sin(index + colsNumbers.back() + 1));
                        diagonal_element += abs(matrix.back());
                        counter++;
                    }
                    matrix[matrix.size() - counter] = 1.1 * diagonal_element;
                }
        rowsOrigins.emplace_back(matrix.size());
        colsNumbers.emplace_back(-1);
    }

    void Revert()
    {
        for (int i = 0; i < rows; i++)
            for (int j = rowsOrigins[i]; j < rowsOrigins[i + 1]; j++)
                if (colsNumbers[j] != i)
                    matrix[j] = 0;
                else
                    matrix[j] = 1 / matrix[j];
    }

    void Print() const
    {
        int counter = 0;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++)
                if (colsNumbers[counter] == j)
                    cout << setprecision(4) << matrix[counter++] << '\t';
                else
                    cout << setprecision(4) << 0.0 << '\t';
            cout << endl;
        }
    }


    vector<double> GetRow(int row) const
    {
        vector<double> res;
        int counter = 0;
        for (int i = 0; i < cols; i++)
            if (colsNumbers[rowsOrigins[row] + counter] == i)
                res.emplace_back(matrix[rowsOrigins[row] + counter++]);
            else
                res.emplace_back(0);
        return res;
    }

    vector<double> GetSparsedRow(int row) const
    {
        return vector<double>(matrix.begin() + rowsOrigins[row], matrix.begin() + rowsOrigins[row + 1]);
    }
};

vector<double> multiplication(const vector<double> &vec, double a)
{
    vector<double> res(vec.size());

    #pragma omp parallel for
    for (size_t i = 0; i < vec.size(); i++)
        res[i] = a * vec[i];

    return res;
}


double multiplication(const vector<double> &vec1, const vector<double> &vec2)
{
    auto start = omp_get_wtime();
    double res = 0;

//    #pragma omp parallel for
//    for (size_t i = 0; i < vec1.size(); i++)
//        res += vec1[i] * vec2[i];
    scalars_amount++;

    scalar_time += omp_get_wtime() - start;
    return res;
}

double multiplication(const SparsedMatrix &matrix, int row, const vector<double> &vec2)
{
    auto start = omp_get_wtime();
    double res = 0;
    vector<double> vec1 = matrix.GetSparsedRow(row);

    #pragma omp parallel for
    for (size_t i = 0; i < vec1.size(); i++) {
        res += vec1[i] * vec2[matrix.colsNumbers[matrix.rowsOrigins[row] + i]];
    }

    sparsed_scalar_time += omp_get_wtime() - start;
    return res;
}

vector<double> multiplication(const SparsedMatrix &matrix, const vector<double> &vec)
{
    auto start = omp_get_wtime();
    vector<double> res(vec.size(), 0);

//    #pragma omp parallel for
//    for (int i = 0; i < matrix.rows; i++) {
//        res[i] = multiplication(matrix, i, vec);
//    }
    multiplications_amount++;

    multiplication_time += omp_get_wtime() - start;
    return res;
}

vector<double> sum(const vector<double> &vec1, const vector<double> &vec2, double a, double b)
{
    auto start = omp_get_wtime();
    vector<double> res(vec1.size());
    vector<double> vec1_multiplied = multiplication(vec1, a);
    vector<double> vec2_multiplied = multiplication(vec2, b);

//    #pragma omp parallel for
//    for (size_t i = 0; i < vec1.size(); i++)
//        res[i] = vec1_multiplied[i] + vec2_multiplied[i];
    sums_amount++;

    sum_time += omp_get_wtime() - start;
    return res;
}

struct result
{
    vector<double> res;
    double tol;
    int maxit;
    int nit;

    result(double tol, int maxit, int nit = 0) : tol(tol), maxit(maxit), nit(nit) {}
};

void solver(const SparsedMatrix &matrix, const vector<double> &right_part, result &data)
{
    SparsedMatrix DD(matrix);
    DD.Revert();
    vector<double> temp_right_part1 = right_part;
    vector<double> temp_right_part2 = right_part;
    vector<double> PP(right_part), PP2(right_part.size(), 0), TT(right_part.size(), 0),
                   SS(right_part.size(), 0), SS2(right_part.size(), 0), VV(right_part.size(), 0);
    double initres = sqrt(multiplication(right_part, right_part));
    double eps = max(numeric_limits<double>::min(), data.tol * initres);
    double res = initres;
    double alphai = 1,
           alphai_1 = 1,
           betai_1  = 1,
           Rhoi_1 = 1,
           Rhoi_2 = 1,
           wi = 1,
           wi_1 = 1;
    int i = 0;
    data.res = vector<double>(right_part.size(), 0);
    for (i = 0; i < data.maxit; i++) {
//        printf("It %d: res = %e; tol = %e\n", i, res, res / initres);
        if (res < eps) break;
        if (res > initres / numeric_limits<double>::min()) {
            data.nit = -1;
            return;
        }
        Rhoi_1 = initres * initres;
        if (i) Rhoi_1 = multiplication(temp_right_part2, temp_right_part1);
        if (abs(Rhoi_1) < numeric_limits<double>::min()) {
            data.nit = -1;
            return;
        }
        if (i) {
            betai_1 = (Rhoi_1 *  alphai_1) / (Rhoi_2 * wi_1);
            PP = sum(PP, temp_right_part1, betai_1, 1.0); //p = r + betai_1 * (p - w1 * v)
            PP = sum(PP, VV, 1.0, -wi_1 * betai_1);
        }
        PP2 = multiplication(DD, PP);
        VV = multiplication(matrix, PP2);
        alphai = multiplication(temp_right_part2, VV);
        if (abs(alphai) < numeric_limits<double>::min()) {
            data.nit = -3;
            return;
        }
        alphai = Rhoi_1 / alphai;
        SS = temp_right_part1; //s = r - alphai * v
        SS = sum(SS, VV, 1.0, -alphai);
        SS2 = multiplication(DD, SS);
        TT = multiplication(matrix, SS2);
        wi = multiplication(TT, TT);
        if (abs(wi) < numeric_limits<double>::min()) {
            data.nit = -4;
            return;
        }
        wi = multiplication(TT, SS) / wi;
        if (abs(wi) < numeric_limits<double>::min()) {
            data.nit = -5;
            return;
        }
        //x = x + alphai * p2 + wi * s2
        data.res = sum(data.res, PP2, 1, alphai);
        data.res = sum(data.res, SS2, 1, wi);
        temp_right_part1 = SS; //r = s - wi * t
        temp_right_part1 = sum(temp_right_part1, TT, 1.0, -wi);
        alphai_1 = alphai;
        Rhoi_2 = Rhoi_1;
        wi_1 = wi;
        res = sqrt(multiplication(temp_right_part1, temp_right_part1));
    }
    printf("Solver_BiCGSTAB: outres: %g\n", res);

    data.nit = i;
}

void test(int a, int b, int c, int threads)
{
    omp_set_num_threads(threads);

    SparsedMatrix matrix = SparsedMatrix(a, b, c);
////    matrix.Print();
////    cout << endl << endl;
////    for (int i = 0; i < 8; i++) {
////        auto row = matrix.GetRow(i);
////        for (size_t j = 0; j < row.size(); j++)
////            cout << row[j] << '\t';
////        cout << endl;
////    }
    vector<double> right_part(a * b * c);
    for (size_t i = 0; i < right_part.size(); i++)
        right_part[i] = sin(i);
    result data(numeric_limits<double>::min(), 1000);

    cout << "Start" << endl;
    #pragma omp parallel
    {
        #pragma omp single
        cout << "Threads: " << omp_get_num_threads() << endl;
    }
    cout << a << ' ' << b << ' ' << c << ' ' << threads << endl;
    auto start = omp_get_wtime();
    solver(matrix, right_part, data);
    auto time = omp_get_wtime() - start;
//    cout << "Time: " << time << endl;
//    cout << "Sum time: " << sum_time << endl;
//    cout << "Scalar time: " << scalar_time << endl;
//    cout << "Multiplication time: " << multiplication_time << endl
//         << "Which includes sparsed scalar time: " << sparsed_scalar_time << endl;
    cout << "Iterations: " << data.nit << endl << endl;
    sum_time = scalar_time = sparsed_scalar_time = multiplication_time = 0;

    const double axpyflop = 2000 * right_part.size() * 3 * 1E-9;
    const double dotflop = 1700 * right_part.size() * 2 * 1E-9;
    const double spmvflop = 1500 * matrix.matrix.size() * 2 * 1E-9;
//    cout << "Sums amount: " << sums_amount << endl;
//    cout << "Scalars amount: " << scalars_amount << endl;
//    cout << "Multiplications amount: " << multiplications_amount << endl;
    cout << "Total flops: " << axpyflop + dotflop + spmvflop << endl;
    cout << "Sum flops: " << axpyflop << endl;
    cout << "Scalar flops: " << dotflop << endl;
    cout << "Multiplication flops: " << spmvflop << endl << endl << endl;
    sums_amount = scalars_amount = multiplications_amount = 0;

////    cout << data.nit << endl;
////    for (auto i : data.res)
////        cout << i << ' ';
////    cout << endl;
//
//    for (int i = 0; i < right_part.size(); i++) {
//        auto row = matrix.GetSparsedRow(i);
//        double left = 0;
//        for (size_t j = 0; j < row.size(); j++)
//            left += row[j] * data.res[matrix.colsNumbers[matrix.rowsOrigins[i] + j]];
//        if (right_part[i] - left > 0.0001)
//            cout << "PROBLEMS" << endl;
//    }
}

int main()
{
//    srand(unsigned(time(nullptr)));

//    omp_set_num_threads(1);
//
//    int sum = 0;
//    double start = omp_get_wtime();
//    #pragma omp parallel for reduction(+:sum)
//    for (size_t i = 0; i < 10000000000; i++)
//        sum += 59 * 46;
//    cout << omp_get_wtime() - start << endl;
//
//    omp_set_num_threads(2);
//
//    start = omp_get_wtime();
//    #pragma omp parallel for reduction(+:sum)
//    for (size_t i = 0; i < 10000000000; i++)
//        sum += 59 * 46;
//    cout << omp_get_wtime() - start << endl;

    int a = 10, b = 10, c = 10;

    for (int i = 1; i < 9; i *= 2)
        test(a, b, c, i);

    a = 100;
    for (int i = 1; i < 9; i *= 2)
        test(a, b, c, i);

    b = 100;
    for (int i = 1; i < 9; i *= 2)
        test(a, b, c, i);

    c = 100;
    for (int i = 1; i < 9; i *= 2)
        test(a, b, c, i);

    return 0;
};