#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include <mpi.h>
#include <omp.h>

using namespace std;

double sum_time, scalar_time, multiplication_time;

class Halo;

class SparsedMatrix
{
    vector<double> data;
    vector<size_t> rowsOrigins, colsNumbers;
    vector<size_t> glob2loc, loc2glob, procMap;
    size_t globalSize, localSize, simpleHaloSubSize;

public:
    SparsedMatrix(int proc_amount, int proc_rank, size_t Nx, size_t Ny, size_t Nz, size_t xParts, size_t yParts, size_t zParts)
    {
        size_t xPart = size_t(ceil(double(Nx) / double(xParts)));
        size_t yPart = size_t(ceil(double(Ny) / double(yParts)));
        size_t zPart = size_t(ceil(double(Nz) / double(zParts)));
        size_t startK = proc_rank / (xParts * yParts) * zPart;
        size_t startJ = (proc_rank % (xParts * yParts)) / xParts * yPart;
        size_t startI = (proc_rank - startK / zPart * (xParts * yParts) - startJ / yPart * xParts) * xPart;

        for (int k = 0; k < Nz; k++)
            for (int j = 0; j < Ny; j++)
                for (int i = 0; i < Nx; i++)
                    procMap.push_back(i / xPart + (j / yPart) * xParts + (k / zPart) * xParts * yParts);

        if (startI / xPart == xParts - 1) xPart = Nx - startI;
        if (startJ / yPart == yParts - 1) yPart = Ny - startJ;
        if (startK / zPart == zParts - 1) zPart = Nz - startK;

        globalSize = Nx * Ny * Nz;
        localSize = xPart * yPart * zPart;
        simpleHaloSubSize = 2 * (xPart * yPart + xPart * zPart + yPart * zPart);

        for (size_t k = startK; k < zPart + startK; k++)
            for (size_t j = startJ; j < yPart + startJ; j++)
                for (size_t i = startI; i < xPart + startI; i++) {
                    double diagonal_element = 0;
                    size_t index = i + Nx * j + Nx * Ny * k;
                    loc2glob.push_back(index);
                    rowsOrigins.push_back(data.size());
                    if (k > 0) {
                        colsNumbers.push_back(index - Nx * Ny);
                        data.push_back(sin(index + colsNumbers.back() + 1));
                        diagonal_element += abs(data.back());
                    }
                    if (j > 0) {
                        colsNumbers.push_back(index - Nx);
                        data.push_back(sin(index + colsNumbers.back() + 1));
                        diagonal_element += abs(data.back());
                    }
                    if (i > 0) {
                        colsNumbers.push_back(index - 1);
                        data.push_back(sin(index + colsNumbers.back() + 1));
                        diagonal_element += abs(data.back());
                    }
                    size_t diagonal_index = data.size();
                    colsNumbers.push_back(index);
                    data.push_back(0);
                    if (i < Nx - 1) {
                        colsNumbers.push_back(index + 1);
                        data.push_back(sin(index + colsNumbers.back() + 1));
                        diagonal_element += abs(data.back());
                    }
                    if (j < Ny - 1) {
                        colsNumbers.push_back(index + Nx);
                        data.push_back(sin(index + colsNumbers.back() + 1));
                        diagonal_element += abs(data.back());
                    }
                    if (k < Nz - 1) {
                        colsNumbers.push_back(index + Nx * Ny);
                        data.push_back(sin(index + colsNumbers.back() + 1));
                        diagonal_element += abs(data.back());
                    }

                    data[diagonal_index] = 1.1 * diagonal_element;
                }
        rowsOrigins.push_back(data.size());
        colsNumbers.push_back(size_t(-1));

        glob2loc = vector<size_t>(globalSize, size_t(-1));
        for (size_t i = 0; i < loc2glob.size(); i++)
            glob2loc[loc2glob[i]] = i;
    }

    void UpdateColsNumbers()
    {
        for (size_t i = 0; i < data.size(); i++)
            colsNumbers[i] = glob2loc[colsNumbers[i]];
    }

    void Revert()
    {
        for (int i = 0; i < localSize; i++)
            for (size_t j = rowsOrigins[i]; j < rowsOrigins[i + 1]; j++)
                if (colsNumbers[j] != i)
                    data[j] = 0;
                else
                    data[j] = 1 / data[j];
    }

    vector<double> multiplication(const SparsedMatrix &matrix, const vector<double> &vec) const
    {
        vector<double> result(vec.size(), 0);
        double start = MPI_Wtime();

        #pragma omp parallel for
        for (size_t i = 0; i < matrix.getLocalSize(); i++)
            for (size_t j = matrix.rowsOrigins[i]; j < matrix.rowsOrigins[i + 1]; j++)
                result[i] += matrix.data[j] * vec[matrix.colsNumbers[j]];

        multiplication_time += MPI_Wtime() - start;
        return result;
    }

    size_t getLocalSize() const { return localSize; }
    size_t getRowOrigin(size_t i) const { return rowsOrigins[i]; }
    size_t getColNumber(size_t i) const { return colsNumbers[i]; }
    size_t getGlobalFromLocal(size_t i) const { return loc2glob[i]; }
    size_t getLocalFromGlobal(size_t i) const { return glob2loc[i]; }
    size_t getProcRank(size_t i) const { return procMap[i]; }
    size_t getSimpleHaloSubSize() const { return simpleHaloSubSize; }

    void setLocal(size_t i, size_t j) { glob2loc[i] = j; }
};

class Halo
{
    vector<vector<pair<size_t, size_t> > > halo;
    vector<int> haloProc;
    int globalHaloSize;
    int proc_amount, proc_rank;

public:
    Halo(int proc_amount, int proc_rank, SparsedMatrix &matrix) : globalHaloSize(0), proc_amount(proc_amount), proc_rank(proc_rank)
    {
        vector<vector<pair<size_t, size_t> > > simpleHalo(proc_amount, vector<pair<size_t , size_t> >(matrix.getSimpleHaloSubSize(), make_pair(-1, -1)));
        for (size_t i = 0; i < matrix.getLocalSize(); i++)
            for (size_t k = matrix.getRowOrigin(i); k < matrix.getRowOrigin(i + 1); k++)
                if (matrix.getLocalFromGlobal(matrix.getColNumber(k)) == -1)
                {
                    for (size_t j = 0; j < simpleHalo.size(); j++)
                        if (simpleHalo[matrix.getProcRank(matrix.getColNumber(k))][j].first == -1)
                        {
                            simpleHalo[matrix.getProcRank(matrix.getColNumber(k))][j].first = matrix.getColNumber(k);
                            break;
                        }

                    for (size_t j = 0; j < simpleHalo.size(); j++)
                        if (simpleHalo[matrix.getProcRank(matrix.getColNumber(k))][j].second == -1)
                        {
                            simpleHalo[matrix.getProcRank(matrix.getColNumber(k))][j].second = matrix.getGlobalFromLocal(
                                    i);
                            break;
                        }
                }

        for (int i = 0; i < proc_amount; i++)
        {
            if (simpleHalo[i][0].first != -1)
            {
                int counter = 0;
                for (size_t j = 0; j < simpleHalo.size(); j++)
                    if (simpleHalo[i][j].first != -1)
                    {
                        matrix.setLocal(simpleHalo[i][j].first, matrix.getLocalSize() + globalHaloSize++);
                        counter++;
                    }
                halo.push_back(vector<pair<size_t, size_t > >(size_t(counter)));
                haloProc.push_back(i);
                for (int j = 0; j < counter; j++)
                    halo.back()[j].first = matrix.getLocalFromGlobal(simpleHalo[i][j].first);
            }
            if (simpleHalo[i][0].second != -1)
            {
                int counter = 0;
                for (size_t j = 0; j < simpleHalo.size(); j++)
                    if (simpleHalo[i][j].second != -1)
                        counter++;
                for (int j = 0; j < counter; j++)
                    halo.back()[j].second = matrix.getLocalFromGlobal(simpleHalo[i][j].second);
            }
        }
    }

    void sendRecv(const SparsedMatrix &matrix, vector<double> &vec) const
    {
        vector<MPI_Request> send_request(halo.size()), recv_request(halo.size());
        vector<vector<double> > send(halo.size()), recv(halo.size());
        for (size_t i = 0; i < halo.size(); i++)
        {
            recv[i] = vector<double>(halo[i].size());
            for (size_t j = 0; j < halo[i].size(); j++)
                send[i].push_back(vec[halo[i][j].second]);
            MPI_Isend(send[i].data(), int(halo[i].size()), MPI_DOUBLE, haloProc[i], 0, MPI_COMM_WORLD, &(send_request[i]));
            MPI_Irecv(recv[i].data(), int(halo[i].size()), MPI_DOUBLE, haloProc[i], MPI_ANY_TAG, MPI_COMM_WORLD, &(recv_request[i]));
        }
        MPI_Waitall(int(halo.size()), send_request.data(), MPI_STATUS_IGNORE);
        MPI_Waitall(int(halo.size()), recv_request.data(), MPI_STATUS_IGNORE);
        for (size_t i = 0; i < halo.size(); i++)
            for (size_t j = 0; j < halo[i].size(); j++)
                vec[halo[i][j].first] = recv[i][j];
    }

    int getGlobalHaloSize() const { return globalHaloSize; }
    int getProcAmount() const { return proc_amount; }
    int getProcRank() const { return proc_rank; }
};

double multiplication(const vector<double> &vec1, const vector<double> &vec2, size_t size)
{
    double local_result = 0, result = 0;
    double start = MPI_Wtime();

    #pragma omp parallel for reduction(+:local_result)
    for (size_t  i = 0; i < size; i++)
        local_result += vec1[i] * vec2[i];

    scalar_time += MPI_Wtime() - start;
    MPI_Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return result;
}

vector<double> sum(const vector<double> &vec1, const vector<double> &vec2, double a, double b, size_t size)
{
    vector<double> result(vec1.size());
    double start = MPI_Wtime();

    #pragma omp parallel for
    for (size_t i = 0; i < size; i++)
        result[i] = vec1[i] * a + vec2[i] * b;

    sum_time += MPI_Wtime() - start;
    return result;
}

struct result
{
    vector<double> res;
    double tol;
    int maxit;
    int nit;

    result(double tol, int maxit, int nit = 0) : tol(tol), maxit(maxit), nit(nit) {}
};

void solver(const SparsedMatrix &matrix, const vector<double> &right_part, const Halo &halo, result &data)
{
    SparsedMatrix reverted_matrix(matrix);
    reverted_matrix.Revert();
    vector<double> right_part_copy = right_part;
    vector<double> PP(right_part), PP2(right_part.size(), 0), TT(right_part.size(), 0),
                   VV(right_part.size(), 0), SS(right_part.size(), 0), SS2(right_part.size(), 0);
    double mineps = 1e-15;
    double initres = sqrt(multiplication(right_part_copy, right_part_copy, matrix.getLocalSize()));
    double res = initres;
    double eps = max(mineps, data.tol * initres);
    double alphai, alphai_1 = 1, betai_1, Rhoi_1, Rhoi_2 = 1, wi, wi_1 = 1;
    int i = 0;
    data.res = vector<double>(right_part.size(), 0);
    for (i = 0; i < data.maxit; i++)
    {
        if (res < eps) break;
        if (res > initres / mineps)
        {
            data.nit = -1;
            return;
        }
        Rhoi_1 = initres * initres;
        if (i) Rhoi_1 = multiplication(right_part, right_part_copy, matrix.getLocalSize());
        if (abs(Rhoi_1) < numeric_limits<double>::min())
        {
            data.nit = -1;
            return;
        }
        if (i)
        {
            betai_1 = (Rhoi_1 * alphai_1) / (Rhoi_2  * wi_1);
            PP = sum(PP, right_part_copy, betai_1, 1, matrix.getLocalSize());
            PP = sum(PP, VV, 1, -wi_1 * betai_1, matrix.getLocalSize());
        }
        PP2 = matrix.multiplication(reverted_matrix, PP);
        halo.sendRecv(matrix, PP2);
        VV = matrix.multiplication(matrix, PP2);
        alphai = multiplication(right_part, VV, matrix.getLocalSize());
        if (fabs(alphai) < numeric_limits<double>::min())
        {
            data.nit = -3;
            return;
        }
        alphai = Rhoi_1 / alphai;
        SS = right_part_copy;
        SS = sum(SS, VV, 1, -alphai, matrix.getLocalSize());
        SS2 = matrix.multiplication(reverted_matrix, SS);
        halo.sendRecv(matrix, SS2);
        TT = matrix.multiplication(matrix, SS2);
        wi = multiplication(TT, TT, matrix.getLocalSize());
        if (abs(wi) < numeric_limits<double>::min())
        {
            data.nit = -4;
            return;
        }
        wi = multiplication(TT, SS, matrix.getLocalSize()) / wi;
        if (abs(wi) < numeric_limits<double>::min())
        {
            data.nit = -5;
            return;
        }
        data.res = sum(data.res, PP2, 1, alphai, matrix.getLocalSize());
        data.res = sum(data.res, SS2, 1, wi, matrix.getLocalSize());
        right_part_copy = sum(SS, TT, 1.0, -wi, matrix.getLocalSize());
        alphai_1 = alphai;
        Rhoi_2 = Rhoi_1;
        wi_1 = wi;
        res = sqrt(multiplication(right_part_copy, right_part_copy, matrix.getLocalSize()));
    }
    if (halo.getProcRank() == 0)
        cout << "Solver_BiCGSTAB: outres: " << res << endl;

    data.nit = i;
}

void test(int proc_amount, int proc_rank, size_t xParts, size_t yParts, size_t zParts, size_t a, size_t b, size_t c, int threads)
{
    omp_set_num_threads(threads);

    SparsedMatrix matrix(proc_amount, proc_rank, a, b, c, xParts, yParts, zParts);
    Halo halo(proc_amount, proc_rank, matrix);
    matrix.UpdateColsNumbers();
    vector<double> right_part(matrix.getLocalSize() + halo.getGlobalHaloSize());
    for (int i = 0; i < matrix.getLocalSize(); i++)
        right_part[i] = sin(matrix.getGlobalFromLocal(i));
    result data(numeric_limits<double>::min(), 1000);

    if (halo.getProcRank() == 0)
    {
        cout << "Start" << endl;
        cout << "Proc amount: " << proc_amount << endl;
        #pragma omp parallel
        {
            #pragma omp single
            cout << "Threads: " << omp_get_num_threads() << endl;
        }
        cout << a << ' ' << b << ' ' << c << ' ' << threads << endl;
    }

    double start = MPI_Wtime();
    solver(matrix, right_part, halo, data);
    MPI_Barrier(MPI_COMM_WORLD);
    double time = MPI_Wtime() - start;

    if (halo.getProcRank() == 0)
    {
        cout << "Time: " << time << endl;
        cout << "Sum time: " << sum_time << endl;
        cout << "Scalar time: " << scalar_time << endl;
        cout << "Multiplication time: " << multiplication_time << endl;
        cout << "Iterations: " << data.nit << endl;
        cout << endl;
    }
    sum_time = scalar_time = multiplication_time = 0;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int proc_amount = 0, proc_rank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &proc_amount);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

    size_t a = size_t(atoi(argv[1])), b = size_t(atoi(argv[2])), c = size_t(atoi(argv[3]));
    size_t xParts = size_t(atoi(argv[4])), yParts = size_t(atoi(argv[5])), zParts = size_t(atoi(argv[6]));
    bool openmp = false;
    if (argc == 8 && string(argv[7]) == "openmp")
        openmp = true;

    if (openmp)
        for (int i = 1; i < 5; i *= 2)
            test(proc_amount, proc_rank, xParts, yParts, zParts, a, b, c, i);
    else
        test(proc_amount, proc_rank, xParts, yParts, zParts, a, b, c, 1);

    MPI_Finalize();
    
    return 0;
}
