#ifndef DATAFACTORY_H
#define DATAFACTORY_H

#include <iostream>
#include <fstream>
#include "qmm.h"
#include "qmmwqi.h"

template<class V, class T>
class StateEquivalence
{
public:
    QuantumMealyMachine<V> machine;
    Matrix<V, T> rho1, rho2;
    AlmostEqual<V> equal;
    int limit;
    StateEquivalence(const QuantumMealyMachine<V> &machine, const Matrix<V, T> &rho1, const Matrix<V, T> &rho2, const AlmostEqual<V> &equal, int limit = -1) : machine(machine), rho1(rho1), rho2(rho2), equal(equal), limit(limit)
    {
        assert(rho1.row == rho1.col && rho1.row == machine.dim);
        assert(rho2.row == rho2.col && rho2.row == machine.dim);
    }
    void solve() const
    {
        if (limit == -1)
        {
            auto ret = machine.check(rho1, rho2, equal);
            if (ret.first)
                puts("rho1 and rho2 are equivalent.");
            else
            {
                puts("rho1 and rho2 are not equivalent.");
                std::cout << ret.second;
            }
        }
        else
        {
            auto ret = machine.check_measure_limit(rho1, rho2, equal, limit);
            if (ret.first)
                printf("rho1 and rho2 are equivalent within %d measurement(s).\n", limit);
            else
            {
                printf("rho1 and rho2 are not equivalent within %d measurement(s).\n", limit);
                std::cout << ret.second;
            }
        }
    }
    void solve_regular_measurement() const
    {
        if (limit == -1)
        {
            auto ret = machine.check_regular_measurement(rho1, rho2, equal);
            if (ret.first)
                puts("rho1 and rho2 are equivalent.");
            else
            {
                puts("rho1 and rho2 are not equivalent.");
                std::cout << ret.second;
            }
        }
    }
    void solve_slow() const
    {
        if (limit == -1)
        {
            auto ret = machine.check_slow(rho1, rho2, equal);
            if (ret.first)
                puts("rho1 and rho2 are equivalent.");
            else
            {
                puts("rho1 and rho2 are not equivalent.");
                std::cout << ret.second;
            }
        }
        else
        {
            auto ret = machine.check_measure_limit_slow(rho1, rho2, equal, limit);
            if (ret.first)
                printf("rho1 and rho2 are equivalent within %d measurement(s).\n", limit);
            else
            {
                printf("rho1 and rho2 are not equivalent within %d measurement(s).\n", limit);
                std::cout << ret.second;
            }
        }
    }
    const static StateEquivalence input(const char *filename)
    {
        FILE *fin = fopen(filename, "rb");

        int dim, sigma, gamma, limit;
        V eps;

        {
            int size = sizeof(int)*4 + sizeof(V);
            int off = 0;
            char buf[size];
            fread(buf, 1, size, fin);
            dim = *(int*)(buf+off); off += sizeof(int);
            sigma = *(int*)(buf+off); off += sizeof(int);
            gamma = *(int*)(buf+off); off += sizeof(int);
            limit = *(int*)(buf+off); off += sizeof(int);
            eps = *(V*)(buf+off); off += sizeof(V);
            assert(off == size);
        }

        int size = sizeof(T)*(sigma+gamma+2)*dim*dim;
        int off = 0;
        char buf[size];
        fread(buf, 1, size, fin);

        std::vector<Matrix<V, T> > unitary;
        for (int a = 0; a < sigma; ++ a)
        {
            Matrix<V, T> mat = Matrix<V, T>::zeros(dim);
            for (int i = 0; i < dim; ++ i)
                for (int j = 0; j < dim; ++ j)
                {
                    mat[i][j] = *(T*)(buf+off); off += sizeof(T);
                }
            unitary.push_back(mat);
        }

        std::vector<Matrix<V, T> > measure;

        for (int b = 0; b < gamma; ++ b)
        {
            Matrix<V, T> mat = Matrix<V, T>::zeros(dim);
            for (int i = 0; i < dim; ++ i)
                for (int j = 0; j < dim; ++ j)
                {
                    mat[i][j] = *(T*)(buf+off); off += sizeof(T);
                }
            measure.push_back(mat);
        }

        Matrix<V, T> rho1 = Matrix<V, T>::zeros(dim);
        Matrix<V, T> rho2 = Matrix<V, T>::zeros(dim);
        for (int i = 0; i < dim; ++ i)
            for (int j = 0; j < dim; ++ j)
            {
                rho1[i][j] = *(T*)(buf+off); off += sizeof(T);
            }

        for (int i = 0; i < dim; ++ i)
            for (int j = 0; j < dim; ++ j)
            {
                rho2[i][j] = *(T*)(buf+off); off += sizeof(T);
            }

        assert(off == size);

        fclose(fin);

        QuantumMealyMachine<V> machine(sigma, gamma, dim, unitary, measure);
        AlmostEqual<V> equal(eps);

        return StateEquivalence(machine, rho1, rho2, equal, limit);
    }
    const static StateEquivalence raw_input(const char *filename)
    {

        int dim, sigma, gamma, limit;
        V eps;

        std::ifstream fin(filename);

        fin >> dim >> sigma >> gamma >> limit >> eps;

        std::vector<Matrix<V, T> > unitary;
        for (int a = 0; a < sigma; ++ a)
        {
            Matrix<V, T> mat = Matrix<V, T>::zeros(dim);
            for (int i = 0; i < dim; ++ i)
                for (int j = 0; j < dim; ++ j)
                    fin >> mat[i][j];
            unitary.push_back(mat);
        }

        std::vector<Matrix<V, T> > measure;

        for (int b = 0; b < gamma; ++ b)
        {
            Matrix<V, T> mat = Matrix<V, T>::zeros(dim);
            for (int i = 0; i < dim; ++ i)
                for (int j = 0; j < dim; ++ j)
                    fin >> mat[i][j];
            measure.push_back(mat);
        }

        Matrix<V, T> rho1 = Matrix<V, T>::zeros(dim);
        Matrix<V, T> rho2 = Matrix<V, T>::zeros(dim);
        for (int i = 0; i < dim; ++ i)
            for (int j = 0; j < dim; ++ j)
                fin >> rho1[i][j];

        for (int i = 0; i < dim; ++ i)
            for (int j = 0; j < dim; ++ j)
                fin >> rho2[i][j];

        fin.close();

        QuantumMealyMachine<V> machine(sigma, gamma, dim, unitary, measure);
        AlmostEqual<V> equal(eps);

        return StateEquivalence(machine, rho1, rho2, equal, limit);
    }
    void save(const char *filename) const
    {
        int size = sizeof(int)*4 + sizeof(V) + sizeof(T)*(machine.sigma+machine.gamma+2)*machine.dim*machine.dim;
        int off = 0;
        char buf[size];

        *(int*)(buf+off) = machine.dim; off += sizeof(int);
        *(int*)(buf+off) = machine.sigma; off += sizeof(int);
        *(int*)(buf+off) = machine.gamma; off += sizeof(int);
        *(int*)(buf+off) = limit; off += sizeof(int);

        *(V*)(buf+off) = equal.eps; off += sizeof(V);

        for (int a = 0; a < machine.sigma; ++ a)
        {
            for (int i = 0; i < machine.dim; ++ i)
                for (int j = 0; j < machine.dim; ++ j)
                {
                    *(T*)(buf+off) = machine.unitary[a][i][j]; off += sizeof(T);
                }
        }

        for (int b = 0; b < machine.gamma; ++ b)
        {
            for (int i = 0; i < machine.dim; ++ i)
                for (int j = 0; j < machine.dim; ++ j)
                {
                    *(T*)(buf+off) = machine.measure[b][i][j]; off += sizeof(T);
                }
        }

        for (int i = 0; i < machine.dim; ++ i)
            for (int j = 0; j < machine.dim; ++ j)
            {
                *(T*)(buf+off) = rho1[i][j]; off += sizeof(T);
            }

        for (int i = 0; i < machine.dim; ++ i)
            for (int j = 0; j < machine.dim; ++ j)
            {
                *(T*)(buf+off) = rho2[i][j]; off += sizeof(T);
            }

        assert(off == size);

        FILE *fout = fopen(filename, "wb");
        fwrite(buf, 1, size, fout);
        fclose(fout);
    }
    void raw_save(const char *filename) const
    {
        std::ofstream fout(filename);
        fout << std::setprecision(12);

        fout << machine.dim << " " << machine.sigma << " " << machine.gamma << " " << limit << " " << equal.eps << std::endl;

        for (int a = 0; a < machine.sigma; ++ a)
        {
            for (int i = 0; i < machine.dim; ++ i)
                for (int j = 0; j < machine.dim; ++ j)
                    fout << machine.unitary[a][i][j] << " ";
            fout << std::endl;
        }

        for (int b = 0; b < machine.gamma; ++ b)
        {
            for (int i = 0; i < machine.dim; ++ i)
                for (int j = 0; j < machine.dim; ++ j)
                    fout << machine.measure[b][i][j] << " ";
            fout << std::endl;
        }

        for (int i = 0; i < machine.dim; ++ i)
            for (int j = 0; j < machine.dim; ++ j)
                fout << rho1[i][j] << " ";
        fout << std::endl;

        for (int i = 0; i < machine.dim; ++ i)
            for (int j = 0; j < machine.dim; ++ j)
                fout << rho2[i][j] << " ";
        fout << std::endl;

        fout.close();
    }
};

namespace DataFactory
{
    typedef double Real;
    typedef Complex<Real> complex;
    typedef Matrix<Real, complex> matrix;
    std::vector<StateEquivalence<double, complex > > query;

    const Real pi = 3.1415926535897932384626433;
    const complex imag(0, 1);
    const matrix I = matrix::eye(2);
    const matrix O = matrix::zeros(2);
    const matrix X = {{0, 1}, {1, 0}};
    const matrix Y = {{0, -imag}, {imag, 0}};
    const matrix Z = {{1, 0}, {0, -1}};
    const matrix H = {{sqrt(2.0)/2, sqrt(2.0)/2}, {sqrt(2.0)/2, -sqrt(2.0)/2}};
    const matrix T = {{1, 0}, {0, exp(imag*pi/4)}};
    const matrix CNOT = direct_sum(I, X);
    const matrix SWAP = {{1, 0, 0, 0},
                         {0, 0, 1, 0},
                         {0, 1, 0, 0},
                         {0, 0, 0, 1}};
    const matrix Toffoli = direct_sum(I, I, I, X);

    void init_data()
    {
        // 1
        AlmostEqual<Real> equal(1e-8);
        std::vector<matrix> unitary{I, X};
        std::vector<matrix> measure{{{sqrt(2.0)/2, 0}, {0, 0}}, {{sqrt(2.0)/2, 0}, {0, 1}}};
        matrix rho1{{0.5, 0.5}, {0.5, 0.5}};
        matrix rho2{{0.5, 0}, {0, 0.5}};
        QuantumMealyMachine<Real> machine(2, 2, 2, unitary, measure);
        query.emplace_back(machine, rho1, rho2, equal);

        // 2
        unitary = {tensor_product(H, I), CNOT};
        measure = {{{1, 0, 0, 0},
                    {0, 1, 0, 0},
                    {0, 0, 0, 0},
                    {0, 0, 0, 0}},
                   {{0, 0, 0, 0},
                    {0, 0, 0, 0},
                    {0, 0, 1, 0},
                    {0, 0, 0, 1}}};
        machine = {2, 2, 4, unitary, measure};
        rho1 = {{1, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0}};
        rho2 = {{0, 0, 0, 0},
                {0, 1, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0}};
        query.emplace_back(machine, rho1, rho2, equal);

        // 3
        rho2 = {{0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 1}};
        query.emplace_back(machine, rho1, rho2, equal);

        // 4
        rho2 = rho1;
        rho1 = matrix(4, 4, [](int, int){ return 0.25; });
        query.emplace_back(machine, rho1, rho2, equal);

        // 5
        unitary = {tensor_product(H, I), tensor_product(I, H), CNOT};
        machine = {3, 2, 4, unitary, measure};
        rho1 = {{1, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0}};
        rho2 = {{0, 0, 0, 0},
                {0, 1, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0}};
        query.emplace_back(machine, rho1, rho2, equal);

        // 6
        rho2 = {{0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 1}};
        query.emplace_back(machine, rho1, rho2, equal);

        // 7
        rho2 = rho1;
        rho1 = matrix(4, 4, [](int, int){ return 0.25; });
        query.emplace_back(machine, rho1, rho2, equal);

        // 8
        unitary = {tensor_product(H, I), SWAP};
        measure = {direct_sum(I, O), direct_sum(O, I)};
        machine = {2, 2, 4, unitary, measure};
        rho1 = {{0.5, 0, 0, 0.5},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0.5, 0, 0, 0.5}};
        rho2 = {{0.5, 0, 0, -0.5},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {-0.5, 0, 0, 0.5}};
        query.emplace_back(machine, rho1, rho2, equal);

        // 9
        query.emplace_back(machine, rho1, rho2, equal, 1);
        // 10
        query.emplace_back(machine, rho1, rho2, equal, 2);

        // 11
        auto H1 = tensor_product(H, I, I);
        auto H2 = tensor_product(I, H, I);
        auto H3 = tensor_product(I, I, H);
        unitary = {H1, H2, H3, Toffoli};
        measure = {direct_sum(I, I, O, O), direct_sum(O, O, I, I)};
        machine = {4, 2, 8, unitary, measure};
        matrix state1{{1}, {0}, {0}, {0}, {0}, {0}, {0}, {0}};
        matrix state2{{0}, {1}, {0}, {0}, {0}, {0}, {0}, {0}};
        rho1 = state1*state1.hermitian();
        rho2 = state2*state2.hermitian();
        query.emplace_back(machine, rho1, rho2, equal);

        // 12
        state1 = {{0}, {1}, {0}, {0}, {0}, {0}, {0}, {0}};
        state2 = {{0}, {0}, {1}, {0}, {0}, {0}, {0}, {0}};
        rho1 = state1*state1.hermitian();
        rho2 = state2*state2.hermitian();
        query.emplace_back(machine, rho1, rho2, equal);

        // 13
        state1 = {{1}, {0}, {0}, {0}, {0}, {0}, {0}, {0}};
        state2 = {{0}, {0}, {1}, {0}, {0}, {0}, {0}, {0}};
        rho1 = state1*state1.hermitian();
        rho2 = state2*state2.hermitian();
        query.emplace_back(machine, rho1, rho2, equal);

        // 14
        state1 = {{0}, {0}, {1}, {0}, {0}, {0}, {0}, {0}};
        state2 = {{0}, {0}, {0}, {1}, {0}, {0}, {0}, {0}};
        rho1 = state1*state1.hermitian();
        rho2 = state2*state2.hermitian();
        query.emplace_back(machine, rho1, rho2, equal);

        // 15
        auto C2H3 = tensor_product(I, direct_sum(I, H));
        unitary = {H1, H2, C2H3, Toffoli};
        machine = {4, 2, 8, unitary, measure};
        state1 = {{1}, {0}, {0}, {0}, {0}, {0}, {0}, {0}};
        state2 = {{0}, {1}, {0}, {0}, {0}, {0}, {0}, {0}};
        rho1 = state1*state1.hermitian();
        rho2 = state2*state2.hermitian();
        query.emplace_back(machine, rho1, rho2, equal);

        // 16
        H1 = tensor_product(H, I, I, I);
        H2 = tensor_product(I, H, I, I);
        H3 = tensor_product(I, I, H, I);
        auto H4 = tensor_product(I, I, I, H);
        auto C3NOT = direct_sum(matrix::eye(14), X);
        unitary = {H1, H2, H3, H4, C3NOT};
        auto M0 = direct_sum(matrix::eye(8), matrix::zeros(8));
        auto M1 = direct_sum(matrix::zeros(8), matrix::eye(8));
        measure = {M0, M1};
        machine = {5, 2, 16, unitary, measure};
        state1 = matrix::zeros(16, 1);
        state1[0][0] = 1;
        state2 = matrix::zeros(16, 1);
        state2[1][0] = 1;
        rho1 = state1*state1.hermitian();
        rho2 = state2*state2.hermitian();
        query.emplace_back(machine, rho1, rho2, equal);

        // 17
        H1 = tensor_product(H, I, I, I, I);
        H2 = tensor_product(I, H, I, I, I);
        H3 = tensor_product(I, I, H, I, I);
        H4 = tensor_product(I, I, I, H, I);
        auto H5 = tensor_product(I, I, I, I, H);
        auto C4NOT = direct_sum(matrix::eye(30), X);
        unitary = {H1, H2, H3, H4, H5, C4NOT};
        M0 = direct_sum(matrix::eye(16), matrix::zeros(16));
        M1 = direct_sum(matrix::zeros(16), matrix::eye(16));
        measure = {M0, M1};
        machine = {6, 2, 32, unitary, measure};
        state1 = matrix::zeros(32, 1);
        state1[0][0] = 1;
        state2 = matrix::zeros(32, 1);
        state2[1][0] = 1;
        rho1 = state1*state1.hermitian();
        rho2 = state2*state2.hermitian();
        query.emplace_back(machine, rho1, rho2, equal);

        // 18
        matrix FT(4, 4, [](int i, int j)
        {
            return exp(imag*pi*i*j/2)/2;
        });
        matrix BPR = CNOT*tensor_product(H, I);
        unitary = {FT, BPR};
        measure = {{{1, 0, 0, 0},
                    {0, 1, 0, 0},
                    {0, 0, 0, 0},
                    {0, 0, 0, 0}},
                   {{0, 0, 0, 0},
                    {0, 0, 0, 0},
                    {0, 0, 1, 0},
                    {0, 0, 0, 1}}};
        machine = {2, 2, 4, unitary, measure};
        rho1 = {{1, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0}};
        rho2 = {{0, 0, 0, 0},
                {0, 1, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0}};
        query.emplace_back(machine, rho1, rho2, equal);

        // 19
        measure = {{{1, 0, 0, 0},
                    {0, 0, 0, 0},
                    {0, 0, 1, 0},
                    {0, 0, 0, 0}},
                   {{0, 0, 0, 0},
                    {0, 1, 0, 0},
                    {0, 0, 0, 0},
                    {0, 0, 0, 1}}};
        machine = {2, 2, 4, unitary, measure};
        rho1 = {{1, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0}};
        rho2 = {{0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 1, 0},
                {0, 0, 0, 0}};
        query.emplace_back(machine, rho1, rho2, equal);

        // 20
        unitary = {FT};
        machine = {1, 2, 4, unitary, measure};
        query.emplace_back(machine, rho1, rho2, equal);

        // 21
        unitary = {FT, tensor_product(T, I)};
        machine = {2, 2, 4, unitary, measure};
        query.emplace_back(machine, rho1, rho2, equal);

        // 22
        unitary = {FT, tensor_product(I, T)};
        machine = {2, 2, 4, unitary, measure};
        query.emplace_back(machine, rho1, rho2, equal);

        // 23
        FT = matrix(8, 8, [](int i, int j)
        {
            return exp(imag*pi*i*j/4)/(2*sqrt(2.));
        });
        unitary = {FT, tensor_product(T, I, I)};
        measure = {direct_sum(O, I, O, I), direct_sum(I, O, I, O)};
        state1 = matrix::zeros(8, 1);
        state1[0][0] = 1;
        state2 = matrix::zeros(8, 1);
        state2[4][0] = 1;
        rho1 = state1*state1.hermitian();
        rho2 = state2*state2.hermitian();
        machine = {2, 2, 8, unitary, measure};
        query.emplace_back(machine, rho1, rho2, equal);

        // 24
        query.emplace_back(machine, rho1, rho2, equal, 1);

        // 25
        query.emplace_back(machine, rho1, rho2, equal, 2);

        // 26
        FT = matrix(16, 16, [](int i, int j)
        {
            return exp(imag*pi*i*j/8)/(4);
        });
        unitary = {FT, tensor_product(T, I, I, I)};
        measure = {direct_sum(O, I, O, I, O, I, O, I), direct_sum(I, O, I, O, I, O, I, O)};
        state1 = matrix::zeros(16, 1);
        state1[0][0] = 1;
        state2 = matrix::zeros(16, 1);
        state2[8][0] = 1;
        rho1 = state1*state1.hermitian();
        rho2 = state2*state2.hermitian();
        machine = {2, 2, 16, unitary, measure};
        query.emplace_back(machine, rho1, rho2, equal);

        // 27
        query.emplace_back(machine, rho1, rho2, equal, 1);

        // 28
        query.emplace_back(machine, rho1, rho2, equal, 2);

        // 29
        unitary = {FT, tensor_product(I, T, I, I), tensor_product(I, I, T, I), tensor_product(I, I, I, T)};
        machine = {4, 2, 16, unitary, measure};
        query.emplace_back(machine, rho1, rho2, equal);

        // 30
        unitary = {tensor_product(H, I), SWAP};
        measure = {direct_sum(I, O), direct_sum(O, I)};
        machine = {2, 2, 4, unitary, measure};
        rho1 = matrix{{1, 1, 1, 1},
                {1, 1, 1, 1},
                {1, 1, 1, 1},
                {1, 1, 1, 1}}/4;
        rho2 = matrix{{1, -1, 1, -1},
                {-1, 1, -1, 1},
                {1, -1, 1, -1},
                {-1, 1, -1, 1}}/4;
        query.emplace_back(machine, rho1, rho2, equal);

        // 31 -- from 26
        FT = matrix(16, 16, [](int i, int j)
        {
            return exp(imag*pi*i*j/8)/(4);
        });
        unitary = {FT, tensor_product(T, I, I, I)};
        measure = {direct_sum(O, I, O, I, O, I, O, I), direct_sum(I, O, I, O, I, O, I, O)};
        state1 = matrix::zeros(16, 1);
        state1[0][0] = 1;
        state2 = matrix::zeros(16, 1);
        state2[4][0] = 1;
        rho1 = state1*state1.hermitian();
        rho2 = state2*state2.hermitian();
        machine = {2, 2, 16, unitary, measure};
        query.emplace_back(machine, rho1, rho2, equal);

        // 32
        state1 = matrix::zeros(16, 1);
        state1[2][0] = 1;
        state2 = matrix::zeros(16, 1);
        state2[6][0] = 1;
        rho1 = state1*state1.hermitian();
        rho2 = state2*state2.hermitian();
        machine = {2, 2, 16, unitary, measure};
        query.emplace_back(machine, rho1, rho2, equal);

        // 33 -- from 23
        FT = matrix(8, 8, [](int i, int j)
        {
            return exp(imag*pi*i*j/4)/(2*sqrt(2.));
        });
        unitary = {FT, tensor_product(T, I, I)};
        measure = {direct_sum(O, I, O, I), direct_sum(I, O, I, O)};
        state1 = matrix::zeros(8, 1);
        state1[0][0] = 1;
        state2 = matrix::zeros(8, 1);
        state2[5][0] = 1;
        rho1 = state1*state1.hermitian();
        rho2 = state2*state2.hermitian();
        machine = {2, 2, 8, unitary, measure};
        query.emplace_back(machine, rho1, rho2, equal);

    }
    void solve()
    {
        int test = 0;
        for (const auto &prob : query)
        {
            if (++ test != 17)
            {
                printf("test%d:\n", test);
                /*
                char filename[111];
                sprintf(filename, "E:/liouzhou_101/qt/data/test%03d.rawdata", test);
                prob.raw_save(filename);

                StateEquivalence<double, complex > pp = StateEquivalence<double, complex >::raw_input(filename);

                std::cout << "test" << test << ":" << std::endl;
                long long stamp = clock();
                prob.solve();
                printf("time: %lldms\n", clock()-stamp);
                pp.solve();
                */
                prob.solve_regular_measurement();
            }
        }
    }

    void init_data_for_qmmwqi()
    {
        puts("qmmwqi tests1:");
        AlmostEqual<Real> equal(1e-8);
        int dim_in = 2;
        int dim_s = 8;
        matrix C = H;
        matrix shift0
        {
            {0, 0, 0, 1},
            {1, 0, 0, 0},
            {0, 1, 0, 0},
            {0, 0, 1, 0}
        };
        matrix shift1
        {
            {0, 1, 0, 0},
            {0, 0, 1, 0},
            {0, 0, 0, 1},
            {1, 0, 0, 0}
        };
        matrix S = direct_sum(shift0, shift1);
        matrix detection = matrix::Toffoli(4, 3, 4, 1);
        //std::cout << detection;
        matrix unitary = detection*tensor_product(I, S)*tensor_product(I, C, I, I);
        std::vector<matrix> measure = {tensor_product(matrix{{1, 0}, {0, 0}}, I, I, I), tensor_product(matrix{{0, 0}, {0, 1}}, I, I, I)};
        QuantumMealyMachineWithQuantumInput<Real> machine(dim_in, dim_s, unitary, measure);
        matrix rho000 = {{1, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0}};
        matrix rho001 = {{0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 1, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0}};
        matrix rho010 = {{0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 1, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0}};
        matrix rho011 = {{0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 1, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0}};
        matrix rho100 = {{0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 1, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0}};
        matrix rho110 = {{0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 0, 0, 0, 1, 0},
                         {0, 0, 0, 0, 0, 0, 0, 0}};

        std::vector<matrix> basis = {
            {{1, 0}, {0, 0}},
            {{0, 0}, {0, 1}},
            {{0.5, 0.5}, {0.5, 0.5}},
            {{0.5, -0.5*imag}, {0.5*imag, 0.5}}
        };
/*
        puts("[000-111] vs [000-111]");
        for (int a1 = 0; a1 < 2; ++ a1) for (int a2 = 0; a2 < 2; ++ a2) for (int a3 = 0; a3 < 2; ++ a3)
        {
            matrix rho1 = tensor_product(basis[a1], basis[a2], basis[a3]);
            for (int b1 = 0; b1 < 2; ++ b1) for (int b2 = 0; b2 < 2; ++ b2) for (int b3 = 0; b3 < 2; ++ b3)
            {
                matrix rho2 = tensor_product(basis[b1], basis[b2], basis[b3]);
                long long stamp = clock();
                printf("(%d,", check(machine, rho1, machine, rho2, basis, equal).first);
                printf("%lld) ", clock()-stamp);
                //return;
            }
            puts("");
        }

        puts("[i00-i11] vs [i00-i11]");
        for (int a2 = 0; a2 < 2; ++ a2) for (int a3 = 0; a3 < 2; ++ a3)
        {
            matrix rho1 = tensor_product(basis[3], basis[a2], basis[a3]);
            for (int b2 = 0; b2 < 2; ++ b2) for (int b3 = 0; b3 < 2; ++ b3)
            {
                matrix rho2 = tensor_product(basis[3], basis[b2], basis[b3]);
                printf("%d ", check_slow(machine, rho1, machine, rho2, basis, equal).first);
            }
            puts("");
        }

        return;
*/
        puts("qmmwqi tests2:");
        QuantumMealyMachineWithQuantumInput<Real> machine1(machine);
        C = {{sqrt(2.0)/2, sqrt(2.0)/2*imag}, {sqrt(2.0)/2*imag, sqrt(2.0)/2}};
        unitary = detection*tensor_product(I, S)*tensor_product(I, C, I, I);
        measure = {tensor_product(matrix{{1, 0}, {0, 0}}, I, I, I), tensor_product(matrix{{0, 0}, {0, 1}}, I, I, I)};
        QuantumMealyMachineWithQuantumInput<Real> machine2(dim_in, dim_s, unitary, measure);
/*
        puts("H[000-111] vs Y[000-111]");
        for (int a1 = 0; a1 < 2; ++ a1) for (int a2 = 0; a2 < 2; ++ a2) for (int a3 = 0; a3 < 2; ++ a3)
        {
            matrix rho1 = tensor_product(basis[a1], basis[a2], basis[a3]);
            for (int b1 = 0; b1 < 2; ++ b1) for (int b2 = 0; b2 < 2; ++ b2) for (int b3 = 0; b3 < 2; ++ b3)
            {
                matrix rho2 = tensor_product(basis[b1], basis[b2], basis[b3]);
                printf("%d ", check(machine1, rho1, machine2, rho2, basis, equal).first);
            }
            puts("");
        }
        return;
*/
        puts("qmmwqi tests3:");
        shift0 = matrix(8, 8, [](int i, int j)
        {
            return i%8 == (j+1)%8;
        });
        shift1 = matrix(8, 8, [](int i, int j)
        {
            return (i+1)%8 == j%8;
        });
        S = direct_sum(shift0, shift1);
        detection = matrix(32, 32, [](int i, int j)
        {
            if ((j&7) == 7) j ^= 1<<4;
            return i == j;
        });
        C = H;
        unitary = detection*tensor_product(I, S)*tensor_product(I, C, I, I, I);
        measure = {tensor_product(matrix{{1, 0}, {0, 0}}, I, I, I, I), tensor_product(matrix{{0, 0}, {0, 1}}, I, I, I, I)};
        machine1 = {2, 16, unitary, measure};

        C = {{sqrt(2.0)/2, sqrt(2.0)/2*imag}, {sqrt(2.0)/2*imag, sqrt(2.0)/2}};
        unitary = detection*tensor_product(I, S)*tensor_product(I, C, I, I, I);
        machine2 = {2, 16, unitary, measure};

        /*
        {
            matrix rho1 = tensor_product(basis[0], basis[0], basis[0], basis[0]);
            matrix rho2 = tensor_product(basis[1], basis[1], basis[1], basis[0]);
            long long stamp = clock();
            printf("%d\n", check(machine1, rho1, machine1, rho2, basis, equal).first);
            printf("%lld\n", clock()-stamp);
        }
        */

        puts("qmmwqi tests4: from test26, paper Case 9");
        std::vector<matrix> state = {matrix{{1, 0}, {0, 0}}, matrix{{0, 0}, {0, 1}}};
        basis.clear();
        basis.push_back(tensor_product(state[0], state[0]));
        basis.push_back(tensor_product(state[0], state[1]));
        matrix FT = matrix(16, 16, [](int i, int j)
        {
            return exp(imag*pi*i*j/8)/(4);
        });
        measure.clear();
        for (int a = 0; a <= 1; ++ a)
            for (int b = 0; b <= 1; ++ b)
            {
                measure.push_back(tensor_product(state[a], state[b], I, I, I, I));
            }
        matrix T1 = tensor_product(T, I, I, I);
        matrix U = matrix::CNOT(6, 5, 1)*tensor_product(I, direct_sum(FT, T1));
        machine = {4, 16, U, measure};
        {
            matrix rho1 = tensor_product(state[0], state[0], state[0], state[0]);
            matrix rho2 = tensor_product(state[1], state[0], state[0], state[0]);
            long long stamp = clock();
            auto ret = check_slow(machine, rho1, machine, rho2, basis, equal);
            printf("%d\n", ret.first);
            if (!ret.first) std::cout << ret.second << std::endl;
            printf("%lldus\n", clock()-stamp);
        }
        puts("qmmwqi tests4: from test31, paper Case 10");
        {
            matrix rho1 = tensor_product(state[0], state[0], state[0], state[0]);
            matrix rho2 = tensor_product(state[0], state[1], state[0], state[0]);
            long long stamp = clock();
            auto ret = check_slow(machine, rho1, machine, rho2, basis, equal);
            printf("%d\n", ret.first);
            if (!ret.first) std::cout << ret.second << std::endl;
            printf("%lldus\n", clock()-stamp);
        }
/*
        puts("qmmwqi tests5: from test23");
        basis.clear();
        basis.push_back(tensor_product(state[0], state[0]));
        basis.push_back(tensor_product(state[0], state[1]));
        FT = matrix(8, 8, [](int i, int j)
        {
            return exp(imag*pi*i*j/4)/(2*sqrt(2.));
        });
        measure.clear();
        for (int a = 0; a <= 1; ++ a)
            for (int b = 0; b <= 1; ++ b)
            {
                measure.push_back(tensor_product(state[a], state[b], I, I, I));
            }
        T1 = tensor_product(T, I, I);
        U = matrix::CNOT(5, 4, 1)*tensor_product(I, direct_sum(FT, T1));
        machine = {4, 8, U, measure};
        {
            matrix rho1 = tensor_product(state[0], state[0], state[0]);
            matrix rho2 = tensor_product(state[1], state[0], state[0]);
            long long stamp = clock();
            auto ret = check(machine, rho1, machine, rho2, basis, equal);
            printf("%d\n", ret.first);
            if (!ret.first) std::cout << ret.second << std::endl;
            printf("%lldus\n", clock()-stamp);
        }

        puts("qmmwqi tests6: from test18");
        basis.clear();
        basis.push_back(tensor_product(state[0], state[0]));
        basis.push_back(tensor_product(state[0], state[1]));
        FT = matrix(4, 4, [](int i, int j)
        {
            return exp(imag*pi*i*j/2)/2;
        });
        matrix BPR = CNOT*tensor_product(H, I);
        measure.clear();
        for (int a = 0; a <= 1; ++ a)
            for (int b = 0; b <= 1; ++ b)
            {
                measure.push_back(tensor_product(state[a], state[b], I, I));
            }
        U = matrix::CNOT(4, 3, 1)*tensor_product(I, direct_sum(FT, BPR));
        machine = {4, 4, U, measure};
        {
            matrix rho1 = tensor_product(state[0], state[0]);
            matrix rho2 = tensor_product(state[0], state[1]);
            long long stamp = clock();
            auto ret = check(machine, rho1, machine, rho2, basis, equal);
            printf("%d\n", ret.first);
            if (!ret.first) std::cout << ret.second << std::endl;
            printf("%lldus\n", clock()-stamp);
        }
*/
        puts("qmmwqi tests7: from test11, paper Case 8");
        basis.clear();
        for (int a = 0; a <= 1; ++ a)
            for (int b = 0; b <= 1; ++ b)
                basis.push_back(tensor_product(state[0], state[a], state[b]));
        matrix H1 = tensor_product(H, I, I);
        matrix H2 = tensor_product(I, H, I);
        matrix H3 = tensor_product(I, I, H);
        U = matrix::CNOT(6, 4, 1)*tensor_product(I, direct_sum(H1, H2, H3, Toffoli));

        measure.clear();
        for (int a = 0; a <= 1; ++ a)
            for (int b = 0; b <= 1; ++ b)
                for (int c = 0; c <= 1; ++ c)
                {
                    measure.push_back(tensor_product(state[a], state[b], state[c], I, I, I));
                }
        machine = {8, 8, U, measure};
        {
            matrix rho1 = tensor_product(state[1], state[1], state[0]);
            matrix rho2 = tensor_product(state[1], state[0], state[1]);
            long long stamp = clock();
            auto ret = check_slow(machine, rho1, machine, rho2, basis, equal);
            printf("%d\n", ret.first);
            if (!ret.first) std::cout << ret.second << std::endl;
            printf("%lldus\n", clock()-stamp);
        }
        puts("qmmwqi tests7: from test11, paper Case 7");
        {
            matrix rho1 = tensor_product(state[0], state[0], state[0]);
            matrix rho2 = tensor_product(state[0], state[0], state[1]);
            long long stamp = clock();
            auto ret = check_slow(machine, rho1, machine, rho2, basis, equal);
            printf("%d\n", ret.first);
            if (!ret.first) std::cout << ret.second << std::endl;
            printf("%lldus\n", clock()-stamp);
        }
    }
}

#endif // DATAFACTORY_H
