#ifndef QMM_H
#define QMM_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <queue>

#include "matrix.h"
#include "orthogonal.h"
#include "almost.h"

template<class T>
std::vector<T> operator + (std::vector<T> a, const T &b)
{
    a.push_back(b);
    return a;
}

template<class V>
// V = float, double
class QuantumMealyMachine
{
public:
    typedef Matrix<V, Complex<V> > matrix;

    class Configuration
    {
    public:
        std::vector<int> a, S, b;
        matrix rho;
        Configuration(const std::vector<int> &a, const std::vector<int> &S, const std::vector<int> &b, const matrix &rho) : a(a), S(S), b(b), rho(rho)
        {
        }
        ~Configuration()
        {
        }
        friend std::ostream &operator << (std::ostream &out, const Configuration &config)
        {
            out << "a = ";
            for (auto x : config.a) out << x << " ";
            out << std::endl;
            out << "S = ";
            for (auto x : config.S) out << x << " ";
            out << std::endl;
            out << "b = ";
            for (auto x : config.b) out << x << " ";
            out << std::endl;
            //out << config.rho;
            return out;
        }
    };

    int sigma, gamma, dim;
    std::vector<matrix> unitary;
    std::vector<matrix> measure;
    QuantumMealyMachine(int sigma, int gamma, int dim, const std::vector<matrix> &unitary, const std::vector<matrix> &measure) : sigma(sigma), gamma(gamma), dim(dim), unitary(unitary), measure(measure)
    {
        assert(sigma == unitary.size());
        assert(gamma == measure.size());
        for (int i = 0; i < sigma; ++ i)
            assert(unitary[i].is_dimension(dim));
        for (int i = 0; i < gamma; ++ i)
            assert(measure[i].is_dimension(dim));
    }
    QuantumMealyMachine(const QuantumMealyMachine &m) : sigma(m.sigma), gamma(m.gamma), dim(m.dim), unitary(m.unitary), measure(m.measure)
    {
    }

    ~QuantumMealyMachine()
    {
    }

    std::pair<bool, Configuration> check(const matrix &rho1, const matrix &rho2, AlmostEqual<V> equal) const
    {
        matrix rho = rho1-rho2;
        std::queue<Configuration> Q;
        Q.push(Configuration({}, {}, {}, rho));
        Orthogonal<V, matrix> ortho(equal);
        while (!Q.empty())
        {
            Configuration config = Q.front();
            Q.pop();
            if (ortho.is_independent(config.rho))
            {
                //std::cout << config << config.rho.trace().x << std::endl;
                if (!equal(config.rho.trace().x, 0))
                    return std::pair<bool, Configuration>(false, config);
                ortho.add(config.rho);
                for (int a = 0; a < sigma; ++ a)
                    Q.emplace(config.a+a, config.S, config.b, unitary[a]*config.rho*unitary[a].hermitian());
                int len = config.a.size();
                for (int b = 0; b < gamma; ++ b)
                    Q.emplace(config.a, config.S+len, config.b+b, measure[b]*config.rho*measure[b].hermitian());
            }
        }
        return std::pair<bool, Configuration>(true, Configuration({}, {}, {}, rho));
    }

    std::pair<bool, Configuration> check_regular_measurement(const matrix &rho1, const matrix &rho2, AlmostEqual<V> equal) const
    {
        matrix rho = rho1-rho2;
        std::queue<Configuration> Q;
        Q.push(Configuration({}, {}, {}, rho));
        Orthogonal<V, matrix> ortho(equal);
        while (!Q.empty())
        {
            Configuration config = Q.front();
            Q.pop();
            if (ortho.is_independent(config.rho))
            {
                //std::cout << config << " " << config.rho.trace().x << std::endl;
                if (!equal(config.rho.trace().x, 0))
                    return std::pair<bool, Configuration>(false, config);
                ortho.add(config.rho);
                int len = config.a.size();
                for (int a = 0; a < sigma; ++ a)
                {
                    for (int b = 0; b < gamma; ++ b)
                        Q.emplace(config.a+a, config.S+len, config.b+b, measure[b]*unitary[a]*config.rho*unitary[a].hermitian()*measure[b].hermitian());
                }
            }
        }
        return std::pair<bool, Configuration>(true, Configuration({}, {}, {}, rho));
    }

    std::pair<bool, Configuration> check_slow(const matrix &rho1, const matrix &rho2, AlmostEqual<V> equal) const
    {
        matrix rho = rho1-rho2;
        std::queue<Configuration> Q;
        Q.push(Configuration({}, {}, {}, rho));
        Orthogonal<V, matrix> ortho(equal);
        while (!Q.empty())
        {
            Configuration config = Q.front();
            Q.pop();
            if (ortho.is_independent_slow(config.rho))
            {
                //std::cout << config << config.rho.trace().x << std::endl;
                if (!equal(config.rho.trace().x, 0))
                    return std::pair<bool, Configuration>(false, config);
                ortho.add_slow(config.rho);
                for (int a = 0; a < sigma; ++ a)
                    Q.emplace(config.a+a, config.S, config.b, unitary[a]*config.rho*unitary[a].hermitian());
                int len = config.a.size();
                for (int b = 0; b < gamma; ++ b)
                    Q.emplace(config.a, config.S+len, config.b+b, measure[b]*config.rho*measure[b].hermitian());
            }
        }
        return std::pair<bool, Configuration>(true, Configuration({}, {}, {}, rho));
    }

    std::pair<bool, Configuration> check_measure_limit(const matrix &rho1, const matrix &rho2, AlmostEqual<V> equal, int limit) const
    {
        if (limit >= dim*dim-1)
            return check(rho1, rho2, equal);
        matrix rho = rho1-rho2;
        std::vector<Orthogonal<V, matrix> > ortho(limit+1, equal);
        std::queue<Configuration> Q;
        Q.push(Configuration({}, {}, {}, rho));
        while (!Q.empty())
        {
            Configuration config = Q.front();
            Q.pop();
            int i = config.S.size();
            if (ortho[i].is_independent(config.rho))
            {
                //std::cout << config;
                if (!equal(config.rho.trace().x, 0))
                    return std::pair<bool, Configuration>(false, config);
                for (int a = 0; a < sigma; ++ a)
                    Q.emplace(config.a+a, config.S, config.b, unitary[a]*config.rho*unitary[a].hermitian());
                if (i < limit)
                {
                    int pos = config.a.size();
                    for (int b = 0; b < gamma; ++ b)
                        Q.emplace(config.a, config.S+pos, config.b+b, measure[b]*config.rho*measure[b].hermitian());
                }
                for (; i <= limit && ortho[i].is_independent(config.rho); ++ i)
                    ortho[i].add(config.rho);
            }
        }
        return std::pair<bool, Configuration>(true, Configuration({}, {}, {}, rho));
    }

    std::pair<bool, Configuration> check_measure_limit_slow(const matrix &rho1, const matrix &rho2, AlmostEqual<V> equal, int limit) const
    {
        if (limit >= dim*dim-1)
            return check(rho1, rho2, equal);
        matrix rho = rho1-rho2;
        std::vector<Orthogonal<V, matrix> > ortho(limit+1, equal);
        std::queue<Configuration> Q;
        Q.push(Configuration({}, {}, {}, rho));
        while (!Q.empty())
        {
            Configuration config = Q.front();
            Q.pop();
            int i = config.S.size();
            if (ortho[i].is_independent_slow(config.rho))
            {
                //std::cout << config;
                if (!equal(config.rho.trace().x, 0))
                    return std::pair<bool, Configuration>(false, config);
                for (int a = 0; a < sigma; ++ a)
                    Q.emplace(config.a+a, config.S, config.b, unitary[a]*config.rho*unitary[a].hermitian());
                if (i < limit)
                {
                    int pos = config.a.size();
                    for (int b = 0; b < gamma; ++ b)
                        Q.emplace(config.a, config.S+pos, config.b+b, measure[b]*config.rho*measure[b].hermitian());
                }
                for (; i <= limit && ortho[i].is_independent_slow(config.rho); ++ i)
                    ortho[i].add_slow(config.rho);
            }
        }
        return std::pair<bool, Configuration>(true, Configuration({}, {}, {}, rho));
    }
};

#endif // QMM_H
