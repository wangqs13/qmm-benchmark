#include <bits/stdc++.h>

#include "complex.h"
#include "matrix.h"
#include "almost.h"
#include "orthogonal.h"
#include "qmm.h"
#include "qmmwqi.h"
#include "datafactory.h"

// using namespace std;

using std::cout;
using std::endl;

typedef double Real;
typedef Complex<Real> complex;
typedef Matrix<Real, complex> matrix;

const complex I(0, 1);
const Real pi = 3.1415926535897932384626433;

int main()
{
/*
    freopen("try.in", "r", stdin);
    freopen("try.out", "w", stdout);
*/

    std::cout << std::setprecision(10);
    DataFactory::init_data();
    //DataFactory::solve();

    DataFactory::init_data_for_qmmwqi();

    return 0;
}
