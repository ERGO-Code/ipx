// Copyright (c) 2023 ERGO-Code. See license.txt for license.

#ifndef IPX_SOLUTION_H_
#define IPX_SOLUTION_H_

#include <vector>
#include "ipx_internal.h"

namespace ipx {

// Interior point corresponding to user model.
struct InteriorSolution {
    Vector x;                   // size num_var primal solution
    Vector xl;                  // size num_var lower bound slack
    Vector xu;                  // size num_var upper bound slack
    Vector slack;               // size num_constr row slack
    Vector y;                   // size num_constr row dual
    Vector zl;                  // size num_var lower bound dual
    Vector zu;                  // size num_var upper bound dual

    InteriorSolution() : InteriorSolution(0, 0) {};
    explicit InteriorSolution(Int num_var, int num_constr)
        : x(num_var), xl(num_var), xu(num_var), slack(num_constr),
        y(num_constr), zl(num_var), zu(num_var)
        {}
};

// Basic point corresponding to user model.
struct BasicSolution {
    Vector x;                   // size num_var primal solution
    Vector slack;               // size num_constr row slack
    Vector y;                   // size num_constr row dual
    Vector z;                   // size num_var bound dual
    std::vector<Int> vbasis;    // size num_var basic statuses of columns
    std::vector<Int> cbasis;    // size num_constr basic statuses of rows

    BasicSolution() : BasicSolution(0, 0) {};
    explicit BasicSolution(Int num_var, int num_constr)
        : x(num_var), slack(num_constr), y(num_constr), z(num_var),
        vbasis(num_var), cbasis(num_constr)
        {}
};

}  // namespace ipx

#endif  // IPX_SOLUTION_H_
