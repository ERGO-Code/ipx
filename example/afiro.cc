// Copyright (c) 2018 ERGO-Code. See license.txt for license.
//
// Example for using IPX from its C++ interface. The program solves the Netlib
// problem afiro.

#include <cmath>
#include <iostream>
#include "lp_solver.h"

using Int = ipxint;

constexpr Int num_var = 12;
constexpr Int num_constr = 9;
const double obj[] = { -0.2194, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.32,
                       -0.5564, 0.6, -0.48 };
const double lb[num_var] = { 0.0 };
const double ub[] = { 80.0, 283.303, 283.303, 312.813, 349.187, INFINITY,
                      INFINITY, INFINITY, 57.201, 500.0, 500.501, 357.501};
// Constraint matrix in CSC format with 0-based indexing.
const Int Ap[] = { 0, 2, 6, 10, 14, 18, 20, 22, 24, 26, 28, 30, 32 };
const Int Ai[] = { 0, 5,
                   1, 6, 7, 8,
                   2, 6, 7, 8,
                   3, 6, 7, 8,
                   4, 6, 7, 8,
                   1, 2,
                   2, 3,
                   2, 4,
                   0, 6,
                   0, 5,
                   2, 5,
                   5, 7 };
const double Ax[] = { -1.0, 0.301,
                      1.0, -1.0, 0.301, 1.06,
                      1.0, -1.0, 0.313, 1.06,
                      1.0, -1.0, 0.313, 0.96,
                      1.0, -1.0, 0.326, 0.86,
                      - 1.0, 0.99078,
                      1.00922, -1.0,
                      1.01802, -1.0,
                      1.4, 1.0,
                      0.109, -1.0,
                      -0.419111, 1.0,
                      1.4, -1.0 };
const double rhs[] = { 0.0, 80.0, 0.0, 0.0, 0.0, 0.0, 0.0, 44.0, 300.0 };
const char constr_type[] = { '<', '<', '=', '<', '<', '=', '<', '<', '<' };

int main() {
    ipx::LpSolver lps;

    // To change parameters from their defaults, create an ipx::Parameters
    // object (which is initialized to default values), make any changes and
    // pass to the LP solver. See the reference documentation for available
    // parameters.
    ipx::Parameters parameters;
    // parameters.crossover = 0;   // turns off crossover
    // parameters.debug = 1;       // sets first debugging level (more output)
    lps.SetParameters(parameters);

    // Load the LP model into IPX.
    Int errflag = lps.LoadModel(num_var, obj, lb, ub, num_constr, Ap, Ai, Ax,
                                rhs, constr_type);
    if (errflag) {
        std::cout << " invalid model (errflag = " << errflag << ")\n";
        return 1;
    }

    // Solve the LP.
    Int status = lps.Solve();
    if (status != IPX_STATUS_solved) {
        // no solution
        // (time/iter limit, numerical failure, out of memory)
        std::cout << " status: " << status << ','
                  << " errflag: " << lps.GetInfo().errflag << '\n';
        return 2;
    }

    // Get solver and solution information.
    ipx::Info info = lps.GetInfo();
    std::cout << "Solver performed " << info.iter << " IPM iterations and "
              << info.updates_crossover << " crossover pivots in "
              << info.time_total << " seconds.\n";

    // Get the interior solution (available if IPM was started).
    double x[num_var], xl[num_var], xu[num_var], slack[num_constr];
    double y[num_constr], zl[num_var], zu[num_var];
    lps.GetInteriorSolution(x, xl, xu, slack, y, zl, zu);

    // Get the basic solution (available if crossover terminated without error).
    double xbasic[num_var], sbasic[num_constr];
    double ybasic[num_constr], zbasic[num_var];
    Int cbasis[num_constr], vbasis[num_var];
    lps.GetBasicSolution(xbasic, sbasic, ybasic, zbasic, cbasis, vbasis);

    return 0;
}
