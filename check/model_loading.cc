// Copyright (c) 2023 ERGO-Code. See license.txt for license.

#include <cmath>
#include "catch.hpp"
#include "lp_solver.h"

using Int = ipxint;

TEST_CASE("argument null", "[model_loading]") {
    ipx::LpSolver lps;

    const Int num_var = 1;
    const Int num_constr = 0;
    const double obj[] = {0};
    const double lb[] = {0};
    const double ub[] = {1};
    const Int Ap[] = {0, 0};
    const Int Ai[] = {};
    const double Ax[] = {};
    const double rhs[] = {};
    const char constr_type[] = {};

    SECTION("model ok") {
        Int errflag = lps.LoadModel(num_var, obj, lb, ub, num_constr, Ap, Ai,
                                    Ax, rhs, constr_type);
        REQUIRE(errflag == 0);
    }
    SECTION("obj missing") {
        Int errflag = lps.LoadModel(num_var, nullptr, lb, ub, num_constr, Ap,
                                    Ai, Ax, rhs, constr_type);
        REQUIRE(errflag == IPX_ERROR_argument_null);
    }
    SECTION("lb missing") {
        Int errflag = lps.LoadModel(num_var, obj, nullptr, ub, num_constr, Ap,
                                    Ai, Ax, rhs, constr_type);
        REQUIRE(errflag == IPX_ERROR_argument_null);
    }
    SECTION("ub missing") {
        Int errflag = lps.LoadModel(num_var, obj, lb, nullptr, num_constr, Ap,
                                    Ai, Ax, rhs, constr_type);
        REQUIRE(errflag == IPX_ERROR_argument_null);
    }
    SECTION("Ap missing") {
        Int errflag = lps.LoadModel(num_var, obj, lb, ub, num_constr, nullptr,
                                    Ai, Ax, rhs, constr_type);
        REQUIRE(errflag == IPX_ERROR_argument_null);
    }
    SECTION("Ai missing") {
        Int errflag = lps.LoadModel(num_var, obj, lb, ub, num_constr, Ap,
                                    nullptr, Ax, rhs, constr_type);
        REQUIRE(errflag == IPX_ERROR_argument_null);
    }
    SECTION("Ax missing") {
        Int errflag = lps.LoadModel(num_var, obj, lb, ub, num_constr, Ap, Ai,
                                    nullptr, rhs, constr_type);
        REQUIRE(errflag == IPX_ERROR_argument_null);
    }
    SECTION("rhs missing") {
        Int errflag = lps.LoadModel(num_var, obj, lb, ub, num_constr, Ap, Ai,
                                    Ax, nullptr, constr_type);
        REQUIRE(errflag == IPX_ERROR_argument_null);
    }
    SECTION("constr_type missing") {
        Int errflag = lps.LoadModel(num_var, obj, lb, ub, num_constr, Ap, Ai,
                                    Ax, rhs, nullptr);
        REQUIRE(errflag == IPX_ERROR_argument_null);
    }
}

TEST_CASE("invalid dimension", "[model_loading]") {
    ipx::LpSolver lps;

    const Int num_var = 1;
    const Int num_constr = 0;
    const double obj[] = {0};
    const double lb[] = {0};
    const double ub[] = {1};
    const Int Ap[] = {0, 0};
    const Int Ai[] = {};
    const double Ax[] = {};
    const double rhs[] = {};
    const char constr_type[] = {};

    SECTION("model ok") {
        Int errflag = lps.LoadModel(num_var, obj, lb, ub, num_constr, Ap, Ai,
                                    Ax, rhs, constr_type);
        REQUIRE(errflag == 0);
    }
    SECTION("num variables is zero") {
        const Int num_var = 0;
        const Int num_constr = 0;
        Int errflag = lps.LoadModel(num_var, obj, lb, ub, num_constr, Ap, Ai,
                                    Ax, rhs, constr_type);
        REQUIRE(errflag == IPX_ERROR_invalid_dimension);
    }
    SECTION("num variables is negative") {
        const Int num_var = -1;
        const Int num_constr = 0;
        Int errflag = lps.LoadModel(num_var, obj, lb, ub, num_constr, Ap, Ai,
                                    Ax, rhs, constr_type);
        REQUIRE(errflag == IPX_ERROR_invalid_dimension);
    }
    SECTION("num constraints is negative") {
        const Int num_var = 1;
        const Int num_constr = -1;
        Int errflag = lps.LoadModel(num_var, obj, lb, ub, num_constr, Ap, Ai,
                                    Ax, rhs, constr_type);
        REQUIRE(errflag == IPX_ERROR_invalid_dimension);
    }
}

TEST_CASE("invalid matrix", "[model_loading]") {
    ipx::LpSolver lps;

    const Int num_var = 2;
    const Int num_constr = 2;
    const double obj[] = {0, 0};
    const double lb[] = {0, 0};
    const double ub[] = {1, 1};
    const Int Ap[] = {0, 1, 2};
    const Int Ai[] = {0, 1};
    const double Ax[] = {1, 2};
    const double rhs[] = {0, 0};
    const char constr_type[] = {'=', '='};

    SECTION("model ok") {
        Int errflag = lps.LoadModel(num_var, obj, lb, ub, num_constr, Ap, Ai,
                                    Ax, rhs, constr_type);
        REQUIRE(errflag == 0);
    }
    SECTION("column pointers do not start at zero") {
        const Int Ap[] = {1, 2, 2};
        Int errflag = lps.LoadModel(num_var, obj, lb, ub, num_constr, Ap, Ai,
                                    Ax, rhs, constr_type);
        REQUIRE(errflag == IPX_ERROR_invalid_matrix);
    }
    SECTION("column pointers are decreasing") {
        const Int Ap[] = {0, 1, 0};
        Int errflag = lps.LoadModel(num_var, obj, lb, ub, num_constr, Ap, Ai,
                                    Ax, rhs, constr_type);
        REQUIRE(errflag == IPX_ERROR_invalid_matrix);
    }
    SECTION("out of bound index") {
        const Int Ai[] = {1, 2};
        Int errflag = lps.LoadModel(num_var, obj, lb, ub, num_constr, Ap, Ai,
                                    Ax, rhs, constr_type);
        REQUIRE(errflag == IPX_ERROR_invalid_matrix);
    }
    SECTION("duplicate index") {
        const Int Ap[] = {0, 2, 2};
        const Int Ai[] = {0, 0};
        Int errflag = lps.LoadModel(num_var, obj, lb, ub, num_constr, Ap, Ai,
                                    Ax, rhs, constr_type);
        REQUIRE(errflag == IPX_ERROR_invalid_matrix);
    }
    SECTION("jumbled indices") {
        const Int Ap[] = {0, 2, 2};
        const Int Ai[] = {1, 0};
        Int errflag = lps.LoadModel(num_var, obj, lb, ub, num_constr, Ap, Ai,
                                    Ax, rhs, constr_type);
        REQUIRE(errflag == 0);
    }
    SECTION("non-finite matrix entry") {
        const double Ax[] = {1, INFINITY};
        Int errflag = lps.LoadModel(num_var, obj, lb, ub, num_constr, Ap, Ai,
                                    Ax, rhs, constr_type);
        REQUIRE(errflag == IPX_ERROR_invalid_matrix);
    }
}

TEST_CASE("invalid vector", "[model_loading]") {
    ipx::LpSolver lps;

    const Int num_var = 1;
    const Int num_constr = 1;
    const double obj[] = {0};
    const double lb[] = {0};
    const double ub[] = {1};
    const Int Ap[] = {0, 0};
    const Int Ai[] = {};
    const double Ax[] = {};
    const double rhs[] = {0};
    const char constr_type[] = {'<'};

    SECTION("model ok") {
        Int errflag = lps.LoadModel(num_var, obj, lb, ub, num_constr, Ap, Ai,
                                    Ax, rhs, constr_type);
        REQUIRE(errflag == 0);
    }
    SECTION("non-finite objective") {
        const double obj[] = {INFINITY};
        Int errflag = lps.LoadModel(num_var, obj, lb, ub, num_constr, Ap, Ai,
                                    Ax, rhs, constr_type);
        REQUIRE(errflag == IPX_ERROR_invalid_vector);
    }
    SECTION("non-finite rhs") {
        const double rhs[] = {INFINITY};
        Int errflag = lps.LoadModel(num_var, obj, lb, ub, num_constr, Ap, Ai,
                                    Ax, rhs, constr_type);
        REQUIRE(errflag == IPX_ERROR_invalid_vector);
    }
    SECTION("lb positive infinity") {
        const double lb[] = {INFINITY};
        Int errflag = lps.LoadModel(num_var, obj, lb, ub, num_constr, Ap, Ai,
                                    Ax, rhs, constr_type);
        REQUIRE(errflag == IPX_ERROR_invalid_vector);
    }
    SECTION("ub negative infinity") {
        const double ub[] = {-INFINITY};
        Int errflag = lps.LoadModel(num_var, obj, lb, ub, num_constr, Ap, Ai,
                                    Ax, rhs, constr_type);
        REQUIRE(errflag == IPX_ERROR_invalid_vector);
    }
    SECTION("lb greater than ub") {
        const double lb[] = {1};
        const double ub[] = {0};
        Int errflag = lps.LoadModel(num_var, obj, lb, ub, num_constr, Ap, Ai,
                                    Ax, rhs, constr_type);
        REQUIRE(errflag == IPX_ERROR_invalid_vector);
    }
    SECTION("unknown constr_type") {
        const char constr_type[] = {'E'};
        Int errflag = lps.LoadModel(num_var, obj, lb, ub, num_constr, Ap, Ai,
                                    Ax, rhs, constr_type);
        REQUIRE(errflag == IPX_ERROR_invalid_vector);
    }
}
