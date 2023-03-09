// Copyright (c) 2023 ERGO-Code. See license.txt for license.

#include <cmath>
#include "catch.hpp"
#include "lp_solver.h"
#include "test_model.h"

namespace {
    constexpr bool kShowOutput = false;

    void set_parameters(ipx::LpSolver& solver) {
        ipx::Parameters parameters = solver.GetParameters();
        parameters.display = kShowOutput;
        parameters.debug = 3; // show more debug output when display is true
        parameters.crossover = 1;
        solver.SetParameters(parameters);
    }

    void check_sign_conditions(const ipx::InteriorSolution& sol,
                               const TestModel& model) {
        const int num_var = model.obj.size();
        const int num_constr = model.rhs.size();
        for (int j = 0; j != num_var; ++j) {
            if (std::isfinite(model.lb[j])) {
                REQUIRE(std::isfinite(sol.xl[j]));
                REQUIRE(std::isfinite(sol.zl[j]));
                REQUIRE(sol.zl[j] >= 0.0);
            } else {
                REQUIRE(sol.xl[j] == INFINITY);
                REQUIRE(sol.zl[j] == 0.0);
            }
            if (std::isfinite(model.ub[j])) {
                REQUIRE(std::isfinite(sol.xu[j]));
                REQUIRE(std::isfinite(sol.zu[j]));
                REQUIRE(sol.zu[j] >= 0.0);
            } else {
                REQUIRE(sol.xu[j] == INFINITY);
                REQUIRE(sol.zu[j] == 0.0);
            }
        }
        for (int i = 0; i != num_constr; ++i) {
            switch (model.constr_type[i]) {
            case '=':
                REQUIRE(sol.slack[i] == 0.0);
                break;
            case '<':
                REQUIRE(sol.slack[i] >= 0.0);
                REQUIRE(sol.y[i] <= 0.0);
                break;
            case '>':
                REQUIRE(sol.slack[i] <= 0.0);
                REQUIRE(sol.y[i] >= 0.0);
                break;
            default:
                REQUIRE(false);
            }
        }
    }

    void check_basis_consistency(const ipx::BasicSolution& sol,
                                 const TestModel& model) {
        const int num_var = model.obj.size();
        const int num_constr = model.rhs.size();
        for (int j = 0; j != num_var; ++j) {
            switch (sol.vbasis[j]) {
            case IPX_basic:
                REQUIRE(sol.z[j] == 0.0);
                break;
            case IPX_nonbasic_lb:
                REQUIRE(std::isfinite(model.lb[j]));
                REQUIRE(sol.x[j] == model.lb[j]);
                break;
            case IPX_nonbasic_ub:
                REQUIRE(std::isfinite(model.ub[j]));
                REQUIRE(sol.x[j] == model.ub[j]);
                break;
            case IPX_superbasic:
                REQUIRE(std::isinf(model.lb[j]));
                REQUIRE(std::isinf(model.ub[j]));
                REQUIRE(sol.z[j] == 0.0);
                break;
            default:
                REQUIRE(false);
            }
        }
        for (int i = 0; i != num_constr; ++i) {
            switch (sol.cbasis[i]) {
            case IPX_basic:
                REQUIRE(sol.y[i] == 0.0);
                break;
            case IPX_nonbasic:
                REQUIRE(sol.slack[i] == 0.0);
                break;
            default:
                REQUIRE(false);
            }
        }
    }

    void check_interior_solution(const ipx::LpSolver& solver,
                                 const TestModel& model) {
        ipx::InteriorSolution sol;
        ipxint errflag = solver.GetInteriorSolution(sol);
        REQUIRE(errflag == 0);
        check_sign_conditions(sol, model);
    }

    void check_basic_solution(const ipx::LpSolver& solver,
                              const TestModel& model) {
        ipx::BasicSolution sol;
        ipxint errflag = solver.GetBasicSolution(sol);
        REQUIRE(errflag == 0);
        check_basis_consistency(sol, model);
    }

    void check_solution(const ipx::LpSolver& solver, const TestModel& model) {
        ipx::Info info = solver.GetInfo();
        if (info.status_ipm != IPX_STATUS_not_run)
            check_interior_solution(solver, model);
        if (info.status_crossover == IPX_STATUS_optimal ||
            info.status_crossover == IPX_STATUS_imprecise)
            check_basic_solution(solver, model);
    }

    void solve(ipx::LpSolver& solver, const TestModel& model,
               ipxint required_status_ipm, ipxint required_status_crossover) {
        ipxint errflag = load(solver, model);
        REQUIRE(errflag == 0);

        SECTION("no dualize") {
            ipx::Parameters parameters = solver.GetParameters();
            parameters.dualize = 0;
            solver.SetParameters(parameters);
            solver.Solve();
            ipx::Info info = solver.GetInfo();
            REQUIRE(info.status_ipm == required_status_ipm);
            REQUIRE(info.status_crossover == required_status_crossover);
            check_solution(solver, model);
        }
        SECTION("dualize") {
            ipx::Parameters parameters = solver.GetParameters();
            parameters.dualize = 1;
            solver.SetParameters(parameters);
            solver.Solve();
            ipx::Info info = solver.GetInfo();
            REQUIRE(info.status_ipm == required_status_ipm);
            REQUIRE(info.status_crossover == required_status_crossover);
            check_solution(solver, model);
        }
    }
} // namespace

TEST_CASE("no constraints", "[solver]") {
    ipx::LpSolver solver;
    set_parameters(solver);

    TestModel model = TestModel(2, 0, '=');

    SECTION("boxed variables") {
        solve(solver, model, IPX_STATUS_optimal, IPX_STATUS_optimal);
    }
    SECTION("lower bounded variables") {
        std::fill(model.ub.begin(), model.ub.end(), INFINITY);
        solve(solver, model, IPX_STATUS_optimal, IPX_STATUS_optimal);
    }
    SECTION("upper bounded variables") {
        std::fill(model.lb.begin(), model.lb.end(), -INFINITY);
        solve(solver, model, IPX_STATUS_optimal, IPX_STATUS_optimal);
    }
    SECTION("free variables") {
        std::fill(model.lb.begin(), model.lb.end(), -INFINITY);
        std::fill(model.ub.begin(), model.ub.end(), INFINITY);
        solve(solver, model, IPX_STATUS_optimal, IPX_STATUS_optimal);
    }
    SECTION("fixed variables") {
        std::fill(model.ub.begin(), model.ub.end(), 0.0);
        solve(solver, model, IPX_STATUS_optimal, IPX_STATUS_optimal);
    }
    SECTION("mixed variables") {
        TestModel model = TestModel(5, 0, '=');
        model.ub = {1.0, INFINITY, -1.0, INFINITY, 1.0};
        model.lb = {0.0, 1.0, -INFINITY, -INFINITY, 1.0};
        solve(solver, model, IPX_STATUS_optimal, IPX_STATUS_optimal);
    }
}

TEST_CASE("single constraint", "[solver]") {
    ipx::LpSolver solver;
    set_parameters(solver);

    TestModel model = TestModel(0, 1, '=')
        .add_column(1.0, {0}, {1.0}, 0.0, 1.0)            // boxed
        .add_column(1.0, {0}, {2.0}, 1.0, INFINITY)       // lower bounded
        .add_column(-1.0, {0}, {3.0}, -INFINITY, -1.0)    // upper bounded
        .add_column(0.0, {0}, {4.0}, -INFINITY, INFINITY) // free
        .add_column(-1.0, {0}, {5.0}, 1.0, 1.0);          // fixed
    model.rhs[0] = 0.5;

    for (char ctype : {'=', '>', '<'}) {
        model.constr_type[0] = ctype;
        DYNAMIC_SECTION(ctype << "-constraint") {
            solve(solver, model, IPX_STATUS_optimal, IPX_STATUS_optimal);
        }
    }
}

TEST_CASE("dependent equality constraints", "[solver]") {
    ipx::LpSolver solver;
    set_parameters(solver);

    // We only want to test that basis construction detects the dependencies and
    // possibly infeasibility. Skip initial iterations.
    ipx::Parameters parameters = solver.GetParameters();
    parameters.switchiter = 0;
    solver.SetParameters(parameters);

    TestModel model = TestModel(0, 2, '=')
        .add_column(0, {0, 1}, {1.0, 2.0}, 0, 1)
        .add_column(0, {0, 1}, {1.0, 2.0}, 0, 1);

    SECTION("consistent rhs") {
        solve(solver, model, IPX_STATUS_optimal, IPX_STATUS_optimal);
    }
    SECTION("inconsistent rhs") {
        model.rhs[0] = 1;
        solve(solver, model, IPX_STATUS_primal_infeas, IPX_STATUS_not_run);
    }
}

TEST_CASE("dependent free variables", "[solver]") {
    ipx::LpSolver solver;
    set_parameters(solver);

    // We only want to test that basis construction detects the dependencies and
    // possibly infeasibility. Skip initial iterations.
    ipx::Parameters parameters = solver.GetParameters();
    parameters.switchiter = 0;
    solver.SetParameters(parameters);

    TestModel model = TestModel(0, 2, '=')
        .add_column(1, {0, 1}, {1.0, 1.0}, -INFINITY, INFINITY)
        .add_column(2, {0, 1}, {2.0, 2.0}, -INFINITY, INFINITY);

    SECTION("consistent obj") {
        solve(solver, model, IPX_STATUS_optimal, IPX_STATUS_optimal);
    }
    SECTION("inconsistent obj") {
        model.obj[0] = -1;
        solve(solver, model, IPX_STATUS_dual_infeas, IPX_STATUS_not_run);
    }
}
