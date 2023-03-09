// Copyright (c) 2023 ERGO-Code. See license.txt for license.

#ifndef IPX_TEST_MODEL_H_
#define IPX_TEST_MODEL_H_

#include <initializer_list>
#include <vector>
#include "lp_solver.h"

// Helper class for unit tests.
struct TestModel {
    std::vector<double> obj;
    std::vector<double> lb;
    std::vector<double> ub;
    std::vector<ipxint> Ap;
    std::vector<ipxint> Ai;
    std::vector<double> Ax;
    std::vector<double> rhs;
    std::vector<char> constr_type;

    // Creates model with empty matrix, zero objective and right-hand side, and
    // bounds 0 <= x <= 1.
    TestModel(ipxint num_var, ipxint num_constr, char constr_type) {
        obj.assign(num_var, 0);
        lb.assign(num_var, 0);
        ub.assign(num_var, 1);
        Ap.assign(num_var + 1, 0);
        rhs.assign(num_constr, 0);
        this->constr_type.assign(num_constr, constr_type);
    }

    // Appends column to the model.
    TestModel& add_column(double obj,
                          std::initializer_list<ipxint> indices,
                          std::initializer_list<double> values,
                          double lb, double ub) {
        this->obj.push_back(obj);
        this->lb.push_back(lb);
        this->ub.push_back(ub);
        Ai.insert(Ai.end(), indices.begin(), indices.end());
        Ax.insert(Ax.end(), values.begin(), values.end());
        Ap.push_back(Ai.size());
        return *this;
    }
};

ipxint load(ipx::LpSolver& solver, const TestModel& model) {
    ipxint num_var = static_cast<ipxint>(model.obj.size());
    ipxint num_constr = static_cast<ipxint>(model.rhs.size());
    return solver.LoadModel(
        num_var, model.obj.data(), model.lb.data(), model.ub.data(), num_constr,
        model.Ap.data(), model.Ai.data(), model.Ax.data(), model.rhs.data(),
        model.constr_type.data());
}

#endif  // IPX_TEST_MODEL_H_
