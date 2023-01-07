// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "model.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include "utils.h"

namespace ipx {

void Model::GetInfo(Info *info) const {
    info->num_rows_solver = num_rows_;
    info->num_cols_solver = num_cols_ + num_rows_; // including slack columns
    info->num_entries_solver = AI_.entries();
    info->dense_cols = num_dense_cols();
}

void Model::clear() {
    dualized_ = false;
    num_rows_ = 0;
    num_cols_ = 0;
    num_dense_cols_ = 0;
    nz_dense_ = 0;
    AI_.clear();
    AIt_.clear();
    b_.resize(0);
    c_.resize(0);
    lb_.resize(0);
    ub_.resize(0);
    norm_bounds_ = 0.0;
    norm_c_ = 0.0;
}

void Model::FindDenseColumns() {
    num_dense_cols_ = 0;
    nz_dense_ = rows() + 1;

    std::vector<Int> colcount(num_cols_);
    for (Int j = 0; j < num_cols_; j++)
        colcount[j] = AI_.end(j)-AI_.begin(j);
    std::sort(colcount.begin(), colcount.end());

    for (Int j = 1; j < num_cols_; j++) {
        if (colcount[j] > std::max(40l, 10l*colcount[j-1])) {
            // j is the first dense column
            num_dense_cols_ = num_cols_ - j;
            nz_dense_ = colcount[j];
            break;
        }
    }

    if (num_dense_cols_ > 1000) {
        num_dense_cols_ = 0;
        nz_dense_ = rows() + 1;
    }
}

void Model::ComputeNorms() {
    norm_c_ = Infnorm(c_);
    norm_bounds_ = Infnorm(b_);
    for (double x : lb_)
        if (std::isfinite(x))
            norm_bounds_ = std::max(norm_bounds_, std::abs(x));
    for (double x : ub_)
        if (std::isfinite(x))
            norm_bounds_ = std::max(norm_bounds_, std::abs(x));
}

double PrimalInfeasibility(const Model& model, const Vector& x) {
    const Vector& lb = model.lb();
    const Vector& ub = model.ub();
    assert(x.size() == lb.size());

    double infeas = 0.0;
    for (Int j = 0; j < x.size(); j++) {
        infeas = std::max(infeas, lb[j]-x[j]);
        infeas = std::max(infeas, x[j]-ub[j]);
    }
    return infeas;
}

double DualInfeasibility(const Model& model, const Vector& x,
                                const Vector& z) {
    const Vector& lb = model.lb();
    const Vector& ub = model.ub();
    assert(x.size() == lb.size());
    assert(z.size() == lb.size());

    double infeas = 0.0;
    for (Int j = 0; j < x.size(); j++) {
        if (x[j] > lb[j])
            infeas = std::max(infeas, z[j]);
        if (x[j] < ub[j])
            infeas = std::max(infeas, -z[j]);
    }
    return infeas;
}

double PrimalResidual(const Model& model, const Vector& x) {
    const SparseMatrix& AIt = model.AIt();
    const Vector& b = model.b();
    assert(x.size() == AIt.rows());

    double res = 0.0;
    for (Int i = 0; i < b.size(); i++) {
        double r = b[i] - DotColumn(AIt, i, x);
        res = std::max(res, std::abs(r));
    }
    return res;
}

double DualResidual(const Model& model, const Vector& y, const Vector& z) {
    const SparseMatrix& AI = model.AI();
    const Vector& c = model.c();
    assert(y.size() == AI.rows());
    assert(z.size() == AI.cols());

    double res = 0.0;
    for (Int j = 0; j < c.size(); j++) {
        double r = c[j] - z[j] - DotColumn(AI, j, y);
        res = std::max(res, std::abs(r));
    }
    return res;
}

}  // namespace ipx
