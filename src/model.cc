// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "model.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include "utils.h"

namespace ipx {

Int Model::Init(const Control& control, const UserModel& user_model) {
    clear();
    user_model_ = &user_model;
    ComputeUserModelAttributes();
    PresolveModel(control);
    return 0;
}

void Model::GetInfo(Info *info) const {
    info->num_rows_solver = num_rows_;
    info->num_cols_solver = num_cols_ + num_rows_; // including slack columns
    info->num_entries_solver = AI_.entries();
    info->dualized = dualized_;
    info->dense_cols = num_dense_cols();
}

void Model::clear() {
    // clear computational form model
    dualized_ = false;
    negated_vars_.clear();
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

    // clear user model
    user_model_ = nullptr;
    num_constr_ = 0;
    num_eqconstr_ = 0;
    num_var_ = 0;
    num_free_var_ = 0;
    boxed_vars_.clear();

    colscale_.resize(0);
    rowscale_.resize(0);
}

void Model::PresolveStartingPoint(const double* x_user,
                                  const double* slack_user,
                                  const double* y_user,
                                  const double* z_user,
                                  Vector& x_solver,
                                  Vector& y_solver,
                                  Vector& z_solver) const {
    const Int m = rows();
    const Int n = cols();
    assert(x_solver.size() == n+m);
    assert(y_solver.size() == m);
    assert(z_solver.size() == n+m);

    Vector x_temp(num_var_);
    Vector slack_temp(num_constr_);
    Vector y_temp(num_constr_);
    Vector z_temp(num_var_);
    if (x_user)
        std::copy_n(x_user, num_var_, std::begin(x_temp));
    if (slack_user)
        std::copy_n(slack_user, num_constr_, std::begin(slack_temp));
    if (y_user)
        std::copy_n(y_user, num_constr_, std::begin(y_temp));
    if (z_user)
        std::copy_n(z_user, num_var_, std::begin(z_temp));
    PresolveGeneralPoint(x_temp, slack_temp, y_temp, z_temp, x_solver, y_solver,
                         z_solver);
}

Int Model::PresolveIPMStartingPoint(const double* x_user,
                                    const double* xl_user,
                                    const double* xu_user,
                                    const double* slack_user,
                                    const double* y_user,
                                    const double* zl_user,
                                    const double* zu_user,
                                    Vector& x_solver,
                                    Vector& xl_solver,
                                    Vector& xu_solver,
                                    Vector& y_solver,
                                    Vector& zl_solver,
                                    Vector& zu_solver) const {
    Int errflag =
        user_model_->CheckInteriorPoint(x_user, xl_user, xu_user, slack_user,
                                        y_user, zl_user, zu_user);
    if (errflag)
        return errflag;

    if (dualized_)
        return IPX_ERROR_not_implemented;

    // Copy user point into workspace.
    Vector x_temp(x_user, num_var_);
    Vector xl_temp(xl_user, num_var_);
    Vector xu_temp(xu_user, num_var_);
    Vector slack_temp(slack_user, num_constr_);
    Vector y_temp(y_user, num_constr_);
    Vector zl_temp(zl_user, num_var_);
    Vector zu_temp(zu_user, num_var_);

    PresolveInteriorPoint(x_temp, xl_temp, xu_temp, slack_temp, y_temp, zl_temp,
                          zu_temp, x_solver, xl_solver, xu_solver, y_solver,
                          zl_solver, zu_solver);
    return 0;
}


void Model::PostsolveInteriorSolution(const Vector& x_solver,
                                      const Vector& xl_solver,
                                      const Vector& xu_solver,
                                      const Vector& y_solver,
                                      const Vector& zl_solver,
                                      const Vector& zu_solver,
                                      double* x_user,
                                      double* xl_user, double* xu_user,
                                      double* slack_user,
                                      double* y_user,
                                      double* zl_user, double* zu_user) const {
    const Int m = rows();
    const Int n = cols();
    assert(x_solver.size() == n+m);
    assert(xl_solver.size() == n+m);
    assert(xu_solver.size() == n+m);
    assert(y_solver.size() == m);
    assert(zl_solver.size() == n+m);
    assert(zu_solver.size() == n+m);

    Vector x_temp(num_var_);
    Vector xl_temp(num_var_);
    Vector xu_temp(num_var_);
    Vector slack_temp(num_constr_);
    Vector y_temp(num_constr_);
    Vector zl_temp(num_var_);
    Vector zu_temp(num_var_);
    PostsolveInteriorPoint(x_solver, xl_solver, xu_solver, y_solver, zl_solver,
                           zu_solver, x_temp, xl_temp, xu_temp, slack_temp,
                           y_temp, zl_temp, zu_temp);
    if (x_user)
        std::copy(std::begin(x_temp), std::end(x_temp), x_user);
    if (xl_user)
        std::copy(std::begin(xl_temp), std::end(xl_temp), xl_user);
    if (xu_user)
        std::copy(std::begin(xu_temp), std::end(xu_temp), xu_user);
    if (slack_user)
        std::copy(std::begin(slack_temp), std::end(slack_temp), slack_user);
    if (y_user)
        std::copy(std::begin(y_temp), std::end(y_temp), y_user);
    if (zl_user)
        std::copy(std::begin(zl_temp), std::end(zl_temp), zl_user);
    if (zu_user)
        std::copy(std::begin(zu_temp), std::end(zu_temp), zu_user);
}

void Model::EvaluateInteriorSolution(const Vector& x_solver,
                                     const Vector& xl_solver,
                                     const Vector& xu_solver,
                                     const Vector& y_solver,
                                     const Vector& zl_solver,
                                     const Vector& zu_solver,
                                     Info* info) const {
    const Int m = rows();
    const Int n = cols();
    assert(x_solver.size() == n+m);
    assert(xl_solver.size() == n+m);
    assert(xu_solver.size() == n+m);
    assert(y_solver.size() == m);
    assert(zl_solver.size() == n+m);
    assert(zu_solver.size() == n+m);

    // Build solution to user model.
    Vector x(num_var_);
    Vector xl(num_var_);
    Vector xu(num_var_);
    Vector slack(num_constr_);
    Vector y(num_constr_);
    Vector zl(num_var_);
    Vector zu(num_var_);
    PostsolveInteriorPoint(x_solver, xl_solver, xu_solver, y_solver, zl_solver,
                           zu_solver, x, xl, xu, slack, y, zl, zu);

    // Evaluate solution to user model.
    user_model_->EvaluateInteriorPoint(x, xl, xu, slack, y, zl, zu, info);
}

void Model::PostsolveBasicSolution(const Vector& x_solver,
                                   const Vector& y_solver,
                                   const Vector& z_solver,
                                   const std::vector<Int>& basic_status_solver,
                                   double* x_user, double* slack_user,
                                   double* y_user, double* z_user) const {
    const Int m = rows();
    const Int n = cols();
    assert(x_solver.size() == n+m);
    assert(y_solver.size() == m);
    assert(z_solver.size() == n+m);
    assert(basic_status_solver.size() == n+m);

    Vector x_temp(num_var_);
    Vector slack_temp(num_constr_);
    Vector y_temp(num_constr_);
    Vector z_temp(num_var_);
    std::vector<Int> cbasis_temp(num_constr_);
    std::vector<Int> vbasis_temp(num_var_);
    PostsolveGeneralPoint(x_solver, y_solver, z_solver, x_temp, slack_temp,
                          y_temp, z_temp);
    PostsolveBasis(basic_status_solver, cbasis_temp, vbasis_temp);
    CorrectBasicSolution(x_temp, slack_temp, y_temp, z_temp, cbasis_temp,
                         vbasis_temp);
    if (x_user)
        std::copy(std::begin(x_temp), std::end(x_temp), x_user);
    if (slack_user)
        std::copy(std::begin(slack_temp), std::end(slack_temp), slack_user);
    if (y_user)
        std::copy(std::begin(y_temp), std::end(y_temp), y_user);
    if (z_user)
        std::copy(std::begin(z_temp), std::end(z_temp), z_user);
}

void Model::EvaluateBasicSolution(const Vector& x_solver,
                                  const Vector& y_solver,
                                  const Vector& z_solver,
                                  const std::vector<Int>& basic_status_solver,
                                  Info* info) const {
    const Int m = rows();
    const Int n = cols();
    assert(x_solver.size() == n+m);
    assert(y_solver.size() == m);
    assert(z_solver.size() == n+m);
    assert(basic_status_solver.size() == n+m);

    // Build basic solution to user model.
    Vector x(num_var_);
    Vector slack(num_constr_);
    Vector y(num_constr_);
    Vector z(num_var_);
    std::vector<Int> cbasis(num_constr_);
    std::vector<Int> vbasis(num_var_);
    PostsolveGeneralPoint(x_solver, y_solver, z_solver, x, slack, y, z);
    PostsolveBasis(basic_status_solver, cbasis, vbasis);
    CorrectBasicSolution(x, slack, y, z, cbasis, vbasis);

    // Evaluate basic solution to user model.
    user_model_->EvaluateBasicPoint(x, slack, y, z, vbasis, cbasis, info);
}

void Model::PostsolveBasis(const std::vector<Int>& basic_status_solver,
                           Int* cbasis_user, Int* vbasis_user) const {
    const Int m = rows();
    const Int n = cols();
    assert(basic_status_solver.size() == n+m);

    std::vector<Int> cbasis_temp(num_constr_);
    std::vector<Int> vbasis_temp(num_var_);
    PostsolveBasis(basic_status_solver, cbasis_temp, vbasis_temp);
    if (cbasis_user)
        std::copy(std::begin(cbasis_temp), std::end(cbasis_temp), cbasis_user);
    if (vbasis_user)
        std::copy(std::begin(vbasis_temp), std::end(vbasis_temp), vbasis_user);
}

void Model::PresolveModel(const Control& control) {
    control.Log()
        << "Input\n"
        << Textline("Number of variables:") << num_var_ << '\n'
        << Textline("Number of free variables:") << num_free_var_ << '\n'
        << Textline("Number of constraints:") << num_constr_ << '\n'
        << Textline("Number of equality constraints:") << num_eqconstr_ << '\n'
        << Textline("Number of matrix entries:")
        << user_model_->A().entries() << '\n';
    PrintCoefficientRange(control);

    // Make an automatic decision for dualization if not specified by user.
    Int dualize = control.dualize();
    if (dualize < 0)
        dualize = num_constr_ > 2*num_var_;
    if (dualize)
        LoadDual();
    else
        LoadPrimal();

    // Scale AI_ and vectors before building AIt_.
    ScaleModel(control);

    AIt_ = Transpose(AI_);
    assert(AI_.begin(num_cols_ + num_rows_) == AIt_.begin(num_rows_));
    FindDenseColumns();
    norm_c_ = Infnorm(c_);
    norm_bounds_ = Infnorm(b_);
    for (double x : lb_)
        if (std::isfinite(x))
            norm_bounds_ = std::max(norm_bounds_, std::abs(x));
    for (double x : ub_)
        if (std::isfinite(x))
            norm_bounds_ = std::max(norm_bounds_, std::abs(x));
    PrintPreprocessingLog(control);
}

void Model::PresolveGeneralPoint(const Vector& x_user,
                                 const Vector& slack_user,
                                 const Vector& y_user,
                                 const Vector& z_user,
                                 Vector& x_solver,
                                 Vector& y_solver,
                                 Vector& z_solver) const {
    const Int m = rows();
    const Int n = cols();

    if (dualized_) {
        assert(num_var_ == m);
        assert(num_constr_ + boxed_vars_.size() == n);

        // Build dual solver variables from primal user variables.
        y_solver = -x_user;
        for (Int j : negated_vars_)
            y_solver[j] *= -1.0;
        for (Int i = 0; i < num_constr_; i++)
            z_solver[i] = -slack_user[i];
        for (Int k = 0; k < boxed_vars_.size(); k++) {
            Int j = boxed_vars_[k];
            z_solver[num_constr_+k] = c(num_constr_+k) + y_solver[j];
        }
        for (Int i = 0; i < m; i++)
            z_solver[n+i] = c(n+i)-y_solver[i];

        // Build primal solver variables from dual user variables.
        std::copy_n(std::begin(y_user), num_constr_, std::begin(x_solver));
        std::copy_n(std::begin(z_user), num_var_, std::begin(x_solver) + n);
        for (Int j : negated_vars_)
            x_solver[n+j] *= -1.0;
        for (Int k = 0; k < boxed_vars_.size(); k++) {
            Int j = boxed_vars_[k];
            if (x_solver[n+j] < 0.0) {
                // j is a boxed variable and z_user[j] < 0
                x_solver[num_constr_+k] = -x_solver[n+j];
                x_solver[n+j] = 0.0;
            } else {
                x_solver[num_constr_+k] = 0.0;
            }
        }
    }
    else {
        assert(num_constr_ == m);
        assert(num_var_ == n);
        std::copy_n(std::begin(x_user), n, std::begin(x_solver));
        std::copy_n(std::begin(slack_user), m, std::begin(x_solver) + n);
        std::copy_n(std::begin(y_user), m, std::begin(y_solver));
        std::copy_n(std::begin(z_user), n, std::begin(z_solver));
        for (Int i = 0; i < m; i++)
            z_solver[n+i] = c(n+i)-y_solver[i];
    }

    // Apply scaling.
    if (colscale_.size() > 0) {
        for (Int j = 0; j < num_cols_; j++) {
            x_solver[j] /= colscale_[j];
            z_solver[j] *= colscale_[j];
        }
    }
    if (rowscale_.size() > 0) {
        for (Int i = 0; i < num_rows_; i++) {
            y_solver[i] /= rowscale_[i];
            x_solver[num_cols_+i] *= rowscale_[i];
            z_solver[num_cols_+i] /= rowscale_[i];
        }
    }
}

void Model::PresolveInteriorPoint(const Vector& x_user,
                                  const Vector& xl_user,
                                  const Vector& xu_user,
                                  const Vector& slack_user,
                                  const Vector& y_user,
                                  const Vector& zl_user,
                                  const Vector& zu_user,
                                  Vector& x_solver,
                                  Vector& xl_solver,
                                  Vector& xu_solver,
                                  Vector& y_solver,
                                  Vector& zl_solver,
                                  Vector& zu_solver) const {
    const Int m = rows();
    const Int n = cols();

    if (dualized_) {
        // Not implemented.
        assert(false);
    }
    else {
        assert(num_constr_ == m);
        assert(num_var_ == n);

        std::copy_n(std::begin(x_user), num_var_, std::begin(x_solver));
        std::copy_n(std::begin(slack_user), num_constr_,
                    std::begin(x_solver) + n);
        std::copy_n(std::begin(xl_user), num_var_, std::begin(xl_solver));
        std::copy_n(std::begin(xu_user), num_var_, std::begin(xu_solver));
        std::copy_n(std::begin(y_user), num_constr_, std::begin(y_solver));
        std::copy_n(std::begin(zl_user), num_var_, std::begin(zl_solver));
        std::copy_n(std::begin(zu_user), num_var_, std::begin(zu_solver));

        for (Int i = 0; i < m; i++) {
            switch (user_model_->constr_type(i)) {
            case '=':
                assert(lb_[n+i] == 0.0 && ub_[n+i] == 0.0);
                // For a fixed slack variable xl, xu, zl and zu won't be used
                // by the IPM. Just put them to zero.
                xl_solver[n+i] = 0.0;
                xu_solver[n+i] = 0.0;
                zl_solver[n+i] = 0.0;
                zu_solver[n+i] = 0.0;
                break;
            case '<':
                assert(lb_[n+i] == 0.0 && ub_[n+i] == INFINITY);
                xl_solver[n+i] = slack_user[i];
                xu_solver[n+i] = INFINITY;
                zl_solver[n+i] = -y_user[i];
                zu_solver[n+i] = 0.0;
                break;
            case '>':
                assert(lb_[n+i] == -INFINITY && ub_[n+i] == 0.0);
                xl_solver[n+i] = INFINITY;
                xu_solver[n+i] = -slack_user[i];
                zl_solver[n+i] = 0.0;
                zu_solver[n+i] = y_user[i];
                break;
            }
        }
    }

    // Apply scaling.
    if (colscale_.size() > 0) {
        for (Int j = 0; j < num_cols_; j++) {
            x_solver[j] /= colscale_[j];
            xl_solver[j] /= colscale_[j];
            xu_solver[j] /= colscale_[j];
            zl_solver[j] *= colscale_[j];
            zu_solver[j] *= colscale_[j];
        }
    }
    if (rowscale_.size() > 0) {
        for (Int i = 0; i < num_rows_; i++) {
            y_solver[i] /= rowscale_[i];
            x_solver[num_cols_+i] *= rowscale_[i];
            xl_solver[num_cols_+i] *= rowscale_[i];
            xu_solver[num_cols_+i] *= rowscale_[i];
            zl_solver[num_cols_+i] /= rowscale_[i];
            zu_solver[num_cols_+i] /= rowscale_[i];
        }
    }
}

void Model::PostsolveGeneralPoint(const Vector& x_solver,
                                  const Vector& y_solver,
                                  const Vector& z_solver,
                                  Vector& x_user,
                                  Vector& slack_user,
                                  Vector& y_user,
                                  Vector& z_user) const {
    const Int m = rows();
    const Int n = cols();

    if (dualized_) {
        assert(num_var_ == m);
        assert(num_constr_ + boxed_vars_.size() == n);
        for (Int j = 0; j < num_var_; j++) {
            x_user[j] = -y_solver[j] * rowscale(j);
            z_user[j] = x_solver[n+j] / rowscale(j);
        }
        for (Int i = 0; i < num_constr_; i++) {
            slack_user[i] = -z_solver[i] / colscale(i);
            y_user[i] = x_solver[i] * colscale(i);
        }
        Int k = num_constr_;
        for (Int j : boxed_vars_) {
            z_user[j] -= x_solver[k] * colscale(k);
            k++;
        }
        assert(k == n);
        for (Int j : negated_vars_) {
            x_user[j] *= -1.0;
            z_user[j] *= -1.0;
        }
    }
    else {
        assert(num_constr_ == m);
        assert(num_var_ == n);
        if (colscale_.size() > 0) {
            x_user = x_solver[std::slice(0, num_var_, 1)] * colscale_;
            z_user = z_solver[std::slice(0, num_var_, 1)] / colscale_;
        } else {
            x_user = x_solver[std::slice(0, num_var_, 1)];
            z_user = z_solver[std::slice(0, num_var_, 1)];
        }
        if (rowscale_.size() > 0) {
            slack_user = x_solver[std::slice(n, num_constr_, 1)] / rowscale_;
            y_user = y_solver * rowscale_;
        } else {
            slack_user = x_solver[std::slice(n, num_constr_, 1)];
            y_user = y_solver;
        }
    }
}

void Model::PostsolveInteriorPoint(const Vector& x_solver,
                                   const Vector& xl_solver,
                                   const Vector& xu_solver,
                                   const Vector& y_solver,
                                   const Vector& zl_solver,
                                   const Vector& zu_solver,
                                   Vector& x_user,
                                   Vector& xl_user,
                                   Vector& xu_user,
                                   Vector& slack_user,
                                   Vector& y_user,
                                   Vector& zl_user,
                                   Vector& zu_user) const {
    const Int m = rows();
    const Int n = cols();

    if (dualized_) {
        assert(num_var_ == m);
        assert(num_constr_ + boxed_vars_.size() == n);
        if (rowscale_.size() > 0)
            x_user = -y_solver * rowscale_;
        else
            x_user = -y_solver;

        // If the solution from the solver would be exact, we could copy the
        // first num_constr_ entries from x_solver into y_user. However, to
        // satisfy the sign condition on y_user even if the solution is not
        // exact, we have to use the xl_solver and xu_solver entries for
        // inequality constraints.
        for (Int i = 0; i < num_constr_; i++) {
            switch (user_model_->constr_type(i)) {
            case '=':
                y_user[i] = x_solver[i] * colscale(i);
                break;
            case '<':
                y_user[i] = -xu_solver[i] * colscale(i);
                break;
            case '>':
                y_user[i] = xl_solver[i] * colscale(i);
                break;
            }
            assert(std::isfinite(y_user[i]));
        }

        // Dual variables associated with lbuser <= x in the scaled user model
        // are the slack variables from the solver. For an exact solution we
        // would have x_solver[n+1:n+m] == xl_solver[n+1:n+m]. Using xl_solver
        // guarantees that zl_user >= 0 in any case. If variable j has no lower
        // bound in the scaled user model (i.e. is free), then the j-th slack
        // variable was fixed at zero in the solver model, but the IPM solution
        // may not satisfy this. Hence we must set zl_user[j] = 0 explicitly.
        if (rowscale_.size() > 0)
            zl_user = xl_solver[std::slice(n, num_var_, 1)] / rowscale_;
        else
            zl_user = xl_solver[std::slice(n, num_var_, 1)];
        for (Int j = 0; j < num_var_; j++)
            if (!std::isfinite(user_model_->lb(j)))
                zl_user[j] = 0.0;

        // Dual variables associated with x <= ubuser in the scaled user model
        // are the primal variables that were added for boxed variables in the
        // solver model.
        zu_user = 0.0;
        Int k = num_constr_;
        for (Int j : boxed_vars_) {
            zu_user[j] = xl_solver[k] * colscale(k);
            k++;
        }
        assert(k == n);

        // xl in the scaled user model is zl[n+1:n+m] in the solver model or
        // infinity.
        for (Int i = 0; i < m; i++) {
            if (std::isfinite(user_model_->lb(i)))
                xl_user[i] = zl_solver[n+i] * rowscale(i);
            else
                xl_user[i] = INFINITY;
        }

        // xu in the scaled user model are the entries in zl for columns of the
        // negative identity matrix (that were added for boxed variables).
        xu_user = INFINITY;
        k = num_constr_;
        for (Int j : boxed_vars_) {
            xu_user[j] = zl_solver[k] / colscale(k);
            k++;
        }
        assert(k == n);

        for (Int i = 0; i < num_constr_; i++) {
            switch (user_model_->constr_type(i)) {
            case '=':
                slack_user[i] = 0.0;
                break;
            case '<':
                slack_user[i] = zu_solver[i] / colscale(i);
                break;
            case '>':
                slack_user[i] = -zl_solver[i] / colscale(i);
                break;
            }
        }

        for (Int j : negated_vars_) {
            assert(std::isfinite(xl_user[j]));
            assert(std::isinf(xu_user[j]));
            assert(zu_user[j] == 0.0);
            x_user[j] *= -1.0;
            xu_user[j] = xl_user[j];
            xl_user[j] = INFINITY;
            zu_user[j] = zl_user[j];
            zl_user[j] = 0.0;
        }
    }
    else {
        assert(num_constr_ == m);
        assert(num_var_ == n);
        if (colscale_.size() > 0)
            x_user = x_solver[std::slice(0, num_var_, 1)] * colscale_;
        else
            x_user = x_solver[std::slice(0, num_var_, 1)];

        // Instead of copying y_solver into y_user, we use the entries from
        // zl_solver and zu_solver for inequality constraints, so that the sign
        // condition on y_user is satisfied.
        for (Int i = 0; i < m; i++) {
            assert(lb_[n+i] == 0.0 || lb_[n+i] == -INFINITY);
            assert(ub_[n+i] == 0.0 || ub_[n+i] ==  INFINITY);
            assert(lb_[n+i] == 0.0 || ub_[n+i] == 0.0);
            switch (user_model_->constr_type(i)) {
            case '=':
                y_user[i] = y_solver[i] * rowscale(i);
                break;
            case '<':
                y_user[i] = -zl_solver[n+i] * rowscale(i);
                break;
            case '>':
                y_user[i] = zu_solver[n+i] * rowscale(i);
                break;
            }
            assert(std::isfinite(y_user[i]));
        }
        if (colscale_.size() > 0) {
            zl_user = zl_solver[std::slice(0, num_var_, 1)] / colscale_;
            zu_user = zu_solver[std::slice(0, num_var_, 1)] / colscale_;
            xl_user = xl_solver[std::slice(0, num_var_, 1)] * colscale_;
            xu_user = xu_solver[std::slice(0, num_var_, 1)] * colscale_;
        }
        else {
            zl_user = zl_solver[std::slice(0, num_var_, 1)];
            zu_user = zu_solver[std::slice(0, num_var_, 1)];
            xl_user = xl_solver[std::slice(0, num_var_, 1)];
            xu_user = xu_solver[std::slice(0, num_var_, 1)];
        }

        // If the solution would be exact, slack_user were given by the entries
        // of x_solver corresponding to slack columns. To satisfy the sign
        // condition in any case, we build the slack for inequality constraints
        // from xl_solver and xu_solver and set the slack for equality
        // constraints to zero.
        for (Int i = 0; i < m; i++) {
            switch (user_model_->constr_type(i)) {
            case '=':
                slack_user[i] = 0.0;
                break;
            case '<':
                slack_user[i] = xl_solver[n+i] / rowscale(i);
                break;
            case '>':
                slack_user[i] = -xu_solver[n+i] / rowscale(i);
                break;
            }
            assert(std::isfinite(slack_user[i]));
        }
    }
}

void Model::PostsolveBasis(const std::vector<Int>& basic_status_solver,
                           std::vector<Int>& cbasis_user,
                           std::vector<Int>& vbasis_user) const {
    const Int m = rows();
    const Int n = cols();

    if (dualized_) {
        assert(num_var_ == m);
        assert(num_constr_ + boxed_vars_.size() == n);
        for (Int i = 0; i < num_constr_; i++) {
            if (basic_status_solver[i] == IPX_basic)
                cbasis_user[i] = IPX_nonbasic;
            else
                cbasis_user[i] = IPX_basic;
        }
        for (Int j = 0; j < num_var_; j++) {
            // slack cannot be superbasic
            assert(basic_status_solver[n+j] != IPX_superbasic);
            if (basic_status_solver[n+j] == 0)
                vbasis_user[j] = std::isfinite(user_model_->lb(j)) ?
                    IPX_nonbasic_lb : IPX_superbasic;
            else
                vbasis_user[j] = IPX_basic;
        }
        Int k = num_constr_;
        for (Int j : boxed_vars_)
            if (basic_status_solver[k++] == IPX_basic) {
                assert(vbasis_user[j] == IPX_basic);
                vbasis_user[j] = IPX_nonbasic_ub;
            }
        for (Int j : negated_vars_) {
            assert(vbasis_user[j] != IPX_nonbasic_ub);
            if (vbasis_user[j] == IPX_nonbasic_lb)
                vbasis_user[j] = IPX_nonbasic_ub;
        }
    }
    else {
        assert(num_constr_ == m);
        assert(num_var_ == n);
        for (Int i = 0; i < num_constr_; i++) {
            // slack cannot be superbasic
            assert(basic_status_solver[n+i] != IPX_superbasic);
            if (basic_status_solver[n+i] == IPX_basic)
                cbasis_user[i] = IPX_basic;
            else
                cbasis_user[i] = IPX_nonbasic;
        }
        for (Int j = 0; j < num_var_; j++)
            vbasis_user[j] = basic_status_solver[j];
    }
}

void Model::CorrectBasicSolution(Vector& x, Vector& slack, Vector& y, Vector& z,
                                 const std::vector<Int> cbasis,
                                 const std::vector<Int> vbasis) const {
    for (Int j = 0; j < num_var_; j++) {
        if (vbasis[j] == IPX_nonbasic_lb)
            x[j] = user_model_->lb(j);
        if (vbasis[j] == IPX_nonbasic_ub)
            x[j] = user_model_->ub(j);
        if (vbasis[j] == IPX_basic)
            z[j] = 0.0;
    }
    for (Int i = 0; i < num_constr_; i++) {
        if (cbasis[i] == IPX_nonbasic)
            slack[i] = 0.0;
        if (cbasis[i] == IPX_basic)
            y[i] = 0.0;
    }
}

void Model::ComputeUserModelAttributes() {
    num_constr_ = user_model_->num_constr();
    num_eqconstr_ = std::count(user_model_->constr_type().begin(),
                               user_model_->constr_type().end(), '=');
    num_var_ = user_model_->num_var();
    num_free_var_ = 0;
    boxed_vars_.clear();
    for (Int j = 0; j < num_var_; j++) {
        bool has_lb = std::isfinite(user_model_->lb(j));
        bool has_ub = std::isfinite(user_model_->ub(j));
        if (!has_lb && !has_ub)
            num_free_var_++;
        else if (has_lb && has_ub)
            boxed_vars_.push_back(j);
    }
}

void Model::ScaleModel(const Control& control) {
    colscale_.resize(0);
    rowscale_.resize(0);

    // Choose scaling method.
    if (control.scale() >= 1)
        EquilibrateMatrix();

    // Apply scaling to vectors.
    if (colscale_.size() > 0) {
        assert(colscale_.size() == num_cols_);
        for (Int j = 0; j < num_cols_; j++) {
            c_[j] *= colscale_[j];
            lb_[j] /= colscale_[j];
            ub_[j] /= colscale_[j];
        }
    }
    if (rowscale_.size() > 0) {
        assert(rowscale_.size() == num_rows_);
        for (Int i = 0; i < num_rows_; i++) {
            b_[i] *= rowscale_[i];
            c_[num_cols_+i] /= rowscale_[i];
            lb_[num_cols_+i] *= rowscale_[i];
            ub_[num_cols_+i] *= rowscale_[i];
        }
    }
}

void Model::LoadPrimal() {
    const SparseMatrix& A = user_model_->A();
    const Vector& obj = user_model_->obj();
    const std::vector<char>& constr_type = user_model_->constr_type();
    const Vector& rhs = user_model_->rhs();
    const Vector& lbuser = user_model_->lb();
    const Vector& ubuser = user_model_->ub();

    num_rows_ = num_constr_;
    num_cols_ = num_var_;
    dualized_ = false;

    // Copy A and append identity matrix.
    AI_ = A;
    for (Int i = 0; i < num_constr_; i++) {
        AI_.push_back(i, 1.0);
        AI_.add_column();
    }
    assert(AI_.cols() == num_var_+num_constr_);

    // Copy vectors and set bounds on slack variables.
    b_ = rhs;
    c_.resize(num_var_+num_constr_);
    c_ = 0.0;
    std::copy_n(std::begin(obj), num_var_, std::begin(c_));
    lb_.resize(num_rows_+num_cols_);
    std::copy_n(std::begin(lbuser), num_var_, std::begin(lb_));
    ub_.resize(num_rows_+num_cols_);
    std::copy_n(std::begin(ubuser), num_var_, std::begin(ub_));
    for (Int i = 0; i < num_constr_; i++) {
        switch(constr_type[i]) {
        case '=':
            lb_[num_var_+i] = 0.0;
            ub_[num_var_+i] = 0.0;
            break;
        case '<':
            lb_[num_var_+i] = 0.0;
            ub_[num_var_+i] = INFINITY;
            break;
        case '>':
            lb_[num_var_+i] = -INFINITY;
            ub_[num_var_+i] = 0.0;
            break;
        }
    }
}

void Model::LoadDual() {
    const SparseMatrix& A = user_model_->A();
    const Vector& obj = user_model_->obj();
    const std::vector<char>& constr_type = user_model_->constr_type();
    const Vector& rhs = user_model_->rhs();
    const Vector& lbuser = user_model_->lb();
    const Vector& ubuser = user_model_->ub();

    num_rows_ = num_var_;
    num_cols_ = num_constr_ + boxed_vars_.size();
    dualized_ = true;

    // Implicitly negate variables with infinite lbuser_ but finite ubuser_.
    std::vector<bool> negated(num_var_);
    for (Int j = 0; j < num_var_; j++) {
        if (std::isinf(lbuser[j]) && std::isfinite(ubuser[j])) {
            negated[j] = true;
            negated_vars_.push_back(j);
        }
    }

    // Build AI.
    AI_ = Transpose(A);
    if (!negated_vars_.empty()) {
        for (Int pos = 0; pos < AI_.entries(); pos++) {
            if (negated[AI_.index(pos)])
                AI_.value(pos) *= -1.0;
        }
    }
    for (Int j : boxed_vars_) {
        AI_.push_back(j, -1.0);
        AI_.add_column();
    }
    assert(AI_.cols() == num_cols_);
    for (Int i = 0; i < num_rows_; i++) {
        AI_.push_back(i, 1.0);
        AI_.add_column();
    }

    // Build vectors.
    b_ = obj;
    for (Int j : negated_vars_)
        b_[j] *= -1.0;
    c_.resize(num_cols_+num_rows_);
    Int put = 0;
    for (double x : rhs)
        c_[put++] = -x;
    for (Int j : boxed_vars_)
        c_[put++] = ubuser[j];
    assert(put == num_cols_);
    for (Int j = 0; j < num_var_; j++) {
        double lb = negated[j] ? -ubuser[j] : lbuser[j];
        // If lb is -infinity, then the variable will be fixed and we can give
        // it any (finite) cost.
        c_[put++] = std::isfinite(lb) ? -lb : 0.0;
    }
    lb_.resize(num_cols_+num_rows_);
    ub_.resize(num_cols_+num_rows_);
    for (Int i = 0; i < num_constr_; i++)
        switch(constr_type[i]) {
        case '=':
            lb_[i] = -INFINITY;
            ub_[i] = INFINITY;
            break;
        case '<':
            lb_[i] = -INFINITY;
            ub_[i] = 0.0;
            break;
        case '>':
            lb_[i] = 0.0;
            ub_[i] = INFINITY;
            break;
        }
    for (Int j = num_constr_; j < num_cols_; j++) {
        lb_[j] = 0.0;
        ub_[j] = INFINITY;
    }
    for (Int j = 0; j < num_var_; j++) {
        double lb = negated[j] ? -ubuser[j] : lbuser[j];
        lb_[num_cols_+j] = 0.0;
        ub_[num_cols_+j] = std::isfinite(lb) ? INFINITY : 0.0;
    }
}

// Returns a power-of-2 factor s such that s*2^exp becomes closer to the
// interval [2^expmin, 2^expmax].
static double EquilibrationFactor(int expmin, int expmax, int exp) {
    if (exp < expmin)
        return std::ldexp(1.0, (expmin-exp+1)/2);
    if (exp > expmax)
        return std::ldexp(1.0, -((exp-expmax+1)/2));
    return 1.0;
}

void Model::EquilibrateMatrix() {
    const Int m = AI_.rows();
    const Int n = AI_.cols() - m;
    const Int* Ap = AI_.colptr();
    const Int* Ai = AI_.rowidx();
    double* Ax = AI_.values();

    colscale_.resize(0);
    rowscale_.resize(0);

    // The absolute value of each nonzero entry of AI can be written as x*2^exp
    // with x in the range [0.5,1). We consider AI well scaled if each entry is
    // such that expmin <= exp <= expmax for parameters expmin and expmax.
    // For example,
    //
    //  expmin  expmax  allowed range
    //  -----------------------------
    //       0       2    [0.5,    4)
    //      -1       3    [0.25,   8)
    //       1       7    [1.0,  128)
    //
    // If AI is not well scaled, a recursive row and column equilibration is
    // applied. In each iteration, the scaling factors are roughly 1/sqrt(max),
    // where max is the maximum row or column entry. However, the factors are
    // truncated to powers of 2, so that no round-off errors occur.

    constexpr int expmin = 0;
    constexpr int expmax = 3;
    constexpr Int maxround = 10;

    // Quick return if entries are within the target range.
    bool out_of_range = false;
    for (Int p = 0; p < Ap[n]; p++) {
        int exp;
        std::frexp(std::abs(Ax[p]), &exp);
        if (exp < expmin || exp > expmax) {
            out_of_range = true;
            break;
        }
    }
    if (!out_of_range)
        return;

    colscale_.resize(n);
    rowscale_.resize(m);
    colscale_ = 1.0;
    rowscale_ = 1.0;
    Vector colmax(n), rowmax(m);

    for (Int round = 0; round < maxround; round++) {
        // Compute infinity norm of each row and column.
        rowmax = 0.0;
        for (Int j = 0; j < n; j++) {
            colmax[j] = 0.0;
            for (Int p = Ap[j]; p < Ap[j+1]; p++) {
                Int i = Ai[p];
                double xa = std::abs(Ax[p]);
                colmax[j] = std::max(colmax[j], xa);
                rowmax[i] = std::max(rowmax[i], xa);
            }
        }
        // Replace rowmax and colmax entries by scaling factors from this round.
        bool out_of_range = false;
        for (Int i = 0; i < m; i++) {
            int exp;
            std::frexp(rowmax[i], &exp);
            rowmax[i] = EquilibrationFactor(expmin, expmax, exp);
            if (rowmax[i] != 1.0) {
                out_of_range = true;
                rowscale_[i] *= rowmax[i];
            }
        }
        for (Int j = 0; j < n; j++) {
            int exp;
            std::frexp(colmax[j], &exp);
            colmax[j] = EquilibrationFactor(expmin, expmax, exp);
            if (colmax[j] != 1.0) {
                out_of_range = true;
                colscale_[j] *= colmax[j];
            }
        }
        if (!out_of_range)
            break;
        // Rescale A.
        for (Int j = 0; j < n; j++) {
            for (Int p = Ap[j]; p < Ap[j+1]; p++) {
                Ax[p] *= colmax[j];     // column scaling
                Ax[p] *= rowmax[Ai[p]]; // row scaling
            }
        }
    }
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

void Model::PrintCoefficientRange(const Control& control) const {
    const SparseMatrix& A = user_model_->A();
    const Vector& obj = user_model_->obj();
    const Vector& rhs = user_model_->rhs();
    const Vector& lbuser = user_model_->lb();
    const Vector& ubuser = user_model_->ub();

    double amin = INFINITY;
    double amax = 0.0;
    for (Int j = 0; j < A.cols(); j++) {
        for (Int p = A.begin(j); p < A.end(j); p++) {
            double x = A.value(p);
            if (x != 0.0) {
                amin = std::min(amin, std::abs(x));
                amax = std::max(amax, std::abs(x));
            }
        }
    }
    if (amin == INFINITY)       // no nonzero entries in A_
        amin = 0.0;
    control.Log()
        << Textline("Matrix range:")
        << "[" << Scientific(amin, 5, 0) << ", "
        << Scientific(amax, 5, 0) << "]\n";

    double rhsmin = INFINITY;
    double rhsmax = 0.0;
    for (double x : rhs) {
        if (x != 0.0) {
            rhsmin = std::min(rhsmin, std::abs(x));
            rhsmax = std::max(rhsmax, std::abs(x));
        }
    }
    if (rhsmin == INFINITY)     // no nonzero entries in rhs
        rhsmin = 0.0;
    control.Log()
        << Textline("RHS range:")
        << "[" << Scientific(rhsmin, 5, 0) << ", "
        << Scientific(rhsmax, 5, 0) << "]\n";

    double objmin = INFINITY;
    double objmax = 0.0;
    for (double x : obj) {
        if (x != 0.0) {
            objmin = std::min(objmin, std::abs(x));
            objmax = std::max(objmax, std::abs(x));
        }
    }
    if (objmin == INFINITY)     // no nonzero entries in obj
        objmin = 0.0;
    control.Log()
        << Textline("Objective range:")
        << "[" << Scientific(objmin, 5, 0) << ", "
        << Scientific(objmax, 5, 0) << "]\n";

    double boundmin = INFINITY;
    double boundmax = 0.0;
    for (double x : lbuser) {
        if (x != 0.0 && std::isfinite(x)) {
            boundmin = std::min(boundmin, std::abs(x));
            boundmax = std::max(boundmax, std::abs(x));
        }
    }
    for (double x : ubuser) {
        if (x != 0.0 && std::isfinite(x)) {
            boundmin = std::min(boundmin, std::abs(x));
            boundmax = std::max(boundmax, std::abs(x));
        }
    }
    if (boundmin == INFINITY)   // no finite nonzeros entries in bounds
        boundmin = 0.0;
    control.Log()
        << Textline("Bounds range:")
        << "[" << Scientific(boundmin, 5, 0) << ", "
        << Scientific(boundmax, 5, 0) << "]\n";
}

void Model::PrintPreprocessingLog(const Control& control) const {
    // Find the minimum and maximum scaling factor.
    double minscale = INFINITY;
    double maxscale = 0.0;
    if (colscale_.size() > 0) {
        auto minmax = std::minmax_element(std::begin(colscale_),
                                          std::end(colscale_));
        minscale = std::min(minscale, *minmax.first);
        maxscale = std::max(maxscale, *minmax.second);
    }
    if (rowscale_.size() > 0) {
        auto minmax = std::minmax_element(std::begin(rowscale_),
                                          std::end(rowscale_));
        minscale = std::min(minscale, *minmax.first);
        maxscale = std::max(maxscale, *minmax.second);
    }
    if (minscale == INFINITY)
        minscale = 1.0;
    if (maxscale == 0.0)
        maxscale = 1.0;

    control.Log()
        << "Preprocessing\n"
        << Textline("Dualized model:") << (dualized() ? "yes" : "no") << '\n'
        << Textline("Number of dense columns:") << num_dense_cols() << '\n';
    if (control.scale() > 0) {
        control.Log()
            << Textline("Range of scaling factors:") << "["
            << Scientific(minscale, 8, 2) << ", "
            << Scientific(maxscale, 8, 2) << "]\n";
    }
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
