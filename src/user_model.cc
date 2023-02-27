// Copyright (c) 2023 ERGO-Code. See license.txt for license.

#include "user_model.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include "utils.h"

namespace ipx {

Int UserModel::Load(const Control& control, Int num_constr, Int num_var,
                    const Int* Ap, const Int* Ai, const double* Ax,
                    const double* rhs, const char* constr_type,
                    const double* obj, const double* lbuser,
                    const double* ubuser) {
    (void)control;
    clear();
    Int errflag = CopyInput(num_constr, num_var, Ap, Ai, Ax, rhs, constr_type,
                            obj, lbuser, ubuser);
    if (errflag)
        return errflag;
    ComputeNorms();
    empty_ = false;
    return 0;
}

void UserModel::GetInfo(Info *info) const {
    info->num_var = num_var_;
    info->num_constr = num_constr_;
    info->num_entries = A_.entries();
}

void UserModel::clear() {
    empty_ = true;
    num_var_ = 0;
    num_constr_ = 0;
    obj_.resize(0);
    constr_type_.clear();
    rhs_.resize(0);
    lb_.resize(0);
    ub_.resize(0);
    A_.clear();
    norm_obj_ = 0.0;
    norm_rhs_ = 0.0;
    norm_bounds_ = 0.0;
}

Int UserModel::CheckInteriorPoint(const double* x, const double* xl,
                                  const double* xu, const double* slack,
                                  const double* y, const double* zl,
                                  const double* zu) const {
    if (!x || !xl || !xu || !slack || !y || !zl || !zu)
        return IPX_ERROR_argument_null;

    // Check that the point is compatible with bounds.
    for (Int j = 0; j < num_var_; j++) {
        if (!std::isfinite(x[j]))
            return IPX_ERROR_invalid_vector;
    }
    for (Int j = 0; j < num_var_; j++) {
        if (!(xl[j] >= 0.0) ||
            (lb_[j] == -INFINITY && xl[j] != INFINITY) ||
            (lb_[j] != -INFINITY && xl[j] == INFINITY))
            return IPX_ERROR_invalid_vector;
    }
    for (Int j = 0; j < num_var_; j++) {
        if (!(xu[j] >= 0.0) ||
            (ub_[j] == INFINITY && xu[j] != INFINITY) ||
            (ub_[j] != INFINITY && xu[j] == INFINITY))
            return IPX_ERROR_invalid_vector;
    }
    for (Int i = 0; i < num_constr_; i++) {
        if (!std::isfinite(slack[i]) ||
            (constr_type_[i] == '=' && !(slack[i] == 0.0)) ||
            (constr_type_[i] == '<' && !(slack[i] >= 0.0)) ||
            (constr_type_[i] == '>' && !(slack[i] <= 0.0)))
            return IPX_ERROR_invalid_vector;
    }
    for (Int i = 0; i < num_constr_; i++) {
        if (!std::isfinite(y[i]) ||
            (constr_type_[i] == '<' && !(y[i] <= 0.0)) ||
            (constr_type_[i] == '>' && !(y[i] >= 0.0)))
            return IPX_ERROR_invalid_vector;
    }
    for (Int j = 0; j < num_var_; j++) {
        if (!(zl[j] >= 0.0 && zl[j] < INFINITY) ||
            (lb_[j] == -INFINITY && zl[j] != 0.0))
            return IPX_ERROR_invalid_vector;
    }
    for (Int j = 0; j < num_var_; j++) {
        if (!(zu[j] >= 0.0 && zu[j] < INFINITY) ||
            (ub_[j] == INFINITY && zu[j] != 0.0))
            return IPX_ERROR_invalid_vector;
    }

    return 0;
}

void UserModel::EvaluateInteriorPoint(const InteriorSolution& point,
                                      Info* info) const {
    const Vector& x = point.x;
    const Vector& xl = point.xl;
    const Vector& xu = point.xu;
    const Vector& slack = point.slack;
    const Vector& y = point.y;
    const Vector& zl = point.zl;
    const Vector& zu = point.zu;

    // rb = rhs-slack-A*x
    // Add rhs at the end to avoid loosing digits when x is huge.
    Vector rb(num_constr_);
    MultiplyAdd(A_, x, -1.0, rb, 'N');
    rb -= slack;
    rb += rhs_;

    // rc = obj-zl+zu-A'y
    // Add obj at the end to avoid loosing digits when y, z are huge.
    Vector rc(num_var_);
    MultiplyAdd(A_, y, -1.0, rc, 'T');
    rc -= zl - zu;
    rc += obj_;

    double presidual = Infnorm(rb);
    double dresidual = Infnorm(rc);

    // Update primal residual with bound residuals.
    for (Int j = 0; j < num_var_; j++) {
        if (std::isfinite(lb_[j]))
            presidual = std::max(presidual, std::abs(lb_[j] - x[j] + xl[j]));
        else
            assert(xl[j] == INFINITY);
        if (std::isfinite(ub_[j]))
            presidual = std::max(presidual, std::abs(ub_[j] - x[j] - xu[j]));
        else
            assert(xu[j] == INFINITY);
    }

    double pobjective = Dot(obj_, x);
    double dobjective = Dot(rhs_, y);
    for (Int j = 0; j < num_var_; j++) {
        if (std::isfinite(lb_[j]))
            dobjective += lb_[j] * zl[j];
        if (std::isfinite(ub_[j]))
            dobjective -= ub_[j] * zu[j];
    }
    assert(std::isfinite(dobjective));
    double objective_gap = (pobjective-dobjective) /
        (1.0 + 0.5*std::abs(pobjective+dobjective));

    double complementarity = 0.0;
    for (Int j = 0; j < num_var_; j++) {
        if (std::isfinite(lb_[j]))
            complementarity += xl[j]*zl[j];
        if (std::isfinite(ub_[j]))
            complementarity += xu[j]*zu[j];
    }
    for (Int i = 0; i < num_constr_; i++)
        complementarity -= y[i]*slack[i];

    info->abs_presidual = presidual;
    info->abs_dresidual = dresidual;
    info->rel_presidual = presidual/(1.0+std::max(norm_rhs_, norm_bounds_));
    info->rel_dresidual = dresidual/(1.0+norm_obj_);
    info->pobjval = pobjective;
    info->dobjval = dobjective;
    info->rel_objgap = objective_gap;
    info->complementarity = complementarity;
    info->normx = Infnorm(x);
    info->normy = Infnorm(y);
    info->normz = std::max(Infnorm(zl), Infnorm(zu));
}

void UserModel::EvaluateBasicPoint(const BasicSolution& point,
                                   Info* info) const {
    const Vector& x = point.x;
    const Vector& slack = point.slack;
    const Vector& y = point.y;
    const Vector& z = point.z;
    const std::vector<Int>& vbasis = point.vbasis;
    double primal_infeas = 0.0;
    double dual_infeas = 0.0;

    for (Int j = 0; j < num_var_; j++) {
        primal_infeas = std::max(primal_infeas, lb_[j] - x[j]);
        primal_infeas = std::max(primal_infeas, x[j] - ub_[j]);
        if (vbasis[j] != IPX_nonbasic_lb)
            dual_infeas = std::max(dual_infeas, z[j]);
        if (vbasis[j] != IPX_nonbasic_ub)
            dual_infeas = std::max(dual_infeas, -z[j]);
    }

    for (Int i = 0; i < num_constr_; i++) {
        switch (constr_type_[i]) {
        case '<':
            primal_infeas = std::max(primal_infeas, -slack[i]);
            dual_infeas = std::max(dual_infeas, y[i]);
            break;
        case '>':
            primal_infeas = std::max(primal_infeas, slack[i]);
            dual_infeas = std::max(dual_infeas, -y[i]);
            break;
        case '=':
            primal_infeas = std::max(primal_infeas, std::abs(slack[i]));
            break;
        }
    }

    info->primal_infeas = primal_infeas;
    info->dual_infeas = dual_infeas;
    info->objval = Dot(obj_, x);
}

namespace {

// Checks if the vectors are valid LP data vectors. Returns 0 if OK and a
// negative value if a vector is invalid.
int CheckVectors(Int m, Int n, const double* rhs,const char* constr_type,
                 const double* obj, const double* lb, const double* ub) {
    for (Int i = 0; i < m; i++)
        if (!std::isfinite(rhs[i]))
            return -1;
    for (Int j = 0; j < n; j++)
        if (!std::isfinite(obj[j]))
            return -2;
    for (Int j = 0; j < n; j++) {
        if (!std::isfinite(lb[j]) && lb[j] != -INFINITY)
            return -3;
        if (!std::isfinite(ub[j]) && ub[j] != INFINITY)
            return -3;
        if (lb[j] > ub[j])
            return -3;
    }
    for (Int i = 0; i < m; i++)
        if (constr_type[i] != '=' && constr_type[i] != '<' &&
            constr_type[i] != '>')
            return -4;
    return 0;
}

// Checks if A is a valid m-by-n matrix in CSC format. Returns 0 if OK and a
// negative value otherwise.
Int CheckMatrix(Int m, Int n, const Int *Ap, const Int *Ai, const double *Ax) {
    if (Ap[0] != 0)
        return -5;
    for (Int j = 0; j < n; j++)
        if (Ap[j] > Ap[j+1])
            return -5;
    for (Int p = 0; p < Ap[n]; p++)
        if (!std::isfinite(Ax[p]))
            return -6;
    // Test for out of bound indices and duplicates.
    std::vector<Int> marked(m, -1);
    for (Int j = 0; j < n; j++) {
        for (Int p = Ap[j]; p < Ap[j+1]; p++) {
            Int i = Ai[p];
            if (i < 0 || i >= m)
                return -7;
            if (marked[i] == j)
                return -8;
            marked[i] = j;
        }
    }
    return 0;
}

} // namespace

Int UserModel::CopyInput(Int num_constr, Int num_var, const Int* Ap,
                         const Int* Ai, const double* Ax, const double* rhs,
                         const char* constr_type, const double* obj,
                         const double* lbuser, const double* ubuser) {
    if (num_constr < 0 || num_var <= 0) {
        return IPX_ERROR_invalid_dimension;
    }
    if (!Ap) {
        return IPX_ERROR_argument_null;
    }
    if (num_var > 0 && !(obj && lbuser && ubuser)) {
        return IPX_ERROR_argument_null;
    }
    if (num_constr > 0 && !(rhs && constr_type)) {
        return IPX_ERROR_argument_null;
    }
    Int num_entries = Ap[num_var];
    if (num_entries > 0 && !(Ai && Ax)) {
        return IPX_ERROR_argument_null;
    }
    if (CheckVectors(num_constr, num_var, rhs, constr_type, obj, lbuser, ubuser)
        != 0) {
        return IPX_ERROR_invalid_vector;
    }
    if (CheckMatrix(num_constr, num_var, Ap, Ai, Ax) != 0) {
        return IPX_ERROR_invalid_matrix;
    }
    num_constr_ = num_constr;
    num_var_ = num_var;
    constr_type_ = std::vector<char>(constr_type, constr_type+num_constr);
    obj_ = Vector(obj, num_var);
    rhs_ = Vector(rhs, num_constr);
    lb_ = Vector(lbuser, num_var);
    ub_ = Vector(ubuser, num_var);
    A_.LoadFromArrays(num_constr, num_var, Ap, Ap+1, Ai, Ax);
    return 0;
}

void UserModel::ComputeNorms() {
    norm_obj_ = Infnorm(obj_);
    norm_rhs_ = Infnorm(rhs_);
    norm_bounds_ = 0.0;
    for (double x : lb_)
        if (std::isfinite(x))
            norm_bounds_ = std::max(norm_bounds_, std::abs(x));
    for (double x : ub_)
        if (std::isfinite(x))
            norm_bounds_ = std::max(norm_bounds_, std::abs(x));
}

}  // namespace ipx
