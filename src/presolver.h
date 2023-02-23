// Copyright (c) 2018-2023 ERGO-Code. See license.txt for license.

#ifndef IPX_PRESOLVER_H_
#define IPX_PRESOLVER_H_

#include <vector>
#include "control.h"
#include "iterate.h"
#include "model.h"
#include "simplex_iterate.h"
#include "solution.h"
#include "sparse_matrix.h"
#include "user_model.h"

namespace ipx {

// Presolver converts between user model and solver model.
//
// Implemented operations:
//
// - dualization
// - scaling

class Presolver {
public:
    // Constructor stores references to the user model and solver model.
    Presolver(const UserModel& user_model, Model& model);

    // Discards presolve information and deallocates memory.
    void clear();

    // Builds new solver model by presolving the user model.
    Int PresolveModel(const Control& control);

    // Writes statistics of presolve to @info.
    void GetInfo(Info* info) const;

    // Transforms point from user model to solver model.
    void PresolveStartingPoint(const BasicSolution& user_point,
                               Vector& x_solver, Vector& y_solver,
                               Vector& z_solver) const;

    // Transforms interior point from user model to solver model. Currently only
    // implemented for the case that the model was not dualized in presolve.
    // Returns:
    //  0
    //  IPX_ERROR_not_implemented if the model was dualized in presolve
    Int PresolveIPMStartingPoint(const InteriorSolution& user_point,
                                 Vector& x_solver,
                                 Vector& xl_solver,
                                 Vector& xu_solver,
                                 Vector& y_solver,
                                 Vector& zl_solver,
                                 Vector& zu_solver) const;

    // Given an IPM iterate, recovers the solution to the user model (see the
    // reference documentation).
    void PostsolveInteriorSolution(const Iterate& iterate,
                                   InteriorSolution& user_point) const;

    // Given a basic solution to the solver model, recovers the basic solution
    // to the user model.
    void PostsolveBasicSolution(const SimplexIterate& iterate,
                                const std::vector<Int>& basic_status_solver,
                                BasicSolution& user_point) const;

    // Given a basic status for each variable in the solver model, recovers the
    // basic statuses for constraints and variables in the user model. Each
    // of the pointer arguments can be NULL.
    void PostsolveBasis(const std::vector<Int>& basic_status_solver,
                        Int* cbasis, Int* vbasis) const;

private:
    // Computes quantities associated with user_model_.
    void ComputeUserModelAttributes();

    // Builds solver model without dualization. In Julia notation:
    // num_rows = nc
    // num_cols = nv
    // AI       = [A eye(nc)]
    // b        = rhs
    // c        = [obj    ; zeros(nc)                      ]
    // lb       = [lbuser ; constr_type_ .== '>' ? -Inf : 0]
    // ub       = [ubuser ; constr_type_ .== '<' ? +Inf : 0]
    // dualized = false
    // Here nc = num_constr and nv = num_var.
    void LoadPrimal();

    // Builds solver model with dualization. In Julia notation:
    // num_rows = nv
    // num_cols = nc + nb
    // AI       = [A' -eye(nv)[:,jboxed] eye(nv)]
    // b        = obj
    // c        = [-rhs                          ; ubuser[jb]  ; -lbuser     ]
    // lb       = [constr_type .== '>' ? 0 : -Inf; zeros(nb)   ; zeros(nv)   ]
    // ub       = [constr_type .== '<' ? 0 : +Inf; Inf*ones(nb); Inf*ones(nv)]
    // dualized = true
    // Here nc = num_constr, nv = num_var, nb is the number of boxed variables
    // and jboxed are their indices. Variables with infinite lbuser but finite
    // ubuser are implicitly scaled by -1. Their indices are stored in
    // flipped_vars_. If a variable j of the user model is a free variable, then
    // the j-th slack variable of the solver model gets a zero upper bound (i.e.
    // it is fixed at zero) and its objective coefficient is set to zero.
    void LoadDual();

    // Scales model_ according to parameter control.scale(). The scaling factors
    // are stored in colscale_ and rowscale_. If all factors are 1.0 (either
    // because scaling was turned off or because the algorithm did nothing),
    // rowscale_ and colscale_ have size 0.
    void ScaleModel(const Control& control);

    // Prints user model attributes to control.Log().
    void PrintUserModelAttributes(const Control& control) const;

    // Prints presolve summary to control.Log().
    void PrintPresolveLog(const Control& control) const;

    // Translates arbitrary primal-dual point from user model to solver model.
    // No sign conditions are assumed for the user point.
    void PresolveGeneralPoint(const Vector& x_user,
                              const Vector& slack_user,
                              const Vector& y_user,
                              const Vector& z_user,
                              Vector& x_solver,
                              Vector& y_solver,
                              Vector& z_solver) const;

    // Translates interior point from user model to solver model. The user point
    // must satisfy the sign conditions imposed by the user model. Currently
    // only implemented for dualized_ == false.
    void PresolveInteriorPoint(const Vector& x_user,
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
                               Vector& zu_solver) const;

    // Translates arbitrary primal-dual point from solver model to user model.
    // No sign conditions are assumed.
    void PostsolveGeneralPoint(const Vector& x_solver,
                               const Vector& y_solver,
                               const Vector& z_solver,
                               Vector& x_user,
                               Vector& slack_user,
                               Vector& y_user,
                               Vector& z_user) const;

    // Translates interior point from solver model to user model. The solver
    // point must satisfy the sign conditions imposed by the solver model.
    void PostsolveInteriorPoint(const Vector& x_solver,
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
                                Vector& zu_user) const;

    // Translates basic statuses from solver model to user model.
    void PostsolveBasis(const std::vector<Int>& basic_status_solver,
                        std::vector<Int>& cbasis_user,
                        std::vector<Int>& vbasis_user) const;

    // Adjusts user model primal-dual point to be consistent with basis:
    // - For nonbasic variables sets the primal variable to its bound.
    // - For basic variables sets the dual variable to zero.
    void CorrectBasicSolution(Vector& x, Vector& slack, Vector& y, Vector& z,
                              const std::vector<Int> cbasis,
                              const std::vector<Int> vbasis) const;

    // Recursively equilibrates model_.AI_ in infinity norm using the algorithm
    // from [1]. The scaling factors are truncated to powers of 2. Terminates
    // when the entries of AI_ are within the range [0.5,8). Preserves the
    // rightmost identity matrix in AI_.
    // [1] P. A. Knight, D. Ruiz, B. Ucar, "A symmetry preserving algorithm for
    //     matrix scaling", SIAM J. Matrix Anal., 35(3), 2014.
    void EquilibrateMatrix();

    double colscale(Int j) const {
        return colscale_.size() > 0 ? colscale_[j] : 1.0;
    }
    double rowscale(Int i) const {
        return rowscale_.size() > 0 ? rowscale_[i] : 1.0;
    }

    // User model and attributes.
    const UserModel& user_model_;
    Int num_constr_{0};           // # constraints
    Int num_eqconstr_{0};         // # equality constraints
    Int num_var_{0};              // # variables
    Int num_free_var_{0};         // # free variables
    std::vector<Int> boxed_vars_; // indices of boxed variables

    // Solver model.
    Model& model_;

    // Presolve information.
    bool dualized_{false};
    std::vector<Int> flipped_vars_; // user variables flipped for dualization

    // Scaling factors. Empty when no scaling was applied.
    Vector colscale_;
    Vector rowscale_;
};

}  // namespace ipx

#endif  // IPX_PRESOLVER_H_
