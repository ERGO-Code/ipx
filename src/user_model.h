// Copyright (c) 2023 ERGO-Code. See license.txt for license.

#ifndef IPX_USER_MODEL_H_
#define IPX_USER_MODEL_H_

#include <vector>
#include "control.h"
#include "solution.h"
#include "sparse_matrix.h"

namespace ipx {

// UserModel holds an LP model in the form given by the user,
//
//   minimize   obj'x                                                    (1)
//   subject to A*x {=,<,>} rhs, lb <= x <= ub.

class UserModel {
public:
    // Constructs an empty model.
    UserModel() = default;

    // Initializes a UserModel object from the form (1).
    // @num_constr: number of rows of A
    // @num_var: number of columns of A
    // @Ap, @Ai, @Ax: matrix A in CSC format, 0-based indexing
    // @rhs: array of size num_constr
    // @constr_type: array of size num_constr with entries '=', '<' or '>'
    // @obj: array of size num_var
    // @lb: array of size num_var, entries can be -INFINITY
    // @ub: array of size num_var, entries can be +INFINITY
    // If the input is invalid an error code is returned and the object becomes
    // empty.
    // Returns:
    //  0
    //  IPX_ERROR_argument_null
    //  IPX_ERROR_invalid_dimension
    //  IPX_ERROR_invalid_matrix
    //  IPX_ERROR_invalid_vector
    Int Load(const Control& control, Int num_constr, Int num_var,
             const Int* Ap, const Int* Ai, const double* Ax,
             const double* rhs, const char* constr_type, const double* obj,
             const double* lb, const double* ub);

    // Writes statistics of input data to @info.
    void GetInfo(Info* info) const;

    // Returns true if the model is empty.
    bool empty() const { return empty_; }

    // Deallocates all memory; the object becomes empty.
    void clear();

    // Returns the number of variables.
    Int num_var() const { return num_var_; }

    // Returns the number of constraints.
    Int num_constr() const { return num_constr_; }

    // Returns a reference to the matrix A in CSC format.
    const SparseMatrix& A() const { return A_; }

    // Returns a reference to a model vector.
    const Vector& obj() const { return obj_; }
    const std::vector<char>& constr_type() const { return constr_type_; }
    const Vector& rhs() const { return rhs_; }
    const Vector& lb() const { return lb_; }
    const Vector& ub() const { return ub_; }

    // Returns an entry of a model vector.
    double obj(Int j) const { return obj_[j]; }
    char constr_type(Int i) const { return constr_type_[i]; }
    double rhs(Int i) const { return rhs_[i]; }
    double lb(Int j) const { return lb_[j]; }
    double ub(Int j) const { return ub_[j]; }

    // Checks whether the vectors hold a valid interior point to the model.
    // Returns:
    //  0 if point is valid
    //  IPX_ERROR_argument_null
    //  IPX_ERROR_invalid_vector
    Int CheckInteriorPoint(const double* x, const double* xl, const double* xu,
                           const double* slack, const double* y,
                           const double* zl, const double* zu) const;

    // Evaluates the residuals, objective, etc. of an interior point. The
    // following info members are set: abs_presidual, abs_dresidual,
    // rel_presidual, rel_dresidual, pobjval, dobjval, rel_objgap,
    // complementarity, normx, normy, normz.
    void EvaluateInteriorPoint(const InteriorSolution& point,
                               Info* info) const;

    // Evaluates the infeasibilities, objective, etc. of a basic point. The
    // following info members are set: primal_infeas, dual_infeas, objval
    void EvaluateBasicPoint(const BasicSolution& point, Info* info) const;

private:
    // Checks that input is valid and copies into the member variables. If the
    // input is invalid, an error code is returned and the object remains
    // unchanged.
    // Returns:
    //  0
    //  IPX_ERROR_argument_null
    //  IPX_ERROR_invalid_dimension
    //  IPX_ERROR_invalid_matrix
    //  IPX_ERROR_invalid_vector
    Int CopyInput(Int num_constr, Int num_var, const Int* Ap, const Int* Ai,
                  const double* Ax, const double* rhs, const char* constr_type,
                  const double* obj, const double* lb,
                  const double* ub);

    // Computes norms of model data.
    void ComputeNorms();

    bool empty_{true};
    Int num_var_{0};
    Int num_constr_{0};
    Vector obj_;
    std::vector<char> constr_type_;
    Vector rhs_;
    Vector lb_;
    Vector ub_;
    SparseMatrix A_;
    double norm_obj_{0.0};
    double norm_rhs_{0.0};
    double norm_bounds_{0.0};
};

}  // namespace ipx

#endif  // IPX_USER_MODEL_H_
