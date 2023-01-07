// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#ifndef IPX_MODEL_H_
#define IPX_MODEL_H_

#include <vector>
#include "control.h"
#include "sparse_matrix.h"
#include "user_model.h"

namespace ipx {

// Model holds an LP model the computational form
//
//   minimize   c'x
//   subject to AI*x = b,                              (dual: y)
//              x-xl = lb, xl >= 0,                    (dual: zl >= 0)
//              x+xu = ub, xu >= 0.                    (dual: zu >= 0)
//
// The matrix AI has m >= 0 rows and n+m columns, where n > 0 is the number of
// "structural" variables. The last m columns of AI form the identity matrix.
// The last m components of c do not need to be zero (can happen when the model
// was dualized by presolver). Entries of -lb and ub can be infinity.
//
// A Model object cannot be modified other than discarding the data and loading
// a new model through a Presolver object.

class Model {
public:
    // Constructs an empty model.
    Model() = default;

    // Writes statistics of solver model to @info.
    void GetInfo(Info* info) const;

    // Returns true if the model is empty.
    bool empty() const { return cols() == 0; }

    // Deallocates all memory; the object becomes empty.
    void clear();

    // Returns the number of rows of AI.
    Int rows() const { return num_rows_; }

    // Returns the number of structural columns of AI (i.e. without the
    // rightmost identity matrix).
    Int cols() const { return num_cols_; }

    // Returns the number of columns classified as dense.
    Int num_dense_cols() const { return num_dense_cols_; }

    // Returns true if column j is classified as dense (0 <= j < n+m).
    bool IsDenseColumn(Int j) const {
        return AI_.entries(j) >= nz_dense_;
    }

    // Returns true if the user model was dualized by presolver.
    bool dualized() const { return dualized_; }

    // Returns a reference to the matrix AI in CSC and CSR format.
    const SparseMatrix& AI() const { return AI_; }
    const SparseMatrix& AIt() const { return AIt_; }

    // Returns a reference to a model vector.
    const Vector& b() const { return b_; }
    const Vector& c() const { return c_; }
    const Vector& lb() const { return lb_; }
    const Vector& ub() const { return ub_; }

    // Returns an entry of a model vector.
    double b(Int i) const { return b_[i]; }
    double c(Int j) const { return c_[j]; }
    double lb(Int j) const { return lb_[j]; }
    double ub(Int j) const { return ub_[j]; }

    // Returns the infinity norm of [b; lb; ub], ignoring infinite entries.
    double norm_bounds() const { return norm_bounds_; }

    // Returns the infinity norm of c.
    double norm_c() const { return norm_c_; }

private:
    // Presolver can access data members to initialize them.
    friend class Presolver;

    // Initializes num_dense_cols_ and nz_dense_. We classify the maximum #
    // columns as "dense" which have more than 40 nonzeros and more than 10
    // times the # nonzeros than any column that is not "dense". If this yields
    // more than 1000 dense columns, then no columns are classified as dense.
    void FindDenseColumns();

    // Computes norms of model data.
    void ComputeNorms();

    bool dualized_{false};        // model was dualized by presolver?
    Int num_rows_{0};             // # rows of AI
    Int num_cols_{0};             // # structural columns of AI
    Int num_dense_cols_{0};       // # columns classified as dense
    Int nz_dense_{0};             // minimum # nonzeros in a dense column
    SparseMatrix AI_;             // matrix AI columnwise
    SparseMatrix AIt_;            // matrix AI rowwise
    Vector b_;
    Vector c_;
    Vector lb_;
    Vector ub_;
    double norm_bounds_{0.0};     // infinity norm of [b;lb;ub]
    double norm_c_{0.0};          // infinity norm of c
};

// Returns the maximum violation of lb <= x <= ub.
double PrimalInfeasibility(const Model& model, const Vector& x);

// Returns the maximum violation of the dual feasibility condition
//   z[j] <= 0 if x[j] > lb[j],
//   z[j] >= 0 if x[j] < ub[j].
// Note that dual feasibility implies complementarity, i.e.
//   x[j] == lb[j] || x[j] == ub[j] || z[j] == 0.
double DualInfeasibility(const Model& model, const Vector& x, const Vector& z);

// Returns the maximum violation of Ax=b.
double PrimalResidual(const Model& model, const Vector& x);

// Returns the maximum violation of A'y+z=c.
double DualResidual(const Model& model, const Vector& y, const Vector& z);

}  // namespace ipx

#endif  // IPX_MODEL_H_
