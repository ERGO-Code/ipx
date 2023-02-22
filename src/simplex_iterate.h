// Copyright (c) 2023 ERGO-Code. See license.txt for license.

#ifndef IPX_SIMPLEX_ITERATE_H_
#define IPX_SIMPLEX_ITERATE_H_

#include "ipx_internal.h"
#include "model.h"

namespace ipx {

// SimplexIterate stores a complementary point x[n+m], y[m], z[n+m], where m, n
// are the number of rows and structural columns of the model. Complementary
// means that for each 0 <= j < n+m either x[j] is at its lower or upper bound
// or z[j] is zero.

class SimplexIterate {
    const Model& model_;

public:
    Vector x, y, z;

    // Constructs a SimplexIterate object associated with a model. @model must
    // be valid as long as the object is used; no data is copied.
    explicit SimplexIterate(const Model& model)
        : model_{model},
        x(model.cols() + model.rows()),
        y(model.rows()),
        z(model.cols() + model.rows())
    {}
};

}  // namespace ipx

#endif  // IPX_SIMPLEX_ITERATE_H_
