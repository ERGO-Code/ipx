// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "indexed_vector.h"

namespace ipx {

IndexedVector::IndexedVector(Int dim) : elements_(dim), pattern_(dim) { }

bool IndexedVector::sparse() const {
    return nnz() >= 0 && nnz() <= kHypersparseThreshold * dim();
}

void IndexedVector::set_to_zero() {
    if (sparse()) {
        for (Int p = 0; p < nnz_; p++)
            elements_[pattern_[p]] = 0.0;
    } else {
        elements_ = 0.0;
    }
    nnz_ = 0;
}

}  // namespace ipx
