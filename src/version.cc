// Copyright (c) 2023 ERGO-Code. See license.txt for license.

#include "version.h"
#include <ostream>

namespace ipx {

std::ostream& operator<<(std::ostream& os, Version const& version) {
    os  << version.major() << '.'
        << version.minor() << '.'
        << version.patch();
    return os;
}

}  // namespace ipx
