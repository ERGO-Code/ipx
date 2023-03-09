// Copyright (c) 2023 ERGO-Code. See license.txt for license.

#ifndef IPX_SRC_VERSION_H_
#define IPX_SRC_VERSION_H_

#include "ipx_internal.h"
#include "ipx_version.h"
#include <iosfwd>

namespace ipx {

struct Version {
    static constexpr ipxint major() { return IPX_VERSION_MAJOR; }
    static constexpr ipxint minor() { return IPX_VERSION_MINOR; }
    static constexpr ipxint patch() { return IPX_VERSION_PATCH; }
};

std::ostream& operator<<(std::ostream& os, Version const& version);

}  // namespace ipx

#endif  // IPX_SRC_VERSION_H_
