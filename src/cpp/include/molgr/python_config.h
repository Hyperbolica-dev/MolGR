#pragma once

#include "molgr/config.h"

#include <pybind11/pybind11.h>

namespace molgr::config
{
    MolGRConfig FromPython(pybind11::handle config);
}
