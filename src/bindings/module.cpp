#include <pybind11/pybind11.h>
#include "mylib/algorithm.h" // 包含你的头文件

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() = "Refactored OpenBabel Extension";
    
    // 将 mylib::calculate_atom_count 暴露为 get_atom_count
    m.def("get_atom_count", &mylib::calculate_atom_count, "Calculate atom count");
}