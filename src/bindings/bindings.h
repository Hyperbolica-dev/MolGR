/*
 * @Author: TMJ
 * @Date: 2026-01-01 23:12:22
 * @LastEditors: TMJ
 * @LastEditTime: 2026-01-01 23:35:25
 * @Description: 请填写简介
 */
#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace molgr {
namespace bind {

    // 声明各个子模块的绑定入口函数
    void bind_consts(py::module_ &m);
    void bind_utils(py::module_ &m);
    void bind_metal(py::module_ &m);
    void bind_scoring(py::module_ &m);
    void bind_reconstruct(py::module_ &m);

    // 如果你有专门的测试函数，也可以单独放一个
    // void bind_testing(py::module_ &m); 

} // namespace bind
} // namespace molgr