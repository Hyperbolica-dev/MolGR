#include "bindings.h"
#include "molgr/consts.h"

void molgr::bind::bind_consts(py::module_ &m)
{
    // 绑定 GetPossibleMetalRadicals
    m.def("get_possible_metal_radicals", &molgr::GetPossibleMetalRadicals,
          R"pbdoc(
            Get possible radical electron counts for a metal given its valence.
            
            Args:
                metal (str): The chemical symbol (e.g., "Fe").
                valence (int): The oxidation state.
            
            Returns:
                set[int]: A set of possible unpaired electron counts.
        )pbdoc",
          py::arg("metal"), py::arg("valence"));

    // 如果有其他常量导出，写在这里
}