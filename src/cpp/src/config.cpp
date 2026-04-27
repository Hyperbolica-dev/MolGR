#include "molgr/python_config.h"

#include <mutex>

#include <pybind11/stl.h>

namespace py = pybind11;

namespace molgr::config
{
    namespace
    {
        std::mutex &ConfigMutex()
        {
            static std::mutex mutex;
            return mutex;
        }

        MolGRConfig &DefaultConfigStorage()
        {
            static MolGRConfig config;
            return config;
        }

        py::object AttrOrNone(py::handle object, const char *name)
        {
            if (object.is_none() || !py::hasattr(object, name))
            {
                return py::none();
            }
            return py::reinterpret_borrow<py::object>(object.attr(name));
        }

        template <typename T>
        T CastAttrOr(py::handle object, const char *name, T fallback)
        {
            py::object value = AttrOrNone(object, name);
            if (value.is_none())
            {
                return fallback;
            }
            return value.cast<T>();
        }

        std::optional<int> CastOptionalIntAttrOr(
            py::handle object,
            const char *name,
            std::optional<int> fallback)
        {
            py::object value = AttrOrNone(object, name);
            if (value.is_none())
            {
                return fallback;
            }
            return value.cast<int>();
        }

        std::vector<std::string> CastStringVectorAttrOr(
            py::handle object,
            const char *name,
            std::vector<std::string> fallback)
        {
            py::object value = AttrOrNone(object, name);
            if (value.is_none())
            {
                return fallback;
            }
            return value.cast<std::vector<std::string>>();
        }

        std::vector<double> CastDoubleVectorAttrOr(
            py::handle object,
            const char *name,
            std::vector<double> fallback)
        {
            py::object value = AttrOrNone(object, name);
            if (value.is_none())
            {
                return fallback;
            }
            return value.cast<std::vector<double>>();
        }

        py::object ResolvePythonConfig(py::handle config)
        {
            if (!config.is_none())
            {
                return py::reinterpret_borrow<py::object>(config);
            }
            py::module_ config_module = py::module_::import("molgr.config");
            return config_module.attr("get_config")();
        }
    }

    const MolGRConfig &GetDefaultConfig()
    {
        std::lock_guard<std::mutex> lock(ConfigMutex());
        return DefaultConfigStorage();
    }

    void SetDefaultConfig(const MolGRConfig &config)
    {
        std::lock_guard<std::mutex> lock(ConfigMutex());
        DefaultConfigStorage() = config;
    }

    MolGRConfig FromPython(py::handle config)
    {
        py::object resolved = ResolvePythonConfig(config);
        MolGRConfig out;

        py::object force_field = AttrOrNone(resolved, "force_field");
        if (!force_field.is_none())
        {
            out.force_field.auto_force_fields_metal_free = CastStringVectorAttrOr(
                force_field,
                "auto_force_fields_metal_free",
                out.force_field.auto_force_fields_metal_free);
            out.force_field.auto_force_fields_with_metals = CastStringVectorAttrOr(
                force_field,
                "auto_force_fields_with_metals",
                out.force_field.auto_force_fields_with_metals);
            out.force_field.organic_force_field = CastAttrOr<std::string>(
                force_field,
                "organic_force_field",
                out.force_field.organic_force_field);
            out.force_field.selection_force_field = CastAttrOr<std::string>(
                force_field,
                "selection_force_field",
                out.force_field.selection_force_field);
            out.force_field.combined_force_field = CastAttrOr<std::string>(
                force_field,
                "combined_force_field",
                out.force_field.combined_force_field);
        }

        py::object resonance = AttrOrNone(resolved, "resonance");
        if (!resonance.is_none())
        {
            out.resonance.max_depth = CastAttrOr<int>(resonance, "max_depth", out.resonance.max_depth);
            out.resonance.limited_discrepancy_max_discrepancy = CastAttrOr<int>(
                resonance,
                "limited_discrepancy_max_discrepancy",
                out.resonance.limited_discrepancy_max_discrepancy);
            out.resonance.traversal_score = CastAttrOr<std::string>(
                resonance,
                "traversal_score",
                out.resonance.traversal_score);
        }

        py::object cpp_backend = AttrOrNone(resolved, "cpp_backend");
        if (!cpp_backend.is_none())
        {
            out.cpp_backend.max_threads = CastOptionalIntAttrOr(
                cpp_backend,
                "max_threads",
                out.cpp_backend.max_threads);
            out.cpp_backend.enable_target_bucket_parallelism = CastAttrOr<bool>(
                cpp_backend,
                "enable_target_bucket_parallelism",
                out.cpp_backend.enable_target_bucket_parallelism);
            out.cpp_backend.enable_candidate_scoring_parallelism = CastAttrOr<bool>(
                cpp_backend,
                "enable_candidate_scoring_parallelism",
                out.cpp_backend.enable_candidate_scoring_parallelism);
            out.cpp_backend.enable_resonance_candidate_parallelism = CastAttrOr<bool>(
                cpp_backend,
                "enable_resonance_candidate_parallelism",
                out.cpp_backend.enable_resonance_candidate_parallelism);
            out.cpp_backend.enable_uff_atom_typing_cache = CastAttrOr<bool>(
                cpp_backend,
                "enable_uff_atom_typing_cache",
                out.cpp_backend.enable_uff_atom_typing_cache);
            out.cpp_backend.resonance_candidate_parallel_threshold = CastAttrOr<int>(
                cpp_backend,
                "resonance_candidate_parallel_threshold",
                out.cpp_backend.resonance_candidate_parallel_threshold);
            out.cpp_backend.candidate_score_parallel_threshold = CastAttrOr<int>(
                cpp_backend,
                "candidate_score_parallel_threshold",
                out.cpp_backend.candidate_score_parallel_threshold);
        }

        py::object metal_scoring = AttrOrNone(resolved, "metal_scoring");
        if (!metal_scoring.is_none())
        {
            out.metal_scoring.organic_score_bucket_relative_ratio = CastAttrOr<double>(
                metal_scoring,
                "organic_score_bucket_relative_ratio",
                out.metal_scoring.organic_score_bucket_relative_ratio);
            out.metal_scoring.organic_force_field_hard_max_ratio = CastAttrOr<double>(
                metal_scoring,
                "organic_force_field_hard_max_ratio",
                out.metal_scoring.organic_force_field_hard_max_ratio);
            out.metal_scoring.open_shell_multimetal_state_penalty_window = CastAttrOr<double>(
                metal_scoring,
                "open_shell_multimetal_state_penalty_window",
                out.metal_scoring.open_shell_multimetal_state_penalty_window);
            out.metal_scoring.open_shell_multimetal_min_state_options = CastAttrOr<int>(
                metal_scoring,
                "open_shell_multimetal_min_state_options",
                out.metal_scoring.open_shell_multimetal_min_state_options);
            out.metal_scoring.same_element_multimetal_unify_threshold = CastAttrOr<int>(
                metal_scoring,
                "same_element_multimetal_unify_threshold",
                out.metal_scoring.same_element_multimetal_unify_threshold);
            out.metal_scoring.max_mixed_valence_spread = CastOptionalIntAttrOr(
                metal_scoring,
                "max_mixed_valence_spread",
                out.metal_scoring.max_mixed_valence_spread);
            out.metal_scoring.max_assignments_per_target = CastAttrOr<int>(
                metal_scoring,
                "max_assignments_per_target",
                out.metal_scoring.max_assignments_per_target);
            out.metal_scoring.selection_weight_values = CastDoubleVectorAttrOr(
                metal_scoring,
                "selection_weight_values",
                out.metal_scoring.selection_weight_values);
            out.metal_scoring.selection_scale_values = CastDoubleVectorAttrOr(
                metal_scoring,
                "selection_scale_values",
                out.metal_scoring.selection_scale_values);
            out.metal_scoring.visible_coordination_reward_weight = CastAttrOr<double>(
                metal_scoring,
                "visible_coordination_reward_weight",
                out.metal_scoring.visible_coordination_reward_weight);
        }

        py::object metal_radical = AttrOrNone(resolved, "metal_radical_inference");
        if (!metal_radical.is_none())
        {
            out.metal_radical_inference.coordination_cutoff_angstrom = CastAttrOr<double>(
                metal_radical,
                "coordination_cutoff_angstrom",
                out.metal_radical_inference.coordination_cutoff_angstrom);
            out.metal_radical_inference.max_considered_donors = CastAttrOr<int>(
                metal_radical,
                "max_considered_donors",
                out.metal_radical_inference.max_considered_donors);
            out.metal_radical_inference.strong_field_threshold = CastAttrOr<double>(
                metal_radical,
                "strong_field_threshold",
                out.metal_radical_inference.strong_field_threshold);
            out.metal_radical_inference.weak_field_threshold = CastAttrOr<double>(
                metal_radical,
                "weak_field_threshold",
                out.metal_radical_inference.weak_field_threshold);
        }

        return out;
    }

    void SetDefaultConfigFromPython(py::handle config)
    {
        SetDefaultConfig(FromPython(config));
    }

    py::dict DefaultConfigSummary()
    {
        const MolGRConfig &config = GetDefaultConfig();
        py::dict out;
        out["resonance_max_depth"] = config.resonance.max_depth;
        out["resonance_max_discrepancy"] = config.resonance.limited_discrepancy_max_discrepancy;
        out["resonance_traversal_score"] = config.resonance.traversal_score;
        out["max_threads"] = config.cpp_backend.max_threads.has_value()
                                 ? py::cast(*config.cpp_backend.max_threads)
                                 : py::none();
        out["enable_uff_atom_typing_cache"] = config.cpp_backend.enable_uff_atom_typing_cache;
        out["max_mixed_valence_spread"] = config.metal_scoring.max_mixed_valence_spread.has_value()
                                              ? py::cast(*config.metal_scoring.max_mixed_valence_spread)
                                              : py::none();
        out["max_assignments_per_target"] = config.metal_scoring.max_assignments_per_target;
        return out;
    }
}
