#pragma once

#include <optional>
#include <string>
#include <vector>

namespace molgr::config
{
    struct ForceFieldConfig
    {
        std::vector<std::string> auto_force_fields_metal_free{"uff"};
        std::vector<std::string> auto_force_fields_with_metals{"uff"};
        std::string organic_force_field = "auto";
        std::string selection_force_field = "auto";
        std::string combined_force_field = "uff";
    };

    struct ResonanceConfig
    {
        int max_depth = 2;
        int limited_discrepancy_max_discrepancy = 1;
        std::string traversal_score = "direct_gain";
    };

    struct CppBackendConfig
    {
        std::optional<int> max_threads;
        bool enable_target_bucket_parallelism = true;
        bool enable_candidate_scoring_parallelism = false;
        bool enable_uff_atom_typing_cache = true;
        int candidate_score_parallel_threshold = 32;
    };

    struct MetalScoringConfig
    {
        double open_shell_multimetal_state_penalty_window = 10.0;
        int open_shell_multimetal_min_state_options = 6;
        int same_element_multimetal_unify_threshold = 3;
        std::optional<int> max_mixed_valence_spread = 3;
        int max_assignments_per_target = 64;
        double metal_coordination_extra_tolerance_angstrom = 0.35;
        double metal_access_radius_scale = 1.0;
        double metal_access_clearance_angstrom = 0.0;
    };

    struct MetalRadicalInferenceConfig
    {
        double coordination_cutoff_angstrom = 3.2;
        int max_considered_donors = 6;
        double square_planar_planarity_tolerance_angstrom = 0.45;
        double trigonal_planar_planarity_tolerance_angstrom = 0.35;
        double linear_angle_min_degrees = 150.0;
        double strong_field_threshold = 1.10;
        double weak_field_threshold = 0.75;
    };

    struct MolGRConfig
    {
        ForceFieldConfig force_field;
        ResonanceConfig resonance;
        CppBackendConfig cpp_backend;
        MetalScoringConfig metal_scoring;
        MetalRadicalInferenceConfig metal_radical_inference;
    };

    const MolGRConfig &GetDefaultConfig();
    void SetDefaultConfig(const MolGRConfig &config);
}
