#pragma once

#include "molgr/config.h"
#include "molgr/types.h"

#include <openbabel/mol.h>

#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

namespace molgr
{
    namespace state
    {
        using MetadataValue = std::variant<int, double, bool, std::string>;
        using MetadataMap = std::unordered_map<std::string, MetadataValue>;
        using CombinedOmolBuilder = std::function<std::shared_ptr<OpenBabel::OBMol>(
            const OpenBabel::OBMol &,
            const std::vector<molgr::MetalAtomPosition> &)>;

        struct PreheatedNoMetalScoreBundle
        {
            std::string post_reinsertion_base_key;
            std::string force_field_config_key;
            std::string force_field_requested;
            std::string force_field_resolved_force_field;
            double organic_core_score = 0.0;
            double base_symmetry_penalty = 0.0;
            molgr::ChargedAtomSnapshotList charged_atom_snapshots;
        };

        struct ReconstructionState
        {
            std::shared_ptr<OpenBabel::OBMol> omol;
            int given_charge = 0;
            int total_charge = 0;
            int total_radical_electrons = 0;
            std::vector<std::string> phase_history;
            mutable MetadataMap metadata;
            int omol_revision = 0;

            mutable std::optional<std::pair<int, std::string>> organic_score_key_cache;
            mutable std::optional<std::pair<int, double>> organic_core_score_cache;
            mutable std::optional<std::pair<int, double>> full_score_cache;
            mutable std::optional<std::pair<int, double>> post_reinsertion_base_symmetry_penalty_cache;
            mutable std::optional<std::pair<int, molgr::ChargedAtomSnapshotList>>
                post_reinsertion_charged_atom_snapshots_cache;
            mutable std::shared_ptr<const PreheatedNoMetalScoreBundle> preheated_score_bundle;

            ReconstructionState() = default;

            ReconstructionState(
                std::shared_ptr<OpenBabel::OBMol> omol_,
                int given_charge_,
                int total_charge_,
                int total_radical_electrons_,
                std::vector<std::string> phase_history_ = {},
                MetadataMap metadata_ = {},
                int omol_revision_ = 0);

            const OpenBabel::OBMol &Mol() const;
            OpenBabel::OBMol &MutableMol();

            void InvalidateOmolDerivedCache();
            PreheatedNoMetalScoreBundle BuildPreheatedScoreBundle(
                const molgr::config::MolGRConfig &config = molgr::config::GetDefaultConfig()) const;
            void PreheatScoreBundle(
                const molgr::config::MolGRConfig &config = molgr::config::GetDefaultConfig());
            const PreheatedNoMetalScoreBundle *PreheatedScoreBundle() const;
            std::string PostReinsertionBaseKey() const;
            double OrganicCoreScore(
                const molgr::config::MolGRConfig &config = molgr::config::GetDefaultConfig()) const;
            std::pair<double, molgr::ChargedAtomSnapshotList> PostReinsertionBaseComponents() const;
            double FullScore(
                const molgr::config::MolGRConfig &config = molgr::config::GetDefaultConfig()) const;
        };

        class OmolStateMachine
        {
        public:
            std::shared_ptr<OpenBabel::OBMol> omol;
            int given_charge = 0;
            std::vector<std::string> phase_history;
            mutable MetadataMap metadata;
            int omol_revision = 0;

            mutable std::optional<std::pair<int, std::string>> organic_score_key_cache;
            mutable std::optional<std::pair<int, double>> organic_core_score_cache;
            mutable std::optional<std::pair<int, double>> full_score_cache;
            mutable std::optional<std::pair<int, double>> post_reinsertion_base_symmetry_penalty_cache;
            mutable std::optional<std::pair<int, molgr::ChargedAtomSnapshotList>>
                post_reinsertion_charged_atom_snapshots_cache;
            mutable std::shared_ptr<const PreheatedNoMetalScoreBundle> preheated_score_bundle;

            OmolStateMachine() = default;

            OmolStateMachine(
                std::shared_ptr<OpenBabel::OBMol> omol_,
                int given_charge_ = 0,
                std::vector<std::string> phase_history_ = {},
                MetadataMap metadata_ = {},
                int omol_revision_ = 0);

            static OmolStateMachine FromReconstructionState(const ReconstructionState &state);

            OpenBabel::OBMol &EnsureUniqueMol();
            void Annotate(const std::optional<std::string> &phase);
            void SetGivenCharge(const std::optional<std::string> &phase, int next_given_charge);

            template <typename Stage, typename... Args>
            bool RunOmolStage(const std::optional<std::string> &phase, Stage &&stage, Args &&...args)
            {
                OpenBabel::OBMol &mol = EnsureUniqueMol();
                const bool hit = std::invoke(
                    std::forward<Stage>(stage),
                    mol,
                    std::forward<Args>(args)...);
                if (hit)
                {
                    ++omol_revision;
                    InvalidateOmolDerivedCache();
                }
                Annotate(phase);
                return hit;
            }

            template <typename Stage, typename... Args>
            bool RunOmolChargeStage(
                const std::optional<std::string> &phase,
                Stage &&stage,
                Args &&...args)
            {
                OpenBabel::OBMol &mol = EnsureUniqueMol();
                const int before_charge = given_charge;
                const bool molecule_hit = std::invoke(
                    std::forward<Stage>(stage),
                    mol,
                    given_charge,
                    std::forward<Args>(args)...);
                const bool charge_hit = before_charge != given_charge;
                if (molecule_hit)
                {
                    ++omol_revision;
                    InvalidateOmolDerivedCache();
                }
                Annotate(phase);
                return molecule_hit || charge_hit;
            }

            void InvalidateOmolDerivedCache();

            OmolStateMachine Branch(
                const std::optional<std::string> &phase,
                std::shared_ptr<OpenBabel::OBMol> next_omol = nullptr,
                std::optional<int> next_given_charge = std::nullopt,
                MetadataMap branch_metadata = {}) const;

            ReconstructionState Freeze(int total_charge, int total_radical_electrons) const;
            ReconstructionState FreezeLike(const ReconstructionState &state) const;
        };

        struct MetalPreparationState
        {
            std::string no_metal_xyz_block;
            std::vector<std::vector<molgr::MetalAtomPosition>> available_valence_radical_states;
            int total_charge = 0;
            int total_radical_electrons = 0;
            std::vector<std::string> phase_history;
            mutable MetadataMap metadata;
        };

        struct MetalCandidateState
        {
            std::vector<molgr::MetalAtomPosition> metal_states;
            int no_metal_charge_target = 0;
            int no_metal_radical_target = 0;
            std::vector<std::string> phase_history;
            mutable MetadataMap metadata;
            std::shared_ptr<ReconstructionState> no_metal_state;
            std::shared_ptr<OpenBabel::OBMol> combined_omol;
            mutable std::optional<double> score;

            mutable std::optional<std::string> metal_state_key_cache;
            mutable std::optional<std::pair<std::string, std::string>> combined_score_key_cache;
            mutable std::optional<std::pair<std::string, double>> combined_score_cache;
            mutable std::optional<std::tuple<const ReconstructionState *, int, std::string>>
                combined_omol_dependency_cache;

            std::string MetalStateKey() const;
            std::pair<std::string, std::string> CombinedScoreKey() const;
            std::shared_ptr<OpenBabel::OBMol> MaterializeCombinedOmol(
                const CombinedOmolBuilder &builder);
            double CombinedScore(
                const molgr::config::MolGRConfig &config = molgr::config::GetDefaultConfig()) const;
        };

        class MetalCandidateStateMachine
        {
        public:
            std::vector<molgr::MetalAtomPosition> metal_states;
            int no_metal_charge_target = 0;
            int no_metal_radical_target = 0;
            std::vector<std::string> phase_history;
            MetadataMap metadata;
            std::shared_ptr<ReconstructionState> no_metal_state;
            std::shared_ptr<OpenBabel::OBMol> combined_omol;
            std::optional<double> score;

            std::optional<std::string> metal_state_key_cache;
            std::optional<std::pair<std::string, std::string>> combined_score_key_cache;
            std::optional<std::pair<std::string, double>> combined_score_cache;
            std::optional<std::tuple<const ReconstructionState *, int, std::string>>
                combined_omol_dependency_cache;

            MetalCandidateStateMachine() = default;

            explicit MetalCandidateStateMachine(
                std::vector<molgr::MetalAtomPosition> metal_states_,
                int no_metal_charge_target_,
                int no_metal_radical_target_,
                std::vector<std::string> phase_history_ = {},
                MetadataMap metadata_ = {});

            static MetalCandidateStateMachine FromCandidateState(const MetalCandidateState &state);

            void Annotate(const std::optional<std::string> &phase);
            void SetNoMetalState(
                const std::optional<std::string> &phase,
                std::shared_ptr<ReconstructionState> next_no_metal_state);
            MetalCandidateState Freeze() const;
        };
    }
}
