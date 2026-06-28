#include "molgr/utils/no_metals/selection.h"

#include "molgr/utils/organic_topology.h"

#include <openbabel/obiter.h>

#include <cmath>

namespace molgr
{
    namespace no_metals
    {
        namespace selection
        {
            namespace
            {
                int FormalChargeAbsoluteSum(const OpenBabel::OBMol &mol)
                {
                    int formal_charge_absolute_sum = 0;
                    for (int atom_idx = 1; atom_idx <= mol.NumAtoms(); ++atom_idx)
                    {
                        const OpenBabel::OBAtom *atom = mol.GetAtom(atom_idx);
                        if (atom == nullptr)
                        {
                            continue;
                        }
                        formal_charge_absolute_sum += std::abs(atom->GetFormalCharge());
                    }
                    return formal_charge_absolute_sum;
                }
            }

            molgr::organic_topology::OrganicTopologyMetrics AnnotateNoMetalCandidateTopology(
                molgr::state::ReconstructionState &candidate,
                const molgr::config::MolGRConfig &config)
            {
                const auto metrics = molgr::organic_topology::ComputeOrganicTopologyMetrics(
                    candidate.Mol(),
                    config.organic_topology);
                candidate.metadata["organic_aromatic_atom_count"] = metrics.aromatic_atom_count;
                candidate.metadata["organic_aromatic_ring_count"] = metrics.aromatic_ring_count;
                candidate.metadata["organic_aromatic_stability_score"] = metrics.aromatic_stability_score;
                candidate.metadata["organic_conjugated_atom_count"] = metrics.conjugated_atom_count;
                candidate.metadata["organic_conjugated_bond_count"] = metrics.conjugated_bond_count;
                candidate.metadata["organic_max_conjugated_component_size"] =
                    metrics.max_conjugated_component_size;
                return metrics;
            }

            double ScoreReconstructionCandidate(
                molgr::state::ReconstructionState &candidate,
                const molgr::config::MolGRConfig &config)
            {
                const double score = candidate.FullScore(config);
                candidate.metadata["score"] = score;
                return score;
            }

            NoMetalTopologySelectionKey NoMetalCandidateTopologySelectionKey(
                molgr::state::ReconstructionState &candidate,
                const molgr::config::MolGRConfig &config)
            {
                const auto metrics = AnnotateNoMetalCandidateTopology(candidate, config);
                const int formal_charge_absolute_sum = FormalChargeAbsoluteSum(candidate.Mol());
                const double conjugation_charge_penalty =
                    static_cast<double>(formal_charge_absolute_sum) / 2.0;
                const double adjusted_max_conjugated_component_size =
                    static_cast<double>(metrics.max_conjugated_component_size) -
                    conjugation_charge_penalty;
                const double adjusted_conjugated_atom_count =
                    static_cast<double>(metrics.conjugated_atom_count) - conjugation_charge_penalty;
                const double adjusted_conjugated_bond_count =
                    static_cast<double>(metrics.conjugated_bond_count) - conjugation_charge_penalty;
                candidate.metadata["organic_formal_charge_absolute_sum"] =
                    formal_charge_absolute_sum;
                candidate.metadata["organic_conjugation_charge_penalty"] =
                    conjugation_charge_penalty;
                candidate.metadata["organic_adjusted_max_conjugated_component_size"] =
                    adjusted_max_conjugated_component_size;
                candidate.metadata["organic_adjusted_conjugated_atom_count"] =
                    adjusted_conjugated_atom_count;
                candidate.metadata["organic_adjusted_conjugated_bond_count"] =
                    adjusted_conjugated_bond_count;
                return {
                    -metrics.aromatic_stability_score,
                    -metrics.aromatic_atom_count,
                    -adjusted_max_conjugated_component_size,
                    -adjusted_conjugated_atom_count,
                    -adjusted_conjugated_bond_count,
                };
            }

            std::tuple<double, int, double, double, double, double> NoMetalCandidateSelectionKey(
                molgr::state::ReconstructionState &candidate,
                const molgr::config::MolGRConfig &config)
            {
                const auto topology_key = NoMetalCandidateTopologySelectionKey(candidate, config);
                const double score = ScoreReconstructionCandidate(candidate, config);
                return {
                    std::get<0>(topology_key),
                    std::get<1>(topology_key),
                    std::get<2>(topology_key),
                    std::get<3>(topology_key),
                    std::get<4>(topology_key),
                    score,
                };
            }
        }
    }
}
