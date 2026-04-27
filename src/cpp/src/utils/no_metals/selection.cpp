#include "molgr/utils/no_metals/selection.h"

#include "molgr/utils/organic_topology.h"

namespace molgr
{
    namespace no_metals
    {
        namespace selection
        {
            molgr::organic_topology::OrganicTopologyMetrics AnnotateNoMetalCandidateTopology(
                molgr::state::ReconstructionState &candidate)
            {
                const auto metrics = molgr::organic_topology::ComputeOrganicTopologyMetrics(candidate.Mol());
                candidate.metadata["organic_aromatic_atom_count"] = metrics.aromatic_atom_count;
                candidate.metadata["organic_aromatic_ring_count"] = metrics.aromatic_ring_count;
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
                molgr::state::ReconstructionState &candidate)
            {
                const auto metrics = AnnotateNoMetalCandidateTopology(candidate);
                return {
                    -metrics.aromatic_atom_count,
                    -metrics.max_conjugated_component_size,
                    -metrics.conjugated_atom_count,
                    -metrics.conjugated_bond_count,
                };
            }

            std::tuple<int, int, int, int, double> NoMetalCandidateSelectionKey(
                molgr::state::ReconstructionState &candidate,
                const molgr::config::MolGRConfig &config)
            {
                const auto topology_key = NoMetalCandidateTopologySelectionKey(candidate);
                const double score = ScoreReconstructionCandidate(candidate, config);
                return {
                    std::get<0>(topology_key),
                    std::get<1>(topology_key),
                    std::get<2>(topology_key),
                    std::get<3>(topology_key),
                    score,
                };
            }
        }
    }
}
