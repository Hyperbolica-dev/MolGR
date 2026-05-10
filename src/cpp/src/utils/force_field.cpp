#include "molgr/utils/force_field.h"

#include "molgr/utils/lru_cache.h"
#include "molgr/utils/perf.h"
#include "molgr/vendor/forcefielduff.h"

#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdint>
#include <cmath>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace
{
    constexpr double kCoordinateScale = 1000000.0;
    constexpr std::size_t kDefaultForceFieldCacheMaxSize = 4096;
    using Clock = std::chrono::steady_clock;

    OpenBabel::OBMol &MutableMol(const OpenBabel::OBMol &mol)
    {
        return const_cast<OpenBabel::OBMol &>(mol);
    }

    std::int64_t QuantizedCoordinate(double value)
    {
        return static_cast<std::int64_t>(std::llround(value * kCoordinateScale));
    }

    std::string Lowercase(std::string value)
    {
        std::transform(
            value.begin(),
            value.end(),
            value.begin(),
            [](unsigned char ch)
            {
                return static_cast<char>(std::tolower(ch));
            });
        return value;
    }

    std::string Uppercase(std::string value)
    {
        std::transform(
            value.begin(),
            value.end(),
            value.begin(),
            [](unsigned char ch)
            {
                return static_cast<char>(std::toupper(ch));
            });
        return value;
    }

    std::string Trim(const std::string &value)
    {
        const auto begin = std::find_if_not(
            value.begin(),
            value.end(),
            [](unsigned char ch)
            {
                return std::isspace(ch) != 0;
            });
        const auto end = std::find_if_not(
            value.rbegin(),
            value.rend(),
            [](unsigned char ch)
            {
                return std::isspace(ch) != 0;
            }).base();
        if (begin >= end)
        {
            return {};
        }
        return std::string(begin, end);
    }

    double ForceFieldEnergyToKjMol(double raw_energy, const std::string &raw_unit)
    {
        const std::string unit = Lowercase(Trim(raw_unit));
        if (unit == "kj/mol")
        {
            return raw_energy;
        }
        if (unit == "kcal/mol")
        {
            return raw_energy * 4.184;
        }
        throw std::runtime_error("Unsupported force-field energy unit: " + raw_unit);
    }

    int CountHeavyAtoms(OpenBabel::OBMol &mol)
    {
        int count = 0;
        FOR_ATOMS_OF_MOL(atom_iter, mol)
        {
            if (atom_iter->GetAtomicNum() != 1)
            {
                ++count;
            }
        }
        return count;
    }

    std::unique_ptr<OpenBabel::MolgrForceFieldUFF> MakeForceFieldInstance(
        const std::string &force_field)
    {
        if (force_field != "uff")
        {
            return nullptr;
        }
        return std::make_unique<OpenBabel::MolgrForceFieldUFF>("MolGR-UFF", false);
    }

    void AppendValue(std::string &out, int value);
    void AppendValue(std::string &out, std::int64_t value);
    void AppendValue(std::string &out, bool value);

    struct ReusableForceField
    {
        std::unique_ptr<OpenBabel::MolgrForceFieldUFF> instance;
        std::string last_exact_setup_key;
        std::string last_openbabel_setup_key;
    };

    std::string BuildOpenBabelSetupKey(const OpenBabel::OBMol &mol)
    {
        OpenBabel::OBMol &mutable_mol = MutableMol(mol);
        std::string key;
        key.reserve(
            32 +
            static_cast<std::size_t>(mutable_mol.NumAtoms()) * 12 +
            static_cast<std::size_t>(mutable_mol.NumBonds()) * 16);

        key += "A";
        AppendValue(key, static_cast<int>(mutable_mol.NumAtoms()));
        key.push_back(':');
        FOR_ATOMS_OF_MOL(atom_iter, mutable_mol)
        {
            AppendValue(key, static_cast<int>(atom_iter->GetAtomicNum()));
            key.push_back(',');
            AppendValue(key, static_cast<int>(atom_iter->GetExplicitDegree()));
            key.push_back(';');
        }

        key += "|B";
        AppendValue(key, static_cast<int>(mutable_mol.NumBonds()));
        key.push_back(':');
        FOR_BONDS_OF_MOL(bond_iter, mutable_mol)
        {
            AppendValue(key, static_cast<int>(bond_iter->GetIdx()));
            key.push_back(',');
            AppendValue(key, static_cast<int>(bond_iter->GetBondOrder()));
            key.push_back(',');
            AppendValue(key, static_cast<int>(bond_iter->GetBeginAtom()->GetAtomicNum()));
            key.push_back(',');
            AppendValue(key, static_cast<int>(bond_iter->GetEndAtom()->GetAtomicNum()));
            key.push_back(';');
        }
        return key;
    }

    std::string BuildExactForceFieldSetupKey(const OpenBabel::OBMol &mol)
    {
        OpenBabel::OBMol &mutable_mol = MutableMol(mol);
        std::string key;
        key.reserve(
            32 +
            static_cast<std::size_t>(mutable_mol.NumAtoms()) * 24 +
            static_cast<std::size_t>(mutable_mol.NumBonds()) * 24);

        key += "A";
        AppendValue(key, static_cast<int>(mutable_mol.NumAtoms()));
        key.push_back(':');
        FOR_ATOMS_OF_MOL(atom_iter, mutable_mol)
        {
            AppendValue(key, static_cast<int>(atom_iter->GetAtomicNum()));
            key.push_back(',');
            AppendValue(key, atom_iter->GetFormalCharge());
            key.push_back(',');
            AppendValue(key, atom_iter->GetSpinMultiplicity());
            key.push_back(',');
            AppendValue(key, static_cast<int>(atom_iter->GetExplicitDegree()));
            key.push_back(',');
            AppendValue(key, atom_iter->IsAromatic());
            key.push_back(';');
        }

        std::vector<std::tuple<int, int, int, bool>> bond_keys;
        bond_keys.reserve(static_cast<std::size_t>(mutable_mol.NumBonds()));
        FOR_BONDS_OF_MOL(bond_iter, mutable_mol)
        {
            int begin_idx = bond_iter->GetBeginAtom()->GetIdx();
            int end_idx = bond_iter->GetEndAtom()->GetIdx();
            if (begin_idx > end_idx)
            {
                std::swap(begin_idx, end_idx);
            }
            bond_keys.emplace_back(
                begin_idx,
                end_idx,
                static_cast<int>(bond_iter->GetBondOrder()),
                bond_iter->IsAromatic());
        }
        std::sort(bond_keys.begin(), bond_keys.end());

        key += "|B";
        AppendValue(key, static_cast<int>(bond_keys.size()));
        key.push_back(':');
        for (const auto &bond_key : bond_keys)
        {
            AppendValue(key, std::get<0>(bond_key));
            key.push_back(',');
            AppendValue(key, std::get<1>(bond_key));
            key.push_back(',');
            AppendValue(key, std::get<2>(bond_key));
            key.push_back(',');
            AppendValue(key, std::get<3>(bond_key));
            key.push_back(';');
        }
        return key;
    }

    ReusableForceField &ThreadLocalForceField(
        const std::string &force_field)
    {
        thread_local auto *force_fields = new std::unordered_map<std::string, ReusableForceField>();
        ReusableForceField &entry = (*force_fields)[force_field];
        if (entry.instance == nullptr)
        {
            entry.instance = MakeForceFieldInstance(force_field);
            if (entry.instance != nullptr)
            {
                entry.instance->SetLogLevel(OBFF_LOGLVL_NONE);
            }
        }
        return entry;
    }

    void PrepareReusableForceFieldForSetup(
        ReusableForceField &entry,
        const std::string &exact_setup_key,
        const std::string &openbabel_setup_key)
    {
        if (entry.instance == nullptr ||
            entry.last_exact_setup_key.empty() ||
            entry.last_exact_setup_key == exact_setup_key ||
            entry.last_openbabel_setup_key != openbabel_setup_key)
        {
            return;
        }

        // OpenBabel::IsSetupNeeded ignores some graph and charge details. If those
        // details changed without changing OpenBabel's coarse signature, reset the
        // reusable instance with an empty molecule so the real Setup rebuilds terms.
        OpenBabel::OBMol empty_mol;
        entry.instance->Setup(empty_mol);
        entry.last_exact_setup_key.clear();
        entry.last_openbabel_setup_key.clear();
    }

    void AppendSeparator(std::string &out)
    {
        out.push_back('|');
    }

    void AppendValue(std::string &out, int value)
    {
        out += std::to_string(value);
    }

    void AppendValue(std::string &out, std::int64_t value)
    {
        out += std::to_string(value);
    }

    void AppendValue(std::string &out, bool value)
    {
        out.push_back(value ? '1' : '0');
    }

    std::string BuildForceFieldEvaluationCacheKey(const OpenBabel::OBMol &mol)
    {
        std::string key = molgr::scoring::BuildScoreKey(mol);
        key += "|uff";
        return key;
    }

    molgr::utils::StringLruCache<molgr::scoring::ForceFieldEvaluation> &
    ForceFieldEvaluationCache()
    {
        static molgr::utils::StringLruCache<molgr::scoring::ForceFieldEvaluation> cache(
            kDefaultForceFieldCacheMaxSize);
        return cache;
    }

    molgr::scoring::ForceFieldEvaluation EvaluateForceFieldUncached(
        const OpenBabel::OBMol &mol,
        const molgr::config::MolGRConfig &config)
    {
        auto *timing_reducer = molgr::pipeline::perf::GetActiveRunTimingReducer();
        const auto total_started = Clock::now();
        const bool contains_metals = molgr::scoring::ContainsMetalAtoms(mol);

        const auto prepare_started = Clock::now();
        OpenBabel::OBMol working_mol(mol);
        working_mol.SetAromaticPerceived(false);
        const int atom_count = static_cast<int>(working_mol.NumAtoms());
        const int heavy_atom_count = CountHeavyAtoms(working_mol);
        if (timing_reducer != nullptr)
        {
            timing_reducer->AddForceFieldPrepareMs(
                std::chrono::duration<double, std::milli>(Clock::now() - prepare_started).count());
        }

        const auto setup_key_started = Clock::now();
        const std::string exact_setup_key = BuildExactForceFieldSetupKey(working_mol);
        const std::string openbabel_setup_key = BuildOpenBabelSetupKey(working_mol);
        if (timing_reducer != nullptr)
        {
            timing_reducer->AddForceFieldSetupKeyMs(
                std::chrono::duration<double, std::milli>(Clock::now() - setup_key_started).count());
        }

        const std::string candidate_force_field_name = "uff";
        ReusableForceField &reusable_force_field = ThreadLocalForceField(candidate_force_field_name);
        PrepareReusableForceFieldForSetup(
            reusable_force_field,
            exact_setup_key,
            openbabel_setup_key);

        OpenBabel::MolgrForceFieldUFF *force_field_ptr = reusable_force_field.instance.get();
        if (force_field_ptr == nullptr)
        {
            throw std::runtime_error("Could not evaluate fixed UFF force-field energy: unavailable");
        }
        force_field_ptr->SetLogLevel(OBFF_LOGLVL_NONE);
        force_field_ptr->ConfigureAtomTypingCache(
            config.cpp_backend.enable_uff_atom_typing_cache,
            exact_setup_key);
        const auto setup_started = Clock::now();
        const bool setup_ok = force_field_ptr->Setup(working_mol);
        if (timing_reducer != nullptr)
        {
            timing_reducer->AddForceFieldSetupMs(
                std::chrono::duration<double, std::milli>(Clock::now() - setup_started).count());
        }
        if (reusable_force_field.instance != nullptr)
        {
            reusable_force_field.last_exact_setup_key = exact_setup_key;
            reusable_force_field.last_openbabel_setup_key = openbabel_setup_key;
        }
        if (!setup_ok)
        {
            throw std::runtime_error("Could not evaluate fixed UFF force-field energy: setup_failed");
        }

        const auto energy_started = Clock::now();
        const double raw_energy = force_field_ptr->Energy();
        if (timing_reducer != nullptr)
        {
            timing_reducer->AddForceFieldEnergyMs(
                std::chrono::duration<double, std::milli>(Clock::now() - energy_started).count());
            timing_reducer->AddForceFieldCalls(1.0);
            timing_reducer->AddForceFieldTotalMs(
                std::chrono::duration<double, std::milli>(Clock::now() - total_started).count());
        }
        const std::string raw_unit = force_field_ptr->GetUnit();
        return molgr::scoring::ForceFieldEvaluation{
            raw_energy,
            raw_unit,
            ForceFieldEnergyToKjMol(raw_energy, raw_unit),
            atom_count,
            heavy_atom_count,
            contains_metals,
        };
    }
}

namespace molgr
{
    namespace scoring
    {
        std::string BuildScoreKey(const OpenBabel::OBMol &mol)
        {
            OpenBabel::OBMol &mutable_mol = MutableMol(mol);
            std::string key;
            key.reserve(
                32 +
                static_cast<std::size_t>(mutable_mol.NumAtoms()) * 56 +
                static_cast<std::size_t>(mutable_mol.NumBonds()) * 24);

            key += "A";
            AppendValue(key, static_cast<int>(mutable_mol.NumAtoms()));
            key.push_back(':');
            FOR_ATOMS_OF_MOL(atom_iter, mutable_mol)
            {
                const OpenBabel::OBAtom &atom = *atom_iter;
                AppendValue(key, static_cast<int>(atom.GetAtomicNum()));
                key.push_back(',');
                AppendValue(key, atom.GetFormalCharge());
                key.push_back(',');
                AppendValue(key, atom.GetSpinMultiplicity());
                key.push_back(',');
                AppendValue(key, QuantizedCoordinate(atom.GetX()));
                key.push_back(',');
                AppendValue(key, QuantizedCoordinate(atom.GetY()));
                key.push_back(',');
                AppendValue(key, QuantizedCoordinate(atom.GetZ()));
                key.push_back(',');
                AppendValue(key, atom.IsAromatic());
                key.push_back(';');
            }

            std::vector<std::tuple<int, int, int, bool>> bond_keys;
            bond_keys.reserve(static_cast<std::size_t>(mutable_mol.NumBonds()));
            FOR_BONDS_OF_MOL(bond_iter, mutable_mol)
            {
                int begin_idx = bond_iter->GetBeginAtom()->GetIdx();
                int end_idx = bond_iter->GetEndAtom()->GetIdx();
                if (begin_idx > end_idx)
                {
                    std::swap(begin_idx, end_idx);
                }
                bond_keys.emplace_back(
                    begin_idx,
                    end_idx,
                    bond_iter->GetBondOrder(),
                    bond_iter->IsAromatic());
            }
            std::sort(bond_keys.begin(), bond_keys.end());

            AppendSeparator(key);
            key += "B";
            AppendValue(key, static_cast<int>(bond_keys.size()));
            key.push_back(':');
            for (const auto &bond_key : bond_keys)
            {
                AppendValue(key, std::get<0>(bond_key));
                key.push_back(',');
                AppendValue(key, std::get<1>(bond_key));
                key.push_back(',');
                AppendValue(key, std::get<2>(bond_key));
                key.push_back(',');
                AppendValue(key, std::get<3>(bond_key));
                key.push_back(';');
            }
            return key;
        }

        std::string BuildMetalStateKey(const std::vector<molgr::MetalAtomPosition> &metal_states)
        {
            std::ostringstream oss;
            for (const auto &metal_state : metal_states)
            {
                oss << metal_state.idx << ','
                    << metal_state.element_idx << ','
                    << metal_state.valence << ','
                    << metal_state.radical_num << ','
                    << QuantizedCoordinate(metal_state.position_x) << ','
                    << QuantizedCoordinate(metal_state.position_y) << ','
                    << QuantizedCoordinate(metal_state.position_z) << ';';
            }
            return oss.str();
        }

        bool ContainsMetalAtoms(const OpenBabel::OBMol &mol)
        {
            OpenBabel::OBMol &mutable_mol = MutableMol(mol);
            FOR_ATOMS_OF_MOL(atom_iter, mutable_mol)
            {
                if (atom_iter->IsMetal())
                {
                    return true;
                }
            }
            return false;
        }

        OpenBabel::OBMol StripMetalAtoms(const OpenBabel::OBMol &mol)
        {
            OpenBabel::OBMol stripped(mol);
            std::vector<unsigned int> metal_indices;
            FOR_ATOMS_OF_MOL(atom_iter, stripped)
            {
                if (atom_iter->IsMetal())
                {
                    metal_indices.push_back(atom_iter->GetIdx());
                }
            }
            if (metal_indices.empty())
            {
                return stripped;
            }

            stripped.BeginModify();
            for (auto it = metal_indices.rbegin(); it != metal_indices.rend(); ++it)
            {
                OpenBabel::OBAtom *atom = stripped.GetAtom(*it);
                if (atom != nullptr)
                {
                    stripped.DeleteAtom(atom);
                }
            }
            stripped.EndModify();
            return stripped;
        }

        ForceFieldEvaluation EvaluateForceField(
            const OpenBabel::OBMol &mol,
            const molgr::config::MolGRConfig &config)
        {
            auto *timing_reducer = molgr::pipeline::perf::GetActiveRunTimingReducer();
            const auto cache_key_started = Clock::now();
            const std::string cache_key =
                BuildForceFieldEvaluationCacheKey(mol);
            if (timing_reducer != nullptr)
            {
                timing_reducer->AddForceFieldCacheKeyMs(
                    std::chrono::duration<double, std::milli>(Clock::now() - cache_key_started).count());
            }
            ForceFieldEvaluation cached;
            if (ForceFieldEvaluationCache().Get(cache_key, cached))
            {
                return cached;
            }
            const ForceFieldEvaluation evaluation =
                EvaluateForceFieldUncached(mol, config);
            ForceFieldEvaluationCache().Put(cache_key, evaluation);
            return evaluation;
        }

        std::tuple<std::size_t, std::size_t, std::size_t> ForceFieldEvaluationCacheInfo()
        {
            return ForceFieldEvaluationCache().Info();
        }

        void ForceFieldEvaluationCacheClear()
        {
            ForceFieldEvaluationCache().Clear();
        }

        ForceFieldEvaluation OrganicForceFieldEvaluation(
            const OpenBabel::OBMol &mol,
            const molgr::config::MolGRConfig &config)
        {
            if (ContainsMetalAtoms(mol))
            {
                throw std::runtime_error(
                    "Organic force-field evaluation only supports metal-free molecules.");
            }
            return EvaluateForceField(mol, config);
        }

        ForceFieldEvaluation CombinedForceFieldEvaluation(
            const OpenBabel::OBMol &mol,
            const molgr::config::MolGRConfig &config)
        {
            return EvaluateForceField(mol, config);
        }

        ForceFieldEvaluation SelectionForceFieldEvaluation(
            const OpenBabel::OBMol &mol,
            const molgr::config::MolGRConfig &config)
        {
            if (ContainsMetalAtoms(mol))
            {
                const OpenBabel::OBMol stripped = StripMetalAtoms(mol);
                return EvaluateForceField(stripped, config);
            }
            return EvaluateForceField(mol, config);
        }

        double OrganicForceFieldEnergy(
            const OpenBabel::OBMol &mol,
            const molgr::config::MolGRConfig &config)
        {
            return OrganicForceFieldEvaluation(mol, config).energy_kj_mol;
        }

        double CombinedForceFieldEnergy(
            const OpenBabel::OBMol &mol,
            const molgr::config::MolGRConfig &config)
        {
            return CombinedForceFieldEvaluation(mol, config).energy_kj_mol;
        }

        double SelectionForceFieldEnergy(
            const OpenBabel::OBMol &mol,
            const molgr::config::MolGRConfig &config)
        {
            return SelectionForceFieldEvaluation(mol, config).energy_kj_mol;
        }
    }
}
