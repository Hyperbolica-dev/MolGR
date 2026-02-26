#pragma once

#include <map>
#include <set>
#include <string>
#include <vector>

namespace molgr
{
    struct ElementInfo
    {
        std::string symbol;
        std::string name;
        int atomic_number;
        double atomic_mass;
        int num_outer_electrons;
        int default_valence;
    };

    struct LennardJonesParams
    {
        double epsilon;
        double sigma;
        double cutoff;
    };

    struct FDSP
    {
        int f, d, s, p;
    };

    extern const std::map<int, ElementInfo> kNonMetalDict;
    extern const std::set<int> kHeteroatoms;
    extern const std::map<std::string, LennardJonesParams> kLJParams;
    extern const std::map<std::string, FDSP> kMetalFDSP;
    extern const std::vector<std::vector<int>> kDElectronsSpin;
    extern const std::map<std::string, std::vector<int>> kMetalValencePrior;
    extern const std::map<std::string, std::vector<int>> kMetalValenceMinor;

    std::set<int> GetPossibleMetalRadicals(const std::string &metal_symbol, int valence);
}
