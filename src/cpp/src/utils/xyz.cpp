#include "molgr/utils/xyz.h"

#include "molgr/vendor/openbabel_threading.h"

#include <openbabel/atom.h>
#include <openbabel/elements.h>

#include <algorithm>
#include <cctype>
#include <locale>
#include <sstream>
#include <string>
#include <vector>

namespace
{
    struct ParsedXyzAtom
    {
        int atomic_num = 0;
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
    };

    bool ParseAtomCountLine(const std::string &line, std::size_t *atom_count)
    {
        std::istringstream stream(line);
        stream.imbue(std::locale::classic());

        std::size_t count = 0;
        stream >> count;
        if (!stream)
        {
            return false;
        }

        std::string trailing;
        if (stream >> trailing)
        {
            return false;
        }

        *atom_count = count;
        return count > 0;
    }

    bool ParseAtomicNumber(std::string token, int *atomic_num)
    {
        if (token.empty())
        {
            return false;
        }

        const bool numeric_token = std::all_of(
            token.begin(),
            token.end(),
            [](unsigned char value)
            {
                return std::isdigit(value) != 0;
            });
        if (numeric_token)
        {
            std::istringstream stream(token);
            int parsed_atomic_num = 0;
            stream.imbue(std::locale::classic());
            stream >> parsed_atomic_num;
            if (!stream || parsed_atomic_num <= 0)
            {
                return false;
            }
            *atomic_num = parsed_atomic_num;
            return true;
        }

        token[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(token[0])));
        for (std::size_t idx = 1; idx < token.size(); ++idx)
        {
            token[idx] = static_cast<char>(std::tolower(static_cast<unsigned char>(token[idx])));
        }

        const int parsed_atomic_num = OpenBabel::OBElements::GetAtomicNum(token.c_str());
        if (parsed_atomic_num <= 0)
        {
            return false;
        }
        *atomic_num = parsed_atomic_num;
        return true;
    }

    bool ParseAtomLine(const std::string &line, ParsedXyzAtom *parsed_atom)
    {
        std::istringstream stream(line);
        stream.imbue(std::locale::classic());

        std::string atom_token;
        stream >> atom_token >> parsed_atom->x >> parsed_atom->y >> parsed_atom->z;
        if (!stream)
        {
            return false;
        }

        return ParseAtomicNumber(atom_token, &parsed_atom->atomic_num);
    }
}

namespace molgr
{
    namespace utils
    {
        bool ReadXyzBlockToMol(const std::string &xyz_block, OpenBabel::OBMol *mol)
        {
            if (mol == nullptr)
            {
                return false;
            }

            std::istringstream input(xyz_block);
            input.imbue(std::locale::classic());

            std::string atom_count_line;
            if (!std::getline(input, atom_count_line))
            {
                return false;
            }

            std::size_t atom_count = 0;
            if (!ParseAtomCountLine(atom_count_line, &atom_count))
            {
                return false;
            }

            std::string title_line;
            if (!std::getline(input, title_line))
            {
                return false;
            }

            std::vector<ParsedXyzAtom> atoms;
            atoms.reserve(atom_count);
            std::string atom_line;
            for (std::size_t idx = 0; idx < atom_count; ++idx)
            {
                if (!std::getline(input, atom_line))
                {
                    return false;
                }

                ParsedXyzAtom parsed_atom;
                if (!ParseAtomLine(atom_line, &parsed_atom))
                {
                    return false;
                }
                atoms.push_back(parsed_atom);
            }

            mol->Clear();
            mol->BeginModify();
            mol->ReserveAtoms(static_cast<unsigned int>(atom_count));
            mol->SetDimension(3);
            mol->SetTitle(title_line.c_str());
            for (const ParsedXyzAtom &parsed_atom : atoms)
            {
                OpenBabel::OBAtom *atom = mol->NewAtom();
                atom->SetAtomicNum(parsed_atom.atomic_num);
                atom->SetVector(parsed_atom.x, parsed_atom.y, parsed_atom.z);
            }
            molgr::vendor::openbabel_threading::ConnectTheDotsAndPerceiveBondOrders(*mol);
            mol->EndModify();
            return true;
        }
    }
}
