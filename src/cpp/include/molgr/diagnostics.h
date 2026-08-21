#pragma once

#include <map>
#include <string>

namespace molgr
{
    namespace diagnostics
    {
        struct ReconstructionDiagnostics
        {
            std::string code;
            std::string stage;
            std::string message;
            std::map<std::string, long long> counts;
            std::map<std::string, std::string> details;

            void Reset()
            {
                code.clear();
                stage.clear();
                message.clear();
                counts.clear();
                details.clear();
            }

            void Count(const std::string &name, long long amount = 1)
            {
                counts[name] += amount;
            }

            void Fail(
                const std::string &code_,
                const std::string &stage_,
                const std::string &message_)
            {
                code = code_;
                stage = stage_;
                message = message_;
            }
        };
    }
}
