/**
 * @file logger.cpp
 * @brief Implementation of logging control.
 */

#include "molgr/utils/logger.h"
#include "molgr/process_guard.h"

#include <iostream>
#include <mutex>
#include <set>
#include <string>

namespace molgr
{

    // 默认日志级别：WARN (生产环境通常只看警告和错误)
    LogLevel g_current_log_level = LogLevel::WARN;

    void SetLogLevel(LogLevel level)
    {
        g_current_log_level = level;
    }

    void WarnUnsafeOpenBabelUse(const char *api, const char *detail)
    {
        EnsureCurrentProcess(api);
        static std::mutex warning_mutex;
        static std::set<std::string> emitted;
        const std::string key = api == nullptr ? "unknown" : api;
        if (g_current_log_level > LogLevel::WARN)
        {
            return;
        }
        std::lock_guard<std::mutex> guard(warning_mutex);
        if (!emitted.insert(key).second)
        {
            return;
        }
        std::cerr << "[WARN] [UNSAFE_OPENBABEL] " << key << ": "
                  << (detail == nullptr ? "raw Open Babel object usage requires thread confinement and manual lifetime management."
                                         : detail)
                  << std::endl;
    }

} // namespace molgr
