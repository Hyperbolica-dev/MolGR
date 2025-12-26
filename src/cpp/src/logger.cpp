/*
 * @Author: TMJ
 * @Date: 2025-12-26 15:50:56
 * @LastEditors: TMJ
 * @LastEditTime: 2025-12-26 15:51:05
 * @Description: 请填写简介
 */
/**
 * @file logger.cpp
 * @brief Implementation of logging control.
 */

#include "molgr/logger.h"

namespace molgr
{

    // 默认日志级别：WARN (生产环境通常只看警告和错误)
    LogLevel g_current_log_level = LogLevel::WARN;

    void SetLogLevel(LogLevel level)
    {
        g_current_log_level = level;
    }

} // namespace molgr