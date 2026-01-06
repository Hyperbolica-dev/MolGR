/*
 * @Author: TMJ
 * @Date: 2025-12-26 15:50:14
 * @LastEditors: TMJ
 * @LastEditTime: 2026-01-01 21:08:04
 * @Description: 请填写简介
 */
/**
 * @file logger.h
 * @brief Simple logging system for C++ extension.
 * @namespace molgr
 */

#pragma once

#include <iostream>
#include <string>

namespace OpenBabel
{
    class OBMol;
}
namespace molgr
{

    // 定义日志级别，模仿 Python 的 logging 级别
    enum class LogLevel
    {
        DEBUG = 0,
        INFO = 1,
        WARN = 2,
        ERROR = 3,
        OFF = 4 // 关闭所有日志
    };

    // 全局日志级别变量声明
    extern LogLevel g_current_log_level;

    // 设置日志级别的函数
    void SetLogLevel(LogLevel level);

    /**
     * @brief Log the detailed information of an OBMol.
     * @details Dumps atom counts, bonds, SMILES, and detailed atom states (charge/spin).
     * @param mol The molecule to log.
     * @param comment Optional comment to print before the details.
     */
    void LogOmolInfos(OpenBabel::OBMol &mol, const std::string &comment = "");

} // namespace molgr

// =============================================================================
// 日志宏 (Macros)
// =============================================================================
// 使用 do-while(0) 技巧确保宏在 if-else 语句中安全使用
// 检查 g_current_log_level，如果级别不够，后面的流操作完全不执行

#define LOG_DEBUG(msg)                                            \
    do                                                            \
    {                                                             \
        if (molgr::g_current_log_level <= molgr::LogLevel::DEBUG) \
        {                                                         \
            std::cerr << "[DEBUG] " << msg << std::endl;          \
        }                                                         \
    } while (0)

#define LOG_INFO(msg)                                            \
    do                                                           \
    {                                                            \
        if (molgr::g_current_log_level <= molgr::LogLevel::INFO) \
        {                                                        \
            std::cerr << "[INFO]  " << msg << std::endl;         \
        }                                                        \
    } while (0)

#define LOG_WARN(msg)                                            \
    do                                                           \
    {                                                            \
        if (molgr::g_current_log_level <= molgr::LogLevel::WARN) \
        {                                                        \
            std::cerr << "[WARN]  " << msg << std::endl;         \
        }                                                        \
    } while (0)

#define LOG_ERROR(msg)                                            \
    do                                                            \
    {                                                             \
        if (molgr::g_current_log_level <= molgr::LogLevel::ERROR) \
        {                                                         \
            std::cerr << "[ERROR] " << msg << std::endl;          \
        }                                                         \
    } while (0)
