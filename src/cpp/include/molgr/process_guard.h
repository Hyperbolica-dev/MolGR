#pragma once

namespace molgr
{
    // Record the process that initialized the native extension.  Open Babel
    // state must not be reused after a POSIX fork because only the calling
    // thread survives in the child.
    void InitializeProcessGuard();

    bool IsCurrentProcess();

    // Throw when the current process is a forked child of the process that
    // initialized MolGR.  Use a fresh process (spawn/exec) for isolation.
    void EnsureCurrentProcess(const char *api);
}
