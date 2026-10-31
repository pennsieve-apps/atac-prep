"""
Container entry point: fixes the thread count from THREADS (or the cgroup CPU quota) before
numpy, polars or SnapATAC2 are imported, then runs main() once. The __main__ guard matters:
SnapATAC2 calls MACS3 and some metrics in spawned worker processes, which import this module again.
"""
import math
import os
import sys

THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS",
               "NUMEXPR_NUM_THREADS", "POLARS_MAX_THREADS", "RAYON_NUM_THREADS")


def cgroup_cpus() -> int | None:
    """CPUs from the cgroup quota (v2 cpu.max, else v1 cfs), None when unlimited or unknown."""
    try:
        quota, period = open("/sys/fs/cgroup/cpu.max").read().split()[:2]
        if quota != "max":
            return max(1, math.ceil(int(quota) / int(period)))
        return None
    except (OSError, ValueError):
        pass
    try:
        quota = int(open("/sys/fs/cgroup/cpu/cpu.cfs_quota_us").read())
        period = int(open("/sys/fs/cgroup/cpu/cpu.cfs_period_us").read())
        if quota > 0:
            return max(1, math.ceil(quota / period))
    except (OSError, ValueError):
        pass
    return None


def set_threads() -> int:
    raw = os.environ.get("THREADS", "").strip()
    n = int(raw) if raw.isdigit() else 0
    if n <= 0:
        n = cgroup_cpus() or 2
    for var in THREAD_VARS:
        os.environ[var] = str(n)
    os.environ.setdefault("TQDM_DISABLE", "1")       # SnapATAC2's progress bars would flood the log
    return n


if __name__ == "__main__":
    set_threads()
    from processor.main import main

    sys.exit(main())
