"""
Container entry point: fixes the thread count from THREADS (or the cgroup CPU quota) before
numpy, polars or SnapATAC2 are imported, then runs main() once. The __main__ guard matters:
SnapATAC2 calls MACS3 and some metrics in spawned worker processes, which import this module again.
"""
import os
import sys

THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS",
               "NUMEXPR_NUM_THREADS", "POLARS_MAX_THREADS", "RAYON_NUM_THREADS")


def set_threads() -> int:
    from processor.safe_inputs import cpu_threads     # standard library only: numpy is not imported yet

    raw = os.environ.get("THREADS", "").strip()
    n = cpu_threads(int(raw) if raw.isdigit() else 0)
    for var in THREAD_VARS:
        os.environ[var] = str(n)
    os.environ.setdefault("TQDM_DISABLE", "1")       # SnapATAC2's progress bars would flood the log
    return n


if __name__ == "__main__":
    set_threads()
    from processor.main import main

    sys.exit(main())
