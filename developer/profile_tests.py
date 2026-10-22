"""Run pytest under cProfile and print the hottest cumulative call sites.

This helper is intended to be executed inside the tox test environment, for
example with ``tox -e profile-tests``.  It keeps profiling setup out of the
production package and accepts normal pytest arguments.
"""

from __future__ import annotations

import cProfile
import pstats
import runpy
from pathlib import Path


def main() -> int:
    """Profile pytest, write a profile artifact, and display top call sites."""
    profile_path = Path(".tox-profile-tests.prof")
    profiler = cProfile.Profile()
    exit_code = 0

    try:
        profiler.enable()
        runpy.run_module("pytest", run_name="__main__")
    except SystemExit as exc:
        exit_code = int(exc.code) if isinstance(exc.code, int) else 1
    finally:
        profiler.disable()
        profiler.dump_stats(profile_path)

    print(f"\nWrote cProfile data to {profile_path}")
    print("\nTop cumulative call sites:")
    pstats.Stats(profiler).strip_dirs().sort_stats("cumulative").print_stats(40)
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
