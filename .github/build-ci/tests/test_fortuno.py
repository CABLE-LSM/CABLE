"""Pytest wrapper around CABLE's Fortuno unit tests.

Run by build-ci's ``Test`` job (``enable-pytest``), which activates the Spack
environment built from this repository's manifest and then invokes pytest. For
each installed ``cable`` spec that supports tests, this module configures the
source tree with ``-DCABLE_TESTS=ON``, builds the ``cable-tests`` target and
runs the ctest suites declared in CMakeLists.txt.

The rebuild is needed because ``spack install`` discards its stage, leaving no
build directory to run ctest in. The Fortuno version is pinned in
``.github/build-ci/config/packages.yaml``.
"""

from __future__ import annotations

import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pytest

TIMEOUT_S = 3600


@dataclass(frozen=True)
class CableSpec:
    """An installed `cable` spec discovered in the active Spack environment."""

    hash: str
    compiler: str
    variants: str

    @property
    def mpi(self) -> bool:
        return "+mpi" in self.variants

    @property
    def pio(self) -> bool:
        return "+pio" in self.variants

    @property
    def library(self) -> str:
        match = re.search(r"library=(\S+)", self.variants)
        return match.group(1) if match else "none"

    @property
    def build_type(self) -> str:
        match = re.search(r"build_type=(\S+)", self.variants)
        return match.group(1) if match else "Release"

    @property
    def id(self) -> str:
        """A pytest id usable in `-k` expressions, which reject `@`, `+` and `~`."""
        variants = f"{'' if self.mpi else 'no'}mpi_{'' if self.pio else 'no'}pio"
        return re.sub(r"[^0-9a-zA-Z]+", "_", f"{variants}_{self.compiler.split('/')[0]}")


def _spack(*args: str) -> str:
    """Run a spack command against the active environment and return stdout."""
    result = subprocess.run(["spack", *args], capture_output=True, text=True, timeout=600)
    if result.returncode != 0:
        raise RuntimeError(
            f"`spack {' '.join(args)}` failed ({result.returncode}):\n{result.stderr}"
        )
    return result.stdout


def discover_cable_specs() -> list[CableSpec]:
    try:
        output = _spack("find", "--format", "{hash}|{compiler}|{variants}", "cable")
    except (RuntimeError, OSError, subprocess.SubprocessError):
        return []
    return [
        CableSpec(*(field.strip() for field in line.split("|")))
        for line in output.splitlines()
        if line.count("|") == 2
    ]


def cable_source_dir() -> Path:
    """Locate the CABLE checkout holding CMakeLists.txt.

    build-ci's manifest checkout is sparse and may lack the sources; the full
    caller checkout is at CALLER_REPOSITORY_PATH.
    """
    candidates = [
        os.environ.get("CABLE_SOURCE_DIR"),
        os.environ.get("CALLER_REPOSITORY_PATH"),
        str(Path(__file__).resolve().parents[3]),
    ]
    for candidate in candidates:
        if candidate and (Path(candidate) / "CMakeLists.txt").is_file():
            return Path(candidate)
    raise RuntimeError(
        "could not locate the CABLE source tree; set CABLE_SOURCE_DIR to the checkout "
        "containing CMakeLists.txt"
    )


def build_environment(spec: CableSpec) -> dict[str, str]:
    """Capture the Spack build environment for a spec.

    `spack build-env` is --clean by default, so variables exported by the CI job
    do not survive it. Capturing the environment and injecting into it keeps
    those settings intact.
    """
    raw = _spack("build-env", f"/{spec.hash}", "--", "env", "-0")
    env = dict(entry.split("=", 1) for entry in raw.split("\0") if "=" in entry)
    # The CI container runs as root; Open MPI refuses to launch without this.
    env.setdefault("OMPI_ALLOW_RUN_AS_ROOT", "1")
    env.setdefault("OMPI_ALLOW_RUN_AS_ROOT_CONFIRM", "1")
    return env


def fortuno_cmake_dir(spec: CableSpec) -> Path | None:
    """Find the installed Fortuno (pinned in packages.yaml) for the spec's compiler.

    Returns None when nothing matches, leaving CMake to fall back to fetching the
    commit pinned in CMakeLists.txt.
    """
    try:
        output = _spack("find", "--format", "{compiler}|{prefix}", "fortuno")
    except (RuntimeError, OSError, subprocess.SubprocessError):
        return None
    for line in output.splitlines():
        compiler, _, prefix = line.partition("|")
        if compiler.strip() != spec.compiler:
            continue
        for libdir in ("lib64", "lib"):
            candidate = Path(prefix.strip()) / libdir / "cmake" / "Fortuno"
            if candidate.is_dir():
                return candidate
    return None


def _tail(result: subprocess.CompletedProcess, limit: int = 4000) -> str:
    return f"{result.stdout.strip()[-limit:]}\n{result.stderr.strip()[-limit:]}"


CABLE_SPECS = discover_cable_specs()

# CMakeLists.txt ignores CABLE_TESTS for library builds, so manifests that only
# build `library=<application>` specs have no Fortuno tests to run.
TESTABLE_SPECS = [spec for spec in CABLE_SPECS if spec.library == "none"]


def test_cable_specs_discovered():
    """Guard against silently skipping everything when discovery fails."""
    assert CABLE_SPECS, (
        "no installed `cable` specs in the active Spack environment; "
        "is the environment activated and installed?"
    )
    if not TESTABLE_SPECS:
        pytest.skip(
            "manifest builds no `library=none` specs; Fortuno tests do not apply to "
            "library builds"
        )


@pytest.mark.parametrize("spec", TESTABLE_SPECS, ids=lambda spec: spec.id)
def test_fortuno_suites_pass(spec: CableSpec, tmp_path: Path):
    """Configure, build and run the Fortuno suites for one installed spec."""
    env = build_environment(spec)
    source_dir = cable_source_dir()

    configure_args = [
        "cmake", "-S", str(source_dir), "-B", str(tmp_path),
        "-DCABLE_TESTS=ON",
        "-DCABLE_LIBRARY=OFF",
        f"-DCABLE_MPI={'ON' if spec.mpi else 'OFF'}",
        f"-DCMAKE_BUILD_TYPE={spec.build_type}",
        # Mirror the package recipe so the build matches the installed spec.
        f"-DCMAKE_{'REQUIRE' if spec.pio else 'DISABLE'}_FIND_PACKAGE_PIO=ON",
    ]
    fortuno_dir = fortuno_cmake_dir(spec)
    if fortuno_dir is not None:
        configure_args.append(f"-DFortuno_DIR={fortuno_dir}")
    configure = subprocess.run(
        configure_args, env=env, capture_output=True, text=True, timeout=TIMEOUT_S
    )
    assert configure.returncode == 0, f"cmake configure failed:\n{_tail(configure)}"

    build = subprocess.run(
        ["cmake", "--build", str(tmp_path), "--target", "cable-tests",
         "-j", str(os.cpu_count() or 1)],
        env=env, capture_output=True, text=True, timeout=TIMEOUT_S,
    )
    assert build.returncode == 0, f"building cable-tests failed:\n{_tail(build)}"

    ctest = subprocess.run(
        # ctest exits 0 when no test matches, which would otherwise pass silently.
        ["ctest", "--test-dir", str(tmp_path), "--output-on-failure", "--no-tests=error"],
        env=env, capture_output=True, text=True, timeout=TIMEOUT_S,
    )
    assert ctest.returncode == 0, f"ctest exited {ctest.returncode}:\n{_tail(ctest)}"
