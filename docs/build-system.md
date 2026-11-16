# Build system

## Current setup

StructPhy uses Python 3.12. `uv sync --locked` installs the Python environment and builds three native tools from source on a fresh checkout. The tools are installed as private package files under `structphy/_bin/` in the environment, including for editable installs. StructPhy uses its installed tools regardless of the working directory or `PATH`. Set `STRUCTPHY_TOOL_DIR` to explicitly override one or more binaries. The old `STRUCTPHY_CACHE_DIR` variable still works as an alias; a local `.structphy/` directory is no longer searched automatically.

| Tool | Pinned version | Role |
| --- | --- | --- |
| TM-align (`TMalign`) | 2022-04-12 | Aligns pairs of PDB structures to compute distances |
| FastME (`fastme`) | 2.1.6.4 | Builds a tree from each bootstrap distance matrix |
| PHYLIP `consense` | 3.697 | Combines bootstrap trees into a consensus |

The native sources and their license notices ship with the project, so the compiler can build them without fetching source archives. The Python dependency versions are tracked in `uv.lock`. [Native source inventory](../native/vendor/README.md)

uv manages Python, the environment, the lockfile, and build isolation. Hatchling writes the source distribution and wheels. Its [build hook](../native/build_hook.py) compiles the tools in a temporary directory and declares which files to bundle. The checkout layout keeps the build recipe and upstream sources together in `native/` and the runtime code in `structphy/`. The example datasets stay in the checkout. See [Astral’s build backend guidance](https://docs.astral.sh/uv/concepts/build-backend/).

Source builds require a C/C++ toolchain and `make`. On macOS, install the Xcode Command Line Tools with `xcode-select --install`. On Debian or Ubuntu, install `build-essential`. PhyloDM has its own native extension and can also require Cargo and a system linker on platforms without a compatible PhyloDM wheel.

Local builds target the host platform and architecture. macOS wheels use the build host’s macOS version as the compatibility floor. Local Linux wheels carry a generic Linux architecture tag and should be used on a compatible distribution. The [CI workflow](../.github/workflows/ci.yml) checks source builds and clean wheel installs on macOS arm64 and Linux x86-64.

## Proposed release workflow

Ship **prebuilt wheels by default** and keep publishing a source distribution. A compatible wheel installs the bundled tools without compiling them on the user's machine. Source builds remain available for unsupported platforms that can build the tools, local development, and performance experiments.

1. Choose a license for StructPhy itself. The root `LICENSE` is currently empty; the included notices apply to the upstream tools.
2. Set and verify the minimum macOS version for release wheels. For Linux, build against a `manylinux` baseline, check the bundled executables' shared-library requirements, and test installs on a compatible target. A generic `linux_*` tag alone gives no glibc compatibility guarantee. [Python platform compatibility tags](https://packaging.python.org/en/latest/specifications/platform-compatibility-tags/)
3. Verify the full install on each release target. In particular, check that the locked PhyloDM version has a compatible wheel before describing the install as compiler-free.
4. Publish the tested wheels and source distribution for the same release. Keep `uv sync --locked` as the development workflow.

## Performance

The current build uses optimization flags such as `-O3` for TM-align. It does not explicitly target the build machine's CPU. FastME builds with OpenMP disabled because StructPhy already distributes independent tree jobs across processes. [Build settings](../native/build_hook.py)

Before adding CPU tuning, compare the release build and a CPU-tuned source build on the same machine. Use the same compiler and tool versions, input PDBs, bootstrap seed, and thread count. Measure both TM-align runtime and total PDB-to-tree wall time over repeated runs. Check the distance matrices and final trees for changes as well.

If the tuned build gives a meaningful end-to-end speedup, make it an explicit source-build option. Include the tuning settings in uv's build cache key and record the compiler, flags, and CPU in benchmark results. Keep release wheels on a portable CPU baseline. CPU-specific instructions can make a binary unusable on older CPUs. [GCC CPU targeting documentation](https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html)

## Development checks

```sh
uv sync --locked
uv run --locked ruff check .
uv run --locked ruff format --check .
uv run --locked ty check
uv run --locked pytest -q -m "not packaging"
uv run --locked pytest -q -m packaging
uv build --no-sources
```

The packaging test builds a wheel from the source distribution and installs it with its dependencies in a clean environment. The other tests exercise the editable install. After changing native build inputs, uv invalidates the cached build; to explicitly refresh it, use `uv sync --locked --reinstall-package structphy`.

Use `uv add` and `uv remove` to manage dependencies. Use `--group dev` for development tools and `--group research` for notebook tools. Start JupyterLab with `uv run --locked --group research jupyter lab`. To update Python dependencies deliberately, run `uv lock --upgrade`, review the lockfile diff, and rerun the checks.
