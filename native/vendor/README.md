# Native sources

These pinned sources are build inputs. The source distribution includes all three
archives. The wheel includes the executables, the notices in `licenses/`, and
the FastME source and build hook under `.dist-info/extra_metadata/native/`.

- `TMalign-20220412.cpp.gz` contains the [TM-align 2022-04-12 C++ source](https://zhanggroup.org/TM-align/TMalign.cpp).
  The build replaces the `malloc.h` include with `stdlib.h` on macOS.
- `fastme-2.1.6.4.tar.gz` contains the [FastME 2.1.6.4 release](https://gite.lirmm.fr/atgc/FastME)
  with upstream prebuilt binaries omitted. It retains `COPYING` and `AUTHORS`.
  The build uses FastME's OpenMP disable flag because StructPhy distributes
  independent tree jobs across processes.
- `phylip-consense-3.697.tar.gz` contains the five C source and header files
  used to build `consense` from the [PHYLIP 3.697 distribution](https://phylipweb.github.io/phylip/download/phylip-3.697.tar.gz).
  The BSD 2-Clause notices are in the C sources.
