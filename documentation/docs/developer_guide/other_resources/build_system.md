For CABLE offline applications, the underlying build system is [CMake](https://cmake.org) based. CMake is widely used in the C, C++ and Fortran communities and there are many great resources out there on the web, in particular:

- [Mastering CMake](https://cmake.org/cmake/help/book/mastering-cmake/index.html)
- [Effective Modern CMake](https://gist.github.com/mbinna/c61dbb39bca0e4fb7d1f73b0d66a4fd1)
- [Creating a CMake project (with Fortran)](https://fortran-lang.org/learn/building_programs/build_tools/#creating-a-cmake-project)
- [An Introduction to Modern CMake](https://cliutils.gitlab.io/modern-cmake/)

## Adding source files

Source files can be added easily by adding to the list of source files in [CMakeLists.txt][CMakeLists.txt]. We recommend keeping the list sorted in alphabetical order when adding new source files.

## Setting compiler flags

The recommended compiler flags for debug and release builds for each compiler are set in the [CMakeLists.txt][CMakeLists.txt] file via the `CABLE_<COMPILER_ID>_Fortran_FLAGS*` variables. These can be modified if required.

???+ warning
    Compiler flags should be guarded by a check on the compiler ID as they are generally compiler specific. See [`CMAKE_<LANG>_COMPILER_ID`](https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html) for a list of supported compiler IDs.

Alternatively, compiler flags can be specified when invoking `cmake` by specifying `-DCMAKE_Fortran_FLAGS=<flags>` when generating the project, however the per configuration default flags may take precedence (see [`CMAKE_<LANG>_FLAGS`](https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_FLAGS.html) for more details).

## Adding external dependencies (libraries)

When linking against external third-party libraries, we strongly recommend not to hard code library paths, include paths, link flags or compiler flags in the CMakeLists.txt file. This ensures that the build system is portable and insulated from changes to library versions and dependencies.

Instead, CMake should query the required flags for a given library from the library itself. To do this, we recommend using either CMake's [`find_package` mechanism](https://cmake.org/cmake/help/book/mastering-cmake/chapter/Finding%20Packages.html) or a tool such as [`pkg-config`](https://en.wikipedia.org/wiki/Pkg-config) (see the [FindPkgConfig](https://cmake.org/cmake/help/latest/module/FindPkgConfig.html) module for using `pkg-config` with CMake.). The dependencies should provide exported targets which can be used via [`target_link_libraries`](https://cmake.org/cmake/help/latest/command/target_link_libraries.html).

For an example, see how the netcdf-fortran dependency was added in [CMakeLists.txt][CMakeLists.txt].

If these approaches are not supported by the external library, please report this as a bug to its maintainers. If the library is an open-source project, consider sending a patch.

???+ warning
    For most cases, CMake's `find_package` should be used in **Config** mode. `find_package` in **Module** mode should only be used if the library is part of the [CMake module distribution](https://cmake.org/cmake/help/latest/manual/cmake-modules.7.html#manual:cmake-modules(7)).

[CMakeLists.txt]: https://github.com/CABLE-LSM/CABLE/blob/main/CMakeLists.txt
