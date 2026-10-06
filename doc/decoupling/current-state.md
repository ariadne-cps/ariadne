# Current decoupling state

Status date: 2026-10-06  
State described here: post-merge architecture after the Algebra extraction,
with an active Windows portability validation on `fix-windows`

The CSV, SVG and include evidence in this directory remain the historical
2026-09-29 baseline. This file describes the current architecture.

## Repository chain

The validated low-level chain is:

```text
utility
  -> foundation
       -> numeric
            -> interval
                 -> algebra
                      -> function
```

Utility, Foundation, Numeric, Interval and Algebra are standalone repositories.
Function is still owned by Ariadne and is the next extraction target.

Ariadne itself consumes the standalone lower stack through its direct
dependencies rather than carrying local copies of those components.

## Direct dependency discipline

A repository consumes only its immediate lower repository. It must not encode
knowledge of nested dependency paths when the direct dependency can propagate
the requirement.

Shared CMake infrastructure follows the same principle. Repositories that use
Configuration carry `submodules/configuration` directly and bootstrap from:

```text
submodules/configuration/cmake/ProjectOptions.cmake
```

When Configuration is also present through a dependency,
`require_same_dependency_commit` checks that all visible gitlinks point to the
same commit. This avoids nested bootstrap paths such as
`submodules/.../submodules/configuration` while still detecting incompatible
pins.

## Python binding layering

Python bindings mirror the C++ repository chain:

```text
pyariadne-foundation
  -> pyariadne-numeric
       -> pyariadne-interval
            -> pyariadne-algebra
                 -> pyariadne
```

Each layer owns and compiles its bindings and exports an aggregate interface for
the next layer. A consumer sees only the Python binding surface of its direct
lower dependency.

Ariadne therefore consumes `pyariadne-algebra`. Its
`python/bindings/utilities.hpp` builds on `algebra-utilities.hpp` instead of
redefining lower-layer helpers. Algebra-owned Python representations and class
registrations are not duplicated in Function/Ariadne bindings.

## Installation

Standalone components register their public headers with
`ariadne_register_public_headers`. The installation helper
`ariadne_install_dependency_bundle` follows registered
`INTERFACE_LINK_LIBRARIES` recursively and installs the reachable public
header directories once.

Ariadne uses dependency bundles for both lower branches:

```cmake
ariadne_install_dependency_bundle(
    TARGET algebra
    DESTINATION include/ariadne
)

ariadne_install_dependency_bundle(
    TARGET threading
    DESTINATION include/ariadne
)
```

This replaces manual installation lists for Algebra, Interval, Numeric,
Foundation, Utility, Threading and Logging.

The remaining CMake-package-specific cleanup is separate from header bundling:
Ariadne still installs its own package configuration and currently obtains the
GMP/MPFR find modules from Numeric. That source-tree knowledge should eventually
be removed from the package-export path.

## Ariadne integration shape

The former `ariadne-core` and `ariadne-kernel` aggregates are gone. Ariadne
builds one aggregate `ariadne` library from Ariadne-owned components plus the
standalone lower-layer objects.

Ariadne no longer contains:

- a local `source/algebra` implementation;
- Function-owned copies of `Polynomial` or Chebyshev polynomial classes;
- Function-owned Sweeper implementations;
- duplicate Algebra Python bindings;
- Algebra-owned tests, demonstrations or tutorials.

In-tree Ariadne tests and benchmarks link the aggregate `ariadne` target.

## Algebra boundary

Algebra depends directly on Interval and has no dependency on Function. It owns
the algebraic representations and infrastructure moved out of Ariadne,
including:

- `Polynomial`;
- `UnivariateChebyshevPolynomial` and
  `MultivariateChebyshevPolynomial`;
- Sweeper implementation code;
- the corresponding C++ tests, Python bindings and standalone tutorials.

Standalone Algebra CI and Ariadne integration CI are green in the state assumed
by this document.

## Active Windows portability validation

A cross-repository `fix-windows` effort is validating that the decoupled stack
can be built on Windows without changing the C++ architecture for that
platform.

The current build-policy streamlining starts in Configuration and is propagated
bottom-up. The active chain is:

```text
configuration
  -> utility
       -> foundation
            -> numeric
```

and independently:

```text
configuration
  -> utility
       -> logging
            -> threading
```

Each repository modified during this decoupling effort uses a `fix-windows`
branch and an open pull request so its CI can be inspected before the revision
is propagated upward.

The current pull requests for this build-policy pass are:

- Configuration PR #2;
- Utility PR #6;
- Foundation PR #6;
- Logging PR #11;
- Threading PR #14;
- Numeric PR #24.

The first propagation pass reached green CI in Utility, Foundation, Logging,
Threading, Numeric and Interval. Algebra was then rebuilt with its previous MSVC
warning suppressions removed and exposed two real diagnostic families:

- Numeric-owned C4244 narrowing warnings triggered by consumer-side template
  instantiation;
- Algebra-owned C4661 warnings caused by broad explicit class instantiation.

The Numeric C4244 sites were fixed at source in
`413a2a00ae61e96d6734196f7b0faac139d5672e`. The next pass will restart from
the latest Numeric main branch,
`bbb57e946f86d6d0134dcb90b81a86102ccae816`, which contains that fix.

Interval and Algebra must therefore be repointed and revalidated from that
Numeric main baseline before propagation continues to Kernel. Algebra C4661
remains open work.

This Windows work also exposed that warning suppression can make integration
progress appear better than it is. The portability rules are therefore
explicit:

1. **No Windows-specific C++ implementation paths.** The same C++ design and
   semantics must be used across supported platforms.
2. **No warning suppression as a fix.** MSVC `/wd*` flags and equivalent
   mechanisms are not accepted. Warnings promoted by the common policy must be
   understood and fixed.
3. **Centralize compiler/build policy.** Common MSVC options and aggregate
   library policy belong in Configuration, not scattered through individual
   repositories.
4. Platform-specific CMake logic is acceptable only for genuine build-system
   or packaging constraints; it must not hide semantic differences in the C++
   implementation.

The current Configuration experiment centralizes the common MSVC policy,
including `/bigobj`, and makes the aggregate project library static on Windows
while remaining shared on Unix-like systems. Repository-local duplicates of
that policy are being removed as the revision is propagated.

## Next boundary: Function

Function should become the repository immediately above Algebra. The remaining
upward dependencies to resolve are principally:

- **Geometry:** Box/domain primitives are mixed with Point/SetInterface and
  higher-level set functionality.
- **Symbolic:** generic symbolic machinery and concrete
  Expression/Function bridges still cross the boundary.
- **Infrastructure cleanup:** verify that Function has no semantic dependency on
  Logging or Threading; the known direct Logging include is cleanup rather than
  part of the intended package boundary.

The Box/domain issue is structural: Function needs a primitive domain
representation, not the higher-level Geometry package as a whole.

## Validation rule

A repository boundary is ready only when it:

1. compiles its public headers in isolation;
2. builds from only its declared lower dependencies;
3. runs focused C++ tests;
4. builds Python bindings using only the aggregate Python interface of its
   direct dependency;
5. installs cleanly and is consumable by an external project;
6. passes the Ariadne integration build after its revision is propagated.

Successful aggregate Ariadne linkage alone is not proof of a valid standalone
boundary.
