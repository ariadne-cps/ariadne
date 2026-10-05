# Ariadne decoupling work plan

Status: active  
Base state: low-level chain through Algebra extracted and green  
Last updated: 2026-10-05

## Established chain

```text
utility -> foundation -> numeric -> interval -> algebra -> function
```

Utility, Foundation, Numeric, Interval and Algebra are standalone. The completed
Foundation/Algebra extraction history is recorded in Git history and
`current-state.md`; completed migration checklists are intentionally omitted
from this active plan.

## F1 - Remove accidental Function infrastructure dependencies

- [ ] Remove the unused `logging/logging.hpp` include from
      `function/calculus_base.hpp`.
- [ ] Verify that Function has no semantic Logging or Threading dependency.
- [ ] Make Function target dependencies explicit rather than relying on
      Ariadne-wide linkage.

## F2 - Remove Function -> high-level Geometry

- [ ] Identify the minimal Box/domain value types required by Function.
- [ ] Separate those primitives from Point/SetInterface and higher-level
      Geometry responsibilities.
- [ ] Move measurable-function and multifunction set integrations above the
      Function boundary where appropriate.
- [ ] Ensure Function public headers no longer depend on the high-level
      `geometry/` package.

## F3 - Remove Function -> high-level Symbolic

- [ ] Separate generic symbolic machinery genuinely required by Function from
      high-level Expression/Space/Variable integration.
- [ ] Move concrete Expression/Function conversion bridges to the higher layer.
- [ ] Ensure Function no longer depends on high-level Symbolic implementation.

## F4 - Make Function standalone

- [ ] Define the standalone Function target above Algebra using only intended
      lower dependencies.
- [ ] Register Function public headers and use the dependency-bundle
      installation model.
- [ ] Give Function a direct Configuration dependency and verify matching
      Configuration gitlinks with its direct dependency.
- [ ] Add public-header isolation checks.
- [ ] Add focused C++ tests.
- [ ] Layer Python bindings on the aggregate `pyariadne-algebra` interface
      only.
- [ ] Add standalone installation and external-consumer/tutorial checks.
- [ ] Get Unix, Windows, Coverage and Python CI green.

## F5 - Integrate standalone Function into Ariadne

- [ ] Add the standalone `ariadne-cps/function` repository as Ariadne's direct
      lower Function dependency.
- [ ] Remove Ariadne's local Function implementation and Function-owned tests or
      bindings that move with the repository.
- [ ] Make Ariadne consume the aggregate Function C++ and Python interfaces.
- [ ] Propagate compatible dependency pins.
- [ ] Get the full Ariadne integration CI green.

## W1 - Windows portability and build-policy streamlining

This work is active and must be completed before the Function extraction is
allowed to rely on Windows integration results.

- [x] Centralize common Windows/MSVC build policy in
      `configuration/fix-windows`.
- [x] Propagate that Configuration revision through Utility and Foundation.
- [x] Propagate the same Configuration revision through Logging and Threading.
- [x] Update Numeric to consume the streamlined Configuration/Foundation chain
      and remove local `/bigobj` duplication.
- [ ] Make Numeric Windows CI green. Current state: Unix and Coverage pass;
      Windows fails in the Build step before tests.
- [ ] Only after Numeric is green, propagate bottom-up through Interval and
      Algebra.
- [ ] Remove all remaining MSVC warning suppressions from Algebra and any other
      repository encountered during upward propagation.
- [ ] Propagate the validated Algebra and Threading revisions into Kernel.
- [ ] Re-run Kernel Windows build and runtime tests with the unsuppressed common
      warning policy.

Non-negotiable constraints for this work:

- **Do not add Windows-specific C++ code paths.**
- **Do not suppress compiler warnings** (`/wd*`, `-Wno-*`, or equivalent)
  to obtain a green build.
- Fix diagnostics at the owning layer and validate that layer before
  propagating its revision upward.
- Keep platform-specific logic in CMake limited to real compiler, linker or
  packaging requirements and centralize common policy in Configuration.
- Until this decoupling/portability pass is closed, use `fix-windows` in every
  repository updated by the effort and expose the change through a PR so CI is
  visible before upward propagation.

## Cross-cutting packaging cleanup

- [ ] Remove Ariadne's remaining source-tree knowledge of Numeric used only to
      copy `FindGMP.cmake` and `FindMPFR.cmake`; keep the installed
      `AriadneConfig.cmake` contract valid while doing so.
- [ ] Continue using `ariadne_register_public_headers` and
      `ariadne_install_dependency_bundle` instead of manual transitive header
      installation lists.

## Validation discipline

For every boundary change:

1. compile affected public headers in isolation;
2. build the component without sibling include leakage;
3. consume only the immediate lower repository at the C++ boundary;
4. for Python, consume only the aggregate
   `pyariadne-<direct-dependency>` interface and utility header;
5. keep layer utility headers free of helpers already owned below;
6. run focused tests and standalone installation/consumer checks;
7. propagate the validated revision upward;
8. run Ariadne integration CI.

Do not treat successful aggregate Ariadne linkage as proof of a standalone
boundary.
