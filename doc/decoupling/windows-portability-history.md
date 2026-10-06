# Windows portability and Kernel decoupling history

Status date: 2026-10-06

This document records the technical directions taken while making the decoupled
Ariadne stack build and link cleanly on Windows. It is intentionally a history
of both successful and unsuccessful approaches, so that failed experiments are
not repeated and green CI is not obtained by hiding diagnostics.

The non-negotiable rules for this work are:

- **no Windows-specific C++ implementation paths**;
- **no compiler-warning suppression**, including MSVC `/wd*` and equivalent
  mechanisms;
- warnings and linker failures must be fixed at the layer that owns the faulty
  code or build semantics;
- platform-specific CMake is acceptable only for genuine compiler, linker or
  packaging constraints;
- while this work is active, modified repositories use `fix-windows` branches
  and PRs so that the CI of each layer is visible before propagation upward.

## Current restart point

The next pass restarts from the latest Numeric main branch:

```text
ariadne-cps/numeric:main
bbb57e946f86d6d0134dcb90b81a86102ccae816
```

This revision descends from
`413a2a00ae61e96d6734196f7b0faac139d5672e`, which made the Numeric
narrowing conversions reported by Algebra/MSVC explicit.

The intended propagation order is:

```text
numeric:main
  -> interval/fix-windows
       -> algebra/fix-windows
            -> kernel/fix-windows
```

The independent lower branch remains:

```text
configuration/fix-windows
  -> utility/fix-windows
       -> logging/fix-windows
            -> threading/fix-windows
```

At the last validated checkpoint, Utility, Foundation, Logging, Threading,
Numeric and Interval were green. Algebra was the first layer above Interval to
fail under the unsuppressed Windows build.

## Directions that worked and are retained

### Generic/concrete virtual split for Function

MSVC exposed ambiguous/final-overrider problems such as C2250 in the Function
hierarchy. The successful direction was to separate generic virtual operations
from concrete clone/create operations and make the mixin the unique final
overrider.

This removed the C2250 family without compiler-specific branches and is part of
the retained architecture.

### Member-only explicit instantiation

Removing all FunctionMixin explicit instantiations fixed one class of MSVC
errors but caused missing out-of-line symbols on macOS. Reintroducing
class-wide explicit instantiation then produced MSVC C2908 failures.

The successful compromise was **member-only explicit instantiation** of the
required out-of-line `FunctionMixin::_call` members. This supplies the symbols
needed by Unix/macOS without asking MSVC to instantiate unrelated or invalid
members.

Class-wide FunctionMixin explicit instantiation must not be reintroduced.

### Complete types before explicit instantiation

MSVC C2139 failures around TaylorModel explicit instantiation were caused by
instantiating against incomplete types. Including the full Taylor model
definition before the relevant FunctionPatch instantiations fixed that failure.

The retained rule is that explicit instantiation sites must see complete types
for every operation they materialize.

### Static Windows aggregate packaging

The original shared aggregate used automatic Windows symbol export and
eventually failed with:

```text
LNK1189: library limit of 65535 objects exceeded
```

Switching the aggregate to a static library on Windows avoided the oversized
import-library problem and allowed the Python extension to link.

This direction was first implemented locally in Kernel, then generalized in
Configuration so individual projects do not carry their own Windows packaging
policy.

The current Configuration direction is:

- aggregate project libraries are `STATIC` on Windows and `SHARED`
  elsewhere;
- `setup_project_library` honors the configured library kind;
- common MSVC `/bigobj` policy lives in Configuration.

This is the current packaging model being validated. A future DLL design with
explicitly curated exports remains possible, but it is not a substitute for the
current correctness work.

### Explicit static-library API semantics for Numeric

After Numeric became static on Windows, treating its public symbols as
`dllimport`/`dllexport` was no longer coherent.

The retained solution is an explicit static-build state
(`ARIADNE_NUMERIC_STATIC`) that makes `ARIADNE_NUMERIC_API` empty for
static Windows consumers/builds, while preserving the build/import distinction
for a possible future DLL configuration.

Using `ARIADNE_NUMERIC_EMBEDDED` as a general synonym for a static library was
rejected because it described a different condition.

### Removing ODR violations exposed by static linking

Static Numeric exposed duplicate definitions that had previously lived in
separate DLL/executable link units, including:

- `class_name<DoublePrecision>()`;
- `class_name<MultiplePrecision>()`;
- `operator""_dec(const char*, std::size_t)`.

The correct fix was to remove the duplicate test definitions and use the
library-owned definitions. No linker workaround was added.

This established an important rule: if the static build exposes an ODR
violation, fix the duplicate ownership rather than restoring a DLL merely to
hide it.

### Replacing layout-dependent interval casts

Interval had conversions implemented by reinterpreting one interval object as a
different interval instantiation. This produced invalid values on Windows,
including the previously observed corrupted `denorm_min<double>` path.

The retained fix constructs the target interval from the underlying bound
values (`.raw()`) rather than relying on object layout.

This is a portable semantic fix, not a Windows branch.

### Local Taylor error value conversion

Taylor model arithmetic contained local layout-dependent error conversions.
Replacing them with explicit value conversion through the underlying raw value
fixed the affected Taylor arithmetic paths without changing the global Numeric
error API.

This local approach is retained.

### Source-level C4244 cleanup in Numeric

Once Algebra warning suppression was removed, Algebra instantiated Numeric
templates that exposed C4244 narrowing diagnostics in Numeric headers. Numeric
standalone CI had not exercised all of those consumer instantiations.

The correct ownership was Numeric. Commit
`413a2a00ae61e96d6734196f7b0faac139d5672e` made the intentional narrowing
conversions explicit, including:

- `long double -> double`;
- integral values -> `double`;
- `unsigned long long -> unsigned long` with the existing range assertion;
- `double -> float`.

That commit is contained in the current Numeric main baseline.

## Directions that were useful diagnostically but were not solutions

### `/WHOLEARCHIVE`

Forcing whole-archive linking was tested while investigating Windows runtime
failures. It did **not** fix the runtime issue.

It did expose real ODR problems in header-defined non-template functions,
including `exp2(Integer)` and `preimage(...)`. Those functions were corrected
at source by making their header definitions inline.

The ODR fixes are retained; `/WHOLEARCHIVE` is not.

### Global rounding-mode compiler flags

A pass propagated strict rounding compiler flags more broadly
(`/fp:strict` on MSVC and `-frounding-math` elsewhere) to test whether the
Taylor runtime failures were caused by compiler floating-point assumptions.

That experiment did **not** fix the runtime failures and was reverted as a
global policy.

Numeric may still carry target-local rounding options where they are part of
Numeric semantics. Such options are not warning suppression and must not be
confused with global portability policy.

### Static linking as a diagnostic amplifier

Static linking is not itself a fix for C4244, C4661, ODR mistakes or invalid
casts. It has nevertheless been useful because it brings more code into the
same link unit and exposes ownership/ODR defects that a DLL boundary can hide.

For the current work this diagnostic property is considered beneficial.

## Directions that failed and were reverted/rejected

### Warning suppression

MSVC suppressions such as:

```text
/wd4244
/wd4250
/wd4267
/wd4459
/wd4661
/wd4702
```

were introduced during earlier attempts to get Windows builds through
`/WX`. Algebra also carried local suppressions for C4244, C4661 and C4702.

This direction is rejected. It hid real portability problems and made the
integration state look healthier than it was.

The policy is now explicit: **no `/wd*` anywhere in the Windows dependency
chain**. Diagnostics must remain visible and be fixed at source.

### Global Numeric Error API rewrite

An attempt to change the global `numeric::Error` conversion/accumulation API
to avoid Windows layout/conversion problems broke valid Unix code such as
Taylor error accumulation.

The experiment was reverted. The lesson is to avoid changing global Numeric
semantics when the actual defect is a local unsafe conversion.

### Class-wide explicit instantiation as a warning workaround

Broad explicit instantiation was repeatedly shown to force compilers to inspect
or instantiate optional/undefined members. On MSVC this manifested as C2908 and
C4661 families.

Suppressing C4661 is rejected, and class-wide instantiation should not be used
merely to force symbol materialization. Prefer explicit instantiation of the
actual required members/functions.

### Kernel-only Windows build policy

Making Kernel alone static and adding local compiler flags was useful to get
past the first Windows packaging failures, but it left submodules compiling
under inconsistent policies.

That direction has been superseded by centralizing common Windows policy in
`configuration` and propagating the same Configuration revision through the
dependency graph.

## Current Algebra failure and next work

After Interval became green, Algebra was rebuilt with all of its former
warning suppressions removed. The resulting Windows failure separated into two
main families:

1. **C4244 originating in Numeric headers** when Algebra instantiated Numeric
   templates. These were correctly fixed in Numeric by
   `413a2a00ae61e96d6734196f7b0faac139d5672e`, now contained in Numeric main.
2. **C4661 originating in Algebra's broad explicit class instantiations**,
   especially Matrix, Polynomial, Differential and Expansion families. These
   remain Algebra-owned work.

The next iteration should therefore:

1. repoint Interval to the current `numeric:main` baseline;
2. validate Interval CI;
3. repoint Algebra to that validated Interval revision;
4. address Algebra C4661 at source, preferring member/function-level explicit
   instantiation over class-wide instantiation where only a subset is actually
   defined or required;
5. keep `/WX` enabled and do not add suppressions;
6. only after Algebra is green, propagate Algebra and the already validated
   Threading branch into Kernel.

## Packaging direction deliberately deferred

A Windows DLL with **explicit** `dllexport`/ `dllimport` on the curated public
ABI could plausibly avoid the original LNK1189 caused by automatic export of a
very large aggregate.

That is a separate packaging/ABI project. It would not solve C4244, C4661,
ODR violations, unsafe casts or runtime numerical defects; at most it would
change where some of them become visible.

For that reason the current work continues with static Windows aggregates until
the C++ and template-instantiation model is warning-clean and correct. DLL
packaging can then be reconsidered without using it to mask defects.
