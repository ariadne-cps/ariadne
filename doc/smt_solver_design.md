# SMT Solver Design and Development Status

This document records the architectural decisions, semantic contracts, testing
policy, and remaining work for the epsilon-SMT solver developed on
`solvers-smt#830`. It is intended to preserve the reasoning behind the code,
not to duplicate the Git history.

## Objective

The long-term objective is a rigorous bounded-real epsilon-SMT solver with
capabilities comparable in purpose to dReal: Boolean reasoning over nonlinear
real theory atoms, validated numerical pruning and certification, controlled
epsilon relaxation, and parallel box search.

The solver currently exposes three outcomes:

- `UNSAT`: the bounded query has been rigorously excluded.
- `EPSILON_SAT`: a validated witness box satisfies the epsilon-relaxed query.
- `UNKNOWN`: the current bounded search cannot certify either result. This is
  used deliberately for exhausted box budgets and terminal boxes that cannot be
  certified; it must never be silently promoted to `EPSILON_SAT`.

## Main implementation

The implementation is split across:

- `source/solvers/smt_solver.hpp/.cpp`: numerical search, DPLL/CDCL
  integration, sequential and parallel search, statistics and result handling.
- `source/solvers/smt_boolean.hpp/.cpp`: Boolean encoding and canonicalization
  of theory atoms.
- `source/solvers/smt_theory.hpp/.cpp`: normalization of real theory
  predicates into primitive literals.
- `tests/solvers/test_smt_solver.cpp` and
  `tests/solvers/test_smt_boolean.cpp`: semantic and coverage tests.

## Theory semantics

Normalized primitive theory relations are represented relative to zero.
Equality and non-strict inequalities use interval bounds; strict positivity
(`GT_ZERO`) is treated separately where closed interval bounds cannot express
strictness directly.

Epsilon bounds are a relaxation of the original theory bounds. A witness is
reported only after validated epsilon certification.

Expressions are simplified before interval solving. Primitive zero expressions
that disappear under simplification are handled without introducing useless
theory work.

Equivalent theory atoms are canonicalized and shared by the Boolean encoder.
Symmetric equalities and disequalities are canonicalized consistently, while
atoms differing in canonical relation, left-hand side, or right-hand side
remain distinct.

## Numerical box processing

Sequential and parallel conjunction solving share the same box-processing
semantics.

A box is processed in this order:

1. validated reduction using the original (non-epsilon) constraints;
2. epsilon certification of the whole reduced box;
3. deterministic witness candidates;
4. sensitivity-guided splitting;
5. optional nonlinear interior-point candidate search for splittable boxes;
6. `UNKNOWN` for a non-splittable box that remains uncertified.

The interior-point candidate is only a candidate. It is never trusted without
validated epsilon certification.

Candidate optimization is skipped on terminal/non-splittable boxes. Once the
deterministic point checks have failed, an optimizer cannot produce a distinct
point inside a singleton box, and invoking it there introduced unnecessary
numerical failure modes.

## Reduction and splitting

The solver combines hull reduction and box shaving. Direct validated range
checks are used where contractor propagation alone is insufficient to expose
infeasibility.

Splitting is sensitivity-guided when a nonzero validated derivative enclosure
provides useful information. A coordinate is considered inactive only when the
evaluated derivative enclosure is representation-wise exactly `[0,0]`; all
other cases are treated conservatively as active.

If sensitivity information is unavailable, the solver falls back to the
geometrically widest coordinate. Statistics distinguish sensitivity-guided
splits and cases where sensitivity overrides the geometric choice.

## Boolean reasoning and CDCL

Boolean combinations are encoded into clauses and solved together with theory
checks.

The current implementation includes:

- unit propagation;
- decision levels and backtracking;
- Boolean conflict analysis;
- learned clauses;
- theory-generated nogoods;
- nonchronological backjump accounting;
- learned-clause activity;
- conservative learned-clause pruning;
- theory-nogood minimization with a configurable budget.

Theory results are interpreted centrally: `EPSILON_SAT` is consistent and
provides a witness, `UNKNOWN` is provisionally consistent but records theory
uncertainty, and `UNSAT` is inconsistent.

Boolean reasoning is allowed to continue after an individual theory search
returns `UNKNOWN`; a local numerical budget exhaustion must not automatically
terminate Boolean reasoning if another Boolean assignment can still decide the
query.

## Parallel search

Parallel conjunction search uses BetterThreads `DynamicWorkload`.

Tests explicitly distinguish the two BetterThreads modes:

- concurrency zero: processing occurs on the calling thread;
- positive concurrency: processing is delegated to worker threads.

The tests also force real parallel box splitting for both validated constraints
and compiled theory literals. This established that `solve_parallel()` is not
silently falling back to the sequential search.

The parallel task is represented by the named `SmtParallelTask` callable
rather than embedding the complete search logic in an anonymous lambda.

The shared parallel state uses atomic stop/result flags and a mutex for
statistics and witness state. Child boxes produced by `SPLIT` are appended
back to the dynamic workload.

## Public preconditions and internal invariants

Public API preconditions are part of the tested contract. Examples include:

- epsilon must be strictly positive;
- solve domains must be bounded;
- constraint argument dimensions must match the box dimension;
- `RealSpace` dimension must match the box dimension;
- accessing a witness from a non-`EPSILON_SAT` result is rejected.

Internal assertions are not retained merely as defensive decoration. The
coverage policy is:

- if an invalid state can enter through a public boundary, validate it at that
  boundary and test the negative case;
- if an internal state is genuinely representable and must be rejected, keep a
  check and provide a meaningful test path;
- if a state is impossible by construction or is already implied by the
  immediately preceding control flow, remove the redundant assertion or encode
  the invariant structurally.

This policy avoids maintaining branches corresponding to impossible states
while preserving checks for genuinely possible invalid inputs.

## Sequential/parallel consistency

The sequential and parallel search paths share conjunction dispatch and box
processing. Their scheduling differs, but result semantics and statistics are
intended to remain consistent.

Empty domains are `UNSAT`. Empty conjunctions are solved without unnecessary
box processing. Global box-processing limits produce `UNKNOWN` rather than an
unsound satisfiability result.

## Testing policy

New SMT functionality is required to have nontrivial tests. Tests should
exercise behavior through public solver APIs when practical. Test-support
helpers are used for deterministic testing of internal algorithms where an
end-to-end construction would be unstable, nonterminating, or would test an
unrelated numerical component instead.

Tests must not be weakened merely to make a regression pass. A failing test
must be classified as a solver defect, an incorrect expectation, or an
undesirable dependency between test phases before it is changed.

Numerically pathological tests that trigger unrelated optimizer invariants are
not accepted as coverage tests. Likewise, deliberately nonterminating search
constructions are removed rather than retained as stress tests without an
explicit bound.

## Coverage contract

The target for SMT functionality is literal 100% coverage for:

- functions;
- lines;
- branches.

LLVM region coverage is useful diagnostically but is not currently a release
criterion. Regions can reflect source-mapping details of templates, macros and
callables that do not correspond to additional semantic behavior.

A missed function, line, or branch must be resolved in one of two ways:

1. add a meaningful test demonstrating the corresponding supported behavior; or
2. remove/refactor the branch when the corresponding state is impossible,
   redundant, or represented at the wrong abstraction boundary.

Coverage must not be increased by manufacturing impossible internal states.

At the last fully measured checkpoint before the most recent coverage cleanup,
`smt_boolean.cpp` and `smt_theory.cpp` had reached 100% functions, lines and
branches. `smt_solver.cpp` had reached 100% functions and was being reduced
from a small residual set of line/branch misses. The exact current percentage
must always be taken from a fresh coverage run after the latest commit rather
than from this document.

## Coverage-related design work

Several implementation changes were motivated by discovering genuine
structural issues while pursuing complete coverage:

- sequential and parallel conjunction search were unified where semantics were
  duplicated;
- unreachable post-contractor checks were removed after verifying contractor
  contracts;
- impossible epsilon-overlap branches were removed rather than tested
  artificially;
- candidate-search states were simplified so that only semantically meaningful
  combinations are represented;
- repeated signed-literal-to-variable conversion is centralized;
- point-box construction used by witness and epsilon checks is explicit rather than
  lambda-generated; this avoids compiler-generated template control-flow being
  attributed to `smt_solver.cpp` as unsupported semantic branches while preserving
  the same tested solver behavior;
- redundant internal assertions have been removed where their failure state was
  impossible by construction;
- test-only classification scaffolding is removed when production control flow already
  expresses the invariant directly; coverage helpers must not create extra state-space
  branches of their own;
- public precondition failure paths are tested explicitly.

Coverage is therefore being used as an audit of the state space, not only as a
test-count metric.

## Important rejected approaches

The following approaches were tried or considered and deliberately rejected:

- treating epsilon overlap alone as `EPSILON_SAT`: overlap is not
  certification;
- running the interior-point optimizer on terminal singleton boxes: it cannot
  discover a different in-domain point and can introduce numerical failures;
- high-precision terminal fallback without a clearly justified validated
  semantics;
- coverage tests based on unstable optimizer degeneracies;
- unbounded/nonterminating DPLL constructions created solely to hit a rare
  branch;
- preserving branches known to be unreachable merely so they can be
  artificially exercised.

The Git history contains the experimental details; this document records the
resulting design decisions.

## Current open work

The immediate work on `solvers-smt#830` is:

1. finish 100% function/line/branch coverage for `smt_solver.cpp`;
2. audit each remaining branch as a real state-space case or a redundant
   representation;
3. preserve deterministic, nontrivial tests for every retained behavior;
4. once the coverage baseline is complete, continue improving epsilon-solver
   completeness and numerical refinement toward the long-term dReal-like goal;
5. later add CI coverage gates so regressions in functions, lines or branches
   fail automatically.

The current development workflow intentionally uses local compilation and test
execution; CI integration is deferred until the solver and its coverage
contract are stable.

## Maintenance rule

When an architectural or semantic decision changes, update this document in the
same commit as the code change. Do not use this file as a chronological
changelog; Git already provides that. Keep it as the current design rationale,
coverage contract, and handoff state for future development.
