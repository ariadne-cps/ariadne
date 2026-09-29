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

### Primary application target: neural certificate verification

The near-term development target is formal certification of neural-network
bounds in workflows such as FOSSIL-style Lyapunov/barrier verification and
CARe-style residual certification. These workflows translate a trained smooth
network and its symbolic derivatives into bounded quantifier-free nonlinear
real arithmetic queries over compact domains.

The solver roadmap therefore prioritizes, in this order:

1. epsilon-safe SAT/theory interaction with explained theory propagation;
2. strong validated ICP for large smooth composed expressions, including
   monotonicity/Newton-style contraction where applicable;
   The SMT ICP fixed point now invokes Ariadne's validated
   `ConstraintSolver::monotone_reduce` after hull reduction and coordinate
   shaving stall, and records monotone rounds/effective contractions separately.
   Because the contractor divides by a derivative interval, it is invoked on a
   coordinate only when validated derivative range analysis proves that interval
   strictly positive or strictly negative on the current box; coordinates whose
   derivative interval contains zero are skipped. Each box-reduction call performs
   at most one monotone sweep. The monotone/Newton phase is disabled by default:
   experiments on the established transcendental regression `sin(x)=0` over
   `[3,4]` showed that even a bounded monotone sweep can substantially perturb
   the normal witness/search path. It is therefore opt-in through solver
   configuration and intended first for measured smooth-network workloads such as
   FOSSIL/CARe. The underlying `ConstraintSolver::monotone_reduce` is also
   bounded to its declared three Newton steps: its previous width-only
   `do/while` could fail to terminate when a Newton step made no progress. The
   direct constraint-solver regression now calls `monotone_reduce` explicitly
   on linear, smooth, singleton and rigorously infeasible monotone cases. SMT
   regressions exercise both validated-constraint and normalized-theory paths,
   including effective contraction and validated nonmonotone-coordinate skipping;
3. efficient handling of large shared expression DAGs produced by feed-forward
   networks and their derivatives. As a first step, normalized SMT theory literals
   cache validated coordinate derivatives at compilation time when the original
   `RealExpression` is structurally differentiable. The AST is inspected before
   calling `derivative()`; expressions containing non-smooth operators such as
   `abs`, `max` or `min` receive no cached derivative and remain supported by
   the normal hull/shaving ICP path. This is necessary because the symbolic
   function backend aborts rather than throws for derivatives of `max/min`.
   The monotone contractor accepts an already compiled derivative so gating and
   Newton reuse the same object instead of rebuilding it per box. Generic
   `ValidatedConstraint` inputs do not retain the originating expression AST, so
   monotone/Newton contraction is intentionally not applied on that entry path;
   the optimization is currently restricted to compiled SMT theory literals;
4. robust support for polynomial and transcendental activations and dynamics;
5. fast counterexample/witness discovery for CEGIS loops;
6. preserve rigorous answer semantics in DP: UNSAT and EPSILON_SAT are
   certified, while UNKNOWN distinguishes explicit resource exhaustion from
   fixed-precision DP-resolution exhaustion;
7. only after those goals, broader language features such as richer SMT-LIB
   integration, quantifiers, or dedicated ODE solving.

For this target, ODE trajectories are normally not integrated by the SMT core.
FOSSIL-like verification supplies the Lie derivative as the symbolic real
expression grad(C)(x) dot f(x), while CARe-like verification supplies symbolic
time/state derivatives and Hamiltonian residual expressions. The SMT task is
therefore primarily bounded QF_NRA with nonlinear/transcendental terms, Boolean
structure, and large composed expression graphs.


The solver currently exposes three outcomes:

- `UNSAT`: the bounded query has been rigorously excluded.
- `EPSILON_SAT`: a validated witness box satisfies the epsilon-relaxed query.
- `UNKNOWN`: neither proof is available. `unknown_reason()` distinguishes
  `RESOURCE_EXHAUSTED` from `DP_RESOLUTION_EXHAUSTED`; `MIXED` records that
  both causes occurred along different explored branches.

### Double-precision numerical contract

The SMT solver is intentionally a double-precision validated solver. Search
geometry, interval evaluation, contractors, deterministic witness candidates,
sensitivity analysis, sequential/parallel queues and the public witness all use
the existing DP Ariadne types.

Multiple-precision evaluation and an unbounded exact search-cell representation
are not part of the architecture of `solvers-smt#830`. They must not be used as
a fallback for terminal boxes or as a prerequisite for completeness.

#### Completeness audit

The delta-completeness argument for DPLL(ICP) in Gao, Avigad and Clarke relies
on two ingredients that must be kept distinct:

1. pruning is well-defined: it contracts boxes, does not retain a box that is
   already proved inconsistent by the interval extension, and never removes a
   real solution;
2. interval extensions become sufficiently narrow when the search box becomes
   sufficiently small. In the paper this is captured by delta-regular interval
   extensions and by choosing a geometric ICP stopping precision from a uniform
   modulus of continuity.

The current SMT reduction is compatible in spirit with the first ingredient:
hull/shaving contraction is followed by validated range rejection, and UNSAT is
reported only after validated exclusion. The second ingredient does not hold
uniformly for arbitrarily small epsilon under a fixed DP arithmetic model.
Validated DP evaluation has a nonzero representation/rounding floor for some
expressions even on a singleton box, while DP bisection itself eventually
reaches adjacent representable endpoints.

The regression
`sqr(sin(x))+sqr(cos(x))-1 == 0` at the singleton `x=1` with
`epsilon=1e-30` is the concrete witness of this limitation. The exact real
expression is zero, but the validated DP enclosure cannot be made narrow enough
to certify that epsilon. No further DP split exists. Returning
`EPSILON_SAT` would therefore require information not supplied by the DP
interval evaluator, while returning `UNSAT` would be incorrect.

Consequently there is no sound general terminal rule, using only the current
fixed-precision interval oracle, that can eliminate terminal `UNKNOWN` for
every positive epsilon and every supported transcendental expression. Doing so
would require at least one of: precision escalation, stronger symbolic
identities/exact reasoning, or a restriction of the numerical contract. The
first option is explicitly rejected for this solver.

This is also consistent with the practical ICP shape used by dReal: its public
dReal4 evaluator documentation accepts a box as delta-satisfying when a formula
is already valid on the box or when the interval evaluation is narrow enough
relative to the requested precision. Such a width criterion is sound only when
the numerical evaluator can actually reach the requested scale; it does not
remove a fixed-precision floor.

The Ariadne SMT contract is therefore:

- `UNSAT` is exact: validated pruning/range exclusion has ruled out the
  original query.
- `EPSILON_SAT` is certified: the returned witness box is validated against
  the epsilon-relaxed constraints.
- `UNKNOWN(RESOURCE_EXHAUSTED)` means the configured search budget stopped a
  search that could otherwise have continued.
- `UNKNOWN(DP_RESOLUTION_EXHAUSTED)` means a box still required refinement but
  no distinct DP children could be produced and no epsilon witness had been
  certified.
- `UNKNOWN(MIXED)` can arise in Boolean/parallel exploration when unresolved
  branches encountered both causes.

This deliberately differs from dReal4's pragmatic sequential escape hatch,
which returns true when the delta condition is unmet but no active dimension is
bisectable. Ariadne does not promote that implementation limit to
`EPSILON_SAT`: fixed-precision exhaustion remains explicit `UNKNOWN`. The
solver follows the delta-complete branch-and-prune ideas where they are sound in
DP, while preserving a stronger result contract.

#### Relation to the current implementation

After original-domain reduction, `_process_box` currently tries a validated
midpoint witness, deterministic point candidates and, for splittable boxes, an
optional interior-point candidate. The function named `_epsilon_satisfied`
certifies a **point box** obtained from `midpoint_box(domain)`; it is not a
whole-box stopping rule. The design documentation must not describe this as
whole-box epsilon certification.

If all point certifications fail and `Box::split` produces identical children,
the only sound current outcome is `UNKNOWN`. This is not an implementation
branch waiting for a generic SAT rule; it is the observable DP-resolution
boundary.

A dReal-style interval-width stopping test remains a valid optional optimization:
if a validated image of the whole reduced box is contained in the epsilon-relaxed
target, the box itself is an epsilon witness. More generally, after a
well-defined prune step, sufficiently narrow interval images can justify the same
conclusion. However, this cannot solve the hard terminal case in which even
singleton DP evaluation is wider than epsilon. It should therefore be treated as
an early certification optimization, not as the missing completeness theorem.

### Whole-box epsilon stopping and DP exhaustion

The production solver applies epsilon certification to the validated image of
the whole reduced box. For each normalized constraint the complete interval
image must satisfy the epsilon-relaxed target; strict positivity keeps its
strict lower-bound test. The certified reduced box itself is returned as the
witness.

This uses the same interval branch-and-prune principle as dReal, but Ariadne's
answer contract is intentionally stricter. Direct containment in the explicit
epsilon-relaxed target is a sound stopping condition and can be stronger than a
generic diameter threshold. If containment and point-candidate certification
both fail, splitting continues.

If no distinct DP children can be produced, Ariadne returns
`UNKNOWN(DP_RESOLUTION_EXHAUSTED)`, never `EPSILON_SAT` without a validated
witness. The statistics counter `dp_resolution_exhaustions` records these
events. Resource limits instead produce `UNKNOWN(RESOURCE_EXHAUSTED)`.

### UNKNOWN taxonomy

There are two independent inconclusive mechanisms:

1. **Resource exhaustion.** The configured box-processing limit is reached while
   pending search remains.
2. **DP-resolution exhaustion.** A surviving box is not epsilon-certified and
   cannot be refined into distinct DP children.

Sequential and parallel conjunction search accumulate these causes separately.
Boolean/CDCL theory search propagates the reason from unresolved theory
branches. If both causes contribute before the overall query terminates,
`SmtUnknownReason::MIXED` is returned.

This distinction is semantic, not cosmetic: increasing a box budget can address
resource exhaustion but cannot repair DP-resolution exhaustion. Conversely,
loosening epsilon or improving validated expression evaluation/contractors can
reduce DP-resolution exhaustion without changing the resource budget.

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
2. validated epsilon certification of the whole reduced box;
3. deterministic witness candidates;
4. sensitivity-guided splitting;
5. optional nonlinear interior-point candidate search for splittable boxes;
6. `UNKNOWN(DP_RESOLUTION_EXHAUSTED)` for a non-splittable box that remains
   uncertified.

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
- validated whole-domain theory-atom classification (`TRUE`, `FALSE` or
  `UNKNOWN`);
- epsilon-safe whole-domain implication analysis for individual theory atoms is
  available as tested infrastructure, but is not currently injected into the CDCL
  search. Domain-only unit implications proved too strong operationally: they
  collapse Boolean branches that intentionally exercise resource budgets,
  backtracking and theory-conflict learning. Future theory propagation therefore
  needs contextual explanations derived from the current partial theory state;
- nonchronological backjump accounting;
- learned-clause activity;
- conservative learned-clause pruning;
  The current CDCL implementation retains learned clauses once created. Candidate
  ranking and pruning statistics remain scaffolding for a future database
  reduction scheme, but clauses are not deactivated until that scheme has an
  explicit progress guarantee; experiments showed that aggressive deactivation
  can rediscover the same conflicts indefinitely. The configured learned-clause
  limit is therefore currently advisory rather than a hard bound.
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

At the baseline immediately before resuming functional development,
`smt_boolean.cpp` and `smt_theory.cpp` had 100% functions, lines and
branches, while `smt_solver.cpp` had 100% functions and branches. One
compiler-mapped line in the parallel task remained the previously audited line
coverage anomaly. The exact current percentage must always be taken from a
fresh coverage run after the latest functional commit rather than from this
document.

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
- repeated signed-literal-to-variable conversion is centralized, including theory
  nogood reconstruction;
- signed Boolean literal truth values are derived directly from the sign bit
  instead of ternary control flow; the representation invariant is negative = false,
  positive = true;
- parallel worker accounting subtracts the Boolean calling-thread observation
  directly, avoiding a redundant ternary branch in statistics-only code;
- box splittability is represented by whether the two children returned by
  `Box::split` differ; for a degenerate split both children equal the parent, so
  separately comparing each child with the parent represented an impossible
  short-circuit state;
- deterministic witness-candidate tests cover both corner-enumeration cutoffs:
  dimensions whose full corner set is capped and dimensions at or above the
  machine shift width;
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
- exhaustive enum switches retain a `default` defensive branch and tests exercise
  that branch with an explicitly invalid enum value; defaults are not removed merely
  to satisfy branch coverage.

Coverage is therefore being used as an audit of the state space, not only as a
test-count metric.
 Structurally infinite contractor fixed-point loops are written without a Boolean
condition; representing them as `while(true)` creates a compiler-visible false
branch that is not part of the solver state space. Boolean combinations of already-computed side-effect-free state
are non-short-circuit where laziness carries no semantics, so LLVM branch coverage
tracks solver decisions rather than evaluation-order edges. Short-circuiting remains
where it protects optional access or changes observable work.
Optional control state is removed when absence is not representable in the real
solver flow. Sensitivity selection, learned-clause protection and CDCL pivot
eligibility use explicit state whose combinations correspond to semantic decisions
rather than evaluation-order branches. DPLL literal assignment is likewise structural: callers select only unassigned
variables, so reassignment/conflict states are not represented inside the assignment
primitive. UNKNOWN box-processing statistics encode the only currently reachable
terminal-uncertified case directly instead of carrying a redundant Boolean flag.
Parallel child suppression remains semantically meaningful under races, but its
decision is isolated in a deterministic helper so both outcomes can be tested.
Conflict analysis relies on the 1-UIP implication-graph invariant: while more than
one current-level literal remains, a resolvable trail pivot exists. Theory nogoods
likewise contain only variables originating from encoded theory atoms. Iteration over
those structures therefore does not represent a fall-through/end state that the
solver can reach.

## Important rejected approaches

The following approaches were tried or considered and deliberately rejected:

- treating epsilon overlap alone as `EPSILON_SAT`: overlap is not
  certification;
- running the interior-point optimizer on terminal singleton boxes: it cannot
  discover a different in-domain point and can introduce numerical failures;
- multiple-precision terminal fallback or unbounded exact search geometry:
  the solver architecture is deliberately DP-only;
- coverage tests based on unstable optimizer degeneracies;
- unbounded/nonterminating DPLL constructions created solely to hit a rare
  branch;
- preserving branches known to be unreachable merely so they can be
  artificially exercised.
- injecting all domain-fixed theory atoms as unconditional level-zero CDCL
  assignments: although exact-domain classification is sound, the first attempt
  bypassed theory-conflict learning/minimization paths and changed the treatment
  of simplified strict atoms under epsilon semantics. Propagation must therefore
  be integrated as theory implications with reasons rather than raw assignments.
  A second experiment restricted propagation to epsilon-infeasible polarities
  and delayed it until after partial-theory consistency. It remained too strong:
  domain-only unit clauses still collapsed Boolean branches that existing CDCL,
  resource-budget and theory-learning regressions intentionally exercise. The
  epsilon-infeasibility classifier is retained, but production propagation is
  deferred until implications can carry contextual theory explanations.

The Git history contains the experimental details; this document records the
resulting design decisions.

## ConstraintSolver coverage prerequisite

Before moving generic ICP mechanics out of `SmtSolver`, the existing
`ConstraintSolver` is being audited and brought under the same coverage
discipline. The first audit identified file-local dead helpers and a permanently
disabled shaving block in the deprecated list-based `reduce` overload; these
are removed rather than covered artificially. Existing public numerical
operations are exercised directly before any SMT refactoring so that later
behavior changes can be attributed to the refactoring rather than to pre-existing
coverage gaps.

The architectural rule remains that this audit must not force SMT policy into
`ConstraintSolver`: generic contraction, feasibility and splitting mechanics
belong there, while epsilon weakening, epsilon-active literal selection and
SMT result semantics remain in `SmtSolver`.

The second coverage tranche directly exercises both `lyapunov_reduce`
overloads, including contraction, witness preservation and empty detection, and
adds a nonlinear infeasibility case whose natural interval image overlaps the
target so that `feasible` cannot terminate by its initial direct range test.
The audit also found that `ConstraintSolver::feasible` still embedded the
legacy `NonlinearInteriorPointOptimiser` iteration and its old dual/Taylor
fallback, whereas current solver code and dedicated tests use
`NonlinearInfeasibleInteriorPointOptimiser`. The constraint solver now delegates
candidate generation and validated feasibility/infeasibility classification to
the latter, retaining its cheap direct interval-disjointness rejection.
A `true` result from `feasible_candidate` is already backed by
`OptimiserBase::validate_feasibility`, so `ConstraintSolver` no longer
re-certifies it with the weaker pointwise `check_feasibility` routine. The
outer `NearBoundaryOfFeasibleDomainException` catch was likewise removed
because `feasible_candidate` converts that condition to `indeterminate`
internally; only numerical failures that can actually propagate, such as a
singular Newton system, remain translated to `indeterminate` here.
The direct attempts to force the two endpoint-clamp branches in
`monotone_reduce` produced invalid numerical states before those branches
could be reached: the validated arithmetic rejected a negative value where a
positive upper bound is required. Those probes are therefore removed rather
than turning an invalid state into a coverage test. The two clamp branches stay
classified as candidate defensive/dead code until a valid reachable
construction or an invariant proof settles them.

This removes the obsolete embedded optimiser algorithm rather than attempting to
cover it artificially. The separate legacy uses in `paver.cpp` are outside this
refactoring and require their own audit.

Coverage from this tranche is used to determine whether the remaining
dual/Taylor infeasibility path is genuinely reachable with the current
interior-point implementation or should be treated as legacy architecture.

The two endpoint-clamp branches in `monotone_reduce` were subsequently
removed after audit. They were unreachable from valid interval states in testing,
and their semantics were also not contractor-correct: when a Newton image lies
strictly outside the current strip, the validated intersection is empty, not a
singleton at the previous endpoint. The implementation now uses the direct
validated intersection in both lower and upper strip updates. A redundant
post-`box_reduce` empty-domain check in `reduce(vector)` was also removed,
because every mutation in the preceding loop is already followed immediately by
the same empty-domain return.

### ConstraintSolver branch-coverage audit

The post-cleanup LLVM report for `constraint_solver.cpp` has 100% function
coverage and 100% line coverage. Its raw branch summary is 110/140 (78.57%),
but the HTML source view contains 40 explicit source branch sites and every one
has both outcomes exercised. The remaining 30 branch sites are not attached to
an explicit branch annotation in the source view; they originate from expanded
validated/macro/inlined machinery rather than uncovered `if`, loop or
conditional behavior in `constraint_solver.cpp`. They are therefore not to be
chased with artificial inputs. Any future explicit source branch introduced in
this component must still be covered in both directions.

The lifetime audit exposed a real interface defect rather than a coverage
artefact: `ConstraintSolverInterface` is polymorphic but previously lacked a
virtual destructor. Deleting `ConstraintSolver` through an interface pointer
was therefore invalid and correctly triggered
`-Wdelete-abstract-non-virtual-dtor`. The interface now has a defaulted virtual
destructor, the derived destructor is marked `override`, and the direct
lifetime test is retained as a regression test for polymorphic destruction.

### First SMT-to-ConstraintSolver extraction

The first behavior-preserving extraction moves the generic hull/shaving fixed
point used by the validated-constraint SMT path into
`ConstraintSolver::propagate`. The new method takes the existing
`List<ValidatedConstraint>` and works directly on the caller's
`UpperBoxType`; no function/codomain reconstruction or additional dynamic
dispatch is introduced. It preserves the exact ordering used by the SMT solver:
hull propagation plus direct validated range rejection, followed on hull stall
by per-constraint/per-coordinate shaving, repeating until a fixed point or an
empty box.

`ConstraintPropagationStatistics` records hull/shaving rounds and effective
rounds so that `SmtSolver` can preserve its observable search statistics
without owning the numerical propagation loop. The theory-literal path is not
moved in this tranche because strict `GT_ZERO` rejection and derivative-cache
use are SMT-specific concerns that still need a clean generic boundary.

The first coverage run after extraction exposed one dead SMT adapter:
`_original_bounds(ValidatedConstraint)` was no longer called because generic
propagation consumes the constraint bounds directly. It is removed rather than
kept solely for coverage. The same build also exposed three
`-Winconsistent-missing-override` diagnostics in `ConstraintSolver`; the
interface overrides are now marked explicitly with `override`.

### Second SMT-to-ConstraintSolver extraction

The theory-literal reduction path now compiles each normalized literal into a
generic `ConstraintPropagationConstraint`: validated function, closed interval
hull, optional precompiled derivatives, and open/closed endpoint flags. This is
a numerical constraint representation, not an SMT representation. Strict
`GT_ZERO` becomes the generic target `[0,+inf)` with an open lower endpoint.
Hull and shaving operate on the closed interval hull, while direct validated
rejection observes endpoint openness exactly as before.

The precompiled `ConstraintSolver::propagate` overload preserves the former
theory ordering: hull plus direct rejection, shaving on hull stall, then at most
one derivative-assisted monotone pass before fixed point. Derivatives are still
compiled once by the SMT frontend and reused for every box, so the extraction
adds neither per-box derivative construction nor generic virtual dispatch. SMT
configuration still controls whether monotone contraction is enabled.

The compiled SMT theory representation is now the generic numerical propagation
constraint itself. Epsilon weakening remains in `SmtSolver`: the numerical
bounds are shifted by epsilon and open endpoints retain strict comparisons when
certifying a witness. Thus numerical propagation moves down without moving SMT
epsilon-result semantics into `ConstraintSolver`.

The first build of this extraction exposed two ownership-boundary issues. The
SMT header now includes `constraint_solver.hpp` because its compiled-theory
alias names `ConstraintPropagationConstraint` directly. Tests that previously
called the SMT-only monotonicity helper were also adjusted: derivative
compilation/differentiability remains tested in SMT, while the actual monotone
gating behavior is exercised through the generic propagation tests in
`test_constraint_solver`.

### Second-extraction coverage cleanup

Coverage after the precompiled-propagation extraction exposed three useful
boundary facts. The relation-based SMT `_epsilon_bounds` helper had become
production-dead and was only kept alive by an invalid-enum test, so it is
removed. Normalized SMT primitives can have an open lower endpoint
(`GT_ZERO`) but never an open upper endpoint, so epsilon certification no
longer carries an unreachable `strict_upper` branch; upper-end openness
remains a generic propagation feature and is tested directly in
`ConstraintSolver`. Finally, direct generic tests now cover a derivative
vector shorter than the box dimension and a negative monotone derivative,
rather than leaving those generic propagation branches unexercised.

### Epsilon-active splitting

Sensitivity-guided splitting now receives only constraints whose validated image
on the current reduced box has not yet met the epsilon stopping condition.
Constraints that are already epsilon-certified remain part of the conjunction
and continue to be checked after every contraction, but they no longer steer the
split coordinate. This applies uniformly to validated constraints and compiled
theory literals, including strict lower endpoints.

The implementation reuses the same per-constraint epsilon certification used by
whole-box stopping, so active-set selection cannot drift semantically from
`EPSILON_SAT` certification. End-to-end regressions use a deliberately wider
second coordinate controlled only by an already certified constraint/literal;
the unresolved first-coordinate constraint must override the geometric split
that the inactive wide coordinate would otherwise induce.

## Neural benchmark plan

The first external neural workload will use the public Barrier 3 benchmark
family rather than an ad-hoc network. Two distinct public references are kept
separate:

- the original FOSSIL Barr3 benchmark uses the two-state dynamics
  `x_dot=y`, `y_dot=-x-y+x^3/3`, domain `[-3,2.5]x[-2,1]`, and a
  sigmoid `2-10-10-1` barrier network;
- the 2026 L4DC reproducibility package for scalable neural-CBF verification
  reuses Barrier 3 with the same dynamics and publishes a pretrained
  `2-64-64-1` tanh network, together with a FOSSIL/dReal verification script.
  The paper reports that the identically trained Barrier 3 network is verified
  by dReal and uses it as a scaling comparison.

This second artifact matches the eventual Ariadne target substantially better
than a synthetic `2-8-8-1` network. The benchmark progression will therefore
start from a smaller internal representation only if needed for suite runtime,
but the external reference target is the published Barrier 3 `2-64-64-1`
tanh model. The Ariadne benchmark must contain only frozen numeric parameters
and internally constructed `RealExpression` objects; it must not depend at
runtime on PyTorch, ONNX, FOSSIL or another ML framework.

The public pretrained `barr3_cbf.pth` and `barr3_cbf.onnx` files are stored
through Git LFS in the external reproducibility repository. Before adding the
Ariadne fixture, extract the exact matrices and biases from one of those files
and record their source and checksum. Do not substitute newly trained or random
weights while describing the result as the published benchmark.

### First in-suite neural scaling fixture

The first permanent neural regression is a `2-8-8-1` prefix extracted
deterministically from the published Barrier-3 `2-64-64-1` checkpoint. The
extraction keeps the first eight first-layer neurons, the top-left `8x8`
second-layer block, the corresponding first eight biases and output weights,
and the original output bias. It is deliberately labelled a **scaling
fixture**, not the published Barrier-3 certificate: truncating the trained
network changes the represented function and carries no published safety
claim.

The source checkpoint SHA-256 is
`788fd56f21abb0cb5b83206b12d4d8d0ed718896c2140222a8fd129475f6c213`.
The concatenated raw float32 bytes of the extracted tensors have SHA-256
`9431d092372b0bc57ffcf4b651d6bcd2f664c29d79988f7055880f79c4c53aec`.
The C++ fixture stores every parameter as a hexadecimal binary floating-point
literal, so each source float32 is promoted exactly to DP.

Ariadne now has a primitive symbolic `tanh` expression node with validated
numeric support. The Barr3 fixtures use this native operator directly rather
than expanding `tanh(z)` through exponentials. This avoids the interval
dependency and overflow problems of the former quotient representation.

The regression has two purposes. First, it evaluates the frozen network at the
origin and checks a tight validated output interval, catching tensor-layout or
parameter-order mistakes. The test avoids naming local variables `result`
inside `ARIADNE_TEST_ASSERT` scopes because the test macro itself declares a
temporary named `result`; using the same identifier causes a self-initializer
compile error after macro expansion. Second, it asks the SMT solver to solve the
nontrivial equation `B(x,y)=0` over the full Barr3 domain
`[-3,2.5]x[-2,1]` with a one-box budget and candidate search disabled. The
expected result is `UNKNOWN(RESOURCE_EXHAUSTED)` after a genuine
sensitivity-guided split. This exercises expression compilation, cached
derivatives, validated reduction, epsilon checking and neural-expression split
selection without turning the ordinary test suite into a long benchmark.

The ordinary neural regression is intentionally kept cheap. The executable
`tests/solvers/test_smt_neural_benchmarks.cpp` retains only the published Barr3
`2-64-64-1` model: it checks the validated origin value, theory-literal
normalization, and one-box end-to-end SMT processing using the cheap geometric
path with sensitivity, witness probing, shaving and hull reduction disabled.
This preserves a real neural-expression regression while avoiding contractor
performance work in every normal CTest and coverage run.

The earlier `8x8`, `16x16` and `32x32` prefix solves remain useful historical
scaling fixtures, but their repeated one-box solves are no longer part of the
ordinary regression suite. Performance and scaling measurements now belong to
`benchmark_smt_barr3_verification`, which exercises the actual published
verification workload and reports detailed phase statistics.

The second scaling fixture is `2-16-16-1`, extracted by the
same prefix rule: first 16 first-layer neurons, top-left `16x16` second-layer
block, corresponding biases/output weights, and the original output bias. Its
concatenated raw float32 tensor bytes have SHA-256
`2aa09bc348938809ce0f904779d8585561e7ce7017bf36d65e05149cfb2f6ce5`.
Independent evaluation of the promoted checkpoint parameters gives
`B(0,0)=-0.46582943379545416`; the test checks a tight validated enclosure
`[-0.466,-0.465]` before running the same one-box Barr3-domain equation
workload used for the 8x8 fixture. This keeps the timing comparison controlled:
only the network width changes.

The third scaling fixture is `2-32-32-1`, using the same
prefix rule. Its concatenated raw float32 tensor bytes have SHA-256
`641b4799849f0057a11000a90b6c0910acc280dc0bd97777ca3142ff3089fecf`.
Independent evaluation after exact float32-to-DP promotion gives
`B(0,0)=1.1903679778286027`; the regression checks the validated enclosure
`[1.190,1.191]` before running the same one-box Barr3-domain equation.
For this larger fixture the source float32 bit patterns are stored directly as
`uint32_t` and converted with C++20 `std::bit_cast<float>`, then promoted to
DP. This is representation-equivalent to hexadecimal double literals while
keeping the fixture substantially smaller.

The measured Debug+coverage solve times are 0.195 s for `8x8`,
0.767 s for `16x16`, and 3.028 s for `32x32`. The two successive ratios
are about 3.93 and 3.95, closely tracking the fourfold growth of the dense
second-layer connection count. Up through `32x32` there is therefore no
evidence of additional super-quadratic scaling from the current symbolic tanh
expansion or SMT machinery.

The exact published `2-64-64-1` model is now included as the fourth scaling
point. Rather than embedding another very large C++ initializer, the six
float32 tensors are stored internally as the raw little-endian tensor payload
`tests/solvers/data/smt_barr3_full64.bin`. The payload is derived directly
from the verified checkpoint by concatenating `W1,b1,W2,b2,W3,b3` in that
order. Its SHA-256 is
`d1c6646a9354d44092e23de07495c40d2a6575f2235aaa5be0d68b1e9b61bdbf`
and its size is 17668 bytes, corresponding to exactly 4417 float32
parameters. The original checkpoint SHA-256 remains
`788fd56f21abb0cb5b83206b12d4d8d0ed718896c2140222a8fd129475f6c213`.

`smt_barr3_full64.hpp` parses that internal payload directly, reconstructs
each IEEE-754 binary32 value from its bytes, promotes it exactly to DP and
builds the same two-layer tanh `RealExpression`. No PyTorch, ONNX or external
framework is used at test runtime. The data path is supplied only to the
dedicated benchmark executable through a CMake compile definition.

Independent exact float32-to-DP evaluation gives
`B(0,0)=3.809158158082951`; the full-model regression checks a validated
origin enclosure `[3.809,3.810]` before executing the same one-box
`B(x,y)=0` workload as the three prefix fixtures. This preserves the scaling
comparison while making the `64x64` point the actual published network rather
than a truncated surrogate.


### Published Barr3 verification workload

The repository now also contains the standalone executable
`benchmark_smt_barr3_verification`. It reproduces the two actual counterexample
queries used by the public FOSSIL `BarrierAlt` verifier for Barr3, rather than
the synthetic `B(x,y)=0` scaling query.

The source benchmark defines
`XD=[-3,2.5]x[-2,1]` and the unsafe set as the union of a radius-0.4 sphere
centred at `(-1,-1)`, the rectangle `[0.4,0.6]x[0.1,0.5]`, and the rectangle
`[0.4,0.8]x[0.1,0.3]`. The current public `BarrierAlt.get_constraints()`
checks only two counterexample formulae: `XU && B>=0`, and
`XD && B>=0 && Bdot+B<0` with alpha=1. The initial-set condition is present
in the benchmark data but commented out in that verifier implementation, so the
Ariadne reproduction deliberately does not invent an initial-set query.

The Lie derivative is built analytically by forward-propagating
`d/dx tanh(z)=(1-tanh(z)^2) dz/dx` and the analogous y derivative through the
two dense layers, then applying the published dynamics
`x_dot=y`, `y_dot=-x-y+x^3/3`. This avoids asking the generic symbolic
differentiator to expand the already large network expression while remaining
mathematically identical to the network derivative.

FOSSIL's default dReal precision is `1e-5`; the standalone Ariadne benchmark
therefore uses epsilon `1e-5`. Candidate witness search and monotone
contraction are enabled. To avoid accidentally turning every development run
into a potentially long verification, the executable defaults to a one-box
smoke budget. A numeric first argument sets the per-query box-processing limit;
the literal argument `full` removes the box budget. The intended progression
is therefore e.g. `benchmark_smt_barr3_verification 16`, then larger budgets,
and only then `benchmark_smt_barr3_verification full`.

The executable reports expression-construction time separately from each unsafe
component and the Lie-query solve time, plus the full SMT statistics and UNKNOWN
reason. It is intentionally not registered with CTest; the ordinary neural
scaling regressions remain deterministic and bounded.

The first full-model attempt exposed an important benchmarking pitfall: with a
nominal box limit of one, the run was still interruptible only after more than
18 minutes, while constructing the network and analytic Lie expression took
only about 0.002 s in Release. A box budget limits completed calls to
`_process_box`; it does not bound the work performed *inside* one box. The
previous benchmark simultaneously routed the unsafe union through DPLL(T),
enabled nonlinear candidate witness search, and enabled derivative-assisted
monotone contraction. Any of those per-box/theory operations may dominate
before the outer box counter can advance.

The diagnostic benchmark therefore now removes those confounders before any
solver change is considered. The unsafe union is decomposed exactly into its
three components and solved as direct primitive conjunctions; rectangle bounds
are represented as query boxes, while the spherical component uses its tight
bounding box plus the sphere inequality. Candidate witness search and monotone
contraction are disabled. Each query prints a flushed `[start]` marker before
entering the solver. This establishes the baseline cost of compilation,
hull/shaving propagation, epsilon checks and sensitivity splitting on one real
64x64 Barr3 box. Features are to be re-enabled one at a time only after this
baseline is measured.

The Release one-box baseline after this isolation is 1.642 s for the unsafe
sphere, 1.768 s and 1.757 s for the two unsafe rectangles, and 15.962 s for
the Lie query. Network+Lie expression construction remains only 0.002 s.
The Lie query is therefore roughly nine times more expensive per box than a
barrier-only query even before candidate search or monotone contraction.

Inspection of theory compilation identified one source of avoidable work:
`_compile_theory_literals()` eagerly constructed every coordinate derivative
for every theory literal even when monotone contraction was disabled. Those
cached derivatives are consumed only by the monotone-contraction path; hull,
shaving, epsilon checking and sensitivity splitting use the compiled functions
directly. Theory compilation now leaves the derivative cache empty unless
`monotone_reduction_enabled()` is true. A dedicated regression checks that a
two-variable literal produces zero cached derivatives with monotone reduction
off and two with it on. This is a semantics-preserving optimization: it changes
only unused precomputation, not the enabled contractor set.


Because that optimization changed the one-box timings only modestly, the solver
now records diagnostic phase timings in `SmtSearchStatistics`: theory
compilation, original reduction, epsilon box checking, midpoint/end-point
witness probing, sensitivity splitting, and nonlinear candidate search. These
are accumulated across processed boxes and are diagnostic only; no solver
decision depends on wall-clock time. The Barr3 verification runner prints all
of them so the expensive phase can be identified before changing contractor
policy.

A subsequent Release one-box run measured 1.479 s for the unsafe sphere,
1.564 s and 1.566 s for the two unsafe rectangles, and 14.387 s for the Lie
query. The Lie-query phase breakdown was 3.240 s theory compilation, 4.783 s
original reduction, 0.020 s whole-box epsilon checking, 1.476 s deterministic
witness probing, and 4.763 s sensitivity splitting. This rules out epsilon
containment checking as the dominant cost and localizes most of the remaining
one-box expense to propagation and sensitivity analysis.

The diagnostic statistics now split sensitivity cost further into symbolic
derivative construction time and validated derivative-evaluation time, with
counts for both operations. Constraint propagation also counts the validated
function evaluations performed by coordinate shaving. These counters are
observational only and do not alter contractor or split semantics. The next
optimization decision is deliberately deferred until the full Barr3 one-box
benchmark identifies whether sensitivity time is dominated by derivative
construction or derivative evaluation, and how many expensive neural-function
evaluations shaving performs.

The first run with those counters resolves that question. On the Lie query,
sensitivity performs four derivative constructions in 0.560 s but four
validated derivative evaluations in 4.052 s; derivative evaluation, not symbolic
derivative construction, is therefore the dominant sensitivity cost. The same
box performs eight shaving function evaluations while the complete reduction
phase costs 4.940 s. Barrier-only rectangles show the same qualitative pattern:
two derivative builds take about 0.07 s while two derivative evaluations take
about 0.44 s. Caching symbolic derivatives can still amortize construction over
multiple boxes, but it cannot address the dominant per-box cost seen here.

Before changing split policy, the benchmark also reports whether sensitivity
actually overrides the geometrically widest coordinate. If the expensive
validated derivative analysis does not change the selected coordinate on the
target workload, retaining it unconditionally would be pure heuristic overhead;
if it does, its search benefit must be measured against that per-box cost before
replacing it.

The override diagnostic shows `sensitivity-overrides=0` for both unsafe
rectangles and for the Lie query, while the unsafe sphere reports one override.
The sphere has equal coordinate widths, so sensitivity is acting as an expensive
tie-breaker there; on the full Lie domain the geometrically widest coordinate
already matches the sensitivity choice despite about 5 seconds of sensitivity
work per box in the measured Release run.

To measure search-quality tradeoffs before changing the default heuristic,
`SmtSolverConfiguration` now exposes an opt-out for sensitivity-guided
splitting. It defaults to enabled, preserving existing solver behavior. The
standalone Barr3 benchmark accepts an optional second argument,
`sensitivity` or `geometric`, so identical box budgets can be compared with
validated sensitivity analysis enabled or with the solver's geometric box split
only. This switch is diagnostic and does not weaken result semantics: it changes
only the search heuristic.

At an eight-box budget the comparison shows no pruning or epsilon certification
under either policy. The Lie query takes 96.985 s with sensitivity and 56.706 s
with geometric splitting; sensitivity itself accounts for 40.277 s, while both
runs process eight boxes and return the same resource-exhaustion outcome. The
geometric run therefore removes about 42 percent of elapsed time at this budget
without losing observable search progress. This is evidence against paying for
full derivative-based sensitivity unconditionally on the Barr3 workload, but it
is not yet sufficient to change the solver default.

With sensitivity removed, deterministic midpoint/endpoint/corner witness
probing becomes the next large heuristic cost: 13.859 s of the 56.706 s Lie
geometric run, with no witness found. The configuration therefore also exposes
an opt-out for deterministic witness probing, enabled by default to preserve
existing behavior. The benchmark accepts a third argument, `witness` or
`no-witness`, so the pure geometric branch-and-prune baseline can be measured
before increasing the box budget. Disabling probing changes only when heuristic
point candidates are attempted; validated whole-box epsilon certification
remains active.

The eight-box geometric/no-witness run reduces the Lie query further to
42.847 s. Of that, 39.241 s are original reduction; the run performs 64 shaving
function evaluations and still records no pruned box or epsilon certification.
The unsafe rectangles likewise spend about 4.81 s of roughly 5.21 s total in
reduction. With sensitivity and deterministic point probing removed, coordinate
shaving is therefore the dominant measured cost.

A further diagnostic switch disables coordinate shaving while retaining hull
propagation, direct validated range rejection and whole-box epsilon
certification. Shaving remains enabled by default. This switch is implemented
at the generic `ConstraintSolver::propagate` boundary rather than by duplicating
propagation in SMT; both validated-constraint and precompiled-constraint paths
are tested with shaving disabled. Monotone contraction remains independently
configurable and, when enabled, is still reached after a hull stall even if
shaving is disabled.

The no-shaving eight-box run leaves 33.855 s of the 37.417 s Lie query inside
original reduction. Inspection shows that the scalar-function `hull_reduce`
overload constructs a new `ValidatedProcedure` from the unchanged validated
function on every call before invoking `simple_hull_reduce`. The benchmark now
separates this procedure-construction time from the hull contraction itself and
from the subsequent direct validated range rejection. This will determine
whether the next optimization should precompile procedures once per theory
literal rather than rebuilding their common-subexpression-eliminated instruction
DAG for every box.

The next diagnostic separates the measured hull-contraction interval further
into temporary-storage allocation, forward procedure execution and backward
constraint propagation. This instrumentation lives in `simple_hull_reduce`
itself and leaves its arithmetic and traversal order unchanged. It distinguishes
the cost of evaluating the Barr3 procedure from the inverse interval operations
performed by the backward contractor before either path is optimized.

The detailed no-shaving run attributes the Lie hull contraction almost entirely
to forward procedure execution: 26.497 s of 28.036 s, versus 1.511 s for
backward propagation and 0.025 s for temporary allocation. Rebuilding the
procedure costs another 4.783 s, while direct validated rejection costs only
1.367 s. Thus the dominant current cost is validated forward execution through
the `Procedure` representation, not inverse/backward contraction.

Before optimizing that execution engine, the benchmark can now disable hull
reduction while retaining direct validated range rejection. Hull remains enabled
by default. This experiment measures whether the expensive contractor produces
enough pruning/contraction on Barr3 to justify optimizing it, or whether a
cheaper/adaptive contractor schedule should be preferred. Both generic and
precompiled propagation paths retain direct rejection when hull is disabled and
have explicit tests for that behavior.

The eight-box no-hull run makes the tradeoff concrete. The Lie query drops from
37.832 s with hull enabled (and shaving, sensitivity and witness probing already
disabled) to 4.975 s with hull disabled. Both runs process eight boxes, prune
zero Lie boxes, certify zero epsilon boxes and perform eight splits. The sphere
also retains the same one pruned box while dropping from 4.504 s to 0.387 s.
Thus, at this budget, hull contraction adds no observable pruning progress on
Barr3 while dominating runtime by more than an order of magnitude over direct
validated rejection.

This does not yet justify removing hull contraction from the solver default:
contractor value can emerge only after further subdivision, and other workloads
may benefit substantially. The next Barr3 experiment therefore increases the
box budget using the cheap geometric/no-witness/no-shaving/no-hull baseline to
locate where direct validated range rejection starts pruning or epsilon
certification starts succeeding. Hull should then be reintroduced at the same
box frontier to measure whether its extra contractions reduce the remaining
search enough to amortize their cost.

At 64 boxes, the no-hull geometric baseline has entered a useful pruning
regime. The sphere prunes 29/64 boxes, each unsafe rectangle prunes 27/64, and
the Lie query prunes 22/64 while splitting 42. The Lie query takes 18.618 s,
including 3.343 s one-time theory compilation, 10.860 s direct rejection and
4.319 s epsilon checking. Thus plain subdivision plus validated range rejection
is already eliminating about one third of the processed Lie boxes without the
expensive hull contractor.

The standalone benchmark now accepts an optional final query selector, `all`
or `lie`, so higher-box-budget experiments can focus on the dominant Lie query
without repeatedly paying for the three unsafe-set components. The default
remains `all`.

At 256 boxes, the Lie-only cheap baseline processes 256 boxes, prunes 119 and
splits 137, so the pruning fraction rises to about 46.5 percent from about 34.4
percent at 64 boxes. Elapsed time is 67.241 s: 3.364 s one-time theory
compilation, 43.484 s direct validated rejection, and 20.294 s whole-box epsilon
checking. No epsilon-certified box has appeared yet. The increasing pruning
fraction is evidence that geometric subdivision is improving direct interval
decision power rather than merely postponing the same unresolved work.

The next measurement continues this budget curve on the Lie query before
introducing adaptive contractors. Explicit search-depth instrumentation is
deliberately deferred because the parallel workload currently carries boxes
without depth metadata; adding sequential-only depth statistics would make the
public statistics inconsistent. The budget/pruning curve already measures the
relevant search effect without changing search-state representation.

At 512 boxes the Lie-only cheap baseline prunes 248/512 boxes (48.4 percent)
and splits 264, suggesting the pruning fraction is approaching a plateau near
one half. Runtime is 130.651 s: 86.169 s direct rejection and 41.033 s whole-box
epsilon checking after 3.352 s compilation. The epsilon phase therefore becomes
the next avoidable duplicate cost: on this contractor-free path it reevaluates
the same unchanged literal functions immediately after direct rejection.

The compiled-theory box processor now has a fused direct-classification path
used only when hull, shaving and monotone reduction are all disabled. Each
literal image is evaluated once and that same validated image is used both for
original infeasibility (including strict endpoint semantics) and whole-box
epsilon certification. If any box-modifying contractor is enabled, the existing
separate reduction and post-reduction epsilon evaluation remain unchanged. Tests
cover original pruning, epsilon certification, an unresolved split, strict `>0`
endpoint rejection, fast-path gating, and statistics aggregation.

The 512-box Lie rerun confirms the optimization is behavior-preserving on the
Barr3 baseline: `pruned=248` and `split=264` are unchanged, while
`fused-direct=512` confirms the fused path classified every processed box.
Whole-box epsilon-check time falls from 41.033 s to 0 because the epsilon decision
now reuses the direct-rejection image, and total elapsed time falls from
130.651 s to 89.703 s (about 31 percent). The remaining measured runtime is now
dominated by the single validated function evaluation per literal used by direct
classification: 86.248 s after 3.363 s one-time theory compilation.

Before changing the evaluator, the standalone Barr3 benchmark now provides an
`eval` mode that compares the current validated function `apply` against a
precompiled `ValidatedProcedure` evaluated forward-only on the same full Barr3
box. Procedure construction is timed separately. This experiment determines
whether procedure precompilation is a viable replacement for the dominant cheap
classification evaluation or whether the existing function evaluator is already
the better execution path. No solver semantics or default configuration change
is made by this diagnostic.

The evaluator comparison decisively rejects `ValidatedProcedure` as a cheap-path
replacement on Barr3. Across eight evaluations of the full input box, the barrier
term takes 0.163 s with normal validated function `apply` versus 3.399 s with
forward-only `ValidatedProcedure` evaluation, while `lie+barrier` takes 1.214 s
versus 23.935 s. This is roughly a 20x slowdown in both cases, even before
amortizing procedure construction (0.073 s and 0.617 s respectively). Both
methods return the same uninformative `[-inf,+inf]` image on the full box, so the
extra procedure cost buys no additional pruning power. Procedure evaluation is
therefore excluded from the cheap-first path; further performance work should
focus on the existing validated function evaluator and on preserving/sharing
expression structure across box evaluations.

Before further performance work, coverage must be returned to 100 percent for
all newly introduced configuration branches, diagnostics and the fused direct
classification behavior. In particular, coverage should confirm both sides of
the fast-path gate and the strict/non-strict classification branches rather than
merely exercising the successful Barr3 path.

## Current open work

The immediate work on `solvers-smt#830` is:

1. validate the explicit UNKNOWN-reason plumbing across sequential, parallel and
   Boolean search and restore full coverage for the new branches;
2. [completed] splitting focuses sensitivity analysis on constraints whose
   validated box evaluation has not yet met the epsilon stopping condition,
   while retaining Ariadne's sensitivity guidance;
3. investigate rigorous DP stopping certificates based on validated continuity
   or derivative information that can reduce `DP_RESOLUTION_EXHAUSTED` without
   arbitrary-precision arithmetic;
4. preserve deterministic, nontrivial tests for every introduced behavior;
5. after every green functional test run, regenerate coverage and restore 100%
   function and branch coverage for newly introduced functionality;
6. improve DP decision power through symbolic simplification,
   correlation-preserving/shared expression evaluation and stronger contractors;
7. retain the audited compiler-mapped line anomaly unless a semantic source
   change resolves it naturally;
8. later add CI coverage gates so regressions in functions, lines or branches
   fail automatically.


### Coverage note: residual LLVM line-report anomaly

The current coverage run leaves a single source line in the parallel SMT worker
reported as uncovered even though the associated parallel solve regression
passes and the surrounding worker logic is exercised. Repeated attempts to
force that exact source location through additional tests did not make the line
appear covered, while the semantic behaviour of the parallel path remained
verified.

This residual is therefore treated as a coverage-instrumentation/toolchain
anomaly, most likely related to LLVM source-line mapping/instrumentation, rather
than as evidence of an untested SMT feature. The project does not add further
artificial tests solely to satisfy that one reported line. Coverage requirements
continue to apply to functional branches and behaviour introduced by the SMT
solver work; this specific line is documented as an explicit exception until a
toolchain-level explanation or reproducible LLVM fix is available.


The evaluator diagnostics were then extended in three directions. First, the
same ValidatedProcedure tape was evaluated on FloatDPBounds rather than
UpperIntervalType; second, the originating RealExpression was converted to
Formula<EffectiveNumber> and evaluated on both interval backends; third, the
benchmark added a Barr3-only reference evaluator that executes the published
2-64-64-1 network directly on FloatDPBounds, propagating the two input
derivatives forward through the network and evaluating tanh with a sign-stable
interval formula. The direct evaluator is diagnostic only and is not used by
the solver.

On the full initial box, eight Lie-plus-barrier evaluations take about 1.257 s
through the current validated-function apply, 0.813 s through the precompiled
procedure on FloatDPBounds, 1.266 s through Formula on FloatDPBounds, but only
0.014 s through the direct network evaluator. The corresponding barrier times
are 0.162 s, 0.106 s, 0.158 s and 0.014 s respectively. Thus replacing apply by
the existing procedure tape alone yields only a modest improvement; it does not
explain the roughly two-orders-of-magnitude gap exposed by the network-specific
reference path.

The direct evaluator also returns finite enclosures on the initial box:
approximately [-24.4,24.5] for the barrier and [-971.8,970.6] for
Lie-plus-barrier. All three generic expression/function paths instead return an
unbounded interval. The published network fixture currently represents tanh(z)
as (exp(2z)-1)/(exp(2z)+1); on wide interval arguments that representation can
overflow and lose all range information even though the mathematical activation
is bounded. The direct diagnostic uses equivalent sign-stable formulas on
one-sided intervals and the conservative range [-1,1] when an interval crosses
zero. This establishes a real numerical representation problem, but the
observed speedup must not be attributed to stable tanh alone: the direct
evaluator also propagates the barrier gradient compactly in forward mode instead
of evaluating the very large pre-expanded symbolic Lie expression.

The structural diagnostic reinforces that distinction. The Lie-plus-barrier
expression reports 1,856,788 recursive node visits but only 31,568 distinct
expression-node pointers; the barrier reports 263,169 visits but only 14,019
distinct pointers. Procedure conversion already memoizes Formula node pointers,
so the next diagnostic prints the resulting procedure instruction count.
Comparing that count with the distinct-node count will determine whether shared
subexpressions are lost before procedure construction or whether the remaining
cost is primarily generic instruction dispatch and the expanded
activation/derivative representation. No core tanh operator or specialized
Barr3 execution path should be introduced until this distinction is measured.


The procedure-instruction diagnostic then identified a generic DAG-preservation
bug outside the SMT solver. For the Barr3 barrier, the original expression has
263,169 recursive node visits but only 14,019 distinct node pointers, while the
procedure built from the expression-space conversion contained 263,169
instructions. For Lie-plus-barrier the corresponding counts were 1,856,788
visits, 31,568 distinct pointers and 1,856,788 instructions. The procedure
converter itself already caches Formula node pointers; inspection showed that
the public make_formula(Expression<Real>, Space<Real>) overload deliberately
called the internal no-cache conversion, unlike the Map<Identifier,SizeType>
overload, which already uses the node-pointer cache.

The Space<Real> overload now delegates to the cached Map overload. A structural
Procedure regression constructs an expression with an explicitly shared
subgraph and checks that conversion produces one instruction per distinct
expression node rather than one per recursive occurrence, while preserving the
evaluated value. This is a general symbolic-function fix rather than an
SMT-specific optimization. The Barr3 evaluator benchmark should be rerun before
adding a native tanh operator: if procedure instruction counts collapse toward
the distinct-node counts and evaluation time falls accordingly, DAG preservation
must be accounted for separately from the remaining numerical instability of
the exp-ratio tanh representation.


The DAG-preservation change was validated by the symbolic expression, Procedure,
and SMT solver test suites. The Barr3 diagnostic run immediately before the
naming-only refactor also confirms the expected structural effect: the
expression-derived procedure now contains exactly 14,019 instructions for the
barrier and 31,568 for Lie-plus-barrier, matching the respective distinct
expression-node-pointer counts rather than the 263,169 and 1,856,788 recursive
node visits. On FloatDPBounds, eight evaluations take about 0.005 s and 0.015 s
respectively, compared with 0.106 s and 0.741 s for procedures produced through
the sharing-losing function conversion measured before the fix. The compact
Lie-plus-barrier procedure is therefore already within measurement noise of the
Barr3-specific direct reference evaluator (0.013 s for eight combined network
evaluations).

This separates the two previously conflated issues. DAG preservation accounts
for essentially all of the large generic-evaluation performance gap. The
remaining material problem is enclosure quality: the generic expression paths
still return an unbounded interval on the initial Barr3 box, while the stable
direct network evaluation returns finite enclosures. Work can therefore move to
a general stable representation of tanh without using a Barr3-specific execution
path as a solver optimization.


### Native tanh validation stage

A general `Tanh` unary elementary operator has now been added across the
numeric, algebraic, symbolic and Procedure layers. Its validated
`Bounds<F>` implementation evaluates the monotone endpoint image with
sign-stable formulae rather than the dependency-prone
`(exp(2*x)-1)/(exp(2*x)+1)` quotient. Support was also completed for the
rounded and geometric interval types used by differential algebra, for
`Graded` propagation through the identity `y' = x'*(1-y^2)`, and for the
symbolic affine/polynomial classifiers. The focused numeric, expression and
Procedure regressions are green.

The published Barr3 fixture now represents each hidden activation with the
native `tanh` expression. Its explicit forward derivative propagation is
otherwise unchanged, so this experiment isolates activation representation
from the already-resolved expression-DAG sharing issue. The fast Full64 neural
smoke test and the standalone evaluator benchmark must be rerun on this exact
revision before drawing conclusions about enclosure quality or performance.
The key comparison remains the initial full-box `eval` profile: generic
function, Formula and expression-derived Procedure images should be compared
with the Barr3-only direct reference evaluator, with particular attention to
whether the generic barrier and Lie-plus-barrier images become finite.

The first Barr3 run with native `tanh` confirms that activation representation
was the enclosure blocker. On the full initial box, all generic paths now return
finite intervals. The barrier image is approximately
`[-24.323104,24.437737]` through the validated-function and Formula paths and
`[-24.4,24.5]` through FloatDPBounds Procedure evaluation, matching the
Barr3-specific direct reference `[-24.4,24.5]` to the expected backend
rounding granularity. Eight compact expression-derived Procedure evaluations
take about 0.028 s for the barrier, versus 0.034 s for the complete direct
Barr3 evaluator.

For Lie-plus-barrier the generic image is finite but still substantially wider:
approximately `[-3244.6,3239.4]` on FloatDPBounds versus
`[-971.8,970.6]` in the direct reference. Inspection of the two otherwise
equivalent derivative propagations shows that the fixture computes the
activation derivative factor as `1-h*h`, while the direct evaluator uses
`1-sqr(h)`. For interval-valued `h` crossing zero, generic multiplication
loses the self-correlation that `sqr(h)` preserves. The fixture now uses
`1-sqr(h)` in both hidden layers. This algebraically equivalent representation
change should be benchmarked before pursuing stronger contractors or evaluator
changes.

The rebuilt Barr3 run after replacing `1-h*h` by `1-sqr(h)` confirms that
the change is present structurally: Lie-plus-barrier recursive expression visits
fall from 497,492 to 308,820 while the distinct-node-pointer count remains
30,672. On FloatDPBounds the Lie-plus-barrier enclosure contracts from about
`[-3244.6,3239.4]` to about `[-889.3,888.3]`. The native-tanh barrier remains
about `[-24.4,24.5]`.

Absolute timing values from this run must not be compared with the immediately
preceding measurements because they were collected on different machines
(Mac Studio versus MacBook). Only within-run timing ratios and structural or
enclosure data are comparable across those reports. On the Mac Studio run,
eight expression-derived Procedure evaluations take about 0.015 s for
Lie-plus-barrier versus 0.013 s for the Barr3-specific direct evaluator, again
showing that DAG-preserved generic Procedure execution is in the same cost
class as the hand-written forward evaluator.

The direct diagnostic previously used a deliberately conservative tanh helper
that returned `[-1,1]` whenever its input crossed zero. Since the native
validated tanh now computes a monotone endpoint image, that old helper makes the
direct reference artificially wider than the generic path. The diagnostic now
uses the native `tanh(FloatDPBounds)` as well, so the next run compares
execution structure rather than two different activation enclosures.

A subsequent same-machine rerun showed that replacing the direct diagnostic's
cross-zero tanh fallback by native validated tanh did not change its final
interval: the direct Lie-plus-barrier enclosure remained approximately
`[-971.8,970.6]`, while the expression-derived Procedure remained tighter at
about `[-889.3,888.3]`. The earlier attribution of that residual difference
to the direct tanh fallback was therefore incorrect.

The benchmark now reports the barrier-gradient components `db/dx`, `db/dy`
and the unfused Lie term separately for both the compact expression-derived
Procedure and the hand-written direct evaluator. This diagnostic is intended to
identify the first stage at which the otherwise algebraically equivalent
evaluations diverge, before making any further solver or arithmetic change.

The intermediate Barr3 diagnostic now closes the residual enclosure question.
On the full initial domain, the compact expression-derived Procedure and the
hand-written direct evaluator produce the same FloatDPBounds images at every
reported stage: `db/dx=[-60.4,60.9]`, `db/dy=[-57.7,59.5]`,
`lie=[-864.9,863.9]`, and `lie+barrier=[-889.3,888.3]`. The barrier image
also agrees at `[-24.4,24.5]`.

This establishes that, after native tanh, DAG-preserving expression conversion,
and use of `sqr(h)` in the activation derivative, the generic symbolic
Expression -> Formula -> Procedure path is no longer losing numerical quality
relative to the Barr3-specific direct forward evaluator. On the same Mac Studio
run, eight Lie-plus-barrier evaluations take about 0.015 s through the compact
generic Procedure and about 0.015 s through the direct evaluator; this
same-machine comparison is meaningful, whereas absolute timings from earlier
MacBook runs are not directly comparable.

The evaluator/enclosure investigation is therefore complete. Further Barr3 work
should return to actual SMT decision power: rerun the cheap
geometric/no-witness/no-shaving/no-hull Lie query at increasing box budgets on a
single machine, record pruning and split counts, and only then decide whether
stronger contractors or search heuristics are still warranted.

The cheap Barr3 Lie-query budget curve was rerun on a single Mac Studio after
the evaluator fixes, using geometric splitting with witness probing, shaving and
hull reduction disabled. At 64 processed boxes the solver prunes 27 and splits
37; at 256 boxes it prunes 123 and splits 133; at 512 boxes it prunes 251 and
splits 261. No whole-box epsilon certification occurs at any of these budgets,
and the fused direct-classification path handles every processed box.

Compared structurally with the earlier pre-native-tanh curve, pruning improves
from 22/64 to 27/64, from 119/256 to 123/256, and from 248/512 to 251/512.
Thus the stronger activation/derivative enclosure improves early direct
classification, but the gain diminishes under subdivision and the pruning
fraction still approaches roughly one half. The remaining difficulty is
therefore no longer an evaluator-quality defect; it is search/constraint
decision power on the unresolved half of the Lie domain.

Absolute runtimes from older MacBook measurements must not be compared with
these Mac Studio timings. Within the Mac Studio curve, elapsed time scales
approximately linearly with processed boxes: about 2.77 s at 64 boxes, 9.24 s at
256 boxes and 18.91 s at 512 boxes, with almost all post-compilation time spent
in fused direct validated rejection. This makes larger-budget geometric runs
cheap enough to characterize the unresolved frontier before reintroducing any
contractor.

The geometric-only Lie curve continues to 1,024 and 2,048 processed boxes with
508 and 1,020 pruned boxes respectively, i.e. 49.6% and 49.8%. No epsilon box
is certified. Together with the 512-box point, this confirms a stable pruning
plateau near one half; further brute-force geometric budget increases are not
expected to change the qualitative result.

The next controlled experiment re-enables only validated monotone/Newton
contraction while keeping sensitivity splitting, deterministic witness probing,
coordinate shaving and hull reduction disabled. The standalone benchmark now
accepts an optional seventh argument, `monotone` or `no-monotone`, after the
query selector; the default is `no-monotone`, so existing benchmark commands
are unchanged. This isolates whether sign-definite cached derivatives can
contract the unresolved frontier cheaply before reconsidering the much more
expensive hull contractor.

The isolated monotone/Newton experiment is not competitive on Barr3 at the
measured budgets. With all other expensive features disabled, enabling monotone
reduction leaves the 64-box result unchanged at 27 pruned and 37 split boxes,
but increases elapsed time from about 2.77 s to 24.44 s. At 256 boxes it leaves
the result unchanged at 123 pruned and 133 split boxes while increasing elapsed
time from about 9.24 s to 99.10 s. The additional cost appears almost entirely
inside reduction.

The benchmark now prints the existing monotone-effective-reduction counter as
well as monotone rounds. A final small-budget rerun can therefore distinguish
between a contractor that changes boxes without affecting aggregate pruning and
one that performs no effective contraction at all. If the effective count is
zero, monotone reduction should be excluded from the Barr3 cheap-first path.


The final requested 64-box monotone diagnostic closes that experiment
decisively. With geometric splitting and witness probing, shaving and hull
reduction disabled, the Lie query again processes 64 boxes, pruning 27 and
splitting 37, but reports 37 monotone rounds and zero effective monotone
reductions. Elapsed time is 25.347 s, of which 23.514 s is reduction. Thus the
validated sign-definite derivative checks and Newton attempts do not contract a
single processed Barr3 box at this frontier; the extra cost cannot be justified
by latent box contraction that merely failed to change the aggregate pruning
count. Monotone/Newton reduction is therefore excluded from the Barr3
cheap-first path.

The next controlled experiment reintroduces only hull reduction at the same
64-box frontier, keeping geometric splitting, deterministic witness probing,
coordinate shaving and monotone reduction disabled. This directly compares
hull contraction against the established contractor-free 64-box result
(27 pruned, 37 split after the evaluator fixes). The benchmark now also prints
the existing hull-effective-reduction counter, so the run distinguishes useful
box contraction from pure contractor overhead even if aggregate pruning is
unchanged. If hull is again ineffective or its contractions do not materially
alter pruning at this frontier, the next development step should move away from
eager per-box contractors and toward a cheaper adaptive trigger or a different
search heuristic.


The 64-box hull-only comparison rejects eager hull contraction on Barr3. With
geometric splitting and witness probing, shaving and monotone reduction
disabled, enabling hull changes the Lie result only from 27 pruned / 37 split
boxes to 28 pruned / 36 split boxes, while elapsed time rises from about
2.77 s to 44.439 s. The run reports 64 hull rounds, zero effective non-empty
hull contractions, 128 procedure builds costing 5.841 s, and 35.503 s in hull
contraction itself. Of that contraction time, 33.670 s is forward Procedure
execution and 1.799 s is backward propagation. Eager hull therefore costs about
an order of magnitude more than the contractor-free path for one additional
pruned box at this frontier. The zero hull-effective count does not imply that
hull never contributed to rejection: a hull contraction that directly empties
a box returns before the non-empty contraction counter is incremented. It does
show that no surviving box was narrowed for later search. Hull is therefore
excluded from the Barr3 cheap-first path alongside monotone/Newton reduction.

Before introducing another search heuristic, the cheap fused-classification
path now counts literal evaluations. This is diagnostic only. On the current
two-literal Lie conjunction the literals are ordered as barrier non-negativity
followed by the Lie violation. Because the fused path stops at the first
originally infeasible literal and no whole-box epsilon certification has yet
occurred on this workload, a 64-box contractor-free run can attribute pruning
without separate SAT queries: with E literal evaluations, the number of boxes
rejected by the first barrier literal is 2*64-E, and the remaining pruned boxes
were rejected only after evaluating the Lie literal. This identifies which
constraint drives the observed pruning plateau before paying for a new
coordinate-selection heuristic.


The 64-box fused-literal diagnostic attributes the cheap-path pruning completely
to the Lie-violation literal. The contractor-free run processes 64 boxes, prunes
27 and splits 37, with 128 fused literal evaluations. Because the conjunction
contains two literals ordered as barrier non-negativity followed by Lie
violation, evaluating both literals on all 64 boxes means that barrier
non-negativity rejects no processed box at this frontier. All 27 pruned boxes
are therefore rejected only by the second literal. The current pruning plateau
is specifically a Lie-expression decision problem rather than a barrier-domain
classification problem.

The next controlled experiment revisits sensitivity-guided splitting on the
current evaluator, with witness probing, shaving, hull and monotone reduction
disabled. The older sensitivity measurements were collected before the
DAG-preserving expression conversion, native tanh support and the subsequent
Barr3 expression fixes, so their absolute derivative-evaluation cost is no
longer a reliable basis for rejecting sensitivity on the current revision. A
64-box sensitivity run can now be compared directly with the current geometric
baseline (27 pruned, 37 split, about 2.835 s) while reporting sensitivity-guided
splits, overrides, derivative build time and derivative evaluation time. If the
updated sensitivity path remains expensive and rarely overrides the geometric
choice, it should be excluded from the Barr3 cheap-first schedule. If it changes
the split sequence materially at acceptable cost, its pruning gain can then be
measured at the same frontier before designing a new heuristic.


The current-evaluator 64-box sensitivity rerun closes that heuristic for Barr3.
With witness probing, shaving, hull and monotone reduction disabled, sensitivity
splitting produces exactly the same observed search outcome as geometric
splitting: 27 boxes pruned and 37 split, with no epsilon certification. All 37
splits are reported as sensitivity-guided but none overrides the geometrically
widest coordinate. Elapsed time rises from about 2.835 s for the geometric
baseline to 27.649 s. The sensitivity phase alone costs 24.866 s, including 94
symbolic derivative builds in 3.376 s and 94 validated derivative evaluations
in 18.576 s. Thus the post-DAG/native-tanh evaluator improvements do not change
the qualitative conclusion: full derivative-based sensitivity is expensive and
does not alter the Barr3 split sequence at this frontier. It is excluded from
the Barr3 cheap-first path together with eager hull and monotone/Newton
contraction.

At this point the cheap Barr3 path is intentionally minimal: fused direct
classification plus geometric splitting. The remaining pruning plateau is
entirely driven by the Lie-violation literal, and the expensive generic
contractors and derivative-guided split heuristic tested so far do not improve
the split sequence enough to justify their cost. Further work should therefore
target the Lie expression itself or a genuinely cheaper split signal, rather
than recombining the rejected eager mechanisms.


An opt-in interval-lookahead split policy is now available for the next Barr3
experiment. It does not build symbolic derivatives and it does not contract the
current box. For each splittable coordinate it bisects the box, evaluates the
currently epsilon-unresolved constraint functions on both candidate children,
and scores the coordinate by the sum of the resulting interval-image widths.
The coordinate with the smallest score is selected. Geometric splitting remains
the fallback and the production default is unchanged; the standalone Barr3
benchmark accepts `lookahead` as a third split-policy value.

The policy reports guided splits, overrides of the geometrically widest
coordinate, function-evaluation count and evaluation time separately. A
two-dimensional regression verifies that lookahead overrides a wider geometric
coordinate when the constraint depends only on the narrower coordinate. This
experiment is deliberately bounded: if the extra child evaluations do not
improve Barr3 pruning enough to compensate for their cost at 64 boxes, the
heuristic should be discarded rather than generalized further.

The 64-box interval-lookahead experiment is negative. It produces the same
27 pruned / 37 split boxes as geometric splitting and no epsilon certification.
Lookahead guides all 37 splits but overrides the geometrically widest coordinate
only three times. Those three changes do not alter aggregate search progress.
Elapsed time rises from about 2.835 s for geometric splitting to 8.483 s, with
5.701 s spent in the split phase. The heuristic performs 188 extra validated
function evaluations, accounting for 4.443 s of measured lookahead evaluation
time. Interval lookahead is therefore excluded from the Barr3 cheap-first path;
the small number of changed split decisions does not repay its evaluation cost.

Before designing another split heuristic, the benchmark now supports one
zero-semantics-change scheduling diagnostic for the fused two-literal Lie query:
an optional final argument selects `barrier-first` (the historical order) or
`lie-first`. The previous literal-attribution run established that the barrier
literal rejects none of the first 64 processed boxes, while the Lie literal
rejects 27. Therefore, if Lie is evaluated first, those 27 boxes can short-circuit
without evaluating the barrier. Under the same search tree the expected fused
literal-evaluation count is 64 Lie evaluations plus 37 barrier evaluations, i.e.
101 instead of 128. This test measures the available benefit from cheap literal
scheduling before implementing any generic adaptive ordering in the solver.

The 64-box literal-order experiment confirms the expected short-circuiting but
shows that its absolute value is modest. Evaluating the Lie-violation literal
before barrier non-negativity preserves the search outcome exactly at 27 pruned
and 37 split boxes, while fused literal evaluations fall from 128 to 101.
Elapsed time falls from about 2.835 s to 2.738 s and direct-classification time
from about 2.233 s to 2.130 s. Thus literal ordering is a valid cheap-path
optimization, but the barrier evaluation is inexpensive enough that the current
64-box gain is only a few percent.

Before promoting adaptive literal ordering into the solver core, the next
measurement repeats the established 2,048-box geometric frontier with Lie first.
Besides timing, the fused literal-evaluation count can attribute pruning at
depth: with B processed boxes and E evaluations, Lie-first short-circuits
`2*B-E` boxes before the barrier is evaluated. Comparing that number with total
pruned boxes shows whether barrier non-negativity starts contributing on the
deeper unresolved frontier. If it remains negligible, a simple stable
cost/effectiveness ordering may be worthwhile; if it becomes important, any
generic reordering policy needs adaptive evidence rather than a fixed Barr3
ordering.

The 2,048-box Lie-first run shows that barrier non-negativity remains irrelevant
to pruning even on the deeper established frontier. With 2,048 processed boxes
and two literals, a full evaluation of both literals would require 4,096 fused
literal evaluations. The run reports 3,076, so exactly 1,020 boxes short-circuit
after the first Lie-violation literal. This exactly matches the total pruned
count of 1,020. Therefore every pruned box is excluded by the Lie literal and
none requires barrier non-negativity for rejection. The search outcome remains
1,020 pruned / 1,028 split with no epsilon certification.

The benchmark now adds a diagnostic `lie-only` query containing only the Lie
violation literal. This is intentionally not treated as equivalent to the
original counterexample query: removing barrier non-negativity enlarges the
counterexample set, so UNSAT of the Lie-only query would prove the original
query but SAT/UNKNOWN would not decide it. The purpose is purely diagnostic.
If the Lie-only run follows the same pruning/split curve while reducing cost
substantially, then barrier non-negativity is not merely failing to prune first;
it is operationally irrelevant to the observed search frontier and should not
drive further Barr3 heuristic design.

The 2,048-box Lie-only diagnostic confirms that barrier non-negativity is
operationally irrelevant on the established cheap-search frontier. Removing the
barrier literal leaves the search outcome exactly unchanged at 1,020 pruned and
1,028 split boxes, with no epsilon certification. Fused literal evaluations
fall to exactly 2,048, one per processed box. Elapsed time decreases only from
64.911 s for the two-literal Lie-first query to 62.567 s for Lie-only, so the
barrier contributes little runtime and no observed decision power.

The next zero-evaluation-cost diagnostic varies only sequential DFS child order.
The existing stack visits the lower/first split child before the upper/second
child. An opt-in `upper-first` configuration reverses that visit order while
preserving the same geometric split boxes, solver semantics and box budget.
The benchmark accepts a final `lower-first` or `upper-first` argument and keeps
`lower-first` as the default. A test-support helper makes the stack push order
explicit and tests both branches. If the 2,048-box pruning fraction changes
materially under upper-first traversal, the apparent plateau is partly a
budget/frontier-order effect; if it remains near one half, search order alone
is not the missing decision mechanism.

The 2,048-box upper-first DFS diagnostic is also negative. Reversing child
visit order changes the Lie-only result only from 1,020 pruned / 1,028 split
to 1,021 pruned / 1,027 split boxes, with no epsilon certification. The pruning
fraction therefore remains essentially one half. The measured elapsed time
drops from 62.567 s to 59.847 s, but the search-progress difference is one box
and is not evidence of improved decision power. Sequential child order is
therefore not pursued as a Barr3 search mechanism.

Work now moves from scheduling heuristics back to enclosure strength. Ariadne
already provides validated multivariate Taylor function models, which preserve
polynomial correlation over a bounded box. The standalone Barr3 benchmark adds
a `taylor` diagnostic that constructs a
`ValidatedScalarMultivariateTaylorFunctionModelDP` for Lie-plus-barrier on the
initial domain using a threshold sweeper of 1e-8. It reports the ordinary
validated interval image and evaluation time, Taylor-model construction time,
and the model range. This is intentionally not integrated into SMT search yet.
The first question is whether the existing Taylor representation yields a
materially tighter finite enclosure at tolerable construction cost. Its generic
`tanh` path currently uses the normed-algebra exp-ratio identity, so this
diagnostic also reveals whether the numerical stability issue previously fixed
for interval evaluation reappears in Taylor arithmetic.

The first Taylor-enclosure diagnostic did not reach a range comparison. During
construction of the Lie-plus-barrier Taylor model it threw
`DivideByZeroException` inside reciprocal evaluation. The reported denominator
Taylor model had average about 6.17 but radius about 149.33, so its polynomial
enclosure spuriously crossed zero. Inspection confirms that validated interval
`tanh` had already been made sign-stable, but `NormedAlgebraOperations<Tanh>`
still used `(exp(2*x)-1)/(exp(2*x)+1)`. Thus this failure is a second instance
of the same numerical representation defect, now exposed in Taylor arithmetic;
it is not evidence that Taylor models themselves lack useful correlation.

`TaylorModel` now overrides `Tanh` for validated models. On a definitely
nonnegative range it uses `(1-exp(-2*x))/(1+exp(-2*x))`; on a definitely
nonpositive range it uses the symmetric `(exp(2*x)-1)/(exp(2*x)+1)`. Before
taking either reciprocal it verifies that the denominator Taylor range remains
strictly positive. If the input crosses zero, or if the modeled denominator
still loses positivity, the implementation returns the rigorous bounded model
`0 +/- 1`. Approximate Taylor models retain the existing normed-algebra path.
A focused Taylor-model regression exercises a wide cross-zero input as well as
positive and negative one-sided inputs. The standalone `taylor` Barr3
diagnostic should now be rerun before deciding whether Taylor enclosure is
worth integrating into SMT classification.

The rerun after stabilizing validated Taylor `tanh` is finite but decisively
negative for Barr3. On the initial domain the ordinary interval evaluation of
Lie-plus-barrier is approximately [-889.218, 888.270] and takes 0.032 s. The
validated Taylor function model takes 16.32 s to construct and its reported
range is approximately [-4789.778, 4785.336], more than five times wider than
the ordinary interval enclosure. Taylor-model construction therefore loses both
on enclosure strength and cost for this workload. It is not integrated into the
SMT classifier.

The Taylor-specific `tanh` patch also triggered `-Wshadow` because its local
`ModelType` typedef duplicated the existing typedef in
`AlgebraOperations<TaylorModel<P,F>>`; the redundant local typedef has been
removed.

One bounded follow-up tests the existing validated affine model before closing
the function-model route entirely. Ariadne's `affine_model(domain,function,dp)`
constructs a Taylor model with `AffineSweeper`, discarding non-affine terms into
the validated remainder as they arise, and converts the result to an affine
model. The standalone benchmark now accepts query `affine` and reports model
construction time and range on the same initial Lie-plus-barrier function. If
this is not substantially cheaper than the full Taylor model or does not beat
the ordinary interval enclosure, affine/Taylor function models should be
excluded from further Barr3 SMT work.

The validated affine-model follow-up is also negative. On the same initial
Lie-plus-barrier domain it takes about 0.622 s to build and reports a range of
approximately [-5097.330, 5099.938]. It is much cheaper than the 16.32 s full
Taylor model, but its enclosure is still roughly 5.7 times wider than the
ordinary interval image. The affine/Taylor function-model route is therefore
closed for Barr3: neither representation improves the enclosure used for
validated rejection.

The next isolated diagnostic is a first-order mean-value enclosure rather than
a function model or contractor. For a box X with midpoint m it computes the
validated enclosure `f(m) + sum_i D_i f(X) * (X_i-m_i)`. This preserves the
first-order dependency between the box displacement and interval gradient
without invoking Newton contraction, sensitivity-guided splitting or Taylor
polynomial construction. The benchmark query `mean-value` reports ordinary
interval evaluation, midpoint evaluation, derivative construction/evaluation
times and the resulting mean-value image on the initial Lie-plus-barrier box.
Only if this enclosure is materially tighter at acceptable cost should it be
considered for per-box SMT classification.

The initial-box mean-value enclosure is decisively worse than direct interval
evaluation. The direct Lie-plus-barrier image is approximately
[-889.218, 888.270] in 0.032 s. The mean-value construction evaluates the
midpoint in about 0.030 s, builds two derivatives in about 0.094 s and evaluates
them over the box in about 0.489 s, but returns approximately
[-48303.826, 48310.506]. The first-order derivative enclosure therefore
amplifies dependency by roughly two orders of magnitude and is excluded from
per-box SMT classification.

Rather than trying another generic enclosure representation, the next
diagnostic decomposes the Barr3 Lie expression itself. The published dynamics
are `dx=y` and `dy=-x-y+x^3/3`, so the benchmark query `lie-components` reports
validated interval images and evaluation times for `db/dx`, `db/dy`, both
dynamics, the products `db/dx*dx` and `db/dy*dy`, their sum, the barrier, and
the final Lie-plus-barrier expression. The goal is to locate the dominant
source of interval width before changing evaluator mathematics again.

The Lie-component decomposition localizes the initial-box interval width very
clearly. `db/dx` is about [-60.365, 60.896] and `db/dy` about
[-57.674, 59.448]. With the published dynamics, `dx=y` is [-2,1] while
`dy=-x-y+x^3/3` is [-12.5,12.5]. Consequently `db/dx*dx` contributes only
about [-121.792,120.729], whereas `db/dy*dy` contributes about
[-743.103,743.103]. Their sum is the Lie image [-864.894,863.832]; the barrier
adds only about [-24.323,24.438]. Thus the dominant source of width is the
`db/dy*dy` product, not the barrier or the x-direction term.

The next diagnostic therefore compares one geometric bisection in each state
coordinate using only the already identified dominant term and the final
Lie-plus-barrier expression. Query `lie-split-profile` evaluates both children
for an x split and for a y split, reporting each child image and the sum of
their widths. This is not yet a search heuristic: it is a single-root
measurement to determine whether the dominant-term structure suggests a cheap
coordinate preference that the ordinary widest-coordinate geometric split is
missing. If the y split materially reduces the dominant and final images, a
specialized cheap split score can be tested; otherwise split selection is not
where the missing correlation should be addressed.

The root split-profile confirms that geometric splitting already chooses the
better coordinate. Splitting x gives a width sum of about 1,912 for the dominant
`db/dy*dy` term and about 2,454 for Lie-plus-barrier, whereas splitting y gives
about 2,704 and 3,156 respectively. Since the initial x interval is also wider
than y, the existing geometric policy already makes this choice. This result
does not justify another split heuristic.

The remaining obvious dependency loss is inside the polynomial dynamics
`dy=-x-y+x^3/3` itself. Direct interval evaluation treats the repeated x
occurrences independently and produces [-12.5,12.5]. Before adding any
polynomial range machinery, the benchmark now tests the algebraically
equivalent Horner-like form `x*(sqr(x)/3-1)-y`. This uses the existing native
`sqr` operation to retain self-correlation in x. Query `lie-dynamics-rewrite`
compares the original and factored dy images, the corresponding `db/dy*dy`
images, and final Lie-plus-barrier images. The main Barr3 query remains
unchanged until this isolated comparison demonstrates a material improvement.

The dynamics rewrite diagnostic is strongly positive. On the initial Barr3
domain, the published polynomial dynamics `-x-y+x^3/3` evaluates to
[-12.5,12.5] in its original syntactic form, while the algebraically equivalent
`x*(sqr(x)/3-1)-y` evaluates to [-7,7]. The dominant product `db/dy*dy` shrinks
from approximately [-743.103,743.103] to [-416.137,416.137]. Consequently the
full Lie-plus-barrier image shrinks from approximately [-889.218,888.270] to
[-562.252,561.305], while evaluation time remains essentially unchanged
(about 0.029 s versus 0.028 s in the isolated diagnostic).

This is not a Barr3-specific relaxation or approximation: the two dynamics
expressions are algebraically identical over the reals. The improvement comes
from retaining the repeated-x correlation through the existing native `sqr`
operator. The published Barr3 fixture therefore now constructs y-dot as
`x*(sqr(x)/3-1)-y`. The original-vs-factored benchmark diagnostic is retained
as evidence for the rewrite. The next search measurement is the 64-box
contractor-free geometric Lie query, which will show whether the tighter root
and descendant enclosures translate into additional exact pruning rather than
only a narrower initial interval.

Deferred direction: dependency-aware algebraic expression optimization.
The Barr3 dynamics rewrite suggests a broader mechanism that should be revisited
only after the current local search/enclosure investigation is exhausted.
Rather than merely minimizing syntactic variable occurrences, such a pass would
generate or recognize algebraically equivalent forms that reduce interval
dependency inflation, for example by preferring native `sqr`, factoring common
terms, Horner-like polynomial forms and other correlation-preserving rewrites.
Candidate forms could eventually be scored on the current box by enclosure
width and evaluation cost. This is deliberately deferred: the current work
continues with direct measurements of the factored Barr3 dynamics before any
general expression optimizer is designed or integrated.

The first search measurement with the factored Barr3 dynamics is neutral at the
64-box frontier. With geometric splitting, Lie-first literal order and all
optional expensive mechanisms disabled, the solver still processes 64 boxes,
prunes 27 and splits 37, with 101 fused literal evaluations and no epsilon
certification. Elapsed time is about 2.713 s versus about 2.738 s for the
previous Lie-first baseline. Thus the much tighter initial Lie-plus-barrier
enclosure does not yet create additional decisions at this shallow frontier.

The next measurement repeats the established 256-box geometric frontier. The
pre-rewrite reference was 123 pruned / 133 split in about 9.24 s. If the
factored dynamics starts changing pruning only at deeper boxes, that run should
reveal it; if the outcome remains identical, the rewrite should be retained for
its enclosure quality but not credited as a search-power improvement.

The factored-dynamics search remains neutral at 256 boxes. The run again
produces exactly 123 pruned and 133 split boxes, with 389 fused literal
evaluations and no epsilon certification. Elapsed time improves modestly from
the pre-rewrite reference of about 9.24 s to about 8.714 s, but there is still
no additional decision power. The rewrite is retained because it is exact and
strictly improves the root enclosure, but current evidence does not credit it
with stronger search pruning.

The dominant remaining dependency is between `db/dy` and the factored y
dynamics, since both ranges can cross zero independently even after the
dynamics range is tightened. The benchmark therefore adds
`lie-correlation-profile`, a benchmark-only replay of the same lower-first
geometric DFS using the Lie-plus-barrier literal. On every split box it records
whether `db/dy` crosses zero, whether y-dot crosses zero, whether both cross
zero, whether both are sign-definite, and average widths of `db/dy`, y-dot and
Lie-plus-barrier. The replay also reports processed/pruned/split counts so its
frontier can be checked against the solver baseline. This diagnostic is meant
to determine whether sign partitioning of the dominant product has realistic
leverage before any new contractor or split policy is implemented.


The 256-box correlation diagnostic reproduces the geometric search frontier:
256 processed, 123 pruned and 133 split. Among the 133 boxes that still require
splitting, the validated range of `db/dy` contains zero in 126 cases (94.7%),
while the factored y dynamics contains zero in only 15 cases (11.3%). Both
contain zero in those same 15 cases, and only 7 split boxes make both factors
sign-definite. The average widths on split boxes are about 5.59 for `db/dy`,
1.32 for y-dot, and 41.26 for Lie-plus-barrier.

This makes y-dot sign partitioning a poor next target: ordinary geometric
splitting already makes the dynamics sign-definite on most unresolved boxes.
The residual ambiguity is concentrated in the neural derivative `db/dy`.
Before changing its enclosure representation, repeat the same diagnostic at the
established 2048-box frontier. If the fraction of unresolved boxes whose
`db/dy` enclosure contains zero remains close to the 256-box value, ordinary
geometric depth is not resolving that derivative efficiently and the next work
should target the derivative representation rather than the dynamics.


The 2048-box correlation profile shows that ordinary geometric splitting does
eventually make the neural derivative more informative, but only gradually.
The replay processes 2048 boxes, pruning 1019 and splitting 1029. Among those
1029 unresolved boxes, `db/dy` still contains zero in 755 cases (about 73.4%),
while the factored y dynamics contains zero in only 63 cases (about 6.1%).
Both contain zero in those same 63 cases; 274 split boxes make both factors
sign-definite. Average widths fall from the 256-box profile to about 3.89 for
`db/dy`, 0.60 for y-dot, and 18.36 for Lie-plus-barrier.

Thus geometric depth is helping, but the derivative enclosure remains the
dominant slow-to-resolve term. The next isolated experiment compares the current
forward-propagated symbolic gradient with an algebraically equivalent
reverse/backprop construction of the same network derivative. The purpose is
not to change differentiation semantics, but to test whether reassociating the
same products and sums yields a tighter natural interval extension before any
new split policy or contractor is designed.


The gradient-reassociation diagnostic is negative. Reassociating the same neural
gradient in a reverse/backprop-like form widens the initial `db/dy` enclosure
from approximately [-57.674,59.448] to [-64.116,66.063] and increases its
evaluation time from about 0.011 s to 0.326 s. The final Lie-plus-barrier
enclosure likewise widens from approximately [-562.252,561.305] to
[-616.745,615.567], with evaluation time rising from about 0.028 s to 0.657 s.
The current forward derivative construction is therefore retained.

The remaining question is whether the frequent zero crossing of the validated
`db/dy` range represents real sign variation or mainly dependency inflation.
Query `lie-gradient-sign-profile` replays the same geometric DFS and, only on
unresolved boxes whose `db/dy` interval contains zero, evaluates `db/dy` at
the midpoint and all four corners. Observing both positive and negative
validated point values proves a real sign change somewhere in the connected
box by continuity. Same-sign samples do not prove sign-definiteness, but a high
same-sign fraction would indicate substantial room for a stronger enclosure.
The first run should use the 256-box frontier to keep this five-point sampling
diagnostic inexpensive.


The 256-box gradient-sign profile gives mixed evidence rather than showing that
the interval zero crossings are mostly real. It reproduces the established
frontier (256 processed, 123 pruned, 133 split) and finds 126 unresolved boxes
whose validated `db/dy` image contains zero. Midpoint-plus-corner sampling
proves both positive and negative values in only 34 of those boxes (27.0%).
The remaining 92 boxes have all five validated samples on one side of zero:
35 positive-only and 57 negative-only, with no point-evaluation ambiguity.
This does not prove those 92 boxes are sign-definite, but it is strong evidence
that dependency inflation may account for a substantial part of the residual
interval ambiguity.

Before introducing a new derivative enclosure, query
`lie-gradient-split-profile` performs one validated candidate bisection in
each state coordinate on every unresolved box whose `db/dy` interval contains
zero. For x and y separately it records how many child boxes become
sign-definite, how often both children become sign-definite, and the sum of
child `db/dy` widths. This directly tests whether a very cheap
derivative-targeted split signal could recover the missing sign information
more efficiently than generic geometric splitting, without evaluating the full
constraint on all candidate children.


The derivative-targeted one-step split profile is negative. On the 126
unresolved boxes whose direct `db/dy` interval contains zero, an x split
produces only 10 sign-definite children in total and makes both children
sign-definite in 5 boxes; a y split produces 13 sign-definite children and
makes both children sign-definite in 6 boxes. The width-score preference is
essentially balanced (x wins 60 boxes, y wins 66), with average child-width
sums about 9.38 for x and 9.17 for y. There is therefore no hidden coordinate
preference that would justify a derivative-specific split heuristic.

One final local enclosure test remains before escalating to a more structural
representation. The earlier mean-value experiment was performed on the complete
Lie-plus-barrier expression over the initial domain and was extremely poor.
Query `lie-gradient-mean-value-profile` instead applies a centered first-order
enclosure only to `db/dy`, and only on the much smaller unresolved frontier
boxes where its direct interval contains zero. The two derivative functions are
built once, then each box uses the midpoint value plus interval derivative
terms. The diagnostic also intersects this enclosure with the ordinary direct
interval, so the combined range can never be weaker than either source alone.
If even this targeted local form fails to remove zero at reasonable cost, the
remaining dependency is unlikely to be recoverable by the generic first-order
machinery already present in Ariadne.


The targeted local mean-value experiment recovers some sign information but is
not competitive enough for integration. At the 256-box frontier, 126 unresolved
boxes have direct `db/dy` ranges containing zero. The local mean-value
enclosure alone is sign-definite in 28 of them, and intersecting it with the
ordinary direct interval remains sign-definite in the same 28 boxes. The direct
range has average width about 5.86; the mean-value range is much worse at about
106.53, while the intersection improves the average to about 4.80. The run
takes about 33.9 s, with 126 midpoint and 252 derivative evaluations. Thus the
intersection extracts a modest amount of extra information, but only by paying
for a very weak first-order enclosure.

The final cheap local dependency test is
`lie-gradient-quadrant-profile`. For each unresolved box whose direct
`db/dy` range contains zero, it bisects x and then y, evaluates the same
validated derivative on the resulting four quadrants, and takes the hull of
those four images. This preserves exact coverage of the parent box and uses
only ordinary interval evaluation. The diagnostic reports whether the
four-quadrant hull becomes sign-definite, whether all four individual
quadrants are sign-definite, whether validated positive and negative quadrants
coexist, and the average width reduction relative to direct evaluation. This
tests whether one small layer of internal domain subdivision can recover
dependency information more cost-effectively than derivative-based
mean-value bounds.


The four-quadrant local subdivision diagnostic tightens the neural derivative
range substantially but resolves relatively few parent boxes. At the 256-box
frontier it again finds 126 unresolved boxes whose direct `db/dy` range
contains zero. Evaluating the derivative on four x/y quadrants and taking the
hull reduces average width from about 5.86 to 3.79, a reduction of roughly 35%.
However, the quadrant hull becomes sign-definite in only 11 boxes, all four
quadrants are individually sign-definite in those same 11 boxes, and 115 parent
boxes remain ambiguous. The run performs 504 extra direct evaluations and takes
about 14.85 s. No box contains both a rigorously positive quadrant and a
rigorously negative quadrant at this resolution.

This result motivates an explicit audit of Ariadne's range-evaluation machinery
before investing in general symbolic expression rewriting. The main findings
are:

* Natural validated interval evaluation through `Formula`/`Procedure` is the
  appropriate tier-zero evaluator for SMT. A `Procedure` compiles the DAG into
  an instruction sequence and evaluates it with outward-rounded interval
  arithmetic. It is cheap and already competitive with direct function
  application, but it inherits the usual dependency problem.
* `Procedure` also implements reverse automatic differentiation. A validated
  gradient over a box is obtained from one forward execution followed by one
  backward sweep, and the DP validated instantiation is already exported. This
  can support centered and monotonicity-aware range bounds without constructing
  and evaluating a separate symbolic derivative function per coordinate.
* The generic function API also provides `gradient_range`,
  `derivative_range`, and `jacobian_range` through differential arithmetic.
  These are useful validated building blocks, but a precompiled `Procedure`
  is the more natural low-overhead representation for repeated SMT box
  evaluation.
* The current validated affine model is not a lightweight affine-arithmetic
  evaluator. `affine_model(domain,function,...)` first constructs a Taylor
  function model using an affine sweeper and then converts it to an affine
  model. Its range is center plus the absolute gradient sum and accumulated
  error. This explains why the Barr3 affine diagnostic can be both more costly
  and much wider than the natural interval extension.
* Validated Taylor models preserve polynomial dependence explicitly and support
  exact domain splitting/restriction without rebuilding the source expression.
  However, their current `range()` keeps special structure only for the
  constant term, linear terms, and pure quadratic terms; mixed quadratic and
  higher-degree terms are accumulated by magnitude into a residual bound.
  Thus a model may retain more correlation than its range routine exploits.
  Strengthening this polynomial range step, for example by intersecting the
  current bound with a validated Horner evaluation of the retained polynomial,
  is a plausible general Ariadne improvement. It does not remove the separate
  cost of constructing the model.
* Taylor patches can be split/restricted algebraically, so a future experiment
  could build a model once and propagate restricted child models down a search
  tree. This is materially different from rebuilding a Taylor model per SMT
  box. The very poor Barr3 root Taylor range means this is not the first
  candidate, but it should not be dismissed solely from the per-root
  construction benchmark.
* Ariadne's Chebyshev polynomial classes currently provide approximate
  polynomial machinery rather than a validated function model with a rigorous
  remainder suitable for SMT exclusion. They are therefore not a drop-in
  validated range evaluator.
* Hull reduction, shaving, monotone/Newton reduction, and interval-Newton-style
  machinery are contractors or equation solvers rather than general scalar
  range evaluators. Their Barr3 cost/effectiveness has already been measured
  separately and should not be used as the baseline range mechanism.

The immediate candidate is therefore a tiered composite evaluator rather than
a wholesale replacement of interval arithmetic. Query
`lie-gradient-composite-profile` compiles the current `db/dy` function to a
`ValidatedProcedure` once. On a box whose natural procedure range still
contains zero, it computes a reverse-AD gradient, a centered mean-value range,
and a monotonicity endpoint range obtained by pinning every coordinate whose
validated derivative is sign-definite to the appropriate endpoint for the
lower and upper evaluations. It intersects all three rigorous enclosures. The
diagnostic reports the separate and combined sign resolutions, average widths,
and time spent in gradient, midpoint, and monotonicity evaluations. This tests
the most promising existing Ariadne machinery before designing a new range
algebra or resuming dependency-aware symbolic rewriting.


The reverse-AD composite diagnostic confirms that a tiered range evaluator can
extract additional information, but also shows that reverse interval AD is not
the right derivative enclosure for this Barr3 frontier. With 256 processed
boxes the replay remains identical to the solver frontier
(`processed=256`, `pruned=123`, `split=133`) and again 126 split boxes
have a natural `db/dy` range crossing zero. The centered reverse-AD bound is
sign-definite on 19 of those 126 boxes. Intersecting it with the natural range
reduces average width from about 5.86 to 5.24, roughly 10.6%. The centered
bound itself is very wide, about 152 on average. No validated reverse-AD
gradient component is sign-definite on any of these boxes, so the monotonicity
endpoint evaluator is never activated. The full diagnostic takes about
13.0 s; the 126 reverse-gradient evaluations account for about 3.33 s and the
126 midpoint evaluations for about 1.01 s.

This is stronger than four-quadrant subdivision in sign resolution
(19 versus 11 parent boxes) at slightly lower total diagnostic time, but it is
weaker than the earlier symbolic-derivative mean-value experiment, which
resolved 28 boxes and produced an average direct/mean-value intersection width
of about 4.80. The difference matters: the centered formula is useful, but the
reverse-AD interval gradient is materially wider than the interval ranges of
the separately differentiated functions. Consequently, spending more effort
on reverse-AD monotonicity is not justified at this point.

The next diagnostic, `lie-gradient-symbolic-procedure-profile`, separates
derivative representation from derivative evaluation cost. It differentiates
`db/dy` symbolically once, compiles the original function and both derivative
functions into `ValidatedProcedure` objects once, and then evaluates those
procedures on the same 256-box frontier. It forms the same centered and
monotonicity endpoint bounds and intersects them with the natural range. If it
recovers the 28-box resolution and approximately 4.80 average intersected
width of the previous symbolic mean-value profile while substantially reducing
runtime, then compiled symbolic derivative procedures are a credible
second-tier range evaluator. If the runtime remains high, the next useful
directions are improving polynomial/Taylor range extraction or returning to
dependency-aware expression rewriting rather than further elaborating this
first-order scheme.


The compiled symbolic-derivative profile reproduces the quality of the earlier
symbolic mean-value bound exactly enough to isolate evaluation cost. At the
256-box frontier it resolves 28 of the 126 natural `db/dy` zero crossings;
the average centered width is about 106.53 and the intersection with the
natural interval has average width about 4.797. Symbolic derivative ranges also
identify 11 sign-definite coordinate instances across 10 boxes, but the
monotonicity endpoint bound is independently sign-definite in only 3 cases and
does not increase the combined resolution beyond the same 28 boxes. Runtime
falls from about 33.9 s for generic function application to about 27.9 s, but
252 derivative-procedure evaluations still consume about 18.0 s. Procedure
construction itself is only about 0.51 s and symbolic derivative construction
about 0.04 s. Therefore the bottleneck is repeated derivative evaluation, not
derivative construction.

A further inspection of `Procedure` corrects an initially pessimistic
assumption: Ariadne already has a multi-output `Vector<Procedure<Y>>`. Its
constructor converts all output formulas using one node cache, so common
Formula DAG nodes shared by multiple outputs become one instruction stream,
and a vector evaluation executes that stream once. Query
`lie-gradient-shared-procedure-profile` therefore constructs the vector
`[db/dy, d(db/dy)/dx, d(db/dy)/dy]`, compiles it as one validated vector
procedure, and uses one shared evaluation per frontier box. It reports the sum
of instruction counts from the three separately compiled scalar procedures and
the instruction count of the shared vector procedure, in addition to the same
centered/monotonicity quality metrics and timing. This directly measures how
much common structure symbolic differentiation preserves and whether Ariadne's
existing multi-output procedure machinery can turn the 28-box centered bound
into a practical second-tier range evaluator without introducing a new range
algebra.


The first shared-procedure performance result (63 of 126 ambiguous boxes made
sign-definite) is invalid and must not be used as evidence for the range
evaluator. The consistency check exposed a pre-existing bug in the validated
fallback overload `join(ValidatedScalarMultivariateFunction,
ValidatedVectorMultivariateFunction)`: it allocated and copied the vector tail
using `f1.result_size()`, i.e. the scalar's result size, rather than
`f2.result_size()`. For a three-output vector
`[db/dy,d(db/dy)/dx,d(db/dy)/dy]`, only the first derivative was copied and
the final derivative component remained the default zero function. The check
therefore reported exact agreement for the natural value and x derivative,
but disagreement for the y derivative on all 126 checked boxes; the shared y
derivative had essentially zero width while the separately compiled derivative
had average width about 48.73.

The validated scalar-vector join fallback now allocates
`f2.result_size()+1` outputs and copies all `f2` components. Regression tests
cover both the three-component validated scalar-vector join and
`Vector<ValidatedProcedure>` evaluation with `FloatDPBounds`, the latter
also guarding the explicit-zero initialization fix made when the shared
procedure experiment first instantiated this path. The shared-procedure
consistency check must be rerun after this fix; only if all three output
mismatch counters are zero should the performance profile be interpreted.


After fixing the validated scalar-vector join fallback, the shared-procedure
consistency check passes exactly on the full 256-box diagnostic frontier. It
checks all 126 unresolved boxes whose natural `db/dy` range contains zero and
reports zero mismatches for the natural value, x derivative and y derivative,
zero negative raw widths, and zero maximum endpoint difference. The average
separate/shared derivative widths also agree exactly: about 41.09 for the x
derivative and 48.73 for the y derivative. This validates the semantics of the
multi-output procedure path after the join fix. The earlier shared-performance
numbers obtained before the fix remain invalid and must not be reused; the
performance profile must now be rerun from the corrected branch.


The corrected shared-procedure performance run preserves the expected
first-order bound quality but provides essentially no cross-output sharing.
At the 256-box frontier it again resolves 28 of 126 natural `db/dy` zero
crossings, with average combined width about 4.797. The three separately
compiled procedures contain 1,705,353 instructions in total, while the
multi-output procedure contains 1,705,349: only four instructions are shared.
The shared evaluation therefore takes about 20.57 s for 133 boxes, compared
with about 19.11 s for the natural plus separate derivative evaluations in the
previous profile, and total diagnostic time rises to about 29.39 s. The
multi-output `Procedure` mechanism is not the limitation; the derivative
functions arrive with almost entirely distinct Formula node identities.

Ariadne already provides structural CSE for
`Vector<Expression<Real>>`, using one eliminator cache across all outputs.
However, the vector `make_formula` conversion previously called the scalar
conversion independently for each component, thereby discarding cross-output
node sharing before `Vector<Procedure>` could exploit it. The vector
conversion now uses one `FormulaSharingCache` across all components, and the
auxiliary-substitution overload delegates to the same shared conversion.
A regression test checks that two structurally identical vector outputs reduce
to a three-instruction procedure after expression CSE.

Query `lie-gradient-expression-cse-profile` is the gate for this direction.
It forms `[db/dy,d(db/dy)/dx,d(db/dy)/dy]` directly as real expressions,
simplifies the two symbolic derivatives, applies vector-wide common
subexpression elimination, and compares the resulting effective multi-output
Procedure size against the same expressions without CSE. This first measures
structural compression only. A validated frontier benchmark is justified only
if the instruction count falls substantially; otherwise first-order symbolic
range evaluation should be considered exhausted and work should move to a
different enclosure representation.


The structural CSE gate is strongly positive. Building the symbolic
`[db/dy,d(db/dy)/dx,d(db/dy)/dy]` expressions yields 917,441 distinct node
pointers and a 917,441-instruction multi-output Procedure without CSE.
Vector-wide common-subexpression elimination leaves the mathematical node
count unchanged but reduces distinct node pointers to about 90,698 and, after
the now-sharing-preserving Expression-to-Formula conversion, reduces the
Procedure to 51,715 instructions. This is a reduction of about 94.4% in the
executed instruction stream. CSE itself is expensive in this prototype
(about 40.8 s), whereas compilation of the CSE DAG takes about 0.016 s versus
0.309 s for the raw DAG. The result therefore separates a potentially very
good repeated evaluator from a currently expensive one-time canonicalization
step.

Query `lie-gradient-expression-cse-frontier` now measures the repeated-use
side directly. It retains the ordinary validated natural interval procedure as
tier zero. Only when that range contains zero does it evaluate the CSE
multi-output effective Procedure on `FloatDPBounds` to obtain the function
and both symbolic derivative ranges in one shared instruction stream. It then
forms the same centered and monotonicity endpoint bounds used by the previous
symbolic profiles. The diagnostic reports the one-time derivative/CSE/build
costs separately from frontier evaluation time. The acceptance criterion is
that the CSE path must reproduce the established 28-of-126 sign resolution
and approximately 4.797 average combined width while reducing repeated
derivative-evaluation cost substantially; the 40 s CSE construction cost is
then a separate optimization/caching problem rather than a range-quality
failure.


The validated CSE frontier run confirms that the compressed symbolic evaluator
preserves the established first-order enclosure quality while making repeated
evaluation cheap. On the same 256-box frontier it resolves 28 of the 126
natural `db/dy` zero crossings, with average centered width about 106.53,
average monotonicity width about 5.82, and average combined width about 4.797,
exactly matching the earlier symbolic-derivative profiles. The 51,715-
instruction CSE procedure evaluates all three outputs on the 126 ambiguous
boxes in about 0.378 s; midpoint evaluations take about 0.126 s and the 20
monotonicity endpoint evaluations about 0.020 s. By comparison, evaluating the
two separately compiled symbolic derivatives previously took about 18.0 s.
The frontier replay itself takes about 9.07 s, so the second-tier range work is
now a small addition to the natural Lie evaluation cost.

The remaining obstacle is preprocessing: symbolic derivative construction
takes about 1.97 s and vector-wide CSE about 45.18 s in this prototype. The
current `CommonSubroutineEliminator` performs structural ordered-set lookup
before child canonicalization, so comparisons repeatedly recurse through large
subtrees. It is now changed to bottom-up interning: original node pointers are
memoized, children are canonicalized first, and only the resulting canonical
candidate is looked up in the structural set. The expression ordering also has
a pointer-identity fast path. These changes preserve the existing structural
CSE semantics while avoiding repeated work on already visited DAG nodes and
making structural comparisons terminate immediately on canonical shared
children. The next measurement must retain the 51,715-instruction result while
reducing the CSE preprocessing time; otherwise a hash-consed key representation
will be required.


The first bottom-up CSE optimisation did not reduce preprocessing cost. It
preserved the same 51,715-instruction result but increased Barr3 CSE time from
about 45.2 s to about 65.5 s. The added node-pointer memo used Ariadne's
`Map`, which is a `std::map`, and most source nodes in this workload are
already pointer-distinct (about 917k distinct pointers for about 1.02M node
occurrences), so the logarithmic memo lookup mostly added work. More
importantly, the canonical cache still ordered full `Expression` objects with
recursive structural comparison.

The CSE cache now uses a canonical local comparator instead. Children are
canonicalised first; for unary, binary and graded internal nodes the cache
ordering compares only the operator and canonical child pointers (plus the
graded payload), so it does not recursively traverse already-canonical
subtrees. Constants and variables retain the existing payload ordering, while
vector getter nodes retain the previous structural fallback because vector
subexpressions are not yet canonicalised by this scalar CSE path. The
pointer-map memo is removed. This preserves the same structural equivalence
criterion for the scalar arithmetic used by Barr3 while making each internal
cache comparison local. The next structural profile must still produce 90,698
distinct CSE node pointers and 51,715 procedure instructions; only the
preprocessing time is expected to change.


The local-key CSE optimisation is successful. It preserves exactly the previous
compressed structure (90,698 distinct CSE node pointers and 51,715 Procedure
instructions) while reducing Barr3 CSE time from about 45.18 s to about 3.23 s.
This is roughly a 14x preprocessing speed-up, and it also reverses the failed
pointer-memo attempt that had increased CSE time to about 65.5 s. Symbolic
derivative construction remains about 1.8 s and CSE compilation about 0.017 s.
At this point the CSE-based first-order range evaluator has both acceptable
one-time compilation cost and very low repeated evaluation cost.

Before integrating it into the SMT solver, query
`lie-gradient-cse-prune-profile` measures solver-relevant decision power
rather than derivative range quality. It replays the same geometric frontier
and considers only boxes that the ordinary `lie+barrier` interval cannot
prune and whose natural `db/dy` range crosses zero. For those boxes it forms
the CSE centered/monotonicity intersection for `db/dy`, recomposes the Lie
range as `db/dx*dx + db/dy*dy + barrier`, and intersects that result with the
ordinary direct Lie range. It reports both the natural recomposed baseline and
the improved recomposed/intersected pruning counts. This isolates whether the
tighter neural derivative actually raises the lower bound of the SMT literal
enough to exclude additional boxes; a range evaluator should not be integrated
based only on narrower intermediate derivatives.


The CSE pruning-impact diagnostic is strongly positive. On the established
256-box geometric frontier, the ordinary direct Lie range prunes 123 boxes and
leaves 133 boxes to split. Among those 133 unresolved boxes, 126 have a natural
`db/dy` interval crossing zero. The CSE-centered/monotonicity evaluator makes
28 of those derivative ranges sign-definite, but the more important result is
the effect on the complete Lie literal: recomposing
`db/dx*dx + db/dy*dy + barrier` with the improved `db/dy` enclosure rejects
65 of the 133 boxes that the natural Lie range could not reject. This is about
48.9% of the previously unresolved frontier boxes, or about 51.6% of the boxes
whose natural `db/dy` range crossed zero.

The natural recomposed Lie range has the same measured average width as the
direct `lie+barrier` range on these boxes (about 43.32), so the decomposition
itself does not introduce a measurable range penalty in this diagnostic. The
improved recomposed range reduces average width to about 39.78, and intersecting
it with the direct Lie range gives the same width and the same 65 rejections.
Thus the extra pruning is attributable to the tighter validated `db/dy`
range rather than to an unrelated change of expression decomposition.

The repeated-use cost remains modest: 126 CSE evaluator calls take about
0.535 s and the 378 component evaluations about 2.03 s. One-time preparation
cost is about 1.79 s for symbolic derivatives, 3.53 s for CSE, and 0.053 s for
Procedure construction. These results make the CSE-based first-order range
evaluator a serious candidate for a second-tier SMT range mechanism. However,
the 65 extra rejections come from replaying the old frontier; integrating the
evaluator would change the search tree, so they must not be extrapolated
directly into a final solver pruning rate. The next step after the planned
rebase is to integrate or emulate the tier in the actual search loop and
measure the resulting frontier and total runtime.


Before integrating the CSE range tier into the SMT search loop, one distinction
must be resolved. The strong 65-box replay gain comes from tightening the
internal neural factor `db/dy` and then recomposing the Lie expression. A
generic SMT literal compiler normally sees only the complete primitive
expression, so a root-level first-order evaluator is not automatically
equivalent to the successful subexpression experiment.

Query `lie-root-cse-range-profile` therefore applies the same validated
CSE/centered/monotonicity construction directly to the complete
`lie+barrier` expression. It replays the established geometric frontier and
counts additional boxes whose improved root range has lower bound at least
zero. If this root evaluator recovers significant pruning, it can be introduced
as a generic optional compiled-literal tier. If it does not, the range-evaluator
architecture must preserve or discover useful internal subexpressions instead
of applying first-order range analysis only at the literal root.


The root-level CSE range profile is stronger than the earlier
subexpression-specific experiment. On the established 256-box frontier, the
natural `lie+barrier` range prunes 123 boxes and leaves 133 unresolved.
Applying a CSE-centered first-order range directly to the complete root
expression proves another 90 of those 133 boxes infeasible, about 67.7% of the
old unresolved frontier. No coordinate derivative of the root is sign-definite
on these boxes, so monotonicity contributes nothing. The centered range itself
is wide on average (about 927), but its one-sided placement is useful; after
intersection with the natural range the average width falls from about 41.26
to 35.73. Repeated evaluation remains cheap: 133 CSE evaluations take about
0.53 s and midpoint evaluations about 0.14 s. One-time preparation is about
4.33 s for symbolic derivatives, 8.49 s for CSE and 0.03 s for Procedure
construction on this larger root expression.

This result means the second-tier evaluator need not depend on Barr3-specific
knowledge of `db/dy`; it can operate generically on a compiled theory literal
root. Before changing the core solver representation, query
`lie-root-cse-search` runs a diagnostic DFS whose frontier is actually changed
by the improved root range. It implements the same cheap setting used by the
Barr3 baseline: natural Lie evaluation first, CSE-centered refinement only
when natural evaluation does not reject, the barrier literal in the selected
literal order, whole-box epsilon certification, geometric splitting and the
selected child order. This provides an end-to-end estimate of the changed
search tree and repeated-evaluation cost without yet adding CSE-specific state
to `CompiledTheoryLiteral`.


The 4,096-box changed-tree comparison closes the root-CSE search experiment.
The post-rebase natural baseline processes all 4,096 allowed boxes, pruning
2,043 and splitting 2,053, with no epsilon certification and about 130.98 s
elapsed. The root-CSE diagnostic processes the same 4,096-box budget, pruning
2,045 and splitting 2,051. Thus the stronger root range avoids only two splits,
exactly the same net split reduction already observed at the 2,048-box gate.
The implied residual frontier remains about 11 boxes for the natural search
versus 7 for the CSE search; the gap does not widen with depth.

The changed pruning attribution is also revealing. At 4,096 boxes the CSE
search reports 478 CSE-Lie prunes and 1,567 barrier prunes, compared with 469
and 552 respectively at the 2,048-box point. The CSE-Lie contribution therefore
barely grows while the barrier absorbs almost all additional pruning deeper in
the changed tree. This supports the interpretation that the stronger root
enclosure mostly moves pruning earlier and changes which literal eventually
rejects descendants, rather than materially reducing total search complexity.

The cost is decisively unfavorable. Natural search takes about 130.98 s. The
CSE changed-tree search alone takes about 165.59 s, roughly 26% slower, and its
one-time derivative/CSE/Procedure preparation adds about 13.57 s more, for an
effective total near 179.2 s, roughly 37% above baseline. Root-level
CSE/centered range evaluation is therefore retained as a useful validated range
technique and diagnostic, but rejected as an eager second-tier Barr3 search
mechanism. The next development direction returns to dependency-aware
algebraic expression optimization, whose goal is to strengthen the natural
interval extension itself without paying a second evaluator on every unresolved
box.


A separate polynomial-model direction is being retained for later rather than
pursued immediately. The most promising architecture is not a pure Bernstein
model, but a validated Chebyshev model used as the approximation/composition
representation together with a Bernstein-form range bounder for the polynomial
part. Chebyshev is the better fit for smooth analytic elementary functions such
as tanh and for repeated composition through the network; Bernstein is the
better fit for certified polynomial range bounds, positivity tests and cheap
subdivision on the low-dimensional input box.

A rigorous Chebyshev implementation would require a genuine model layer rather
than merely changing the numeric type of the existing
`ChebyshevPolynomial`: polynomial plus uniform remainder, explicit roundoff
accounting, sweeping/truncation, validated elementary-function approximation,
composition, scaling/restriction and a validated range routine. The existing
Chebyshev polynomial code can likely serve as the algebraic kernel, but is
currently approximate-only and lacks model semantics. A C0 model is sufficient
for the Barr3 range-evaluation use case; derivative-aware parity with Taylor
would require additional derivative remainder information rather than simply
differentiating a polynomial-plus-uniform-error model.

The preferred eventual vertical slice is therefore: validated C0 Chebyshev
model on [-1,1]^n, validated tanh, polynomial/remainder propagation, then
Chebyshev-to-Bernstein conversion for range certification and subdivision.
This direction remains deferred while work returns to dependency-aware
algebraic expression optimisation, whose immediate goal is to improve the
natural interval extension itself without the cost of a second evaluator.


Dependency-aware algebraic optimisation resumes after closing the eager root-CSE
search experiment. Ariadne's current symbolic `simplify(Expression)` is not an
algebraic optimiser: it primarily performs constant folding and structural
simplification and does not systematically factor, expand, choose Horner forms,
or score equivalent forms by validated interval quality.

The first gate therefore avoids designing a general rewrite engine prematurely.
Query `lie-algebraic-form-profile` compares three exactly equivalent outer
forms of the Barr3 Lie expression on the same baseline DFS frontier: the
current factored-dynamics form, an expanded product form, and a form that
groups the repeated y dependence as
`y*(db_dx-db_dy) + db_dy*x*(sqr(x)/3-1) + barrier`. It reports root widths,
average frontier widths, pruning counts and evaluation times. If one equivalent
form improves validated natural interval bounds materially, it supplies a
concrete rewrite pattern and scoring target for a generic dependency-aware
optimizer. If all outer forms are effectively equivalent, optimisation effort
should move inside the network derivative expression rather than adding broad
symbolic rewriting machinery without evidence.


The 256-box outer-form gate is decisively negative for the two alternative
root rewrites. The established factored-dynamics baseline again prunes 123
boxes and splits 133, with average Lie-plus-barrier width about 24.325 on the
replayed frontier. Both the expanded-product and grouped-y root forms prune
zero boxes. Their root widths increase from about 1123.56 to about 1228.26,
and their average frontier widths increase to about 33.415 and 33.412
respectively. They are also slower to evaluate (about 10.06 s and 10.04 s
versus 7.17 s for the baseline over this diagnostic). This is strong evidence
that dependency-aware optimisation cannot be reduced to syntactically grouping
a repeated variable at the literal root: separating the network-dependent
\`db/dy\` factor from the already successful \`polynomial-y\` dynamics destroys
the useful interval correlation.

The next gate moves the algebraic combination inside the network while retaining
ordinary natural interval evaluation. Query
\`lie-directional-propagation-profile\` propagates the directional derivative
of the network directly along the Barr3 vector field rather than constructing
\`db/dx\` and \`db/dy\` separately and combining them only at the root. It
compares two exactly equivalent first-layer seeds on the unchanged baseline DFS
frontier. The factored seed uses
\`w_x*y + w_y*(polynomial-y)\`; the grouped seed uses
\`y*(w_x-w_y) + w_y*polynomial\`, where the regrouping is performed only while
\`w_x\` and \`w_y\` are constants. Both variants then propagate one directional
quantity through the second layer and output weights. This remains a single
natural interval evaluator per candidate expression; there is no CSE-centered,
mean-value, Taylor, Chebyshev/Bernstein, or other second-tier range mechanism.
The measurement reports root widths, average baseline-frontier widths, pruning
overlap/differences, and evaluation cost. A positive result would identify
early directional combination as a concrete dependency-aware rewrite pattern;
a negative result would push the optimisation search further inside the
network's weighted sums/activation factors rather than back toward a second
evaluator.


The 256-box directional-propagation gate is negative. The established baseline
again prunes 123 boxes and splits 133, with average width about 24.325. Direct
directional propagation from the first layer prunes only 9 boxes with the
factored dynamics seed and 12 with the constant-coefficient grouped seed; none
of those prunes is additional to the baseline. The candidate root widths are
about 2259.27 and 1994.06 versus 1123.56 for the baseline, and their average
frontier widths are about 42.11 and 39.15. Both candidates are individually
faster to evaluate (about 5.15 s each versus 7.22 s for the baseline), but the
loss of enclosure strength is overwhelming. Combining the vector field with
the network at the first layer therefore destroys more useful dependency
structure than it recovers.

The next gate tests the intermediate algebraic location rather than abandoning
directional combination entirely. Query \`lie-layer2-directional-profile\`
retains the baseline's separate \`dx\` and \`dy\` derivative propagation
through the first-layer weighted sums. For each second-layer neuron it then
forms the local directional term
\`a2_i*(dzdx_i*y + dzdy_i*(polynomial-y))\` before multiplying by the output
weight and performing the final output sum. A second candidate groups the
local y dependence as
\`a2_i*(y*(dzdx_i-dzdy_i) + dzdy_i*polynomial)\`.
This tests whether the useful correlation is lost specifically by the global
output aggregation: it combines the two gradient components only after their
separate first-layer structure has been preserved, but before the final
\`sum_i w3_i\`. As in the previous gates, the baseline expression alone
determines the DFS frontier, and each candidate remains a single ordinary
natural interval evaluator.


The 256-box layer-2 directional gate is also negative, but it localises the
damage more precisely. Combining the two gradient components only at each
second-layer neuron recovers much of the enclosure quality lost by first-layer
directional propagation: the factored candidate prunes 89 boxes rather than 9.
It nevertheless remains strictly weaker than the established baseline, which
again prunes 123 boxes. The factored layer-2 form loses 34 baseline prunes,
gains none, widens the root from about 1123.56 to 1525.54 and raises average
frontier width from about 24.325 to 28.823. It is somewhat cheaper to evaluate
(about 6.03 s versus 7.36 s), but the decision-power loss is still material.
The grouped layer-2 form is worse again: 35 prunes, root width about 1494.79,
average width about 31.773 and about 7.75 s evaluation time. Together with the
root and first-layer experiments, this closes the search for a better location
at which to merge `db/dx` and `db/dy`: preserving the two derivatives
separately through the final output aggregation is part of the useful natural
interval structure, not an avoidable source of dependency inflation.

Dependency-aware optimisation therefore moves inside the dominant `db/dy`
expression itself while leaving the outer Lie structure unchanged. The current
forward derivative and the previously tested full reverse/backprop
reassociation are two extreme parenthesisations of the same double sum. The
next gate tests intermediate blockwise reassociations of second-layer neurons,
so that shared first-layer derivative terms are factored only within small
blocks rather than across the whole layer. This directly asks whether there is
a useful granularity between the strong forward form and the weak full reverse
form, while retaining one ordinary natural interval evaluation.


The 256-box blockwise reassociation profile is decisively negative for the
entire tested family. The forward baseline again prunes 123 boxes and splits
133, with average Lie-plus-barrier width about 24.325 and about 5.50 s of
baseline evaluation time. Blocks of 2, 4, 8, 16 and 32 second-layer neurons
prune only 11, 17, 21, 27 and 31 boxes respectively; none gains a single prune
that the baseline misses. Average widths remain much worse than baseline,
decreasing only from about 31.261 at block size 2 to about 30.129 at block size
32. The root enclosure is identical across all five candidates
(approximately [-608.557,607.609], width 1216.17), still wider than the
baseline root width 1123.56.

The computational profile is equally unfavorable. The baseline Procedure has
267,856 instructions, whereas the block candidates contain roughly
3.08--3.11 million instructions each and require about 69--70 s of evaluation
time apiece on the 256-box diagnostic, versus about 5.5 s for the baseline.
The monotone improvement in pruning as block size increases is not evidence for
a useful intermediate granularity: even the largest tested block loses 92
baseline prunes, and the whole family is both much wider and an order of
magnitude larger. This closes the forward-to-reverse reassociation family as a
practical dependency-aware optimisation target.

The next step therefore stops guessing whole-expression parenthesisations and
attributes the width of the successful forward `db/dy` form to its actual
second-layer neuron contributions. Query
`lie-gradient-width-attribution-profile` retains the established baseline DFS
frontier and, on boxes that survive the Lie range test, measures for every
second-layer neuron the validated width of its local
`w3_i*(1-h2_i^2)*dzdy_i` contribution, together with the widths and
zero-crossing frequency of `dzdy_i` and the activation derivative factor
`1-h2_i^2`. It reports the highest-width contributors aggregated across the
frontier. The purpose is diagnostic: if a small subset of neurons dominates the
forward derivative width, the next algebraic rewrite can target their local
structure; if width is diffuse across the layer, local hand rewrites are
unlikely to scale and the optimiser needs a more global scoring/search
mechanism.


The 2048-box width-attribution profile confirms that the 256-box pattern is
stable with depth. Average db/dy width falls from about 5.588 to about 3.887,
but the four leading contributors remain neurons 42, 49, 59 and 57. Their
combined share rises slightly from about 25.3% to 26.2%, while the top twelve
contributors rise from about 44.7% to 51.8%. The activation-derivative factor
1-h2^2 still never touches zero for the reported dominant neurons, and dzdy
crosses zero only on a small minority of the 1029 inspected boxes. Thus the
observed dependency pattern is not a shallow-frontier artifact and does not
reduce to sign ambiguity.

One caution is important: equality between db/dy interval width and the sum of
local term widths is a property of interval addition, not proof that the final
sum introduces no dependency loss. The remaining inflation may arise inside
each product (1-h2^2)*dzdy, across different neuron contributions in the final
sum, or both.

Query lie-gradient-correlation-attribution-profile therefore performs one
validated two-coordinate quadrant subdivision of every baseline-surviving box.
It compares three widths: the direct db/dy interval; the sum of per-neuron
term hulls after evaluating each term on the four common quadrants; and the
hull of the complete db/dy sum evaluated on those same quadrants. Reduction
from the first to the second isolates dependency recoverable within individual
terms. The additional reduction from the second to the third measures
cross-neuron correlation/cancellation preserved by evaluating the whole sum on
a common partition. This diagnostic decides whether the next algebraic search
should target the local activation-derivative product or the structure of the
neuron sum.


The first 256-box run of
`lie-gradient-correlation-attribution-profile` is invalid as a baseline-frontier
measurement and its reported recovery fractions must not be used as evidence.
It produced `pruned=0`, `split=256`, `inspected=256` instead of the
established `123/133` baseline replay. Inspection found that the profile
accidentally reused its diagnostic x-coordinate subdivision to advance the DFS,
thereby replacing the normal geometric `box.split()` traversal with repeated
x splitting. The measured values (about 39.7% local recovery, 12.2% additional
cross-term recovery and 51.8% total recovery) therefore refer to a different
frontier.

The profile is corrected so the baseline replay first computes and retains the
ordinary geometric children from `box.split()`. The x/y four-quadrant
subdivision is now used only to evaluate local term and complete-sum hulls on
the current box; traversal always continues through the retained geometric
children. A valid 256-box rerun must reproduce `processed=256`,
`pruned=123` and `split=133` before any attribution result is interpreted.
