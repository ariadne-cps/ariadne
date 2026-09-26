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

The only uncovered function attributed to `constraint_solver.hpp` is the
defaulted virtual destructor. A direct lifetime test through
`ConstraintSolverInterface*` is added so that the interface destruction path
is exercised explicitly rather than being left as a coverage artefact.

## Current open work

The immediate work on `solvers-smt#830` is:

1. validate the explicit UNKNOWN-reason plumbing across sequential, parallel and
   Boolean search and restore full coverage for the new branches;
2. make splitting focus on constraints whose validated box evaluation has not
   yet met the epsilon stopping condition, mirroring the active-formula idea in
   dReal while retaining Ariadne's sensitivity guidance;
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
