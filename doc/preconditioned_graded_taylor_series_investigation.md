# Preconditioned Graded Taylor Series Integrator: Investigation Log

**Branch:** `solvers-integrator#357`  
**Last updated:** 2026-09-23  
**Current HEAD when this log was created:** `cb1496eb436a1d4ed226554a4f18eaa4da39f29a`  
**Latest analysed investigation HEAD:** `f2fb94ce5a3b2dabd9341e3ab051d2c625dbf916`

## Purpose of this document

This is the persistent research log for the investigation of a preconditioned graded Taylor-series integrator in Ariadne. It exists so that the investigation can be resumed from a fresh chat or development session without reconstructing the reasoning from Git history.

**Before continuing this investigation:**
1. Read this file completely.
2. Check the current branch HEAD and subsequent commits.
3. Do not repeat an experiment listed under "Rejected or currently unsupported hypotheses" unless new evidence justifies it.
4. Preserve the distinction between measured results, conclusions supported by those results, and hypotheses still under investigation.
5. Update this file after every experiment that materially changes the conclusions or next step.

The original motivation was comparison with Taylor-model flowpipe techniques described by Xin Chen / Flow*, but the new Ariadne class should be judged on its own merits and should not be designed as a Flow* clone.

---

## 1. Main question

Can a `PreconditionedGradedTaylorSeriesIntegrator`, based on Ariadne's existing `GradedTaylorSeriesIntegrator` machinery but maintaining an explicit affine local coordinate system,

```
x = c + A y
```

reduce wrapping / accumulated Taylor-model error enough to improve Ariadne's existing integrators at competitive computational cost?

A related but distinct question is whether this explains previously observed differences between Ariadne and Flow*. **That has not been established.** Be careful not to conflate these questions.

---

## 2. Architecture that has survived testing

The useful preconditioned architecture is a persistent two-layer representation

```
X_l(s) = c_l + A_l y_l(s)
```

where:
- `c_l` is the physical centre,
- `A_l` is the local linear coordinate map,
- `y_l(s)` is the Taylor-function representation of the accumulated initial-set dependence.

The key implementation types are:
- `PreconditionedGradedTaylorSeriesIntegrator`
- `PreconditionedTaylorSeriesState`
- `PreconditionedTaylorSeriesStep`

The evolver carries the specialised `PreconditionedTaylorSeriesState` from one workload item to the next instead of reconstructing it from the already-composed physical enclosure at every step.

### Important ordering

For the evolved set, evaluate the local flow at the final time **before** composing with the accumulated local mapping:

```
local flow Phi(y,t)
       |
       | t = h
       v
local endpoint Phi(y,h)
       |
       | y = y_l(s)
       v
X_{l+1}(s)
```

Composing the full space-time flowpipe first and evaluating `t=h` afterwards produced unnecessary mixed terms and larger sweep/remainder errors.

### Important negative result

Do **not** reconstruct a fresh QR preconditioner from the fully composed physical Taylor mapping at every step. That experiment caused rapid/explosive growth of the remainder. Precondition the fresh local transition and preserve the two-layer state across steps.

---

## 3. Naming / conceptual classification

The new class is called:

`PreconditionedGradedTaylorSeriesIntegrator`

This naming is intentional.

It is based primarily on the machinery of `GradedTaylorSeriesIntegrator`: the flow expansion is generated through graded series / graded differentials, not by the Picard contractor used by Ariadne's Picard integrators.

Earlier discussion explicitly compared:
- `PicardIntegrator`
- `TaylorSeriesIntegrator`
- `GradedTaylorSeriesIntegrator`

The graded Taylor-series implementation was retained as the base because it behaved better for the investigation and accurately describes the construction mechanism.

---

## 4. Confirmed results

### 4.1 Persistent two-layer state is essential

Re-preconditioning the globally composed physical mapping caused an artificial wrapping/remainder explosion.

Carrying `(c,A,y)` persistently removed that instability and allowed QR preconditioning to complete the Van der Pol benchmark.

This is an architectural result, not merely tuning.

### 4.2 QR has a real geometric / error effect

A controlled comparison using the same persistent two-layer propagation gave approximately:

| Configuration | Reach steps | Time | Final radius | Final error |
|---|---:|---:|---:|---:|
| Preconditioned IDENTITY, max error 1e-6 | 353 | 6.01 s | 0.239 | 0.03111 |
| Preconditioned QR, max error 1e-6 | 440 | 12.86 s | 0.222 | 0.01182 |

Thus QR materially reduced final accumulated error and slightly reduced the final radius, but at a large cost in step count and runtime.

### 4.3 QR is not currently competitive at equal final accuracy

A later equal-accuracy experiment compared the standard `GradedTaylorSeriesIntegrator` at tighter local tolerances with preconditioned QR at `1e-6`:

| Method | StepMaximumError | Reach steps | Time | Final error | Final radius |
|---|---:|---:|---:|---:|---:|
| GradedTaylorSeries | 1e-6 | 351 | 5.82 s | 0.03617 | 0.244 |
| GradedTaylorSeries | 3e-7 | 373 | 6.37 s | 0.01523 | 0.225 |
| GradedTaylorSeries | 1e-7 | 421 | 7.46 s | 0.00581 | 0.217 |
| GradedTaylorSeries | 3e-8 | 467 | 8.32 s | 0.00389 | 0.216 |
| Preconditioned QR | 1e-6 | 440 | 14.57 s | 0.01182 | 0.222 |

The standard graded integrator can therefore reach comparable or better final accuracy substantially faster. The current QR implementation is promising geometrically but **not yet an efficiency improvement**.

### 4.4 The QR penalty appears immediately after the first rotation

A diagnostic compared IDENTITY and QR at the **same physical state after the first step**.

The physical state had error approximately:

```
5.4532979603e-7
```

For IDENTITY at the second step:

```
h = 0.02
physical error = 6.345026606e-7
=> accepted
```

For QR, on the same physical state:

```
h = 0.02   physical error = 4.366899840e-6   => rejected
h = 0.015  physical error = 1.037166203e-6   => rejected
h = 0.01   physical error = 1.369072124e-7   => accepted
```

So the key efficiency problem is already present in the **single local graded-series calculation after rotation**. It is not primarily caused by later composition with the accumulated global Taylor mapping.

### 4.5 QR worsens the interval problem seen by the graded machinery

At that same second-step diagnostic:

IDENTITY local bounding box:

```
[-0.30155446, 0.34887864] x [-0.28472898, 0.14196793]
```

QR local bounding box:

```
[-0.34362090, 0.41198114] x [-0.33444289, 0.19626299]
```

IDENTITY transformed vector-field range:

```
[1.9404668, 2.3671638] x [-7.0467233, -1.3338121]
```

QR transformed vector-field range:

```
[2.0080381, 3.4843396] x [-7.4608531, -0.28910370]
```

This is evidence of interval wrapping introduced by rotating and then axis-aligning the boxes. The remaining question is where the graded-differential calculation amplifies this into the approximately 6.9x larger physical remainder at `h=0.02`.

---

## 5. Remainder investigation

### 5.1 Polynomial defect of the ordinary graded-series first step

For the first Van der Pol step at `h=0.02`, order 5, the ordinary graded Taylor-series model had approximately:

```
model errors:
component 0 ~ 3.97e-8
component 1 ~ 4.95e-7
```

The defect of the polynomial part,

```
f(p) - dp/dt
```

had ranges of approximately:

```
component 0: [-1.43e-6, 2.33e-6]
component 1: [-1.33e-5, 3.13e-5]
```

A crude multiplication of the defect scale by `h=0.02` gives the same order of magnitude as the Taylor-model remainder. Therefore there was no evidence that the first local GradedSeries remainder was grossly overestimated.

### 5.2 GradedTaylorPicard does refine its remainder, but this did not explain the issue

A first-step comparison at order 5 found approximately:

```
GradedTaylorSeries final dominant error: 4.94915e-7

GradedTaylorPicard:
after construction: 2.66024e-6
after validation:   2.09081e-6
refinement 1:       2.05149e-6
refinement 2/final: 2.04826e-6
```

Thus Ariadne's Picard refinement works, but on this step the refined Picard remainder was still about 4.1 times the GradedSeries remainder.

**Conclusion:** the hypothesis "Flow*/Picard wins because Ariadne GradedSeries lacks Picard remainder refinement" is not supported by this experiment.

---

## 6. Rejected or currently unsupported hypotheses

These should not be retried without new evidence.

### 6.1 "Just add QR to the existing physical mapping"

Rejected.

Applying QR by repeatedly rebuilding the preconditioner from the globally composed physical Taylor mapping caused remainder explosion. Persistent two-layer state is required.

### 6.2 SCALED_QR is the missing ingredient

Not supported.

A variant using

```
x = c + Q D z,  z in [-1,1]^n
```

was tested after the persistent-state architecture was in place.

It did not recover the larger step. In one controlled run both QR and SCALED_QR needed about 440 steps; SCALED_QR was slower and gave essentially the same final enclosure/error. In later variants it could require even more steps.

Keep scaling separate from the core QR design unless new evidence justifies it.

### 6.3 The transformed physical bounding box alone causes the smaller QR step

Not supported as a complete explanation.

An experiment attempted to bound the transformed ODE directly in local coordinates rather than using the interval image of the physical bounding box. This did not recover the desired step behaviour. One intermediate implementation also exposed a domain problem: casting a restricted FunctionPatch to an unrestricted interface does not remove its underlying domain checks.

The latest second-step diagnostic does nevertheless show that the QR local bounding box is wider and the transformed vector-field range is worse. This is relevant, but not yet the complete mechanism.

### 6.4 Moving-centre / Chen Approach III solves the QR step problem

Rejected for the implementation tested.

A moving-centre experiment used the decomposition

```
x = p_c(t) + A y
y' = A^{-1} ( f(p_c(t)+Ay) - p_c'(t) )
```

The resulting runs still used roughly 456 steps and became substantially slower (roughly 28--32 s in the tested configurations).

Important qualification: Chen's Approach III also uses a specific degree-by-degree construction of the deviation polynomial. Passing the transformed non-autonomous field through Ariadne's ordinary graded-series machinery is not a full reproduction of Chen's algorithm. Therefore the experiment rejects that implementation strategy, not Chen's method in general.

### 6.5 The original Ariadne-vs-Flow* comparison is already established apples-to-apples

Not established.

At one point numbers from Chen's Van der Pol benchmark (including widths around 0.6308 / 0.6120 for the cited settings) were being compared with Ariadne quantities such as final radius or Taylor-model error. Those quantities may not be identical metrics.

Do not claim that the new integrator explains a Flow* advantage until:
- the exact Flow* metric is identified,
- fixed/adaptive step semantics are aligned,
- order/cutoff/settings are aligned,
- enclosure width/radius/remainder are compared consistently.

---

## 7. Experimental dead ends and lessons

### API/type lessons encountered

Ariadne distinguishes exact floats, validated bounds, upper intervals, boxes, etc. Several failed compilation attempts came from silently assuming conversions that do not exist.

Examples:
- `Bounds<FloatDP>` is not implicitly an `Interval<UpperBound<FloatDP>>`.
- Products involving validated bounds may return `Bounds<FloatDP>`, not `FloatDP`.
- `Positive<Bounds<...>>` and `Positive<UpperNumber<...>>` do not expose the assumed `.raw()` API.
- `Vector<Graded<Differential<Bounds<FloatDP>>>>` and `Vector<Differential<Bounds<FloatDP>>>` do **not** have an `.error()` member.
- `VectorFunctionPatchElement` is not necessarily accepted where a scalar FunctionPatch is required; use the appropriate `.get(i)` API when needed.
- `graded_series_flow_step` overloads distinguish `ExactIntervalType` from `Interval<StepSizeType>`.

Do not add speculative conversions without first checking the actual Ariadne types and existing usage in the codebase.

### Build/link lesson

Calling the `jacobian_value<FloatDP>` template from the new QR code compiled but failed at link time because that instantiation was not exported. The implementation was changed to extract first-order Taylor coefficients directly instead.

### QR implementation detail

The public orthogonal decomposition was treated carefully:
- columns were explicitly normalized to obtain a true rotation,
- columns of the linear coefficient matrix were sorted by decreasing norm before QR, following the investigated Flow* strategy.

Do not assume an orthogonal-decomposition helper necessarily returns already unit-normalized columns without checking its contract.

---

## 8. Important commits / experimental milestones

The exact branch history may later be rebased, but these SHAs identify the experiments as originally performed.

- `9722ccd68698edb50648ca1e9f08ceffad537284`  
  Avoid unexported `jacobian_value` template in QR preconditioner.

- `7109203e17af3f4fd88f19923ac59cf0d351de4b`  
  Rotation-only QR experiment; removed explicit second scaling that was unstable in the then-current architecture.

- `750c585c97bb466e81d75c582b65b9f534978fd1`  
  **Major architectural milestone:** carry persistent two-layer preconditioned Taylor state across steps.

- `ed5f58abf219beaf8e4257ef7f0842563af89340`  
  Controlled IDENTITY-vs-QR comparison; added QR column ordering.

- `de48ea128ac7214d8bff263d0af4e06928b35002`  
  Fix approximate comparison in QR column ordering. This revision corresponds closely to the useful stable IDENTITY/QR architecture used as the base for later comparisons.

- `404330f8eeea16f38fc2d4cb752c8b5805cbf5c6`  
  SCALED_QR experiment with persistent state.

- `5c9411298480d708e491ee000d41313d5eb2268a`  
  Experiment bounding preconditioned flow directly in local coordinates.

- `a3ee8fba57fb77839b5aa68b037f1c574abc466d`  
  Build transformed vector field without restricted FunctionPatch domain; part of local-bound experiment.

- `af31ef56394f088b8feaa79e2c2c819a744e0f8f`  
  Moving-centre / Approach III-inspired experiment.

- `48bc7087e5eb1541917fdaeb0868c4c76e6f0ab8`  
  API fixes for moving-centre experiment.

- `cb563d583dc6adeb166b36a675e73d1be305575a`  
  Instrument first-step graded Taylor polynomial defect.

- `0909432b5ae82f4bbbf9b26e86f739fd15f195bb`  
  Compare GradedSeries and GradedTaylorPicard remainder mechanisms.

- `83f674bcce733f3cb0ceaf0ef46f65ba836b1187`  
  Equal-final-accuracy benchmark: standard GradedSeries vs preconditioned QR.

- `e3c64f0237ef0eb630dd62a67c891edab46d8b0d`  
  **Current key diagnostic:** compare IDENTITY and QR from exactly the same physical state at the second step.

- `f26f33d5b1afd3df025873f9ea0b02898b8bf554`  
  Initial attempt to trace internal graded differential error; did not compile because raw differential containers have no `.error()`.

- `cb1496eb436a1d4ed226554a4f18eaa4da39f29a`  
  Fix that diagnostic by measuring maximum stored interval-coefficient magnitude instead.

---

## 9. Current diagnostic and immediate next step

**Current question:**

Why does QR produce about 6.9x more physical Taylor-model remainder than IDENTITY at `h=0.02` on the same second-step physical state?

We already know:
- the difference exists before accumulated-state composition;
- QR widens the axis-aligned local bounding box;
- QR widens the transformed vector-field interval range;
- this makes the QR step fail the `1e-6` physical error test at `h=0.02`, while IDENTITY passes.

The current HEAD instruments the fixed-degree `graded_series_flow_step` pipeline:

```
graded_flow_init
    -> graded_flow_iterate
    -> flow_differential
    -> flow_function
```

Because the intermediate `Graded<Differential<Bounds<FloatDP>>>` objects have no explicit remainder/error field, the diagnostic records the maximum magnitude of their stored interval coefficients:

```
fdphic_coeff_mag
fdphib_coeff_mag
dphic_coeff_mag
dphib_coeff_mag
dphi_coeff_mag
tphi_errors
```

### NEXT STEP

Build and run current HEAD:

```bash
git pull
ninja vanderpol
./examples/continuous/vanderpol > graded_internal.txt 2>&1
```

Inspect `graded_internal.txt`.

The goal is to determine the earliest stage at which IDENTITY and QR diverge strongly:
- If centre data remain similar but bounding data (`fdphib` / `dphib`) inflate, the main mechanism is interval evaluation on the rotated axis-aligned bounding box.
- If inflation first appears later in `flow_differential` or `flow_function`, inspect that conversion/truncation logic instead.
- Do not modify the preconditioner again until this diagnostic is understood.

---



### 9.1 Result of coefficient-magnitude diagnostic (2026-09-23)

The diagnostic introduced by `cb1496eb436a1d4ed226554a4f18eaa4da39f29a` localised the QR inflation further.

On the same second-step physical state, at `h=0.02`:

```
IDENTITY:
fdphic_coeff_mag = 973.45623
fdphib_coeff_mag = 3616.2787
dphic_coeff_mag  = 194.69125
dphib_coeff_mag  = 723.25573
dphi_coeff_mag   = 723.25573
tphi dominant error ~ 6.35e-7

QR:
fdphic_coeff_mag = 1023.8441
fdphib_coeff_mag = 26610.952
dphic_coeff_mag  = 204.76882
dphib_coeff_mag  = 5322.1903
dphi_coeff_mag   = 5322.1903
tphi dominant error ~ 4.30e-6
```

The centre-based branch (`c`) changes by only about 5%, whereas the validated bounding-box branch (`b`) inflates by about 7.36x. The final dominant Taylor-model error inflates by about 6.8x.

This is strong evidence that the QR penalty is generated in the **bounding differential computation**, before `flow_differential` and before conversion to the Taylor function model. In `flow_differential`, highest-spatial-degree terms and the highest temporal degree are intentionally taken from `dphib`; consequently the inflated bounding derivatives directly become the certified remainder terms.

The fact that `fdphib` and `dphib` retain the same large magnitude for the QR candidates at `h=0.02`, `0.015`, and `0.01`, while `tphi_errors` shrink strongly with `h`, is consistent with the derivative enclosure being a property of the spatial bounding box and the temporal powers subsequently scaling its contribution.

### Updated NEXT STEP

Find **which temporal iteration first creates the ~7x inflation** in the bounding branch. Instrument `graded_flow_init` / each subsequent `graded_flow_iterate` separately for centre and bounding data. Compare IDENTITY and QR at the second physical step.

If the inflation is already present immediately after `graded_flow_init`, inspect direct interval evaluation of the transformed vector field and its spatial differential on the rotated box.

If it emerges only after one or more `graded_flow_iterate` calls, inspect repeated interval composition / antidifferentiation at that temporal degree.

Do not change the preconditioning strategy before this is identified.

---



### 9.2 Temporal-iteration localisation (2026-09-23)

The per-iteration diagnostic shows that the QR penalty is **not present at the first vector-field evaluation** and then grows rapidly with temporal order.

For the same second-step state:

| temporal iteration | IDENTITY `dphib` magnitude | QR `dphib` magnitude | QR / IDENTITY |
|---:|---:|---:|---:|
| 1 | 9.494 | 9.444 | 0.99x |
| 2 | 23.183 | 29.576 | 1.28x |
| 3 | 75.192 | 147.420 | 1.96x |
| 4 | 282.395 | 914.136 | 3.24x |
| 5 | 723.256 | 5322.190 | 7.36x |

The centre branch remains close throughout; at order 5 it is about `194.69` for IDENTITY and `204.77` for QR.

Also, `fdphib / dphib = k` at temporal order `k` (for example `26610.952 / 5322.1903 = 5`), so the antidifferentiation is not generating the excessive ratio. The inflation is already present in the output of `compute_procedure`.

**Conclusion:** the current bottleneck is repeated validated evaluation of the vector-field Procedure on the bounding graded differential. The next diagnostic must identify the first Procedure instruction / elementary operation that amplifies QR relative to IDENTITY.

### Updated NEXT STEP

Instrument the temporary values produced by `compute_procedure` for the same-state second-step IDENTITY and QR probes, for every temporal iteration. Record:
- Procedure instruction index and operation,
- centre/bounding branch,
- maximum stored interval-coefficient magnitude.

The key target is the bounding branch at iterations 2--5. Determine whether the amplification is associated with `sqr`, multiplication, subtraction, or another operation in the Van der Pol procedure.

Do not modify the QR preconditioner until that operation-level source is known.

---


### 9.3 Procedure-level diagnostic: rotation destroys the sparse expression structure (2026-09-23)

The output from commit `24022ccf8abba99b74dc44c55e22f68d8eb2bc18` gives a more specific explanation than "a multiplication inside `compute_procedure` is bad".

First, the IDENTITY and QR Procedures are **not the same instruction stream**, so instruction numbers cannot be compared one-for-one:

- the IDENTITY transformed field uses 38 instructions (0--37);
- the QR transformed field uses 70 instructions (0--69).

The reason is structural. With IDENTITY coordinates, the Van der Pol field retains its sparse original algebraic form. After the affine QR transformation

```
x = c + A y
y' = A^{-1} f(c + A y)
```

the same polynomial field becomes a denser mixed polynomial in both local variables. The generated Procedure therefore contains longer Horner-like multiplication chains and both output components depend nontrivially on both local variables.

This matters because the validated bounding branch repeatedly feeds interval-valued graded coefficients through those chains. The first temporal iteration is still benign:

```
IDENTITY dphib_1 = 9.4940427
QR       dphib_1 = 9.4436234
```

but the QR Procedure then recursively feeds already enlarged coefficients back through the dense polynomial. In particular, in the QR bounding branch the second local input (`x[1]`, instruction 27) evolves as

```
iteration 1:   1.000
iteration 2:   9.444
iteration 3:  29.576
iteration 4: 147.420
iteration 5: 914.136
```

and is multiplied by the fixed factor represented by instruction 30 (`c[13] ~ 3.7686`), after which it participates in further multiplication chains. By iteration 5:

```
instruction 31: mul(v[29],v[30])  ~ 3.445e3
instruction 36: mul(v[35],v[31])  ~ 4.572e3
instruction 66: mul(v[65],v[31])  ~ 2.482e4
final output                         2.661e4
```

For comparison, the final IDENTITY bounding Procedure output at iteration 5 is about `3.616e3`.

The important conclusion is therefore **not** that a single erroneous primitive operation has been found. Ordinary validated multiplication is doing what it is asked to do. The QR penalty comes from the interaction of:

1. coordinate rotation,
2. algebraic densification of the polynomial vector field,
3. axis-aligned interval bounding of graded coefficients,
4. repeated validated multiplication/composition across temporal orders.

This also explains why the centre branch remains close: point/centre coefficients preserve cancellations and correlations that interval bounding cannot retain.

The result strengthens the hypothesis that the efficiency gap is a **representation/evaluation issue**, not a defect in the QR matrix itself. A rotation can geometrically reduce wrapping of the reachable set while simultaneously making the local polynomial vector field much worse for interval-based graded-differential evaluation.

### Updated NEXT STEP

Do **not** tune the QR matrix yet.

The next experiment should test whether preserving the polynomial/correlated part of the transformed field longer avoids this interval dependency amplification. Two concrete directions are worth separating:

1. **Expression/evaluation experiment:** compare the current generated QR Procedure against an algebraically equivalent representation chosen to minimise interval dependency (e.g. retain factored/Horner structure derived from the original Van der Pol field rather than fully composing/expanding the affine transform). The purpose is to see whether the 70-instruction dense Procedure itself is responsible for most of the growth.

2. **Flow*-style Taylor-model experiment:** evaluate the transformed vector field on Taylor-model/polynomial objects and push only truncation/residual terms into interval remainders, instead of using the bounding graded differential as the principal high-order enclosure. This is closer to the mechanism described by Chen, where polynomial dependence is preserved and remainder validation/refinement is handled separately.

A useful short diagnostic before implementing either full approach is to compute the QR second step with the same physical state and same `h=0.02`, but replace the bounding Procedure evaluation by a hand-factored implementation of the transformed Van der Pol polynomial. If the order-5 `dphib` drops materially below `5322`, expression dependency is confirmed as a major cause. If it does not, the problem lies deeper in the bounding graded-differential representation itself.

Do not interpret this result as evidence that QR preconditioning is intrinsically unsuitable. It shows that **QR plus the current interval-valued graded Procedure evaluator** is a poor combination on this benchmark.

---


### 9.4 Sparse physical-Procedure experiment rejects expression densification as the main cause (2026-09-23)

Commit `1e0c8c4e752da91baea705cdd0f47c8c37f69a1b` evaluated the same QR local step through the **original sparse physical Procedure** instead of the 70-instruction dense transformed Procedure. The affine maps `x=c+A*y` and `A^{-1}` were applied directly to graded differentials, while the production path remained unchanged.

At the first step, where `A=I`, the two paths agree closely, as expected. At the important same-state second-step QR probe with `h=0.02`:

```
dense transformed Procedure:
  dphib order 5       = 5322.1903
  local errors        = [8.3626e-7, 4.2908e-6]

sparse physical Procedure:
  dphib order 5       = 5177.6996
  local errors        = [8.60e-7, 4.41e-6]
```

The per-order sparse-path bounding magnitudes were:

```
iteration 1: local_f ~ 10.078
iteration 2: local_f ~ 62.946
iteration 3: local_f ~ 478.875
iteration 4: local_f ~ 3724.779
iteration 5: local_f ~ 25888.498
```

These are essentially the same explosive progression as the dense transformed Procedure. The final Taylor-model error is even slightly worse in the sparse path. The same conclusion holds at `h=0.015` and `h=0.01`.

**Conclusion:** algebraic densification from 38 to 70 Procedure instructions is **not the main cause** of the QR penalty. The previous hypothesis in section 9.3 was useful but is now rejected by direct experiment. The dominant loss occurs in the interval-valued graded-differential representation/evaluation itself once rotated coordinates are bounded axis-aligned.

A further observation is important: the sparse diagnostic currently obtains the physical zero-order enclosure by mapping the axis-aligned local bounding box back through `A`. But that local box was itself obtained by interval-evaluating `A^{-1}` on the physical flow bound. The round trip

```
physical box -> interval(A^{-1} box) -> interval(A local box)
```

can introduce a second wrapping enlargement before the sparse physical Procedure is evaluated.

### Updated NEXT STEP

Test whether this **box round-trip**, rather than expression structure, accounts for a significant fraction of the remaining QR growth.

For the bounding branch only, initialise the zero-temporal-order physical differentials directly with the already validated `physical_bounding_box`, while retaining the affine gradient with respect to local variables given by `A`. Higher temporal coefficients continue to be transformed by `A` from the local graded state.

This is a diagnostic separation:
- if order-5 growth drops strongly, the dominant mechanism is repeated axis-aligned box conversion around the coordinate transform;
- if growth remains near `5e3`, the deeper limitation is the interval-coefficient graded differential representation under rotated dependencies.

Do not use the sparse-Procedure path as a production optimisation: this experiment showed no accuracy benefit by itself.

---


### 9.5 Direct physical bounding box helps materially, but does not remove high-order growth (2026-09-23)

Commit `8221d5c5ebc2d0b5d54977f7dbbd73402a3edd9c` tested the sparse physical-Procedure path again, but initialized the bounding branch zero-order physical differentials directly from the already validated `physical_bounding_box`, avoiding the interval round trip

```
physical box -> interval(A^{-1} box) -> interval(A local box).
```

For the same second-step QR probe at `h=0.02`, the result improved materially:

```
dense transformed Procedure:
  dphib order 5       = 5322.1903
  dominant local err  = 4.2908e-6

sparse physical Procedure, with box round-trip:
  dphib order 5       = 5177.6996
  dominant local err  = 4.41e-6

sparse physical Procedure, direct physical box:
  dphib order 5       = 3525.2330
  dominant local err  = 2.91e-6
```

The direct-physical-box bounding progression is:

```
iteration 1: local_f ~ 9.397
iteration 2: local_f ~ 51.208
iteration 3: local_f ~ 381.622
iteration 4: local_f ~ 2854.998
iteration 5: local_f ~ 17626.165
```

Thus the box round-trip accounts for a significant part of the QR penalty: order-5 `dphib` falls by about one third relative to the dense production path, and the dominant Taylor-model error falls by about 32%. The improvement persists when the candidate step is reduced:

```
h=0.015: 1.019e-6 -> 6.91e-7
h=0.010: 1.345e-7 -> 9.12e-8
```

However, the high-order coefficient growth is still severe. Even after removing this avoidable box wrapping, the bounding branch grows from roughly `9.4` at order 1 to `1.76e4` in the vector-field coefficient at order 5, corresponding to `dphib ~ 3525`. This is still far above the IDENTITY second-step `dphib ~ 723`.

**Conclusion:** there are two distinct effects.

1. The current QR implementation pays an avoidable penalty by converting a validated physical box to an axis-aligned local box and then, directly or indirectly, losing that correlation again. This should not be part of a final design.
2. Removing that penalty is not sufficient. The remaining factor of roughly 4.9 in order-5 `dphib` versus IDENTITY points to the interval-coefficient graded representation itself as the deeper limitation under rotated dependencies.

The expression-densification hypothesis remains rejected: keeping the sparse physical Procedure only becomes useful once the physical enclosure is also kept in physical coordinates.

### Updated NEXT STEP

Stop refining the dense transformed-Procedure path as the main architecture. The evidence now supports a split representation closer to Flow*:

- keep the polynomial/Taylor dependence in correlated symbolic or Taylor-model form;
- keep the validated physical flow box for range/remainder validation rather than repeatedly rotating it into axis-aligned local boxes;
- apply QR/preconditioning to the polynomial dependence, not to every interval enclosure;
- maintain and refine a separate interval remainder.

The next prototype should therefore be a **single-step polynomial + remainder experiment**, not another QR-matrix tweak. On the Van der Pol second-step state, compute the temporal Taylor polynomial using centre/polynomial coefficients, evaluate its image over the local initial Taylor model, and validate only the missing/truncation part against the physical bounding box. The immediate target is to reproduce the same `h=0.02` step while avoiding the bounding `Graded<Differential<Bounds>>` recursion that produces `dphib ~ 3525--5322`.

---


### 9.6 Clarification: the persistent affine-local state is a factorisation, not extra symbolic information (2026-09-23)

A conceptual correction is important for interpreting the improvements seen so far.

The ordinary `GradedTaylorSeriesIntegrator`, when used by the evolver, does **not** simply discard all dependence on the previous evolved set. It constructs a fresh one-step flow map `Phi_k(x,t)`, and the evolved set is obtained by composing that map with the previous Taylor-model state `X_k(s)`:

```
X_{k+1}(s) = Phi_k(X_k(s), h).
```

Therefore the previous symbolic dependence is present in the standard pipeline as well.

The preconditioned prototype instead stores the current state in the factorised form

```
X_k(s) = c_k + A_k Y_k(s),
```

where:
- `centre` is `c_k`;
- `linear_map` is `A_k`;
- `normalised_mapping` is `Y_k(s)`;
- `local_domain` is only an enclosure of the range of `Y_k`.

These objects do **not** contain more mathematical information than the corresponding full Taylor map. The benefit is representational: the affine frame is kept outside the accumulated Taylor-model mapping.

For a local endpoint transition `E_k(y)`, the implementation first preconditions the fresh local endpoint,

```
E_k(y) = c_{k+1} + A_{k+1} T_k(y),
```

and only afterwards composes its local map with the accumulated previous local map:

```
Y_{k+1}(s) = T_k(Y_k(s)),
X_{k+1}(s) = c_{k+1} + A_{k+1} Y_{k+1}(s).
```

This ordering matters numerically even though it is mathematically equivalent to composing the physical maps directly. It prevents the new affine transformation from repeatedly acting on and being swept together with the already accumulated interval remainder. Earlier experiments showed that re-preconditioning the already-composed physical map causes artificial remainder growth; preconditioning the fresh local endpoint before accumulated composition avoids that failure mode.

**Current limitation:** this factorisation is used primarily **between steps**. During construction of a new local flow, the current `graded_series_flow_step` still receives the axis-aligned `local_domain` / `local_bounding_box`, rather than exploiting the full symbolic structure of `Y_k(s)`. Thus the high-order validated recurrence inside the step remains essentially the same machinery as `GradedTaylorSeriesIntegrator`.

This distinction should be kept explicit:

- current improvement: better factorisation and ordering of inter-step representation/composition;
- not yet achieved: preserving the correlated polynomial dependence of the incoming Taylor model inside the high-order validated flow computation itself.

A Flow*-like polynomial+remainder prototype should target the second point rather than claiming that the current persistent state already preserves symbolic information absent from the standard evolver.

---


### 9.7 Centre-polynomial defect experiment prepared (2026-09-23)

The next experiment targets the key architectural question directly: is the centre/polynomial temporal expansion already a good approximation, with the large error coming mainly from the way the bounding recurrence certifies it?

A diagnostic path now constructs a **centre-only graded Taylor polynomial** by running the ordinary centre branch of the graded recurrence and using it for all retained coefficients, including the highest temporal/spatial terms. This object is not claimed to be a validated flow enclosure by itself.

For this polynomial candidate `P(y,t)`, the diagnostic computes the ODE defect

```
R(y,t) = dP/dt - g(P(y,t))
```

over the same local domain and time interval, together with the initial-condition mismatch `P(y,0)-y`.

Interpretation:

- if the defect and initial mismatch are small while the validated bounding branch still produces `dphib ~ 3.5e3--5.3e3`, then the polynomial approximation is intrinsically good and the main missing piece is a separate a-posteriori remainder validation/refinement mechanism;
- if the defect is already large, then simply separating polynomial and remainder will not be enough and the polynomial construction itself must change.

This is deliberately a diagnostic before implementing a full Flow*-style remainder iteration. A small defect would justify the next step: derive a rigorous remainder enclosure from the defect plus a Lipschitz/Jacobian bound on the validated physical flow box.

---


### 9.8 Centre-polynomial defect is small: separate remainder validation is justified (2026-09-23)

The diagnostic from `e797d6c2164adf6a9e04f6108303471c26d6057a` compiled and ran successfully.

For the important second-step QR state, the centre-only polynomial has essentially negligible construction error (about `2e-11`) and an initial-condition mismatch of the same order. At `h=0.02`:

```
polynomial_errors ~ [1.97e-11, 1.67e-11]

defect R = dP/dt - g(P):
  component 0 ~ [-1.93e-7,  2.02e-7]
  component 1 ~ [-1.09e-6,  1.39e-6]

initial mismatch:
  component 0 ~ +/-1.97e-11
  component 1 ~ +/-1.66e-11
```

When the candidate step is reduced, the defect decreases rapidly:

```
h=0.020: max |R| ~ 1.39e-6
h=0.015: max |R| ~ 3.39e-7
h=0.010: max |R| ~ 5.47e-8
```

The raw scale `h * sup|R|` at `h=0.02` is only about `2.8e-8`, before accounting for dynamical amplification. This is orders of magnitude below the current QR graded-series local remainder (about `4.3e-6`) and also below the improved direct-physical-box diagnostic (about `2.9e-6`).

**Conclusion:** the polynomial approximation itself is not the source of the large QR remainder. The evidence now directly supports the Flow*-like split architecture: retain the centre/polynomial expansion and certify a separate remainder around it, instead of replacing highest-order coefficients by the recursively inflated bounding graded differential.

This does **not** yet prove that a rigorous remainder of order `1e-8` is attainable: the defect must be propagated through a validated Jacobian/Lipschitz bound on a tube containing the exact flow. The next diagnostic therefore computes a conservative infinity-norm Lipschitz bound `L` over the already validated local flow box. With

```
epsilon = sup |dP/dt - g(P)|
e0      = sup |P(y,0)-y|
```

the first a-posteriori enclosure to test is the Gronwall bound

```
|e(h)| <= exp(L h) e0 + (exp(L h)-1)/L * epsilon.
```

If this already beats the current graded bounding remainder by a large factor, implement it as the first rigorous polynomial+remainder prototype. If it is too pessimistic, move directly to componentwise/matrix remainder refinement rather than returning to the old graded bounding recurrence.

---


### 9.9 Scalar a-posteriori Gronwall bound preserves the defect advantage (2026-09-23)

The Lipschitz diagnostic confirms that the simple scalar amplification bound is not large enough to destroy the centre-polynomial advantage.

For the important second-step QR state, the validated infinity-norm Jacobian bound over the local flow box is

```
L = 13.210077
```

and the centre-polynomial diagnostics give:

```
h=0.020:
  epsilon = sup |dP/dt-g(P)| ~ 1.3928304e-6
  e0      ~ 1.9695e-11
  exp(Lh) ~ 1.30239
  Gronwall endpoint bound ~ 3.19e-8

h=0.015:
  epsilon ~ 3.3896547e-7
  e0      ~ 1.8710e-11
  exp(Lh) ~ 1.21915
  Gronwall endpoint bound ~ 5.65e-9

h=0.010:
  epsilon ~ 5.4732665e-8
  e0      ~ 2.4786e-11
  exp(Lh) ~ 1.14122
  Gronwall endpoint bound ~ 6.13e-10
```

These values use

```
e(h) <= exp(Lh)e0 + (exp(Lh)-1)/L * epsilon.
```

At `h=0.02`, the resulting ~`3.2e-8` bound is roughly 17 times smaller than the current dense QR local-flow error (~`5.45e-7`) for this same second-step probe, even before implementing componentwise refinement. The gap grows rapidly as `h` decreases.

The first-step diagnostic has `L=11.445905` and similarly small defect, so the effect is not unique to the rotated second step. Later carried states show `L` around `13.87`, still moderate.

**Conclusion:** a separate a-posteriori remainder is now quantitatively justified. The old bounding graded recurrence is not needed to obtain a competitive local error on this benchmark, provided the defect bound and Jacobian bound are validated on a tube known to contain the exact solution.

### Updated NEXT STEP

Implement a first rigorous prototype which returns the centre polynomial plus a uniform interval remainder derived from the scalar Gronwall estimate.

There is one validation issue that must not be skipped: the Jacobian bound currently uses the pre-existing validated `local_bounding_box`, which itself comes from the ordinary flow bounder, while the defect is evaluated on the polynomial image. For a production replacement, verify that the proposed polynomial-plus-remainder tube is contained in a domain on which both the defect and Jacobian bounds were computed. Initially it is acceptable to reuse the existing validated flow box as the certification domain; this isolates the remainder mechanism without changing the bounder at the same time.

The prototype should:
1. retain the centre-only Taylor polynomial;
2. compute validated componentwise defect magnitudes;
3. compute a validated scalar `L` over the existing local flow box;
4. attach the Gronwall endpoint/tube remainder to the Taylor model;
5. compare the resulting physical-coordinate error and accepted step against the current graded-series path.

Only after this succeeds should the scalar bound be replaced by Flow*-style componentwise/fast remainder refinement.

---


### 9.10 First polynomial+remainder run exposed a containment-test API issue (2026-09-23)

Commit `9993215f9e0cb3857738f9a69fe6ff7f0821a290` ran the first prototype that would attach the Gronwall remainder to the centre polynomial. The numerical ranges show that the centre polynomial is comfortably inside the existing certification box; for example on the second-step QR probe at `h=0.02`:

```
polynomial range:
  [-0.16156445, 0.21679093]
  [-0.14619583, 0.071965785]

local certification box:
  [-0.34362090, 0.41198114]
  [-0.33444289, 0.19626299]
```

Nevertheless the prototype reported `rejected=polynomial_outside_certification_box` for every tested candidate. This is not evidence against the polynomial+remainder method: it is a comparison-type issue. Existing integrator code performs the analogous check as

```
definitely(subset(phi.range(), cast_exact_box(bx)))
```

rather than comparing a validated range directly against an `UpperBoxType`.

The next commit changes only this guard to compare against `cast_exact_box(local_bounding_box)`. No remainder formula or polynomial construction is changed. The experiment must be rerun before drawing conclusions about the attached remainder.

---


### 9.11 Attached Gronwall remainder succeeds; next experiment propagates it (2026-09-23)

After fixing the containment guard, the first actual polynomial+remainder Taylor models were produced.

On the important second-step QR state:

```
h=0.020:
  Gronwall physical error = 3.6709e-8
  dense graded QR error   = 4.3669e-6
  improvement             ~119x

h=0.015:
  Gronwall physical error = 6.3262e-9
  dense graded QR error   = 1.0372e-6
  improvement             ~164x

h=0.010:
  Gronwall physical error = 6.9795e-10
  dense graded QR error   = 1.3691e-7
  improvement             ~196x
```

The first step at `h=0.02` similarly improves from about `5.45e-7` to `3.51e-8`.

Later carried-state diagnostics are also stable: with `L ~ 13.87`, the prototype remains around `3.65e-8` at `h=0.02` while the dense graded path reaches about `9.50e-6`.

This is substantially stronger evidence than the earlier scalar estimate because the remainder has now actually been attached to the Taylor model, transformed back to physical coordinates, and its resulting model error measured.

**Next experiment:** use the certified polynomial+remainder flow for the actual step acceptance and subsequent state propagation, while still computing the dense graded flow in parallel as an A/B diagnostic. This tests whether the local gain survives composition, preconditioning, and multiple carried steps. The dense path remains as a temporary fallback only if the Gronwall candidate cannot be certified.

---


### 9.12 Multi-step propagation preserves the gain (2026-09-23)

The first run in which the certified Gronwall polynomial flow was actually propagated across steps succeeded.

Step 0 at h=0.02:

- Gronwall physical error: 3.5113e-8
- dense graded error: 5.4533e-7
- both methods accept the step
- propagated mapping error: 3.5113e-8
- normalised state error: 3.5402e-8

At the next carried QR state, again at h=0.02:

- Gronwall physical error: 3.6709e-8
- dense graded error: 4.3669e-6
- Gronwall accepts the step
- dense graded QR rejects the step

After composition and preconditioning of that second step:

- flowpipe error: 7.5228e-8
- final mapping error: 7.3449e-8
- final normalised error: 7.2455e-8

Thus the local improvement survives the actual inter-step machinery. More importantly, the new method accepts h=0.02 on the second carried step where the old graded QR path would already reduce the step.

**Next step:** stop paying for the dense graded bounding recurrence in normal execution. Keep it only under diagnostics for A/B comparison. The production candidate should be the centre polynomial plus separately certified remainder; if that candidate cannot be certified, reduce h and retry. Only after this change are runtime comparisons meaningful.

---


### 9.13 Clean production benchmark established (2026-09-23)

After gating the evolver investigation probes on the integrator diagnostics flag, the production-only Gronwall run over Van der Pol from t=0 to t=0.40 is clean:

```
elapsed_seconds = 0.700001
reach_sets       = 21
intermediate_sets= 21
```

With maximum step size 0.02, 21 reach/intermediate sets are consistent with taking the full requested step throughout this run. This is now a usable timing baseline because the dense graded A/B path and the IDENTITY/QR investigation probes are absent.

The next experiment runs the ordinary `GradedTaylorSeriesIntegrator` immediately after the Gronwall integrator in the same executable, with the same:
- Van der Pol initial set and horizon 0.40;
- `StepMaximumError(1e-6)`;
- threshold sweeper `1e-12`;
- spatial/temporal order fixed at 5;
- maximum evolver step 0.02;
- enclosure radius and spatial-error configuration;
- reconditioning disabled.

The output marker `[IntegratorBenchmark]` reports elapsed time and set counts for both methods. This gives the first controlled speed/accepted-step comparison without relying on separate process timings.

---


### 9.14 First controlled speed comparison: accuracy gain does not yet imply per-step speed gain (2026-09-23)

The controlled t=0.40 benchmark with identical order/tolerance/evolver configuration produced:

```
GRONWALL: 0.706001 s, 21 reach sets
GRADED:   0.221001 s, 21 reach sets
```

At a forced maximum step of 0.02 both methods complete with the same number of sets, so the current Gronwall prototype is about 3.2x slower in wall-clock time on this short benchmark despite its much smaller local remainder.

This is not surprising architecturally: the new path currently pays for
- a centre graded recurrence;
- construction and range evaluation of the polynomial ODE defect `dP/dt-g(P)`;
- a first-order differential/Jacobian evaluation on the certification box;
- QR/state factorisation and the specialised carried-state evolver path.

The old graded integrator pays for the centre+bounding recurrence but no separate residual composition/range pass.

The result changes the immediate optimisation question. The main potential speed advantage of the new method is not cheaper work per fixed h=0.02 step; it is the ability to take larger validated steps because its remainder is much smaller. A fair next experiment must therefore remove the artificial common maximum-step bottleneck.

The next benchmark extends the horizon to t=1.0 and compares both methods at max_step=0.02, then runs the Gronwall method at max_step=0.04. This tests whether the accuracy headroom can be converted into fewer steps. If 0.04 is accepted robustly, follow with the ordinary graded method at 0.04 to determine whether its local-error criterion forces reductions.

---


### 9.15 Larger steps convert accuracy headroom into fewer steps, but current cost remains dominant (2026-09-23)

On the t=1.0 benchmark:

```
GRONWALL max_step=0.02: 1.45701 s, 50 reach sets
GRADED   max_step=0.02: 0.499001 s, 51 reach sets
GRONWALL max_step=0.04: 1.32101 s, 28 reach sets
```

Increasing the Gronwall maximum step from 0.02 to 0.04 reduces the number of reach sets from 50 to 28, confirming that the smaller certified remainder can be converted into substantially larger accepted steps. Runtime, however, improves only from 1.46 s to 1.32 s. The current per-step cost of residual construction/range evaluation, Jacobian evaluation, QR factorisation and carried-state composition therefore dominates enough that halving the nominal step count is not yet sufficient to beat the old graded method.

The 28-set result also shows that h=0.04 is not accepted uniformly; the method reduces the step on some states. Adaptive step choice will matter.

**Next experiment:** run the ordinary graded integrator with the same max_step=0.04. This distinguishes two possibilities:
- if graded is forced close to its previous ~0.02 step count, the new method has a real step-size advantage and optimisation should focus on reducing per-step residual/Jacobian cost;
- if graded also takes near-0.04 steps, then the current tolerance/benchmark regime does not expose the accuracy advantage strongly enough for speed comparison, and tighter tolerances should be tested.

---


### 9.16 Gronwall has a real step-size advantage at max_step=0.04, but is still slower (2026-09-23)

The direct larger-step comparison gives:

```
tolerance 1e-6, horizon 1.0:

GRONWALL max_step=0.02: 1.48901 s, 50 reach sets
GRADED   max_step=0.02: 0.515001 s, 51 reach sets

GRONWALL max_step=0.04: 1.35201 s, 28 reach sets
GRADED   max_step=0.04: 0.746001 s, 45 reach sets
```

Thus the new method has a genuine accepted-step advantage: with the same 0.04 ceiling it needs 28 sets versus 45 for the ordinary graded integrator. The graded method is being forced to reduce h much more often by its larger local remainder.

However, the current prototype is still about 1.8x slower overall at max_step=0.04. The per-step overhead is therefore the dominant performance problem now; accuracy/step-size behaviour is already moving in the desired direction.

The next experiment tightens `StepMaximumError` from 1e-6 to 1e-8 at max_step=0.04 for both methods. This directly probes the original scaling question: if the graded remainder hits its accuracy floor sooner while the separate Gronwall remainder continues to scale, the reach-set count ratio should widen substantially. If both methods simply shrink h with similar asymptotics, optimisation rather than accuracy architecture becomes the main task.

---


### 9.17 Benchmark fairness requires a longer horizon (2026-09-23)

A short horizon can understate an important cost of the persistent Taylor-model approach: as evolution proceeds, the carried symbolic map is repeatedly composed with fresh local transitions, so the cost of later steps may differ substantially from the cost of early steps. A benchmark over t=0.4 or t=1.0 therefore does not yet characterise steady multi-step behaviour.

The benchmark has been changed to a longer horizon, t=5.0, and stripped down to the comparisons that now matter:

```
tolerance = 1e-6, max_step = 0.04:
  GRONWALL
  GRADED

tolerance = 1e-8, max_step = 0.04:
  GRONWALL
  GRADED
```

The output marker is `[IntegratorLongBenchmark]` and reports elapsed time and reach/intermediate-set counts.

Reconditioning remains disabled deliberately for this experiment. This exposes the raw long-term cost of the two propagation architectures instead of allowing periodic reconditioning to reset the symbolic complexity. If one method becomes pathological purely because reconditioning is disabled, that itself must be recorded; a later benchmark can then re-enable a matched reconditioning policy for both methods.

The long-horizon result should be interpreted along two axes:
- accepted-step efficiency: number of reach sets;
- accumulated per-step cost: elapsed time divided by the number of sets, and how total runtime scales relative to the earlier t=1.0 benchmark.

---


### 9.18 Long-horizon benchmark confirms both better scaling and increasing per-step cost (2026-09-23)

At horizon t=5.0 and max_step=0.04:

```
tolerance 1e-6:
  GRONWALL: 5.805 s, 134 sets  (~43.3 ms/set)
  GRADED:   3.485 s, 169 sets  (~20.6 ms/set)

tolerance 1e-8:
  GRONWALL: 14.592 s, 215 sets (~67.9 ms/set)
  GRADED:    8.643 s, 401 sets (~21.6 ms/set)
```

The accuracy scaling advantage is robust and increases at tighter tolerance: at 1e-8 the Gronwall path uses about 46% fewer sets. But the long-horizon result also confirms that its average cost per set grows substantially: roughly 43 ms/set at 1e-6 and 68 ms/set at 1e-8, compared with about 21 ms/set for the ordinary graded path.

Relative to the t=1.0 runs, this indicates that the persistent/composed representation and/or the residual certification work becomes increasingly expensive as the carried state evolves. The next task is therefore profiling, not another blind algorithmic change.

A lightweight production profiler has been added around three major local certification phases:
- centre graded polynomial construction;
- residual construction/range evaluation;
- Jacobian/Lipschitz evaluation.

It emits cumulative `[GronwallCostProfile]` lines every 100 candidate calls with diagnostics disabled. Comparing these cumulative costs with total wall time will tell us whether the dominant growth is inside local certification or outside it (composition/preconditioning/evolver state propagation). This distinction determines the next optimisation target.

---


### 9.19 Profiling identifies residual construction/range as the dominant local-certification cost (2026-09-23)

The first production profile gives a clear result.

At the first 100 Gronwall candidate calls:

```
centre polynomial: 0.732 s
residual:          1.385 s
Jacobian:          0.0027 s
```

During the tighter-tolerance run, cumulative totals by 700 calls are:

```
centre polynomial: 4.137 s
residual:          9.362 s
Jacobian:          0.0175 s
```

Thus the Jacobian/Lipschitz calculation is negligible. The dominant measured local-certification cost is constructing and ranging the generic function-patch residual `dP/dt-g(P)`, roughly 2.3 times the centre-polynomial construction cost by 700 calls.

The cumulative profiler is shared across the two Gronwall benchmark instances, so absolute totals at calls 200--700 include both tolerance runs. This does not affect identification of the dominant phase, but a later profiler should be per-integrator if exact per-run attribution is needed.

There is still meaningful unaccounted wall time. The next profile should measure the physical affine reconstruction and then carried-map composition/preconditioning. The optimisation direction suggested by this result is already clear: do not spend effort optimising the Jacobian bound. A high-value Flow*-like change is to avoid constructing/ranging the residual as a generic composed Taylor-function patch, deriving the defect/remainder directly from Taylor coefficients/recurrence or using specialised polynomial arithmetic.

---


### 9.20 Carried-state profiling: residual remains the largest cost, but flowpipe composition is also material (2026-09-23)

The carried-state profile resolves the previous uncertainty.

For the first Gronwall run, by 100 propagated steps:

```
flowpipe composition: 1.273 s
endpoint composition: 0.359 s
state composition:    0.324 s
preconditioning:      0.0058 s
state range:          0.00045 s
```

At 300 propagated steps cumulatively across the two Gronwall runs:

```
flowpipe composition: 3.415 s
endpoint composition: 1.012 s
state composition:    0.943 s
preconditioning:      0.017 s
state range:          0.0014 s
```

At roughly the same stage the local-certification cumulative profile reaches:

```
centre polynomial: 3.863 s
residual:          8.932 s
Jacobian:          0.0168 s
```

So the dominant single measured phase is still residual construction/evaluation, but the three composition paths together are also significant. Preconditioning itself and range extraction are negligible.

The next profile splits the residual into:
- `compose(g,P)`;
- derivative/subtraction assembly;
- initial-condition defect;
- explicit range evaluation;

and separately measures physical affine reconstruction. This is needed before changing the algorithm, because a dominant `compose(g,P)` suggests replacing generic composition with a recurrence/coefficient-level defect, while a dominant `range()` suggests a cheaper specialised enclosure may be sufficient.

---


### 9.21 Residual split pinpoints generic composition as the optimisation target (2026-09-23)

The split residual profile is decisive. At 700 cumulative Gronwall candidate calls:

```
total residual:          8.960 s
  compose(g,P):          8.234 s   (~91.9% of residual)
  initial defect:        0.690 s
  derivative/subtraction:0.036 s

explicit range:          0.024 s
physical reconstruction:0.053 s
Jacobian:                0.017 s
centre polynomial:       3.872 s
```

Therefore neither `range()`, derivative assembly, Jacobian evaluation nor physical affine reconstruction explains the performance gap. The overwhelmingly dominant local-certification operation is the generic function composition `compose(g,centre_polynomial)`.

This materially narrows the Flow*-like optimisation target. The next implementation experiment should avoid generic `FunctionPatch` composition for the residual. Since the centre polynomial is generated by the same graded Taylor recurrence for the ODE, the high-order defect should be recoverable from recurrence/coefficient information without re-evaluating the whole vector field through a generic Taylor-model composition. A specialised residual evaluator that executes the vector-field procedure directly on polynomial coefficients is the preferred direction.

Carried flowpipe composition remains the second major cost (~3.44 s cumulatively at 300 propagated steps), but it is smaller than residual composition and is part of producing the actual reach flowpipe. Optimise the residual kernel first; then re-profile before redesigning persistent-state composition.

---


### 9.22 Recurrence-level residual prototype (2026-09-23)

Profiling identified generic `compose(g,P)` as the dominant residual-certification cost. The next experiment computes an ODE-defect candidate directly from the graded recurrence.

After the final centre iteration, `dphic` contains the retained polynomial state. The prototype executes the vector-field `ValidatedProcedure` once more directly on this graded state, avoiding generic Taylor-function composition, and compares the resulting recurrence-level residual with the existing generic defect. The generic residual remains the production certification path until equivalence/enclosure behaviour has been checked.

The marker `[RecurrenceResidualProfile]` reports cumulative construction time and residual range every 100 calls. In diagnostics mode, `recurrence_defect_range` is printed beside the existing `defect_range`.

Two gates must be passed before replacing the generic path:
1. numerical/enclosure agreement must be understood and justified;
2. the recurrence path must be substantially cheaper than generic `compose(g,P)`.

---


### 9.24 Degree-by-degree procedure evaluation reproduces g(P), but is not yet faster enough (2026-09-24)

The same-step probe now shows exact agreement at printed precision:

```
generic_field_range    = [{2.1514581:2.3727625},{-5.1240475:-2.4433552}]
recurrence_field_range = [{2.1514581:2.3727625},{-5.1240475:-2.4433552}]
```

Thus evaluating the vector-field Procedure degree-by-degree on the graded centre state is semantically consistent with generic `compose(g,P)` at the range level.

However, the straightforward implementation is only moderately cheaper: by 700 calls it costs about 4.85 s versus 8.18 s for generic composition (~1.7x), and because both are currently executed in parallel the benchmark becomes slower. This cost is expected because the prototype rebuilds argument prefixes and calls `compute_procedure` once per degree.

The next correctness gate constructs

```
recurrence_defect = derivative(P) - recurrence_field
```

using the ordinary Taylor-model derivative but the recurrence-evaluated field, and compares its range directly with the existing generic defect on the same step. If these agree, the generic `compose(g,P)` can be replaced without changing the mathematical certification architecture.

After correctness is established, optimise the recurrence field evaluation by retaining/reusing the incremental Procedure state already generated while constructing the centre polynomial, rather than replaying all degrees in a second pass. That is the route expected to recover most of the potential speedup.

---


### 9.25 Recurrence field passes the defect-level correctness gate; switch production certification (2026-09-24)

On the identical first step, using the same Taylor-model derivative, the generic and recurrence-based defects are nearly identical:

```
generic defect:
  x: [-1.3388421e-7,  1.7941606e-7]
  y: [-9.8547251e-7,  1.3952296e-6]

recurrence defect:
  x: [-1.3388421e-7,  1.7941606e-7]
  y: [-9.6209909e-7,  1.3814724e-6]
```

The field ranges themselves are identical at printed precision. The small difference in the second defect component comes from the different Taylor-model construction/sweeping route, not from a macroscopic semantic mismatch. It is tiny relative to the failed earlier order-one recurrence defect and is slightly narrower in this sample.

This passes the practical correctness gate for an experiment using the recurrence-evaluated field in production Gronwall certification. The generic `compose(g,P)` is now removed from the production path and retained only in diagnostics mode for comparison.

Important: the current degree-by-degree recurrence evaluation still replays the Procedure after centre-polynomial construction. Therefore this commit tests the end-to-end benefit of replacing generic composition, not the final intended optimisation. If performance improves but remains insufficient, the next step is to retain the final incremental Procedure state from centre construction so that `g(P)` does not require a replay.

---


### 9.26 Production recurrence field yields a substantial end-to-end speedup (2026-09-24)

Replacing generic `compose(g,P)` in production certification with the recurrence-evaluated field preserves essentially the same integration behaviour and substantially reduces runtime.

At t=5, max_step=0.04:

```
tolerance 1e-6:
  previous Gronwall baseline: ~5.8 s, 134 sets
  recurrence production:      4.903 s, 135 sets
  graded:                     3.432 s, 169 sets

tolerance 1e-8:
  previous Gronwall baseline: ~14.6 s, 215 sets
  recurrence production:      11.680 s, 215 sets
  graded:                      8.728 s, 401 sets
```

The one-set change at 1e-6 is consistent with the slightly different/narrower recurrence defect enclosure seen in the same-step comparison; the tight-tolerance count is unchanged. The local first-step Gronwall error also changes only slightly, from about 3.51e-8 to 3.48e-8.

The generic residual-composition timer is now zero in production, confirming that the intended bottleneck has been removed. The remaining recurrence-field replay itself costs about 4.66 s cumulatively at 700 calls, and centre-polynomial construction about 8.51 s. Flowpipe composition remains the next independent large cost.

The next optimisation removes the replay. At exit from the centre Taylor recurrence, `dphic` is the final P_m while `fdphic/tmpdphic` retain the incremental Procedure state used to generate P_m from the previous degree. Because `compute_procedure` is incremental, one additional update on the retained state should append the missing final degree of `g(P_m)`. This should replace the current full degree-by-degree replay with a single incremental Procedure update. Correctness must again be checked against the generic field/defect diagnostic before trusting the timing.

---


### 9.27 Retaining Procedure state is correct and improves runtime further (2026-09-24)

The retained-state optimisation passes the same-step semantic check: generic and recurrence field ranges remain identical at printed precision, and the defect ranges are unchanged from the validated recurrence-production experiment.

End-to-end at t=5:

```
tolerance 1e-6:
  retained-state Gronwall: 4.734 s, 135 sets
  previous recurrence:     4.903 s, 135 sets
  graded:                  3.472 s, 169 sets

tolerance 1e-8:
  retained-state Gronwall: 10.891 s, 215 sets
  previous recurrence:     11.680 s, 215 sets
  graded:                   8.792 s, 401 sets
```

The optimisation is therefore valid and useful, though smaller than initially hoped. The existing `RecurrenceResidualProfile` still reports ~3.47 s at 700 calls, which means that timer includes not only the single incremental `compute_procedure` update but also conversion of the resulting graded field through `differential(...)` and `flow_function(...)` into a Taylor function patch.

The next measurement splits those two costs:
- retained incremental Procedure update;
- graded-differential -> Taylor-function conversion.

This matters because if conversion dominates, further optimising `compute_procedure` is the wrong target. In that case the defect/remainder should remain in graded/coefficient form longer and avoid materialising a full field Taylor patch.

---


### 9.28 Conversion, not the retained Procedure update, dominates recurrence-field cost (2026-09-24)

The split profile at 700 candidate calls is:

```
retained compute_procedure update: 0.824 s
graded -> Taylor conversion total:  2.587 s
```

Thus about 76% of the remaining recurrence-field materialisation cost is conversion, not vector-field evaluation. The end-to-end benchmark remains stable:

```
1e-6: GRONWALL 4.645 s / 135 sets; GRADED 3.409 s / 169 sets
1e-8: GRONWALL 10.715 s / 215 sets; GRADED 8.907 s / 401 sets
```

This rules out further optimisation of `compute_procedure` as the immediate priority. The next measurement splits conversion into:
- `differential(final_f,n,so,to)`;
- `flow_function(...)` materialisation/sweeping.

If `flow_function` dominates, a promising design is to form/range the defect directly from graded differentials and only materialise the final physical flow polynomial. If `differential` dominates, the graded-to-multivariate coefficient extraction itself needs a specialised residual path.

Flowpipe composition (~3.95 s cumulatively at 350 propagated steps) is now comparable to or larger than any single local-certification subphase and remains the other major optimisation target.

---


### 9.29 flow_function materialisation is the recurrence-field conversion bottleneck (2026-09-24)

At 700 candidate calls:

```
retained compute_procedure: 0.866 s
differential extraction:    0.00145 s
flow_function:              2.717 s
```

So `differential(...)` is effectively free; almost the entire conversion cost is `flow_function`, which constructs a Taylor function model on the widened time domain and then restricts it to the step domain.

This strongly supports keeping the residual in coefficient/graded form rather than materialising `g(P)` as a standalone Taylor patch. Before implementing that larger change, the next profile splits `flow_function` itself into:
- `make_taylor_function_model(...)`;
- `restriction(...)`.

If model construction dominates, the direct graded-defect route is the right target. If restriction dominates, a cheaper domain construction or direct creation on the final domain may recover much of the cost with less architectural change.

The long-run timing fluctuates between runs (this run: 11.196 s vs 9.008 s at 1e-8), so optimisation decisions should be based primarily on cumulative internal timers and set counts, not sub-second wall-clock differences between individual runs.

---


### 9.30 restriction dominates flow_function; test direct final-domain materialisation (2026-09-24)

The `flow_function` profile is decisive. Across all calls, restriction dominates model construction by roughly 4--5x. For example at 2000 calls:

```
make_taylor_function_model: 1.017 s
restriction:                4.671 s
```

and at 3500 calls:

```
make_taylor_function_model: 1.483 s
restriction:                7.210 s
```

For the recurrence-field path specifically, `flow_function` costs about 2.42 s at 700 candidate calls, while the retained Procedure update costs only 0.82 s.

The standard `flow_function` deliberately builds on a widened time domain `[t-h,t+h]` and then restricts to the actual step domain `[t,t+h]`. That construction is appropriate for the original flow-model path, but for the auxiliary residual field it may be unnecessary: the field is used only to form/range the defect on the actual step domain.

The next experiment therefore leaves centre-polynomial construction unchanged but materialises the recurrence field directly on `join(domx,domt,doma)`, bypassing `restriction`. This is a semantic experiment, not yet an assumed-safe optimisation. The existing diagnostics compare its field and defect ranges with the generic path on the identical first step. If those remain valid enclosures and the accepted-step counts remain stable, this cheaper final-domain construction can replace the widened-domain route for residual certification.

---


### 9.31 Direct final-domain materialisation is invalid; revert (2026-09-24)

The direct-domain experiment fails the semantic gate immediately.

On the identical first step:

```
generic_field_range:
  [{2.1514581:2.3727625},{-5.1240475:-2.4433552}]

direct recurrence_field_range:
  [{2.2007242:2.3979110},{-5.0048405:-2.3084512}]
```

and the corresponding defect explodes from roughly 1e-7--1e-6 to order 1e-2--1e-1:

```
generic defect:
  x ~ [-1.34e-7, 1.79e-7]
  y ~ [-9.85e-7, 1.40e-6]

direct recurrence defect:
  x ~ [-5.07e-2,-2.38e-2]
  y ~ [-1.57e-1,-9.75e-2]
```

The Gronwall remainder therefore jumps to about 3.9e-3 on h=0.02 and the step is rejected. The adaptive search is then forced down to h=0.0003125 in many states, explaining the runaway runtime. The widened-domain construction plus restriction is not an incidental implementation detail here; it is part of how the Taylor coefficients are interpreted/scaled by `make_taylor_function_model`.

This experiment is reverted.

The profiling result remains useful: restriction is expensive, but bypassing it by changing the domain is not semantically valid. The next optimisation should instead avoid materialising the auxiliary field Taylor patch entirely. Since `differential(final_f,...)` is essentially free, the promising route is to construct the defect directly at the Differential/graded level and only create the minimal enclosure needed by Gronwall, rather than converting `g(P)` into a full function patch and then subtracting/ranging it.

---


### 9.32 Direct Differential-level defect range experiment (2026-09-24)

The direct final-domain field experiment was invalid because changing the model domain changes the coefficient scaling. The next experiment keeps the correct widened-domain scaling but avoids restricting the auxiliary field patch.

The candidate path is:

```
dphi = Differential representation of P
field = Differential representation of g(P)
defect = derivative(dphi,time) - field
wide_defect = make_taylor_function_model(defect, domx x [t-h,t+h])
range = evaluate(wide_defect, normalised box with time in [0,1])
```

The key observation is that, in the widened time model, the actual forward step `[t,t+h]` is exactly the normalised half-interval `[0,1]`. Therefore the defect can be ranged directly on that sub-box without constructing a restricted Taylor function patch.

This experiment runs in parallel with the current rigorous production path. It reports `direct_differential_defect_range` in the same-step diagnostic and `direct_defect_seconds` in `[RecurrenceResidualProfile]`.

Acceptance gates:
1. its range must safely contain or match the current recurrence defect range on the identical step;
2. the cost must be materially below the current `flow_function` materialisation;
3. production remains unchanged until those checks pass.

---


### 9.33 Differential-level defect ranging passes first gate and is much cheaper; switch production range source (2026-09-24)

The direct Differential-level defect experiment is successful on the identical first step:

```
generic defect:
  x [-1.3388421e-7, 1.7941606e-7]
  y [-9.8547251e-7, 1.3952296e-6]

recurrence-patch defect:
  x [-1.3388421e-7, 1.7941606e-7]
  y [-9.6209909e-7, 1.3814724e-6]

direct Differential defect range:
  x [-1e-9, 1.76e-7]
  y [-1e-8, 1.39e-6]
```

The direct range is not a superset of the generic patch range; it is a different, substantially tighter enclosure obtained by forming the algebraic defect before Taylor-model materialisation and evaluating the resulting validated model directly on the forward half-box. Because all operations remain validated, this is a candidate rigorous enclosure of the same defect, but it must be validated operationally by checking acceptance behaviour and long-run set counts.

Cost is strongly favourable. At 700 calls:

```
current recurrence-field flow_function: 2.482 s
direct Differential defect path:        0.588 s
```

roughly a 4.2x reduction for this subphase. The direct path itself still includes `make_taylor_function_model`; the expensive `restriction` is absent.

The next production experiment uses `direct_defect_range` as the Gronwall forcing bound while retaining the old recurrence-field patch solely for diagnostics/profiling comparison. This deliberately isolates the semantic/step-size effect before deleting the old path. If accepted-step counts improve or remain stable and no enclosure failures appear, the old recurrence-field `flow_function` can then be removed entirely, recovering its ~2.5 s/700-call cost.

---


### 9.34 Rigour caveat: do not use the direct Differential defect in production yet (2026-09-24)

The direct Differential-level defect path is computationally promising, but the previous commit switched it into production too early. A close inspection of the representations shows that numerical agreement and use of interval arithmetic are not by themselves a proof that the resulting range encloses the full ODE defect required by the Gronwall argument.

There are two distinct issues to settle:

1. `final_f` is a finite graded/Differential representation of `g(P)`. We must prove that all terms omitted by the spatial/temporal truncation are either absent for the relevant vector field or are enclosed by an explicit remainder. The generic Taylor-model composition has its own truncation/error machinery; a coefficient-level recurrence cannot silently assume the omitted terms are zero.

2. The existing Taylor-model diagnostic `derivative(centre_polynomial,...)` is not a suitable proof oracle for this question: `TaylorModel::differentiate` calls `clobber()`, discarding the model's uniform error before differentiating. Therefore the close agreement of `direct_differential_defect_range` with the patch-based defect is useful diagnostically but does not constitute a proof of rigour.

For these reasons the production Gronwall forcing range is reverted to the previously used patch-based path. The direct Differential path remains diagnostic only.

The next rigorous route should construct the centre polynomial and its derivative directly from the same validated coefficient representation, and attach an explicit enclosure for the omitted tail of `g(P)`. Only after that remainder is accounted for can the direct coefficient-level defect replace the generic/patch path while preserving a validated Gronwall certificate.

---


### 9.35 Rigorous polynomial-field residual experiment (2026-09-24)

The unresolved issue in the cheap Differential-level residual is truncation of `g(P)`. For polynomial vector fields this can be removed exactly rather than estimated.

A Procedure degree analyser now accepts only operations that preserve polynomiality:
constants/variables, add/subtract, multiply, division by a degree-zero quantity, sign/halve, square, and non-negative integer power. Any transcendental, reciprocal of a non-constant, root, min/max, etc. rejects the exact-polynomial path.

If the vector field has algebraic degree `q` and the retained polynomial `P` has Differential degree `d`, the Procedure is re-evaluated on a copy of `P` padded to degree `q*d`. Differential multiplication can then no longer discard any polynomial composition term. The defect

```
dP/dt - g(P)
```

is formed at that full degree, converted once on the correct widened time domain, and evaluated on the normalised forward half-box. All coefficients remain validated intervals.

For Van der Pol the Procedure is cubic, so with the current centre polynomial degree 10 the exact composition degree is 30. This path is diagnostic only until its output and cost are measured.

This is materially different from the previous cheap truncated Differential path: for an accepted polynomial Procedure there is no unrepresented algebraic tail of `g(P)`. The remaining rigour question is then only whether the polynomial candidate represented by `dphi` is exactly the candidate around which the final Taylor patch/remainder certificate is constructed; its conversion/sweeping error is already carried by the output Taylor model and must not be differentiated as though it were a smooth error function.

---


### 9.36 Full exact polynomial composition is rigorous-looking but computationally non-competitive (2026-09-24)

The degree-complete polynomial diagnostic works and confirms that Van der Pol is recognised as polynomial:

```
exact_polynomial_available=true
exact_polynomial_degree=30
```

On the identical first step:

```
cheap truncated Differential defect:
  x [-1e-9, 1.76e-7]
  y [-1e-8, 1.39e-6]

full degree-30 polynomial defect:
  x [-1e-9, 1.77e-7]
  y [-2e-8, 1.42e-6]
```

So the omitted algebraic tail is small on this example, but non-zero. This confirms that the cheap truncated path cannot simply be declared exact.

Performance is unacceptable for production. At 100 candidate calls the full exact-polynomial diagnostic already costs 10.75 s; at 600 calls it costs 66.93 s. It also inflates the enclosing centre-polynomial timer because the diagnostic currently executes inside the centre helper. This explains the runaway benchmark (23.56 s already for the 1e-6 run, versus ~4.6 s before the exact diagnostic).

Conclusion: full dense degree-q*d expansion is a useful rigour oracle, not a viable integrator kernel. It should not be computed in normal runs.

The implementation is therefore changed so the degree-complete exact-polynomial path runs only with diagnostics enabled. Production returns to the previously measured cost.

The next promising rigorous direction is **not** to materialise all degree-30 coefficients, but to bound only the omitted tail beyond the retained Differential degree. For a polynomial Procedure this can be done instruction-by-instruction using a split representation:
- retained polynomial coefficients up to degree d;
- a scalar/interval tail magnitude for degrees > d.

Addition/subtraction combine tails additively. Multiplication combines retained-retained overflow plus retained-tail and tail-tail bounds. This is analogous to a Taylor-model algebra and can provide the missing rigorous correction to the cheap direct defect without constructing the full high-degree expansion.

---


### 9.37 Course correction: the integrator must remain general, not polynomial-specialised (2026-09-24)

The degree-complete Van der Pol experiment was a **diagnostic oracle only**. It answered one narrow question: whether the cheap truncated Differential residual was missing a non-zero composition tail. It did; the omitted tail was small but non-zero. The degree-30 construction is therefore useful evidence, but it is **not** the intended production architecture.

The integrator target remains general validated ODE dynamics, including non-polynomial vector fields. Van der Pol is only the current benchmark. We should not introduce an algorithmic dependency on polynomial degree or a polynomial-only residual certificate.

The actual problem to solve is more general:

> Given a validated centre polynomial/Taylor representation P and a general Procedure g, compute a rigorous and cheap enclosure of the defect dP/dt - g(P), including the effect of all truncation/remainder terms, without paying for a full generic FunctionPatch composition plus restriction.

The promising direction is therefore a **general validated truncated algebra** for Procedure evaluation, not a polynomial-tail algebra. Conceptually each intermediate Procedure value should carry:
- a retained polynomial/Differential part up to the chosen degree/order;
- an explicit validated remainder enclosure for everything not represented in that retained part.

For algebraic operations, the remainder is propagated by standard Taylor-model rules. For general elementary operations (reciprocal, exp, log, sin, cos, etc.), the remainder must be obtained from validated range/derivative bounds or the same elementary-function remainder machinery already used by Ariadne's Taylor models. This is the general analogue of what Flow*-style Taylor-model evaluation needs: truncate aggressively, but account rigorously for the discarded part.

The polynomial degree-30 diagnostic should remain optional and diagnostic-only as a reference check on polynomial examples. It must not drive API design, step selection, or the production residual algorithm.

Immediate next investigation:
1. inspect Ariadne's existing TaylorModel elementary-operation implementation and Procedure execution machinery to identify reusable validated remainder propagation;
2. prototype a lightweight Procedure evaluator whose value type is a retained Differential plus a scalar/vector validated remainder;
3. compare its defect enclosure and cost against the current recurrence-field Taylor-patch path on Van der Pol;
4. then test on at least one genuinely non-polynomial continuous example before considering the residual kernel successful.

This supersedes the polynomial-tail-only production direction suggested at the end of section 9.36.

---


### 9.38 General TaylorModel residual experiment: no polynomial assumption (2026-09-24)

The next residual experiment deliberately removes the polynomial specialisation.

The candidate P is defined as the polynomial part of the already-built centre Taylor patch. Its uniform enclosure error is clobbered **only to define the differentiable candidate polynomial itself**; that error is not treated as a differentiable function.

The unrestricted vector field `g` is then evaluated directly on the underlying `Vector<ValidatedTaylorModelDP>`:

```
P_models = polynomial models of centre_polynomial
G_models = g(P_models)
R_i      = derivative(P_i,t) - G_models[i]
```

This is a general path: `Function` dispatches TaylorModel arguments through Ariadne's validated TaylorModel elementary algebra. Operations such as multiplication and elementary transcendental functions carry validated remainder/truncation information in the TaylorModel error term. No algebraic degree of `g` is assumed.

This also avoids the expensive generic FunctionPatch `compose(g,P)` and avoids materialising a separate field patch followed by `restriction`.

The experiment is diagnostic-only for now. It reports:
- `general_tm_defect_range` beside the existing defect diagnostics;
- `[GeneralTaylorModelResidualProfile]` cumulative cost every 100 calls.

Correctness gate: for the same candidate polynomial, the resulting TaylorModel must be a validated enclosure of `dP/dt-g(P)`. Before production use we still need to check that the Function/TaylorModel call path indeed preserves the expected validated remainder semantics for all supported elementary operations and that the candidate used for the final flow is exactly the same polynomial whose defect is certified.

The polynomial degree-30 path remains only an optional oracle and is no longer part of the architectural direction.

---


### 9.39 Direct general TaylorModel evaluation is rigorous-looking but too slow as implemented (2026-09-24)

The same-step diagnostic is encouraging:

```
generic patch defect:
  x [-1.3388421e-7, 1.7941606e-7]
  y [-9.8547251e-7, 1.3952296e-6]

general TaylorModel defect:
  x [-1.34e-7, 1.80e-7]
  y [-9.9e-7, 1.40e-6]
```

The ranges agree at the displayed precision. Unlike the Differential-only experiment, this path evaluates the unrestricted general vector field directly in Ariadne's validated TaylorModel algebra, so non-polynomial elementary operations can carry Taylor-model remainder information.

However, performance is unacceptable in its naive form. At 700 candidate calls the general TaylorModel residual costs 8.56 s cumulatively, versus:
- retained Procedure update: 0.85 s;
- recurrence field flow_function: 2.51 s;
- cheap truncated Differential defect: 0.60 s.

The extra experiment raises the 1e-8 Gronwall benchmark to 18.4 s while preserving the same 215 sets. Therefore evaluating `g(candidate_polynomial.models())` from scratch as a full TaylorModel expression on every candidate is not the desired Flow*-like kernel.

This result refines the architectural target. We need the **remainder semantics of TaylorModel**, but the **incremental coefficient evaluation of compute_procedure**. The likely useful design is to augment the existing graded Procedure recurrence with a lightweight validated remainder propagated alongside the retained coefficients, reusing already-computed Procedure temporaries instead of re-evaluating the full function in TaylorModel algebra.

The full general-TaylorModel residual is now diagnostics-only, like the degree-30 polynomial oracle. It remains a useful reference enclosure for future lightweight remainder implementations, including non-polynomial test systems.

Next step: inspect the internal TaylorModel arithmetic, especially multiplication and elementary-series composition, and identify the minimal scalar error state required to reproduce its rigorous truncation bounds alongside the existing Graded<Differential> Procedure evaluation.

---


### 9.40 Lightweight validated remainder prototype: start from generic multiplication semantics (2026-09-24)

Inspection of Ariadne's `TaylorModel` arithmetic shows the key multiplication error identity used after multiplying retained polynomial parts:

```
re += xe*ye + xs*ye + ys*xe
```

where `xe,ye` are uniform input errors and `xs,ys` are l1 coefficient norms. Sweeping additionally moves discarded coefficients into the uniform error. This is precisely the semantic ingredient missing from the cheap Differential residual.

A first lightweight evaluator has therefore been added. Each Procedure temporary carries:

```
(retained Differential, uniform FloatDPUpperBound remainder)
```

It currently implements the generic arithmetic primitives needed by Van der Pol without using polynomial degree information:
- constants and variables;
- +, -, unary +/- and half;
- multiplication;
- square.

For multiplication, the remainder contains:
1. the l1 magnitude of retained-retained coefficient products whose total degree exceeds the retained Differential degree;
2. `||p1|| e2 + ||p2|| e1 + e1 e2`.

Thus this is not the previous "assume truncated terms are zero" Differential path. It explicitly encloses the discarded multiplication tail on the normalised unit box.

Unsupported Procedure operations currently make the lightweight path unavailable rather than silently losing rigour. The next stages will add reciprocal/division and analytic unary operations using the same validated Taylor-series truncation logic already present in `TaylorModel::compose(AnalyticFunction,...)`.

The experiment reports `lightweight_defect_range` and cumulative `lightweight_defect_seconds`. Production certification is unchanged. The important first gate on Van der Pol is whether this cheap enclosure contains/agrees with the full general-TaylorModel oracle and the degree-complete polynomial oracle while costing materially less.

---


### 9.41 Naive scalar-tail Procedure algebra fails: dependency blow-up and excessive cost (2026-09-24)

The first lightweight `(Differential, scalar remainder)` prototype is not promising and is removed.

On the identical first step it produces:

```
general TaylorModel defect:
  x [-1.34e-7, 1.80e-7]
  y [-9.9e-7, 1.40e-6]

naive scalar-tail defect:
  x [-1e-9, 1.76e-7]
  y [-7.1026861e7, 7.1026861e7]
```

The second component explodes by roughly fourteen orders of magnitude. Later candidate calls show scalar remainders from 1e6 to 1e10. This is not a small implementation-tuning issue: independently collapsing the tail of every Procedure temporary to a scalar interval destroys correlations, then multiplication repeatedly feeds those independent errors back through the expression DAG.

Performance is also poor. At 700 calls the lightweight experiment costs 12.63 s cumulatively, versus:
- retained incremental Procedure update: 0.84 s;
- recurrence-field flow_function: 2.46 s;
- cheap truncated Differential defect: 0.59 s.

The 1e-8 benchmark consequently rises to 22.15 s while retaining the same 215 accepted sets. The experiment is therefore removed from the active code, not merely disabled.

Architectural conclusion: reproducing TaylorModel error formulas with a scalar remainder per Procedure temporary is too coarse. Ariadne's full TaylorModel path remains tight because the polynomial part and its error are propagated together and truncation/sweeping occurs at controlled representation boundaries, rather than replacing every discarded dependency by an independent scalar error at every DAG node.

The next investigation should exploit the already-existing incremental `Graded<Differential>` state more directly. In particular, instead of propagating a scalar tail through the whole Procedure DAG, measure whether only the **newly generated final grade(s)** omitted from the retained field are sufficient to bound the residual tail. For analytic operations the graded recurrence already computes coefficients degree-by-degree. A local tail estimate based on the first omitted grade and a validated convergence/majorant bound may preserve dependency information in retained grades and add a scalar enclosure only once at the final output.

This is closer to Flow*-style high-order remainder estimation than a scalar interval algebra threaded through every elementary operation.

---


### 9.42 Project success criteria: three simultaneous objectives (2026-09-24)

The investigation must be judged against **three distinct objectives**, not only local runtime or step count:

1. **Remove or substantially lower the approximation-error plateau as the integration step decreases.**
   The key symptom that motivated this branch is that Ariadne's current Taylor/graded integrators stop gaining accuracy below a residual over-approximation floor as the step size is reduced, whereas Flow* appears to continue scaling better. Any new integrator architecture must demonstrate that reducing the step continues to reduce the over-approximation error over a meaningfully wider regime.

2. **Permit larger validated integration steps.**
   The new method should obtain a useful validated flow enclosure at larger `h` than the current GradedTaylorSeriesIntegrator for comparable approximation quality. The present Gronwall-based prototype already shows evidence in this direction through substantially fewer accepted sets at tight tolerances.

3. **Improve total execution time at equal over-approximation error and equal simulated time horizon.**
   Raw cost per step is not the right final metric. The comparison must fix both the achieved over-approximation error and the simulated time interval (same initial time and same final time) and then compare total runtime. A method that reaches the same error only by simulating a shorter time interval has not achieved a fair performance improvement. A method that takes fewer/larger steps but has more expensive individual steps is successful only if the end-to-end runtime is lower for the same simulated horizon and the same achieved over-approximation error.

These objectives are coupled but must be measured separately. In particular:
- fewer steps alone does not prove better accuracy scaling;
- a lower local residual alone does not prove better end-to-end runtime;
- a faster step alone is irrelevant if the method reaches the same approximation-error plateau;
- wall-clock comparisons at different achieved errors are not fair;
- wall-clock comparisons over different simulated time horizons are not fair.

The benchmark protocol should therefore retain, for each method and tolerance/step regime:
- identical simulated horizon (same initial and final time) across methods;
- accepted/rejected step count and effective step sizes;
- final or maximum over-approximation error using the same metric;
- total runtime over that identical horizon;
- runtime versus achieved error curve;
- error versus step-size curve.

A successful Flow*-like integrator should improve all three axes together: a lower error floor, larger viable steps, and a better runtime/error Pareto curve.

This triple objective is now the governing criterion for future experiments. Optimisations that improve only a local subphase but do not plausibly contribute to at least one of these three goals should be deprioritised.

---


### 9.43 Equal-horizon runtime-vs-achieved-error benchmark (2026-09-24)

A dedicated benchmark now compares the current Gronwall prototype and GradedTaylorSeriesIntegrator over the **same simulated interval [0,5]** and the same maximum configured step 0.04.

The sweep uses step-error tolerances:

```
1e-5, 1e-6, 1e-7, 1e-8, 1e-9
```

For every run it reports:
- wall-clock runtime;
- achieved final enclosure error, defined as the maximum Taylor-model error among all state components of all final enclosures;
- reach/intermediate/final set counts.

The output marker is `[IntegratorAccuracyBenchmark]`.

This is intended to build the first fair runtime-versus-achieved-overapproximation-error curve at fixed time horizon. The old one-step residual probe is removed from this benchmark program so it does not contaminate timing.

The final-state Taylor-model error is a representation-level over-approximation metric, not yet a complete geometric distance to the exact reachable set. It is nevertheless common to both integrators and directly measures the accumulated model remainder that motivated this investigation. If later a stronger common geometric metric is introduced, this benchmark should retain both metrics.

---


### 9.44 First fixed-horizon runtime-vs-achieved-error result (2026-09-24)

The [0,5] benchmark with max_step=0.04 gives:

| tolerance | Gronwall time | Gronwall achieved error | Gronwall sets | Graded time | Graded achieved error | Graded sets |
|---|---:|---:|---:|---:|---:|---:|
| 1e-5 | 4.047 s | 5.3598e-4 | 125 | 2.619 s | 5.6439e-3 | 133 |
| 1e-6 | 4.766 s | 1.9037e-4 | 135 | 3.469 s | 1.3084e-3 | 169 |
| 1e-7 | 7.205 s | 3.7193e-5 | 165 | 5.859 s | 2.0525e-4 | 268 |
| 1e-8 | 11.330 s | 8.2638e-6 | 215 | 8.847 s | 4.6583e-5 | 401 |
| 1e-9 | 20.542 s | 1.8731e-6 | 327 | 15.357 s | 8.8297e-6 | 675 |

Important conclusions:

1. At equal nominal tolerance, Gronwall is slower, but it also delivers about 4.7x--10.5x lower final Taylor-model error. Equal-tolerance wall-clock comparisons therefore substantially understate its efficiency at equal achieved error.

2. The curves overlap enough for a first equal-error comparison. Gronwall at tolerance 1e-6 gives error 1.90e-4 in 4.77 s, slightly better error than Graded at 1e-7 (2.05e-4) in 5.86 s. Thus at approximately 2e-4 achieved error, Gronwall is already about 19% faster while using 135 vs 268 sets.

3. Gronwall at tolerance 1e-8 gives error 8.26e-6 in 11.33 s, slightly better than Graded at 1e-9 (8.83e-6) in 15.36 s. Thus around 9e-6 achieved error, Gronwall is about 26% faster while using 215 vs 675 sets.

4. The achieved error continues to decrease strongly across the tested Gronwall sweep:
5.36e-4 -> 1.90e-4 -> 3.72e-5 -> 8.26e-6 -> 1.87e-6.
There is no plateau visible down to the current 1e-9 tolerance in this final-Taylor-error metric. The Graded curve also continues decreasing over this range, so this benchmark alone does not yet reproduce the original plateau symptom; a dedicated fixed-step/error-vs-h experiment remains necessary for objective 1.

5. Larger-step capability is strongly supported by set counts. The advantage grows with accuracy: 125 vs 133 sets at 1e-5, reaching 327 vs 675 at 1e-9. Since both runs cover exactly [0,5], these counts imply materially larger average accepted steps for Gronwall.

This is the first evidence that the current prototype meets objective 3 at two directly overlapping achieved-error levels, despite being slower at equal nominal tolerance. The next benchmark should target objective 1 explicitly: disable tolerance-driven ambiguity as much as possible and sweep maximum/fixed step size while measuring achieved over-approximation, so that an error plateau as h decreases can be observed directly.

---


### 9.45 Fixed-horizon error-vs-step benchmark for plateau detection (2026-09-24)

The next benchmark targets objective 1 directly. Both integrators cover the identical interval [0,5] and use the same order-5 configuration and sweeper. Instead of sweeping the nominal error tolerance, the local `StepMaximumError` is deliberately loosened to 1e-2 and the configured maximum step is swept through:

```
0.04, 0.02, 0.01, 0.005, 0.0025
```

The purpose is to make `max_step`, rather than the local error threshold, control the effective step whenever the validated flow bound permits it.

Each run reports the existing common final Taylor-model error plus set counts under the marker:

```
[IntegratorStepScalingBenchmark]
```

Interpretation rule:
- first verify from the set count that the actual average step tracks the requested maximum step (for horizon 5, ideal counts are approximately 125, 250, 500, 1000, 2000);
- only then use the resulting `achieved_final_error` as an error-vs-h curve;
- if either integrator takes materially more sets than 5/max_step, its curve at that point is not a clean fixed-step point and must be labelled as adaptivity/bounder-limited.

The key test is whether the Graded error curve flattens as max_step decreases while the Gronwall/preconditioned curve continues decreasing. Runtime is still recorded, but this sweep is primarily a plateau/scaling experiment, not the objective-3 benchmark.

---


### 9.46 Fixed-step scaling result: strong early gain, then a new carried-state floor (2026-09-24)

The max-step sweep is clean: set counts are essentially exactly 5/h for both methods, so the requested max step controls the integration rather than the loose StepMaximumError threshold.

| max h | Gronwall time | Gronwall final error | sets | Graded time | Graded final error | sets |
|---|---:|---:|---:|---:|---:|---:|
| 0.04 | 3.968 s | 5.3598e-4 | 125 | 2.397 s | 1.3749e-2 | 125 |
| 0.02 | 6.467 s | 1.4429e-5 | 251 | 3.473 s | 8.6237e-4 | 251 |
| 0.01 | 10.622 s | 1.8470e-6 | 500 | 5.241 s | 5.6970e-5 | 500 |
| 0.005 | 18.081 s | 2.4749e-6 | 1000 | 8.027 s | 7.6255e-6 | 1000 |
| 0.0025 | 31.193 s | 3.8885e-6 | 2000 | 12.966 s | 6.9174e-6 | 2000 |

This is a crucial result.

1. Gronwall/preconditioned improves dramatically from h=0.04 to h=0.01: roughly 290x lower final model error for a 4x smaller step. Over the same range Graded improves about 241x. The new method is also much more accurate at identical h: about 26x at 0.04, 60x at 0.02, and 31x at 0.01.

2. However, the Gronwall curve reaches its minimum around h=0.01 and then **gets worse** as h is reduced further:
1.847e-6 -> 2.475e-6 -> 3.888e-6.
Therefore objective 1 is not solved end-to-end yet. The new local flow certificate is much tighter, but another error source accumulates with the number of steps and dominates for small h.

3. Graded shows the expected flattening: 5.697e-5 -> 7.626e-6 -> 6.917e-6. It still improves at h=0.005 but essentially plateaus by h=0.0025. The Gronwall floor is lower in this experiment, but its upward turn is especially diagnostic.

4. Because both methods take exactly the intended number of steps, the worsening Gronwall error cannot be blamed on adaptive step selection. The likely culprit is per-step carried-state representation/composition/sweeping error. Halving h doubles the number of endpoint/flowpipe compositions and preconditionings. This matches the earlier observation that the flow certificate itself becomes extremely small at small h while carried-state operations remain nonzero.

5. This redirects the next investigation: do not spend the next effort tightening the local residual. At h <= 0.01 the local residual is no longer the visible limiting factor. Instrument the final state-function error immediately before and after endpoint composition, preconditioning, state composition and sweeping, and determine which carried-state operation injects the approximately per-step floor.

This result also explains why a Flow*-like integrator cannot be obtained solely by improving the one-step flow polynomial/certificate. The representation of the propagated set between steps must preserve the tighter local accuracy.

---


### 9.47 Instrument carried-state error injection below the h=0.01 minimum (2026-09-24)

The next diagnostic targets the newly exposed small-step floor rather than the local residual.

For each accepted preconditioned step it samples the Taylor-model error vectors at these representation boundaries:

```
input_normalised_errors
physical_local_flow_errors
local_endpoint_errors
evolved_physical_errors
local_transition_errors
next_normalised_errors
```

under the marker `[CarriedStateErrorProfile]`. The first five steps and then every 100th step are reported.

Interpretation:
- `physical_local_flow_errors -> local_endpoint_errors` isolates endpoint evaluation of the fresh one-step flow;
- `local_endpoint_errors -> local_transition_errors` shows the effect of preconditioning the fresh endpoint;
- `input_normalised_errors + local_transition_errors -> next_normalised_errors` exposes the accumulation introduced by composing the new local coordinate transition with the carried symbolic state;
- `evolved_physical_errors` gives a physical-coordinate reference after direct endpoint composition with the carried state.

The benchmark is reduced to the Gronwall/preconditioned integrator at max_step 0.01, 0.005 and 0.0025, exactly the regime where the final error changes from 1.85e-6 to 2.47e-6 to 3.89e-6. This avoids rerunning the already-understood Graded comparison and keeps the diagnostic focused.

No production arithmetic or certification decision is changed by this instrumentation.

---


### 9.48 Carried-state probe localises the small-step floor to accumulated state composition (2026-09-24)

The carried-state error probe gives a clear separation of scales.

At h=0.01, final step 500:
- fresh physical local-flow error is about 1.0e-10;
- fresh local-transition error is about 1.5e-10;
- carried normalised-state error is about 2.1e-6 in the dominant component;
- evolved physical error is 1.847e-6.

At h=0.005, final step 1000 of that run (global profile step 1500):
- fresh local-flow/local-transition errors remain about 0.7e-10 / 0.95e-10;
- carried normalised-state error is about 2.77e-6;
- evolved physical error is 2.475e-6.

At h=0.0025, final step 2000 of that run (global profile step 3500):
- fresh local-flow error is only about 4.7e-11;
- fresh local-transition error is about 6.4e-11;
- carried normalised-state error has grown to about 4.35e-6;
- evolved physical error is 3.888e-6.

Endpoint evaluation changes the fresh local-flow error only at approximately machine-level relative scale. Preconditioning increases the fresh one-step error modestly (typically a factor around 1.3--1.5), but it remains O(1e-10), four to five orders of magnitude below the carried error. Therefore neither the Gronwall local certificate nor fresh endpoint/preconditioning error explains the observed O(1e-6) floor.

The dominant error is already present in `state.normalised_mapping()` and is propagated/augmented by:

```
compose(local_transition.normalised_mapping(),
        state.normalised_mapping())
```

This confirms the carried-state composition as the next target. The data do **not** yet distinguish how much of the new error at each composition is unavoidable propagation/amplification of the previous uniform remainder versus newly generated polynomial truncation/sweeping error.

A sparse follow-up diagnostic is added every 100 steps. It repeats the same composition after clobbering only the incoming uniform errors of both operands, preserving their polynomial parts. The resulting `clean_composition_errors` measure the error newly generated by polynomial composition/sweeping/roundoff in that operation, while the normal `actual_composition_errors` include propagation of the incoming remainders.

Marker:

```
[CarriedStateCompositionProbe]
transition_input_errors=...
state_input_errors=...
actual_composition_errors=...
clean_composition_errors=...
```

This diagnostic is not used for certification or returned state; it is an attribution experiment only.

---


### 9.49 Composition probe: new composition error is tiny; accumulated incoming remainder dominates (2026-09-24)

The clean-composition experiment answers the attribution question.

At h=0.01, step 500:
- state input error: 2.055e-6 / 1.959e-7;
- actual composed error: 2.071e-6 / 1.948e-7;
- clean composition error (same polynomial parts, incoming uniform errors clobbered): 1.30e-10 / 1.27e-10.

At h=0.005, the same pattern persists. For example at the end of the run (global step 1500):
- state input error: 2.762e-6 / 2.627e-7;
- actual composition error: 2.773e-6 / 2.620e-7;
- clean composition error: 1.12e-10 / 1.23e-10.

At h=0.0025, final global step 3500:
- state input error: 4.345e-6 / 4.137e-7;
- actual composition error: 4.354e-6 / 4.131e-7;
- clean composition error: 8.92e-11 / 1.13e-10.

Thus the repeated composition is **not creating O(1e-6) fresh truncation error at each step**. Its newly generated error, when incoming uniform remainders are removed, stays around O(1e-10). The dominant O(1e-6) quantity is the already-carried uniform remainder.

The small-step plateau/upturn is therefore mainly the result of repeatedly propagating a finite-order carried-state remainder through many steps, not a catastrophic loss of precision in one composition operation. This is consistent with the expected existence of a spatial-order floor: at fixed spatial degree, the per-step map becomes more accurate with h, but the finite-order state representation still has a non-zero truncation/remainder scale that is transported over T/h steps.

This changes the interpretation of objective 1. We do **not** require the error to tend to zero for fixed spatial order. We require the plateau to occur at the natural finite-order floor and to move downward when the spatial order is increased, rather than being dominated by an avoidable implementation artefact.

A direct spatial-order test is therefore the next diagnostic. At fixed horizon [0,5], fixed max_step=0.0025, fixed temporal order 5 and fixed sweeper, run the Gronwall integrator with spatial orders 4, 5 and 6. If the final error decreases materially with spatial order, the observed small-h floor is primarily the expected representation-order floor. If it barely moves, a carried-state representation issue remains.

Temporary carried-state attribution logging is removed before this benchmark so its timing is not polluted.

---



### 9.50 Spatial-order test rejects the simple "natural order-5 floor" explanation (2026-09-24)

At fixed horizon [0,5], max_step=0.0025 and temporal order 5, varying only spatial order gave:

| spatial order | runtime | achieved final error | sets |
|---:|---:|---:|---:|
| 4 | 26.992 s | 3.58278e-6 | 2000 |
| 5 | 32.110 s | 3.88849e-6 | 2000 |
| 6 | 36.652 s | 3.89376e-6 | 2000 |

The floor did not decrease with spatial order. Orders 5 and 6 were essentially identical, while order 4 was slightly better. This rejected the simple hypothesis that the observed floor was set directly by the configured spatial degree.

Combined with the composition probe, this pointed to representation simplification rather than the degree ceiling: small coefficients were likely being swept into the uniform error before the extra degree could matter.

---

### 9.51 Sweeper threshold is the dominant small-step accuracy floor (2026-09-24)

At fixed horizon [0,5], max_step=0.0025, spatial order 5, temporal order 5 and loose local tolerance 1e-2, varying only the ThresholdSweeper cutoff gives:

| sweep threshold | runtime | achieved final error | sets |
|---:|---:|---:|---:|
| 1e-10 | 21.215 s | 2.06453e-4 | 2000 |
| 1e-12 | 32.214 s | 3.88849e-6 | 2000 |
| 1e-14 | 56.618 s | 6.67444e-8 | 2000 |
| 1e-16 | 103.174 s | 1.35309e-9 | 2000 |

This is decisive: the apparent small-h plateau is controlled overwhelmingly by the sweeper threshold, not by the configured spatial degree. Tightening the threshold by two decades lowers final error by approximately 53x, then 58x, then 49x. Across 1e-10 -> 1e-16 the final error improves by more than five orders of magnitude while the number of steps is unchanged.

Therefore objective 1 must be interpreted relative to both representation degree **and representation simplification threshold**. The previous q=4/5/6 test did not move the floor because all three representations were being simplified at the same 1e-12 cutoff before the extra degree could matter.

The cost tradeoff is equally strong: retaining these small coefficients makes the representation much more expensive. Runtime grows 21.2 -> 32.2 -> 56.6 -> 103.2 s. Thus simply setting the threshold near machine precision is not a production solution, even though it demonstrates that the new flow/carry architecture can reach much lower errors.

Architectural implication: the remaining accuracy/performance problem is now sharply identified as **representation management**. We need to preserve the coefficients that matter to long-horizon dependency without paying the full cost of a globally tiny threshold. Candidate directions include graded/adaptive sweeping, scale-aware sweeping, and retaining terms based on their propagated impact rather than instantaneous coefficient magnitude.

This also strengthens the interpretation of the earlier equal-error benchmark: the current 1e-12 configuration is not a hard mathematical accuracy limit. Its floor is a tunable simplification/cost tradeoff.

---


### 9.52 Preserve low-degree small coefficients: sweeper-policy benchmark (2026-09-24)

The threshold sweep established that coefficient sweeping controls the observed accuracy floor, but globally lowering the threshold is too expensive. The next test asks whether a **degree-based policy** can preserve the symbolic information that matters without retaining every small high-degree term.

At fixed horizon [0,5], max_step=0.0025, spatial order 5, temporal order 5 and loose local tolerance 1e-2, compare:

1. `ThresholdSweeper(1e-12)`: current baseline;
2. `ThresholdSweeper(1e-14)`: high-accuracy/high-cost reference;
3. `GradedSweeper(degree=5)`: keep every coefficient of degree <=5, regardless of magnitude, discard only terms above degree 5;
4. `GradedThresholdSweeper(degree=5, threshold=1e-14)`: explicit degree ceiling plus tight magnitude threshold.

The most informative comparison is (2) versus (3). If `GradedSweeper(5)` approaches the 1e-14 error with runtime much closer to the 1e-12 baseline, then the problem is specifically that absolute threshold sweeping destroys small but structurally useful low-degree coefficients. That would give us a strong production direction without inventing a new sweeper yet.

If `GradedSweeper(5)` is both expensive and not substantially more accurate, then coefficient magnitude alone is not enough to decide what to retain and a more selective propagated-impact/adaptive policy is needed.

Output marker: `[IntegratorSweeperPolicyBenchmark]`.

---


### 9.54 RelativeThresholdSweeper(1e-12) is computationally explosive (partial run, 2026-09-24)

The relative-threshold benchmark was interrupted during the third run, `relative_1e-12`, after about 1900 of the intended 2000 steps. The two absolute baselines completed normally:

| policy | runtime | achieved final error | sets |
|---|---:|---:|---:|
| absolute 1e-12 | 33.079 s | 3.88849e-6 | 2000 |
| absolute 1e-14 | 59.523 s | 6.67444e-8 | 2000 |

The relative 1e-12 run did not finish. By global profile step 5900 (approximately 1900 steps into that run), cumulative carried-state costs imply enormous incremental cost compared with the preceding absolute runs. In particular, cumulative flowpipe composition grew from about 30.0 s at the end of the absolute 1e-14 run to 173.4 s, state composition from about 12.1 s to 52.6 s, and endpoint composition from about 11.0 s to 44.3 s.

This means the relative policy is retaining far more structure than the absolute 1e-12 baseline in this problem. Since its cutoff is
`relative_threshold * (radius(polynomial) + uniform_error)`, the current model scale can make the effective threshold substantially smaller than 1e-12. Repeated composition then grows the retained expansion and causes a severe cost explosion.

Conclusion: `RelativeThresholdSweeper(1e-12)` does not provide the desired accuracy/cost compromise on this benchmark. The concept of scale-aware sweeping remains relevant, but the raw existing relative policy is too permissive at this parameter value. The next useful test should first measure the effective norm/cutoff or use substantially larger relative thresholds (for example 1e-10 and 1e-8) rather than attempting `relative_1e-14`, which would almost certainly be even more expensive.

---


### 9.55 Calibrate relative sweeping with coarser thresholds (2026-09-24)

The interrupted `relative_1e-12` run showed that the raw relative threshold was far too permissive on this benchmark: the carried polynomial representation grew enough to make repeated compositions prohibitively expensive.

The next benchmark therefore keeps the two absolute references and replaces the overly aggressive relative values with:

```
relative_1e-8
relative_1e-10
```

at the same horizon [0,5], max_step=0.0025 and spatial/temporal order 5.

The intent is not to match numeric threshold values across absolute and relative sweepers. Their semantics differ:
- absolute: discard if `|c| < tau_abs`;
- relative: discard if `|c| < tau_rel * (radius(polynomial) + uniform_error)`.

Hence the useful comparison is empirical: find whether a relative threshold can produce an error between the absolute 1e-12 and 1e-14 references at runtime close to or below the absolute 1e-12 baseline.

If relative 1e-8 is already too inaccurate, while relative 1e-10 is still too expensive, the existing relative policy is too sensitive to model scale for this workload and a more controlled adaptive/budget policy will be needed.

---


### 9.56 Existing RelativeThresholdSweeper does not improve the accuracy/cost frontier (2026-09-24)

The calibrated relative-threshold benchmark completed:

| policy | runtime | achieved final error | sets |
|---|---:|---:|---:|
| absolute 1e-12 | 32.172 s | 3.88849e-6 | 2000 |
| absolute 1e-14 | 57.576 s | 6.67444e-8 | 2000 |
| relative 1e-8 | 36.408 s | 4.32786e-4 | 2000 |
| relative 1e-10 | 96.740 s | 7.36170e-6 | 2000 |

Neither relative point is competitive with the absolute baseline:
- relative 1e-8 is about 1.13x slower than absolute 1e-12 while roughly 111x less accurate;
- relative 1e-10 is about 3.0x slower than absolute 1e-12 while roughly 1.9x less accurate;
- absolute 1e-14 is both faster and about 110x more accurate than relative 1e-10.

Therefore the existing `RelativeThresholdSweeper` is not the selective sweeping policy needed for this carried-state workload. Its single model-wide scale `radius(polynomial)+uniform_error` does not discriminate which coefficients preserve important dependency. Tuning the scalar relative threshold merely moves between aggressive information loss and representation explosion, without improving the Pareto frontier established by absolute threshold sweeping.

This closes the immediate sweeper-reuse branch. The next step should not be more scalar-threshold tuning. Instrument the carried expansion to determine which retained terms account for the accuracy gain from absolute 1e-12 to 1e-14: coefficient counts by degree and magnitude band, ideally sampled at representative late steps. That evidence can support a genuinely selective policy (degree-dependent threshold, retained-term budget, or propagated-impact criterion) rather than another global scalar cutoff.

---

## 10. Direction of the project

The preconditioned direction remains worth investigating, but the target is now precise:

> Preserve the proven geometric/error benefit of QR while eliminating or reducing the penalty it creates inside Ariadne's validated graded-series calculation.

Do not currently frame the goal as "make Ariadne behave like Flow*". The useful Ariadne-specific design is:

```
PreconditionedGradedTaylorSeriesIntegrator
    = GradedTaylorSeries machinery
    + persistent affine local coordinates (c,A,y)
    + optional QR preconditioning
```

The current bottleneck is the interaction between QR coordinates and validated interval evaluation in the graded-series step.

---

## 11. Rules for future investigation

1. **Always compare from the same physical state** when diagnosing coordinate effects.
2. **Separate local step error from accumulated global enclosure error.**
3. **Compare efficiency at equal final accuracy**, not merely equal `StepMaximumError`.
4. **Do not infer Flow* equivalence from similarly named quantities.**
5. **Keep diagnostic experiments short** (one or two steps) whenever possible.
6. **Remove or isolate temporary diagnostics before final cleanup.**
7. **Do not optimise for merely passing tests:** preserve mathematically faithful validated behaviour.
8. After each meaningful experiment, record:
   - commit,
   - exact configuration,
   - measured output,
   - conclusion,
   - whether the hypothesis was confirmed/rejected,
   - next question.



### 9.57 Profile the coefficient band that buys accuracy (2026-09-24)

The intermediate absolute-threshold control benchmark showed that the tested degree-selective policy does not beat the Pareto frontier of plain absolute sweeping. In particular, absolute `3e-14` was slightly faster and more accurate than the selective `tight d<10` policy at essentially the same runtime.

The next diagnostic therefore stops tuning degree cutoffs and asks a more basic attribution question:

> Which coefficients are actually preserved by `3e-14` and lost by `1e-12`, and at which spatial degrees do they occur?

The Van der Pol benchmark now runs only the two decisive absolute policies, `1e-12` and `3e-14`. A diagnostic sweeper classifies every coefficient seen by the sweeper into three magnitude bands:

```
|c| < 3e-14
3e-14 <= |c| < 1e-12
|c| >= 1e-12
```

For each total spatial degree it records coefficient count and summed absolute coefficient mass. The middle band is the important one: these terms are discarded by the `1e-12` policy but retained by `3e-14`.

The same classification is collected separately on the trajectory generated by each active threshold, since the two policies can subsequently produce different carried Taylor expansions.

Output markers:

```
[IntegratorSweepBandBenchmark]
[SweepBandProfile]
```

This diagnostic does not alter the certification rule beyond selecting the active absolute threshold itself.

Interpretation:
- if the bridge-band population is concentrated in a small set of degrees, a more targeted degree/magnitude policy may still be justified;
- if it is broad across degrees but concentrated in coefficient mass, a retained-term budget or magnitude-rank policy is more plausible;
- if it is broad both in degree and count, the gain from the tighter threshold likely comes from a long tail of dependencies, and a new sweeper heuristic alone may not reduce cost substantially.

Do not design another production sweeper before inspecting this profile.


### 9.58 Snapshot the actually carried expansion (2026-09-24)

The cumulative sweep-band profile showed that the coefficients between `3e-14` and `1e-12` are not confined to low degree. The tighter run encounters a broad bridge population through approximately degrees 0--16, with especially large counts around degrees 4--8. This explains why the previous degree-selective policy could not beat the plain absolute `3e-14` point.

However, those counts are cumulative sweeper events, not unique coefficients resident in the carried state. They may count the same persistent structure many times during repeated arithmetic and composition.

The next diagnostic therefore removes the profiled sweeper and restores ordinary `ThresholdSweeper` instances for the two decisive policies:

```
absolute 1e-12
absolute 3e-14
```

At carried-state steps 500, 1000, 1500 and 2000, the evolver directly inspects `state.normalised_mapping()` immediately before the local step and `local_step.final_state().normalised_mapping()` immediately after it. For every total degree it records the resident coefficient count and absolute mass in the same three bands:

```
|c| < 3e-14
3e-14 <= |c| < 1e-12
|c| >= 1e-12
```

Output marker:

```
[CarriedExpansionSnapshot]
```

The benchmark summary marker is:

```
[IntegratorCarriedExpansionBenchmark]
```

Interpretation:
- a small, stable resident bridge population together with millions of cumulative sweep events would point to repeated processing of persistent structure as the main cost opportunity;
- substantial resident growth under `3e-14` would instead show that the tighter policy genuinely carries a much larger state representation;
- comparing input/output snapshots at the same milestone reveals how much one local transition changes the resident structure.

No certification rule is changed by this diagnostic.


### 9.59 Separate carried-expansion snapshots from verbose integrator diagnostics (2026-09-24)

The first carried-expansion snapshot run was interrupted before step 500 because it enabled the existing global `diagnostics` flag. That flag activates the full historical diagnostic stack, including per-instruction `[GradedProcedureDiagnostic]` output, producing tens of thousands of lines long before the first requested snapshot.

This was a diagnostic-control mistake, not evidence that the carried-expansion inspection itself is expensive.

A dedicated `carried_expansion_diagnostics` flag is now added to `PreconditionedGradedTaylorSeriesIntegrator`. The Van der Pol probe uses:

```
diagnostics = false
carried_expansion_diagnostics = true
```

The evolver therefore emits only the requested `[CarriedExpansionSnapshot]` records at steps 500, 1000, 1500 and 2000, while all pre-existing verbose diagnostics remain disabled. The numerical algorithm and certification rules are unchanged.


### 9.60 Profile temporary Taylor-product generation (2026-09-24)

The carried-expansion snapshots rule out resident representation size as the main explanation for the cutoff/runtime trade-off. The `3e-14` run carries only tens of additional resident coefficients at the sampled late steps, while the earlier cumulative sweeper experiment observed millions of sweep events.

The next diagnostic therefore instruments the Taylor-model product kernel `_ifma`, where multiplication is performed monomial-by-monomial and the intermediate result is swept after every source monomial.

A diagnostic-only global profile is enabled only around the two Van der Pol benchmark runs. It records:

- number of Taylor-model product calls;
- exact number of coefficient product pairs `|x|*|y|` processed by those calls;
- number of intermediate sweep passes;
- total number of materialised terms immediately before and after those sweeps;
- total number of terms removed by the sweeps;
- maximum intermediate sweep input/output size;
- wall-clock time spent inside the profiled product kernel.

Products executed from `TaylorModel::_compose(x,y)` are tagged separately from all other Taylor-model products, yielding two contexts:

```
general
compose
```

Output marker:

```
[TaylorProductGenerationProfile]
```

This experiment does not change sweeping or certification. Its purpose is to determine whether the tighter cutoff primarily increases (a) the number of coefficient products generated, (b) the size of intermediate merged expansions, or (c) both, and how much of the additional product-kernel time is specifically attributable to composition.


### 9.61 Classify individual Taylor coefficient products before merging (2026-09-24)

The product-generation profile localized most of the `3e-14` runtime increase inside the general Taylor-model product kernel rather than the carried-state composition context. The next diagnostic therefore classifies every coefficient product `x_i*y_j` at the point where it is generated inside `_ifma`, before it is merged with an existing coefficient carrying the same multi-index.

For the `ValidatedTag,FloatDP` path used by this benchmark, the profiler records:

- number of individual coefficient products with `|x_i*y_j| < sweep_threshold`;
- number with `|x_i*y_j| >= sweep_threshold`;
- summed absolute mass of each class.

These counters are reported separately for the existing `general` and `compose` contexts through the same marker:

```
[TaylorProductGenerationProfile]
```

The classification uses ordinary double values only for diagnostics and does not affect the Taylor arithmetic, sweep decisions, error term, or certification.

The decisive question is whether a large majority of generated products are already individually below the active cutoff. If so, a rigorous early-discard product kernel becomes worth prototyping; if not, the main opportunity lies elsewhere in merge/sweep implementation rather than product generation.


### 9.62 Fix individual-product diagnostic double counting (2026-09-24)

The first individual-product classification run exposed a diagnostic bug: the reported
`individual_below_threshold + individual_above_threshold` count exceeded
`product_pairs`.

The cause was that the classification occurred at the top of the merge loop, before
knowing whether the current `yiter` product was actually consumed. In the
`ra < ta` branch only `riter` advances, so the same `x_i*y_j` pair could be
classified repeatedly on successive merge iterations.

The diagnostic now classifies `x_i*y_j` only in branches that consume `yiter`:

- `ra == ta`, where the product is fused into an existing coefficient;
- `ta < ra`, where the product creates a new coefficient;
- the trailing `while(yiter!=y.end())`, once per remaining product.

The Van der Pol reporting path also asserts the invariant

```
individual_below_threshold + individual_above_threshold == product_pairs
```

for both `general` and `compose` contexts. The arithmetic and certification remain
unchanged; this commit only repairs the profiler.


### 9.63 Split below-threshold products by merge role (2026-09-24)

The corrected individual-product profile shows that roughly 83--85% of products in the
general Taylor-model product context are individually below the active cutoff, and that
their count is close to the number of terms later removed by intermediate sweeps.

Before implementing early discard, the profiler now partitions below-threshold products
according to the merge role in which they are consumed:

- `collision`: `ra == ta`, where the product is fused into an already resident coefficient;
- `new_term`: `ta < ra`, where the product creates a new coefficient before the next resident term;
- `trailing`: the final `while(yiter!=y.end())`, where every remaining product necessarily creates a new coefficient.

New output fields on `[TaylorProductGenerationProfile]` are:

```
individual_below_collision
individual_below_new_term
individual_below_trailing
```

The reporting code asserts both diagnostic invariants:

```
below + above == product_pairs
below_collision + below_new_term + below_trailing == below
```

The sum `below_new_term + below_trailing` is the conservative candidate pool for an
early-discard implementation. No arithmetic or certification rule is changed by this diagnostic.


### 9.64 Prototype conservative early discard in Taylor products (2026-09-24)

The merge-role profile shows that, for the `3e-14` general Taylor-product workload,
the overwhelming majority of individually below-threshold products occur in the
`new_term` and `trailing` branches rather than in collisions with existing coefficients.
These are precisely the branches where a coefficient is materialised only to be removed
by the intermediate sweep.

An experimental global switch now enables a conservative early-discard path in `_ifma`
for `ValidatedTag,FloatDP` Taylor models using a concrete `ThresholdSweeper<FloatDP>`.

The prototype deliberately keeps `mul_err(xv,yv,te)`: it therefore preserves the
existing coefficient computation and its roundoff contribution. Only after `tv` has
been computed does it ask the same threshold sweeper whether the newly-created term
would be discarded. If so, it adds `abs(tv)` rigorously to the error term and does not
append the coefficient to the temporary expansion.

The optimisation is applied only in:
- the `ta < ra` new-term branch;
- the trailing `while(yiter!=y.end())` branch.

The `ra == ta` collision branch is intentionally unchanged.

This prototype therefore targets materialisation, subsequent merge traffic, and sweep
traffic without yet attempting to avoid the floating-point multiplication itself.

The Van der Pol benchmark is temporarily reduced to an apples-to-apples comparison at
`3e-14` with all product profiling and carried-expansion snapshots disabled:

```
absolute_3e-14_baseline   early_discard=false
absolute_3e-14_early      early_discard=true
```

Output marker:

```
[IntegratorEarlyDiscardBenchmark]
```

The first acceptance criterion is that the early path remains rigorous and preserves
essentially the same final accuracy while reducing runtime. If successful, a second
stage can investigate a rigorous pre-product cutoff test that avoids `mul_err` itself.


### 9.65 Fix early-discard benchmark cleanup error (2026-09-24)

The first build of the early-discard benchmark failed because the transition from the
previous product-profiling benchmark left three stale lines in
`examples/continuous/vanderpol.cpp`: two calls to `print_product_profile(...)` referring
to the now-removed `product_profile` variable, plus an extra lambda terminator.

This was a benchmark-cleanup mistake only. The stale lines are removed; the intended
two-run comparison remains:

```
absolute_3e-14_baseline   early_discard=false
absolute_3e-14_early      early_discard=true
```

No Taylor arithmetic or early-discard logic is changed by this fix.


### 9.66 Conservative early-discard result (2026-09-24)

The conservative early-discard prototype was benchmarked at absolute cutoff `3e-14`
with all product profiling and carried-expansion snapshots disabled.

Results:

```
baseline:      49.3621 s   final error 1.7496349409095074e-7
early discard: 48.1391 s   final error 1.7496349409095074e-7
```

The observed speedup is about 2.5%, while the reported final error is identical. This
shows that avoiding materialisation of newly-created terms which the same threshold
sweeper would immediately remove is a valid optimisation opportunity, but materialisation
alone is not the dominant cost. The prototype still pays for `mul_err`, multi-index
construction and merge traversal.

Before attempting a more aggressive pre-product cutoff test, the investigation now
revisits the semantics of the existing multiplication kernel. In the current `_ifma`,
`t.sweep()` is executed after each source monomial of `x`. Consequently, partial
contributions to the same final multi-index can be swept into the remainder before later
source monomials contribute to that index. This remains rigorous as an enclosure, but it
is not equivalent to accumulating the complete coefficient for each multi-index and then
applying the cutoff once to the multiplication result.

### 9.67 Compare incremental sweep with final-product sweep (2026-09-24)

A diagnostic switch now allows `_ifma` to use either:

```
incremental_sweep=true   # existing behaviour: sweep after every source monomial
incremental_sweep=false  # accumulate the full product, sweep once at the end
```

The merge arithmetic itself is unchanged. In final-sweep mode the temporary expansion is
carried unswept through all source monomials; after the complete product has been merged,
a single `r.sweep()` is applied. Early discard is disabled for both runs so that the test
isolates only the placement of the cutoff operation.

The Van der Pol comparison uses the same absolute cutoff `3e-14`:

```
absolute_3e-14_incremental
absolute_3e-14_final
```

Output marker:

```
[IntegratorSweepSemanticsBenchmark]
```

This test answers whether the existing incremental cutoff is merely a performance
engineering choice or whether delaying cutoff until complete multi-index aggregation
materially changes the accuracy/runtime frontier.


### 9.68 Final-sweep cutoff frontier (2026-09-24)

The direct semantic comparison at absolute cutoff `3e-14` gave:

```
incremental sweep: 49.1471 s   final error 1.7496349409095074e-7
final sweep:       56.2521 s   final error 8.5378288508794491e-8
```

Thus delaying the cutoff until all contributions to each product multi-index have been
aggregated reduces the final error by about 51% at the same threshold, at a runtime cost
of about 14%.

This confirms that the placement of the sweep is not merely an implementation detail:
the existing incremental sweep loses polynomial information before later contributions to
the same multi-index arrive. The enclosure remains rigorous, but the cutoff is applied to
partial coefficients rather than to the fully accumulated product coefficient.

The next experiment maps the accuracy/runtime frontier of final-sweep multiplication at:

```
1e-12
3e-13
1e-13
3e-14
```

All runs use `incremental_sweep=false` and early discard remains disabled. The goal is
to determine whether final-sweep semantics allow a looser cutoff to match or improve the
accuracy of the current incremental `3e-14` point while reducing runtime.

The decisive comparison is against the established incremental frontier, especially:

```
incremental 1e-12  ~32 s   3.888e-6
incremental 1e-13  ~43 s   5.026e-7
incremental 3e-14  ~49 s   1.750e-7
incremental 1e-14  ~57 s   6.674e-8
```

If a final-sweep point lies below and to the left of one of these points, then the new
sweep semantics improve the Pareto frontier rather than merely trading time for accuracy.

Output marker remains:

```
[IntegratorSweepSemanticsBenchmark]
```


### 9.69 Final-sweep frontier result (2026-09-24)

The four-point final-sweep frontier completed as follows:

```
cutoff   elapsed_s   final_error
1e-12    34.6421     2.1585938662392414e-6
3e-13    40.8661     7.4109181449436055e-7
1e-13    46.6021     2.5932372140202589e-7
3e-14    55.4641     8.5378288508794491e-8
```

Final sweeping consistently improves accuracy at a fixed cutoff, but the naive
implementation does not dominate the established incremental-sweep frontier. The curves
interlace: delaying the sweep preserves more polynomial structure but carrying the
unswept expansion through the legacy repeated merge/swap kernel costs enough runtime to
consume the numerical advantage.

The conclusion is therefore not to return to incremental semantics, but to change the
multiplication algorithm so that complete coefficients are accumulated by multi-index
without repeatedly rebuilding the whole expansion.

### 9.70 Prototype a direct full-product accumulator (2026-09-24)

Inspection of `Expansion<MultiIndex,Coefficient>` shows that Ariadne already provides
the operations needed for a first accumulator prototype:

- reserve sparse storage;
- append index/coefficient pairs;
- sort by multi-index;
- `TaylorModel::unique()`, which combines adjacent equal multi-indices using
  `add_err` and transfers addition roundoff to the model error;
- a final `sweep()`.

The experimental product path therefore constructs one temporary expansion containing:

```
existing r terms
+
all x_alpha * y_beta coefficient products
```

Each coefficient product is still computed with `mul_err`, so multiplication roundoff
is rigorously accumulated. The temporary expansion is then sorted once, equal
multi-indices are combined once, and the configured sweeper is applied once to the fully
aggregated coefficients.

This is deliberately a first structural prototype rather than the final data structure.
It can temporarily hold `|r| + |x||y|` entries before sorting, but avoids the legacy
algorithm's repeated merge/copy/swap of the entire accumulated expansion after every
source monomial. It also implements the intended final-coefficient cutoff semantics.

The benchmark compares legacy final-sweep merging against the accumulator at the two
most informative cutoffs:

```
final_merge_1e-13
accumulator_1e-13
final_merge_3e-14
accumulator_3e-14
```

All runs disable early discard and product profiling. Output marker:

```
[IntegratorProductAccumulatorBenchmark]
```

Acceptance criteria are:
- rigorously completed runs;
- final errors consistent with final-sweep semantics at the same cutoff;
- a material runtime reduction versus the repeated-merge final-sweep kernel.

If this succeeds, the next refinement should replace append-all-plus-sort with a true
multi-index keyed accumulator or a degree/index-addressed workspace, depending on the
observed temporary product sizes and MultiIndex structure.


### 9.71 Direct accumulator result and compression profiling (2026-09-24)

The append-all/sort/unique accumulator prototype completed successfully and preserved
final-sweep accuracy while reducing runtime:

```
cutoff   legacy final-merge   accumulator   final error (accumulator)
1e-13    46.9231 s            43.8901 s     2.5932365467959273e-7
3e-14    ~55.5 s              53.0491 s     8.5378293221320844e-8
```

At `1e-13` the prototype saves about 6.5% versus the repeated-merge final-sweep
implementation, despite deliberately materialising every coefficient product before
sorting. The final error differs only at the level expected from a different floating-point
summation order, while the validated enclosure is preserved.

This is strong evidence that repeated full-expansion merge/copy/swap is a real algorithmic
cost. The next decision is which accumulator structure should replace append-all-plus-sort.

A dedicated accumulator profile now records, for each direct-product call:

- `product_pairs = |x||y|`;
- `temporary_entries = |r| + |x||y|` before sorting;
- `unique_entries` immediately after `sort()+unique()`;
- cumulative and maximum temporary/unique sizes;
- the overall compression ratio `temporary_entries / unique_entries`.

The Van der Pol benchmark is temporarily reduced to one accumulator run at `3e-14`,
where product volume is highest:

```
accumulator_profile_3e-14
```

Output marker:

```
[TaylorProductAccumulatorProfile]
```

Interpretation:
- a large compression ratio favours a keyed or degree-addressed accumulator;
- a ratio near one suggests that contiguous append/sort may already be close to optimal,
  and effort should focus on sorting/allocation rather than hashing/tree structures.


### 9.72 Multi-index enumeration and dense-workspace feasibility (2026-09-24)

The compression profile at `3e-14` reported:

```
calls                 1,256,468
product_pairs         208,224,840
temporary_entries     208,224,840
unique_entries         36,018,650
duplication_ratio           5.78103
max_temporary_entries       7,777
max_unique_entries             325
```

This strongly favours accumulating directly by destination multi-index rather than
materialising duplicate products and sorting them afterward.

Inspection of `MultiIndex` found that Ariadne already provides `operator++()`, which
enumerates successive multi-indices by total degree/composition, but there is no active
public rank/position method: the declarations for `position()` / `number()` are
commented out. Therefore a dense accumulator would need either a small combinatorial
ranking helper or a precomputed index-to-slot table.

Before choosing that representation, the accumulator profiler now also records:

- maximum Taylor-model argument size;
- maximum degree present in `x`;
- maximum degree present in `y`;
- maximum possible product degree `degree(x)+degree(y)`;
- maximum number of dense slots needed to represent all monomials up to that degree,
  computed as `C(argument_size + product_degree, product_degree)`.

These values decide whether a dense combinatorial workspace is genuinely small enough
for the observed workload. If the slot count remains modest, direct rank-to-array
accumulation should avoid hashing, tree nodes, duplicate MultiIndex storage, and sorting.
If it grows too large, the next implementation should instead use a compact sparse
index-to-slot structure.


### 9.73 Dense-workspace feasibility and keyed accumulator prototype (2026-09-24)

The `3e-14` feasibility run reported:

```
max_argument_size   3
max_x_degree       16
max_y_degree       17
max_product_degree 31
max_dense_slots  5456
max_unique_entries 325
```

The complete total-degree monomial space is therefore small, and the actually occupied
set is smaller still. This makes a dense index-to-slot workspace practical.

The profiler parameter named `argument_size` also shadowed Ariadne's namespace-level
`argument_size` attribute generator under `-Wshadow`. It has been renamed
`num_variables`.

A first keyed accumulator is now implemented. Rather than materialising all coefficient
products and sorting them, it uses a collision-free mixed-radix rank for each multi-index:

```
slot = a0 + base*a1 + base^2*a2 + ...
base = maximum_degree + 1
```

Only the slot-to-touched map is dense. Multi-indices and coefficients are stored only
for occupied slots, so the large dense space does not require constructing thousands of
floating-point coefficient objects on every call. Repeated contributions are accumulated
immediately with `add_err`; individual products still use `mul_err`.

At the end of a product, only the occupied indices are sorted into Ariadne's reverse
lexicographic order, an `Expansion` is constructed once, and a single final sweep is
performed. For the observed workload this means sorting at most a few hundred occupied
indices instead of thousands of duplicate coefficient-product entries.

The benchmark now compares at `3e-14`:

```
sort_unique_3e-14
dense_3e-14
```

Both use final-coefficient sweep semantics and disable the earlier incremental
early-discard path. The comparison will show whether direct keyed accumulation recovers
the remaining cost of append/sort/unique.


### 9.74 Fix dense-accumulator MultiIndex construction (2026-09-24)

The first build of the keyed dense accumulator failed because the touched-index vector
attempted to construct a `MultiIndex` directly from a `MultiIndexData const&`.
`MultiIndex` has no such constructor; the available owning constructor takes the
argument size and a pointer to the degree array.

The touched-index insertion now uses:

```
MultiIndex(as, index.begin())
```

via `emplace_back(as,index.begin())`, which creates an owning copy of the current
multi-index. This is a construction/API fix only; the ranking, accumulation, roundoff
accounting and final-sweep semantics are unchanged.


### 9.75 Dense accumulator removes per-index MultiIndex allocation (2026-09-24)

The first keyed dense accumulator produced the first clear Pareto improvement:

```
kernel                 elapsed_s   final_error
sort+unique 3e-14      52.2391     8.5378293221320844e-8
dense keyed 3e-14      45.4911     8.5378297219108396e-8
```

Compared with the established incremental-sweep `3e-14` point (about 49.1 s and
`1.75e-7` final error), the dense final-sweep kernel is both faster and about twice as
accurate.

Inspection of the first dense prototype exposed another avoidable cost: every newly
occupied slot stored an owning `MultiIndex`. Since `MultiIndex` owns a dynamically
allocated degree array, this caused one allocation/copy for every distinct coefficient
encountered across all products.

The mixed-radix slot number already contains the full multi-index and, importantly, its
numeric descending order is exactly Ariadne's reverse lexicographic order because the
last variable is the most significant radix digit. The refined kernel therefore stores
only:

- the dense `slot -> touched-position` map;
- a vector of touched integer slots;
- a vector of coefficients.

Touched slots are sorted as integers in descending order. The corresponding
`MultiIndex` is reconstructed from radix digits only when the final `Expansion` is
emitted. This removes the per-distinct-index `MultiIndex` allocation and replaces
MultiIndex comparisons during sorting with integer comparisons.

The benchmark is reduced to a single `dense_3e-14` run. The previous
`TaylorProductAccumulatorProfile` output is suppressed for the dense branch because
those counters instrument only the append/sort/unique path and would otherwise print
misleading zeros.


### 9.76 Refined dense result and apples-to-apples architecture benchmark (2026-09-24)

Removing per-distinct-index `MultiIndex` allocations and sorting mixed-radix slots
directly improved the dense `3e-14` run from 45.4911 s to 45.0521 s. The final error
remained exactly `8.5378297219108396e-8` in the reported run. The additional gain is
about 1%, so further micro-optimisation of index storage is not currently the highest
priority.

The next benchmark addresses the comparison with Ariadne's original
`GradedTaylorSeriesIntegrator`. A two-way comparison would be misleading because the
new dense product kernel is global Taylor-model machinery and can also benefit the
original graded integrator. The benchmark therefore runs three configurations under the
same Van der Pol setup:

```
graded_legacy_3e-14
graded_dense_3e-14
preconditioned_dense_3e-14
```

All three use:
- absolute Taylor sweep threshold `3e-14`;
- `StepMaximumError(1e-2)`;
- Lipschitz tolerance `0.5`;
- fixed spatial and temporal orders 5;
- evolver maximum step `0.0025`;
- maximum enclosure radius 1;
- maximum spacial error `1e-6`;
- evolver reconditioning disabled;
- evolution time 5.

The legacy graded run uses the original incremental Taylor-product kernel. The graded
dense run uses the same final-coefficient dense accumulator as the preconditioned run.
This separates the benefit of the product-kernel redesign from the benefit or overhead
of the persistent QR-preconditioned architecture.

Each run reports:

```
[IntegratorArchitectureBenchmark]
method
elapsed_seconds
achieved_final_error
final_radius
reach_sets
```

The key comparison for architectural competitiveness is
`graded_dense_3e-14` versus `preconditioned_dense_3e-14`; the legacy point is retained
to quantify how much of the improvement comes purely from the Taylor-product kernel.


### 9.77 Graded-dense equal-accuracy frontier (2026-09-24)

The completed apples-to-apples architecture benchmark gave:

```
method                         elapsed_s   final_error              final_radius  reach_sets
graded_legacy_3e-14            23.2251     5.5511748688431288e-7   0.0403        2000
graded_dense_3e-14             21.0271     3.3536682299918048e-7   0.0403        2000
preconditioned_dense_3e-14     44.7441     8.5378297219108396e-8   0.0403        2000
```

The dense full-coefficient product kernel therefore improves the original
`GradedTaylorSeriesIntegrator` as well as the preconditioned integrator. At the same
`3e-14` cutoff, graded+dense is about 9.5% faster than the legacy graded kernel and
reduces final error by about 40%.

The remaining architectural question is now an equal-accuracy comparison. At `3e-14`,
the preconditioned+dense integrator is about 2.1x slower than graded+dense, but its final
Taylor-model error is about 3.9x smaller. The final enclosure radius is identical in the
reported run, so the observed benefit is primarily in the carried Taylor-model remainder.

The next experiment therefore maps only the `GradedTaylorSeriesIntegrator + dense
accumulator` frontier below `3e-14`:

```
3e-14
1e-14
3e-15
1e-15
```

All other benchmark parameters remain unchanged. The target is the established
preconditioned+dense error:

```
8.5378297219108396e-8
```

If graded+dense reaches this error in less than 44.7441 s, the QR-preconditioned
architecture is still not competitive at equal final Taylor-model error on this workload.
If it requires more time, then the preconditioned architecture has produced a genuine
accuracy/runtime advantage rather than only an accuracy improvement at fixed cutoff.

Output marker remains:

```
[IntegratorArchitectureBenchmark]
```


### 9.78 Benchmark warning cleanup (2026-09-24)

After switching the Van der Pol driver from the three-way architecture comparison to the
graded+dense cutoff frontier, the local `run_preconditioned_dense` lambda remained in
the source but was no longer called, triggering `-Wunused-variable`.

The obsolete lambda has been removed. This is benchmark-driver cleanup only and does not
change the four graded+dense frontier runs or the Taylor-product kernel.


### 9.79 Stop graded-frontier work; return to preconditioned dense optimisation (2026-09-24)

The graded+dense cutoff frontier was completed and is now considered sufficient for the
current investigation:

```
cutoff   elapsed_s   final_error
3e-14    21.2951     3.3536682299918048e-7
1e-14    25.3351     2.7345746418298614e-7
3e-15    31.2151     2.4947609102080698e-7
1e-15    37.8401     2.4157692596632773e-7
```

The graded integrator is clearly approaching an accuracy plateau well above the
preconditioned+dense `3e-14` result (`8.5378297219108396e-8`). Further tightening of
the graded cutoff is therefore not pursued here.

Development focus returns to the QR-preconditioned integrator and specifically to the
dense full-product accumulator. The benchmark driver is reduced again to one stable
preconditioned+dense `3e-14` run so subsequent kernel changes can be compared directly
against the established ~44.7--45.1 s baseline.

The next optimisation target is the remaining per-call dense-workspace overhead:
`slot_to_touched`, `touched_slots`, and `touched_coefficients` are currently
allocated afresh in every `_ifma`. The immediate investigation should determine whether
reusing workspace storage (or moving it to a small reusable helper object) yields a
material gain before changing arithmetic or cutoff semantics.


### 9.80 Reusable dense accumulator workspace (2026-09-24)

The QR-preconditioned dense kernel still allocated and initialised its dense workspace
vectors afresh for every `_ifma` call. A reusable `thread_local` workspace is now
used instead. The slot map grows only when a larger product is encountered; touched-slot
and coefficient vectors retain capacity between calls. At the end of each product only
the slots actually used are reset to the unused sentinel.

The arithmetic and cutoff semantics are unchanged: products still use `mul_err`,
collisions still use `add_err`, and sweeping still happens only after complete
coefficient aggregation.

The benchmark remains a single `preconditioned_dense_3e-14` run. Lightweight workspace
counters report total calls, slot-map growth events, touched/coefficient capacity growth
events, maximum slot count, and maximum touched count.


### 9.81 Fix dense-workspace statistics linkage (2026-09-24)

The first reusable-workspace build failed because the header-only/template workspace code
in `taylor_model.tpl.hpp` referenced the implementation-only global
`g_taylor_model_dense_workspace_stats`, which is defined in `taylor_model.cpp` and is
not visible at template instantiation sites.

The template no longer touches that global directly. Two public recording helpers declared
in `taylor_model.hpp` and implemented in `taylor_model.cpp` now update the statistics:
one for workspace preparation/growth and one for maximum touched count. This fixes linkage
and keeps the profiling state encapsulated in the implementation unit. Dense arithmetic
and workspace reuse are otherwise unchanged.


### 9.82 Pre-rank dense operand slots (2026-09-24)

The reusable-workspace experiment completed at:

```
elapsed_s   44.7151
final_error 8.5378297219108396e-8
calls       1,256,468
slot_resizes 32
capacity_grows 144
max_slot_count 29,791
max_touched_count 325
```

This is indistinguishable from the prior ~44.74--45.05 s baseline. Workspace allocation
and growth are therefore not a material bottleneck.

The next dense-kernel optimisation removes redundant multi-index work from the hot
coefficient-product loop. Previously each pair performed:

```
product_index = alpha + beta
slot = rank(product_index)
```

which constructs the summed `MultiIndex` and then traverses its coordinates again to
compute the mixed-radix slot.

Because the mixed-radix ranking is linear in the exponent vector,

```
rank(alpha + beta) = rank(alpha) + rank(beta)
```

the dense kernel now precomputes one slot per source term of `x` and `y`. The inner
pair loop then computes the destination slot with a single integer addition and updates
the dense accumulator directly. `mul_err`, `add_err`, full coefficient aggregation,
and final sweeping are unchanged.

This targets the ~208 million product-pair iterations observed in the preconditioned
workload rather than per-call allocation overhead.


### 9.83 Dense hot-loop phase profile after operand pre-ranking (2026-09-25)

Pre-ranking the operand multi-indices produced a material improvement:

```
before pre-ranking: 44.7151 s
after pre-ranking:  39.6101 s
final error:        8.5378297219108396e-8 (unchanged)
```

The improvement is about 11.4%, confirming that repeated multi-index construction and
ranking inside the ~208 million coefficient-product iterations was a significant cost.

The next diagnostic instruments the dense kernel at phase granularity, deliberately
avoiding clock reads inside the coefficient-pair loop so the measurement does not
substantially perturb the hot path. It records cumulative time for:

- degree/workspace preparation;
- operand-slot pre-ranking;
- the complete coefficient-pair loop;
- final touched-slot sort, Expansion emission, and sweep.

It also counts product pairs, newly occupied destination slots, and collision updates.
Because every product pair necessarily executes one `mul_err` and every collision
executes one `add_err`, these counters quantify the arithmetic call volumes without
timing every arithmetic operation individually.

Output marker:

```
[TaylorDenseHotLoopProfile]
```

The benchmark remains the single QR-preconditioned dense `3e-14` run.


### 9.84 Dense hot-loop profile result (2026-09-25)

The phase-profiled preconditioned+dense `3e-14` run completed with unchanged final
error and reported:

```
elapsed_s          40.6941
product_pairs      208,224,840
new_slots           36,018,650
collision_slots    172,206,190
prepare_seconds          0.130099
prerank_seconds          0.162433
pair_loop_seconds       22.2758
emit_sweep_seconds       1.03001
```

The profiler adds about 1.08 s versus the immediately preceding unprofiled 39.6101 s run,
so its absolute runtime should not replace the uninstrumented performance baseline.

The result nevertheless localises the dense-kernel cost very clearly. Preparation and
operand pre-ranking are negligible, and final sort/emission/sweep is only about one
second. The coefficient-pair loop accounts for about 22.3 seconds. Of 208.2 million
products, 172.2 million (82.7%) update an already occupied slot and therefore execute
both `mul_err` and `add_err`; only 36.0 million create a new slot.

This confirms that further work should target the arithmetic/update path inside the pair
loop rather than workspace allocation, ranking, sorting, or sweeping.

The benchmark driver also still contained the no-longer-used local `run_graded` lambda
after development focus returned to the preconditioned integrator. It has been removed
to eliminate the `-Wunused-variable` warning; this does not alter the benchmark.


### 9.85 Fuse collision multiply-and-add with validated fma (2026-09-25)

The dense hot-loop profile showed:

```
product_pairs       208,224,840
new_slots            36,018,650
collision_slots     172,206,190
pair_loop_seconds        22.2758
```

Thus 82.7% of coefficient products hit an already occupied destination slot. The previous
collision path performed two validated arithmetic primitives:

```
product = mul_err(x,y,error)
accumulator = add_err(accumulator,product,error)
```

Inspection of `model_utilities.hpp` found Ariadne already provides a validated
`fma_err(x,y,z,error)` primitive. It computes the nearest `x*y+z` result and
accumulates an outward-rounded error contribution while managing the rounding mode only
once for the fused operation.

The dense collision path now uses:

```
accumulator = fma_err(x,y,accumulator,error)
```

New destination slots still use `mul_err`, since there is no prior coefficient to add.
This preserves validated arithmetic while removing the intermediate product and one
separate validated addition on the 172.2 million collision updates.

The temporary phase timers/counters inserted for diagnosis are disabled in the benchmark
path so the next runtime is directly comparable with the uninstrumented 39.6101 s
pre-ranked baseline.


### 9.86 Validated fused collision update result (2026-09-25)

Replacing the dense collision path's separate `mul_err` plus `add_err` with Ariadne's
validated `fma_err` produced another large speedup:

```
pre-ranked mul_err+add_err   39.6101 s   8.5378297219108396e-8
pre-ranked fma_err           34.7261 s   8.5378288508794491e-8
```

Runtime improves by about 12.3% relative to the pre-ranked baseline, while the reported
final error changes only at the 1e-14 absolute scale and is slightly smaller. The run
completes all 2000 reach sets with final radius 0.0403.

The final composition profiles are `flowpipe_compose_seconds=10.2579`,
`endpoint_compose_seconds=4.01251`, and `state_compose_seconds=4.67633`, consistent
with removing one validated arithmetic primitive from the 172.2 million collision
updates measured previously.

The fused path remains rigorous: new slots use `mul_err`, while repeated contributions
use `fma_err`.


### 9.87 Re-profile the dense pair loop after validated fma (2026-09-25)

The previous 22.2758 s pair-loop measurement predates the switch to `fma_err`. To decide
whether further optimisation should remain inside `_ifma`, the benchmark now measures
only the coefficient-pair loop around the fused implementation. This keeps instrumentation
lighter than the earlier four-phase profile and reports:

```
[TaylorDenseFmaPairProfile]
calls
product_pairs
new_slots
collision_slots
pair_loop_seconds
```

The uninstrumented performance baseline remains 34.7261 s. The profiled runtime is used
only to localise the remaining cost, not as the new performance baseline.


### 9.88 Rounding-mode switch cost probe (2026-09-25)

The fused dense pair-loop profile reports 16.0404 s for 208,224,840 coefficient pairs,
down from 22.2758 s before `fma_err`. The uninstrumented fused baseline remains about
34.7 s.

Inspection of Ariadne's rounding implementation is especially relevant on the current
Apple arm64 build. The x86-specific SSE/GCC branches do not apply there, so
`set_builtin_rounding_to_nearest()` and `set_builtin_rounding_upward()` use the C99
`fesetround` path. Both `mul_err` and `fma_err` switch to nearest for the central
coefficient and back to upward for the rigorous error bound. Consequently every one of
the ~208 million product pairs performs two rounding-mode changes.

Before redesigning validated arithmetic, the benchmark now restores the uninstrumented
integration path and runs a separate post-orbit calibration of ten million
nearest/upward pairs. Because this probe runs after the integration stopwatch has
stopped, it does not contaminate the reported integrator runtime.

Output marker:

```
[TaylorRoundingModeProbe]
pairs
switches
elapsed_seconds
seconds_per_switch
```

The result will quantify how much of the remaining 16 s dense pair loop can plausibly be
attributed to repeated `fesetround` calls and whether a batched-rounding redesign is
worth the additional numerical complexity.


### 9.89 Two-pass batched rounding for the dense raw-float product (2026-09-25)

The isolated rounding probe measured:

```
20,000,000 rounding-mode switches = 0.300001 s
15.0 ns per switch
```

At 208,224,840 product pairs, the current `mul_err`/`fma_err` implementation performs
about 416.45 million nearest/upward mode changes, corresponding to roughly 6.25 s at the
measured switch cost. This is large enough to justify an experimental batched-rounding
implementation.

For raw floating-point coefficients only, the dense product now performs two passes.

**Pass 1: centre coefficients under one nearest-rounding phase**

- set rounding to nearest once;
- new slots compute `mul(rounded,x,y)`;
- collisions compute `fma(rounded,x,y,z)`;
- for each collision, save the exact pre-update centre coefficient `z`;
- save one byte per pair identifying new-slot versus collision.

These centre updates are the same operations used internally by `mul_err` and
`fma_err`, in the same pair order.

**Pass 2: rigorous roundoff accounting under one upward-rounding phase**

- set rounding upward once;
- replay the same product-pair order;
- new slots compute the same upward `u` and `ml` values used by `mul_err`;
- collisions use the saved pre-update `z` to compute the same upward `u` and `ml`
  values used by `fma_err`;
- call `acc_err` in the original pair order.

The key correctness point is that `product_roundoff` does not feed back into the centre
coefficient recurrence. Delaying its updates until the second pass therefore preserves
the centre sequence, while saving the pre-update collision coefficient preserves the
inputs needed for the original per-operation error bound. Error contributions are then
accumulated in the same order as before.

Non-raw coefficient types retain the existing per-operation `mul_err`/`fma_err`
fallback. The post-run synthetic rounding probe has been removed so the next runtime is
a clean integrator measurement. The reference fused baseline is about 34.7--34.8 s with
final error `8.5378288508794491e-8`.


### 9.90 Fix batched-rounding build: local error accumulation (2026-09-25)

The first batched-rounding build failed because `acc_err` is not part of the public
Taylor-model template interface. It is a helper defined in the anonymous namespace of
`model_utilities.hpp`, so `taylor_model.tpl.hpp` cannot name it directly.

The batched path now spells out exactly the same upward-rounded accumulation performed
by `acc_err`:

```
error = add(rounded,error,hlf(add(rounded,ml,u)))
```

This is applied in both the new-slot and collision branches of the second (upward)
pass. No arithmetic formula or ordering has changed; this commit only removes the
invalid dependency on an implementation-local helper.


### 9.91 Batched-rounding dense result (2026-09-25)

The two-pass raw-float batched-rounding implementation produced the largest dense-kernel
speedup so far:

```
per-pair validated fma baseline   34.7261--34.8241 s
batched nearest/upward passes     22.4231 s
final error before                8.5378288508794491e-8
final error after                 8.5378223539938479e-8
final radius                      0.0403
reach sets                        2000
```

This is approximately a 35.5% runtime reduction relative to the 34.7261 s fused baseline,
and about a 59.6% reduction relative to the ~55.5 s first final-sweep implementation.
The reported final error is slightly smaller by about 6.50e-14 absolute; there is no
accuracy regression in this benchmark.

The composition timings collapse correspondingly:
`flowpipe_compose_seconds=3.98765`,
`endpoint_compose_seconds=1.69752`, and
`state_compose_seconds=1.60962`.

The improvement is substantially larger than the ~6.25 s estimate from the isolated
`fesetround` probe. The batched implementation therefore benefits not only from reducing
rounding-mode switches from O(product-pairs) to O(_ifma calls), but also from separating
the centre and error-bound arithmetic into uniform rounding phases. This likely improves
the generated hot loops and avoids repeated helper/rounding-control overhead.

Correctness rationale remains the one established in section 9.89: the nearest centre
coefficient recurrence is preserved in pair order; every collision stores its exact
pre-update centre value; the upward pass reconstructs the same `u` and `ml` inputs
used by `mul_err`/`fma_err`; and roundoff contributions are accumulated in the same
pair order. The small final-error difference is consistent with changed execution
context/rounding boundaries but should still be covered by dedicated equivalence tests
before this optimisation is considered production-ready.

The new performance baseline for the experimental QR-preconditioned dense integrator at
cutoff `3e-14` is 22.4231 s.


### 9.92 Dense batched-rounding equivalence test harness (2026-09-25)

Before treating the 22.4231 s batched kernel as production-ready, a direct A/B correctness
path has been added to the Taylor-model tests.

A new runtime switch,
`taylor_model_dense_batched_rounding_enabled()`, selects between two implementations
inside the same dense pre-ranked accumulator:

- **off:** the validated per-pair reference path using `mul_err` for first contributions
  and `fma_err` for collisions;
- **on:** the two-pass nearest/upward batched-rounding path.

This isolates arithmetic scheduling from indexing, sweeping, and accumulator structure.

`tests/function/test_taylor_model.cpp` now multiplies deterministic Taylor models under
both settings and requires exact equality of both the resulting Expansion and Error. The
cases cover:

- dense univariate collision chains with mixed signs and non-dyadic values;
- multivariate collisions and cancellation;
- a wide finite dynamic range;
- nonzero input-model errors.

The test is instantiated by the existing Taylor-model suite for both FloatDP and FloatMP.
The previous global toggle state is restored after the comparison.

This is intentionally stronger than merely checking enclosure overlap or approximate
equality: if the two-pass transformation really reconstructs the same per-operation
centre recurrence and error contributions, the results should be bit-identical for these
deterministic cases. Any failure is evidence that the equivalence argument is incomplete
and must be investigated before enabling batched rounding by default.


### 9.93 Rigorous batched-rounding enclosure oracle after non-bit-identical result (2026-09-25)

The first A/B test rejected bit identity: in the collision-heavy FloatDP case, two centre
coefficients differed by one last-place rounding step while the Error remained identical.
Bit equality is therefore too strong a criterion for the changed arithmetic schedule.

The replacement test now checks the property that matters mathematically. For every
zero-input-error product it constructs an independent outward-rounded coefficient oracle
with `Bounds<F>`: source coefficients are singleton intervals, products and collision
sums are interval operations, and the maximum distance from each exact-coefficient
interval to the batched centre coefficient is summed upward. Since Taylor-model variables
are normalised to [-1,1], that sum is a rigorous sup-norm bound for the centre-polynomial
error. The test requires

```
coefficient_error_bound <= batched.error()
```

and also requires the batched Error not to be smaller than the trusted per-pair Error.

The per-pair comparison remains diagnostic rather than normative. It now reports the
number of differing centre coefficients, maximum absolute difference, and for FloatDP the
maximum ULP distance. Nonzero-input-error cases retain the validated error-budget
comparison; the independent coefficient oracle is intentionally restricted to the
zero-error polynomial core.


### 9.94 Fix rigorous-oracle test compilation (2026-09-25)

The first rigorous-oracle test did not compile for two independent test-code reasons:

1. `abs(Float)` resolves to validated arithmetic and therefore returns `Bounds<Float>`,
   so it cannot be assigned to a raw `Float`. The maximum absolute centre difference is
   diagnostic only, so it is now computed as a `double` from `get_d()`; it is not used
   in any rigorous assertion.
2. `MultiIndex` has equality but no ordering relation suitable for `std::map`. The
   small deterministic oracle now uses a vector of `(MultiIndex, Bounds<F>)` pairs and
   a linear equality lookup. This affects test bookkeeping only; the rigorous Bounds
   arithmetic and the enclosure assertion are unchanged.

The correctness condition remains
`coefficient_error_bound <= batched.error()`, with the coefficient error bound built
entirely using outward-rounded Ariadne arithmetic.


### 9.95 Root cause of the batched/reference mismatch: fused versus non-fused collision arithmetic (2026-09-25)

The failed exact-equivalence test exposed a real implementation mismatch in the
experimental batched path, but not a flaw in the batching idea itself.

The dense reference path calls the `fma_err` helper defined locally in
`taylor_model.tpl.hpp`. Despite its name, that helper does **not** use a hardware/semantic
fused multiply-add. Its centre update is:

```
rv = xv * yv + zv
```

through `Rounded<F>` operators, i.e. a rounded multiplication followed by a rounded
addition. Its upward error-bound calculations use the same separate multiply/add
sequence.

The batched implementation had instead used:

```
fma(rounded,xv,yv,zv)
```

for both the centre and error-bound collision calculations. That is a genuinely fused
operation and can differ by one ulp from the separate multiply-then-add sequence. This
exactly explains the observed FloatDP differences while the accumulated Error happened
to remain equal in the first regression case.

The batched collision path has now been corrected to reproduce the local reference
primitive exactly:

```
add(rounded,mul(rounded,xv,yv),zv)
```

in both nearest and upward phases. New-slot multiplication was already identical.

Because the batched path now executes the same low-level arithmetic operations in the
same product-pair order and accumulates the same error contributions in the same order,
the direct test once again requires exact representation equality of Expansion and Error
between batched and per-pair modes. The previous Bounds oracle is removed: it was testing
a different, stronger coefficient-wise interval property and failed even when batched
and reference results were bit-identical, so it was not a discriminator for this
transformation.

This correction may reduce some of the 22.4231 s performance gain because the earlier
batched version benefited from fused arithmetic in addition to amortised rounding-mode
changes. The next correctness run must pass exactly before re-benchmarking performance.


### 9.96 Exact batched/reference equivalence confirmed (2026-09-25)

After replacing the accidental fused collision arithmetic with the exact
multiply-then-add sequence used by the local `fma_err`, the dedicated A/B suite passes
exactly for all four deterministic products in both FloatDP and FloatMP.

For every case:

```
differing_coefficients = 0
max_abs_difference     = 0
same(batched.expansion(), per_pair.expansion()) = true
same(batched.error(),     per_pair.error())     = true
```

This includes the collision-heavy non-dyadic case that previously exposed two one-ulp
differences, the multivariate cancellation case, the wide-dynamic-range case, and the
nonzero-input-error case.

The experiment therefore isolates the optimisation cleanly: the batched implementation
now changes only the placement/frequency of rounding-mode switches. The coefficient
operations, their order, the error-bound operations, and the order of Error accumulation
match the trusted per-pair path exactly for the regression suite.

The next required measurement is the Van der Pol benchmark with this corrected,
bit-equivalent batching. The previous 22.4231 s result cannot be retained as the final
batched baseline because it also benefited from the accidental fused arithmetic.


### 9.97 Performance of the bit-equivalent batched-rounding kernel (2026-09-25)

The corrected batched implementation, which is bit-identical to the trusted per-pair
arithmetic in the dedicated FloatDP/FloatMP regression suite, completes the Van der Pol
benchmark at:

```
elapsed_seconds       24.0561
achieved_final_error  8.5378288508794491e-8
final_radius          0.0403
reach_sets            2000
```

The final error is exactly the same reported value as the 34.7--34.8 s per-pair
`mul_err/fma_err` reference. This is consistent with the dedicated exact-equivalence
tests and confirms that the previous 22.4231 s result also included an additional
benefit from accidentally fused collision arithmetic.

Against the clean 34.7261 s per-pair baseline, pure rounding-mode batching saves
10.6700 s, about 30.7%. Against the roughly 55.5 s first final-sweep implementation, the
total dense-kernel work has reduced runtime by about 56.6%.

Final cumulative composition costs are:

```
flowpipe_compose_seconds  4.61170
endpoint_compose_seconds  1.95459
state_compose_seconds     1.90666
```

The workspace profile remains:
1,256,468 dense calls, 32 slot-map resizes, 144 coefficient-capacity growth events,
maximum slot count 29,791, and maximum touched count 325.

The remaining global profile is now qualitatively different from the earlier dense
kernel. Gronwall centre construction is the largest named cumulative component at
13.9638 s. Recurrence residual procedure time is 2.4839 s, while the flow-function
profile reports 2.4645 s in model construction and 1.3425 s in restriction. These
profiles may have nested scopes and must not be summed directly against wall time.

The next optimisation phase should therefore re-profile the architecture rather than
continue assuming `_ifma` is dominant. In particular, the 13.9638 s Gronwall centre
path should be decomposed before further dense-kernel micro-optimisation.


### 9.98 Decompose the centre-polynomial recurrence cost (2026-09-25)

With the bit-equivalent batched dense kernel, the Van der Pol baseline is 24.0561 s and
the existing Gronwall profile reports 13.9638 s in centre-polynomial construction. That
outer timer includes the entire `graded_series_centre_polynomial_step`, so it is too
coarse to identify the next optimisation target.

The recurrence helper is now instrumented at phase granularity, with clock reads only at
phase boundaries:

- `graded_flow_init`;
- all temporal `graded_flow_iterate` calls;
- `flow_differential` for the centre polynomial;
- `flow_function` materialisation of the centre polynomial;
- the retained-state final Procedure evaluation used for `g(P_m)`;
- conversion of that result to a Differential;
- `flow_function` materialisation of the recurrence field;
- the existing direct-defect path.

The new cumulative marker is:

```
[CentreRecurrenceCostProfile]
```

This profile should explain most of the 13.96 s outer centre timer and distinguish
recurrence arithmetic from Taylor-patch materialisation. The 24.0561 s run remains the
performance baseline; the instrumented run is diagnostic only.


### 9.99 Compare direct Differential defect against production patch defect (2026-09-25)

The centre-recurrence profile attributes about 6.54 s cumulatively to the final
`g(P_m)` Procedure evaluation, recurrence-field Taylor-patch materialisation, and the
coefficient-level direct-defect path. Before removing any patch-level work, the direct
defect must be shown to be at least as conservative as the range currently used by the
production Gronwall bound.

The production semantics are unchanged. Whenever the centre polynomial lies inside the
certification box, the code now compares, component by component,

```
mag(centre_result.direct_defect_range[i])
```

against

```
mag(defect.range()[i])
```

where the latter is the current production residual bound. The cumulative diagnostic

```
[DirectDefectRangeComparison]
```

reports:

- number of calls and compared components;
- components for which the direct range magnitude is at least the production magnitude;
- counts of strictly larger and strictly smaller direct magnitudes;
- maximum direct/production magnitude ratio;
- maximum absolute magnitude difference.

No Gronwall remainder uses the direct range yet. The purpose of this run is only to
establish whether replacing the patch-level residual range would be semantically safe,
and how much extra conservatism that replacement would introduce.


### 9.100 Fix direct-defect comparison build (2026-09-25)

The first direct-defect comparison build failed because `mag(Bounds<FloatDP>)`
returns `Error<FloatDP>`, not a raw `FloatDP`. The diagnostic conversion now uses
`Error::raw().get_d()` for both production and direct magnitudes.

This affects diagnostic reporting only; the production Gronwall path and all validated
arithmetic remain unchanged.


### 9.101 Decompose the patch residual into polynomial core and Taylor-model Error (2026-09-25)

The direct Differential residual is systematically much smaller than the production
patch residual: all 4000 compared components were smaller, with a maximum observed
direct/production magnitude ratio of about 0.0696. This means it cannot simply replace
the production range without a separate correctness argument.

The next diagnostic isolates where the gap is introduced. For each residual component
the code now materialises the same derivative and recurrence-field Taylor models used by
the production path, records their attached Errors, copies them, calls `clobber()` on
the copies to remove only those Error terms, and subtracts the resulting polynomial
cores. It then records:

- the ratio between the production residual magnitude and the error-free core residual;
- the ratio between the direct Differential residual and that same core residual;
- the maximum sum of the two source Taylor-model Errors;
- the maximum Error attached to the materialised defect model itself;
- the cumulative cost of this diagnostic decomposition.

The production Gronwall remainder is still computed from the original
`defect.range()`; this instrumentation does not alter integration semantics.

If the direct residual tracks the clobbered polynomial core while the production
residual is much larger, the excess is attributable to patch materialisation/sweeping
Errors rather than to a different polynomial defect. That would identify the exact
piece that must be certified directly before patch-level materialisation can be removed.


### 9.102 Fix defect-error decomposition build (2026-09-25)

The first defect-error decomposition build failed because `defect.get(i)` is exposed as
a type-erased `ValidatedScalarMultivariateFunctionPatch`, which has no direct
`.model()` accessor.

The diagnostic already has the two underlying Taylor models used to construct that
component:

```
derivative_model
field_model
```

so the model-level residual is now formed directly as

```
materialised_defect_model = derivative_model - field_model
```

and its attached Error is inspected from that object.

This changes only diagnostic bookkeeping. The production `defect` FunctionPatch and the
Gronwall residual range remain untouched.


### 9.103 Audit coefficient loss before patch-level residual subtraction (2026-09-25)

The previous decomposition showed that removing the attached Taylor-model Errors changes
the production residual magnitude by only about 0.2%, while the direct Differential
residual can be about 14 times smaller. The gap therefore arises before the final model
Error term.

The new diagnostic audits the two residual operands across the
Differential-to-TaylorModel materialisation boundary:

- `derivative(dphi)` versus the derivative Taylor model obtained from the already
  materialised centre polynomial;
- `recurrence_field_differential` versus the materialised recurrence-field Taylor model.

For every source Differential term it counts coefficients that disappear entirely and
coefficients whose stored value changes, accumulating the magnitudes of the missing and
changed contributions. It also compares the coefficient L1 magnitude of the residual
formed before materialisation with the coefficient L1 magnitude of the clobbered
patch-level residual after separate materialisation and subtraction.

The cumulative marker is:

```
[DefectSweepCoefficientAudit]
```

with fields for missing/changed term counts and magnitude sums for both operands,
`pre_subtract_magnitude_sum`, `post_subtract_magnitude_sum`, and diagnostic runtime.

The production Gronwall path is unchanged. The purpose is to determine whether the large
direct/production residual gap is explained by coefficients being dropped or modified
before cancellation can occur.


### 9.104 Fix coefficient audit scope and compare on the widened Taylor domain (2026-09-25)

The first coefficient audit did not compile because it was inserted in the outer
preconditioned-step routine, while the source `dphi` and
`recurrence_field_differential` objects exist only inside
`graded_series_centre_polynomial_step`. It also incorrectly treated Differential
coefficients, which are `Bounds<FloatDP>`, as raw FloatDP values.

The invalid outer-scope audit has been removed. The replacement diagnostic is placed
inside the recurrence helper, where both source Differentials are available, and avoids
comparing differently scaled representations.

On the exact same widened domain used by the direct defect, it now constructs:

```
wide_derivative = make_taylor_function_model(derivative_dphi)
wide_field      = make_taylor_function_model(recurrence_field_differential)
wide_defect     = make_taylor_function_model(
                      derivative_dphi-recurrence_field_differential)
```

It then clobbers only the attached Errors of the two separately materialised operands,
subtracts their stored polynomial cores, and compares that polynomial with the clobbered
directly materialised defect. This directly tests the order-of-operations hypothesis:

```
materialise(A-B)   versus   materialise(A)-materialise(B)
```

without mixing source Differential coordinates with restricted Taylor-patch
coordinates.

The cumulative `[DefectSweepCoefficientAudit]` reports coefficient differences,
coefficients present only on either side, cumulative coefficient L1 magnitudes, maximum
L1 ratios, maximum coefficient difference, and the respective materialisation Error
budgets. The production path remains unchanged.


### 9.105 Fix coefficient-audit build warning and zero construction (2026-09-25)

The widened-domain coefficient audit initially failed to compile because `FloatDP(0)`
selects a private raw-double constructor in this context. Zero comparisons now use the
public integral-plus-precision constructor `FloatDP(0u,dp)`.

The local diagnostic variable `da` also shadowed the existing parameter vector
`da`; it has been renamed to `direct_abs` (and the paired value to
`separate_abs`) to keep the build warning-free.

These changes affect only the diagnostic audit.


### 9.106 Widened-domain coefficient audit result: materialisation commutes with subtraction (2026-09-25)

The widened-domain coefficient audit decisively rejects the hypothesis that the
direct/production residual gap is caused by coefficient sweeping before cancellation.

Across all 2000 recurrence calls:

```
differing_coefficients        = 0
direct_only_coefficients      = 0
separate_only_coefficients    = 0
direct_coefficient_l1         = 2.27459e-8
separate_coefficient_l1       = 2.27459e-8
max_direct_to_separate_ratio  = 1
max_separate_to_direct_ratio  = 1
max_coefficient_difference    = 0
```

Thus, on the widened domain and with attached Errors removed,

```
poly(materialise(A-B))
```

is exactly coefficient-identical to

```
poly(materialise(A)-materialise(B))
```

for this benchmark. There is no lost polynomial cancellation at this boundary.

The attached Error budgets differ: the directly materialised defect reaches
`1.26506e-13`, while the sum of separately materialised source Errors reaches
`4.24978e-13`. This is real but still does not explain the much larger production
range gap observed after `flow_function` restriction.

The remaining structural difference is therefore the restriction from the widened time
domain to the forward time interval. The direct diagnostic evaluates the widened defect
model directly on the forward half of its normalised time coordinate, whereas the
production path separately restricts the centre polynomial and recurrence field before
forming the residual. The next investigation should isolate this restriction step.

This audit is intentionally expensive: it adds two extra Taylor-model materialisations
per recurrence call and accumulates 2.69575 s of diagnostic work, increasing the observed
wall time to 27.3211 s. That runtime is not a performance baseline.

The remaining local variable shadow warning in the audit has also been removed by
renaming `separate_abs` to `separate_coeff_abs`.


### 9.107 Restriction-order audit for the centre defect (2026-09-25)

The widened-domain coefficient audit established exact polynomial identity between
`materialise(A-B)` and `materialise(A)-materialise(B)` on all 2000 calls, so the
direct/production residual gap does not arise from pre-restriction sweeping or lost
cancellation.

That expensive audit has now been removed. The remaining structural difference is
restriction from the widened time domain to the forward time interval. The new diagnostic
therefore compares, on identical widened and forward domains,

```
restriction(materialise(A-B))
```

with

```
restriction(materialise(A)) - restriction(materialise(B))
```

where `A=derivative(dphi)` and `B=recurrence_field_differential`.

The cumulative marker `[DefectRestrictionAudit]` reports coefficient differences after
clobbering attached Errors, coefficients present only on either side, maximum coefficient
difference, direct and separate Error budgets, direct and separate range magnitudes, and
the maximum separate/direct range-magnitude ratio.

The production Gronwall path is unchanged. This isolates whether restriction is the first
operation at which the two residual constructions cease to agree.


### 9.108 Restriction is the first non-commuting transformation (2026-09-25)

The restriction-order audit identifies the first concrete transformation at which the
direct and separately materialised residual constructions diverge.

Before restriction, the widened-domain audit found exact coefficient identity. After
restriction, by 2000 recurrence calls the cumulative diagnostic reports:

```
differing_coefficients                  33626
direct_only_coefficients                    3
separate_only_coefficients               2737
max_coefficient_difference          5.25375e-14
max_direct_error                    6.94198e-13
max_separate_error                  1.18291e-12
max_direct_range_mag                4.22746e-11
max_separate_range_mag              4.26437e-11
max_separate_to_direct_range_ratio  104.674
```

Thus restriction does not commute with subtraction in the current Taylor-model
implementation. Restricting the two operands separately introduces thousands of
coefficients absent from the restricted combined defect and changes many shared
coefficients. The separate path also carries a larger Error budget.

This explains the qualitative source of the large direct/production residual gap:
the two paths are identical before restriction and diverge at restriction. The very large
maximum range ratio occurs on components or steps where the direct residual is extremely
small, so it should not be interpreted as a typical factor; nevertheless it establishes
that separate restriction can strongly inhibit cancellation.

The audit is expensive: it adds extra restriction/materialisation work, accumulates
4.45028 s of diagnostic time, raises dense calls from 1,256,468 to 1,866,748, and
increases wall time to 28.9991 s. This run is diagnostic only; the clean performance
baseline remains 24.0561 s.

The production Gronwall path still uses the separately restricted centre polynomial and
recurrence field. A production optimisation must not simply substitute the smaller direct
range: the next step is to establish a rigorous validated argument for forming the
residual on the widened domain and restricting the residual once, including the
Taylor-model Error semantics of restriction.


### 9.109 Production experiment: subtract on the widened domain, restrict once (2026-09-25)

The restriction audit in 9.108 identified the first non-commuting transformation: the
validated derivative and recurrence-field operands are coefficient-identical before
restriction, but restricting them separately introduces thousands of coefficient
differences and a larger Taylor-model Error budget before the residual subtraction.

The production path is now changed in the narrowest way that preserves the existing
validated operand semantics:

```
wide_derivative = materialise(dP/dt, widened_domain)
wide_field      = materialise(g(P), widened_domain)
wide_defect     = wide_derivative - wide_field
defect          = restriction(wide_defect, forward_domain)
```

This is deliberately **not** the cheaper `materialise(dP/dt-g(P))` shortcut.  The two
operands are still materialised separately with the existing Taylor-model Error
propagation, so no new assumption is made about omitted Procedure terms or generic
non-polynomial dynamics.  The only semantic change is to subtract the two validated
Taylor models before the already-validated restriction operation.

The old production sequence,

```
restriction(materialise(dP/dt))
    - restriction(materialise(g(P)))
```

and its now-obsolete defect decomposition/restriction diagnostics have been removed from
the hot path.  The cheap combined-Differential residual remains diagnostic only.

**Expected effect:** preserve more polynomial cancellation, reduce the residual range and
remove one restriction of a large operand.  This may both tighten the Gronwall remainder
and reduce runtime.  The clean reference remains:

```
elapsed_seconds       24.0561
achieved_final_error  8.5378288508794491e-8
final_radius          0.0403
reach_sets            2000
```

This commit has not yet been benchmarked in the connected environment.  Required next
step:

```bash
ninja vanderpol
./examples/continuous/vanderpol > widened_defect.txt 2>&1
```

Acceptance gates:
1. build and regression tests remain valid;
2. all 2000 reach sets complete;
3. final radius does not regress;
4. final Taylor-model error is no worse than the 24.0561 s reference path;
5. wall time improves materially, or the tighter residual produces a measurable accuracy
   benefit that justifies its cost.

If the new path is slower despite the tighter residual, profile the two widened
materialisations separately before changing the recurrence architecture again.


### 9.110 Widened-domain subtract-before-restrict benchmark result (2026-09-25)

The production experiment from 9.109 completed all 2000 Van der Pol reach sets.

Measured result:

```
method                preconditioned_dense_3e-14
elapsed_seconds       27.3231
achieved_final_error  4.7221144502652554e-8
final_radius          0.0403
reach_sets            2000
```

Reference bit-equivalent batched baseline:

```
elapsed_seconds       24.0561
achieved_final_error  8.5378288508794491e-8
final_radius          0.0403
reach_sets            2000
```

Thus subtracting the two validated operands on the widened domain and restricting the
combined residual once improves the final Taylor-model error by about 44.7%, while
preserving the final radius and reach-set count.  It is, however, about 13.6% slower
(+3.267 s) than the 24.0561 s reference.  This is therefore a real accuracy improvement,
but not yet a runtime optimisation.

The run also changes the dense workload materially:

```
Taylor dense calls: 985,603
```

versus 1,256,468 calls in the 24.0561 s reference profile.  Despite fewer dense calls,
the centre-polynomial path grows to 15.9423 s.  Its final cumulative decomposition is:

```
graded_flow_iterate              6.15741 s
final Procedure g(P)             2.81835 s
validated defect materialisation 2.88531 s
centre polynomial flow_function  1.84397 s
cheap direct defect diagnostic   2.03148 s
```

The 2.03148 s cheap direct-defect path is still diagnostic-only and is not used by the
production Gronwall certificate.  It is therefore now pure benchmark overhead and should
be removed from non-diagnostic runs before judging the widened-domain production design.
Likewise, exact-polynomial diagnostic bookkeeping remains negligible when disabled.

**Conclusion:** keep the subtract-before-restrict construction for the next experiment
because it produces a substantially tighter certified result.  First remove the
non-production direct-defect calculation from clean runs and re-benchmark.  Do not
revert solely from the 27.3231 s timing: roughly two seconds of that run are explicitly
known diagnostic work.

If the cleaned run remains slower than 24.0561 s, decompose the two widened Taylor-model
materialisations and investigate whether the already available centre-polynomial model
can supply the derivative operand without rematerialising it.


### 9.111 Remove direct-defect diagnostic work from clean production runs (2026-09-25)

The widened-domain subtract-before-restrict benchmark in 9.110 completed at 27.3231 s,
but its centre-recurrence profile showed 2.03148 s spent in the cheap
materialise-after-subtraction direct-defect path.  That path is diagnostic only and does
not contribute to the production Gronwall certificate.

It is now executed only when integrator diagnostics are enabled.  Clean benchmark runs
retain the validated widened-domain production defect from 9.109 but skip:

```
direct_defect_differential = dP/dt - g(P)
make_taylor_function_model(direct_defect_differential, ...)
evaluate(..., forward_half_box)
```

No production enclosure semantics are changed by this commit; it only removes known
diagnostic overhead from the benchmark path.

Reference to beat for the tightened widened-domain variant:

```
elapsed_seconds       27.3231
achieved_final_error  4.7221144502652554e-8
final_radius          0.0403
reach_sets            2000
```

The historical behaviour-preserving baseline remains 24.0561 s at
8.5378288508794491e-8 final error.  The immediate objective is runtime: first recover
the known ~2 s diagnostic cost, then continue profiling the widened-domain residual
construction if the cleaned result is still slower than 24.0561 s.


### 9.112 Clean widened-domain benchmark: near parity with tighter error (2026-09-25)

After moving the direct-defect diagnostic behind the diagnostics switch, the clean
Van der Pol run completes with:

```
elapsed_seconds       24.6201
achieved_final_error  4.7221144502652554e-8
final_radius          0.0403
reach_sets            2000
dense_calls           985603
```

Compared with the previous widened-domain run (27.3231 s), this recovers 2.703 s
(~9.9%) while leaving the final error exactly unchanged.  Compared with the historical
bit-equivalent 24.0561 s baseline, the tightened variant is only 0.564 s (~2.35%) slower
while its final Taylor-model error is about 44.7% smaller.

The profile now reports the principal centre costs at 2000 calls as:

```
graded_flow_iterate              5.68528 s
final Procedure g(P)             2.49065 s
validated widened defect         2.67860 s
centre polynomial flow_function  1.72418 s
centre total                    14.5972 s
```

The log still reports `direct_defect_seconds=1.85034`, despite the direct-defect
calculation being guarded by the diagnostics flag.  Inspection shows that the stopwatch
itself was left outside the guard, so this number is timing/accounting overhead rather
than the removed diagnostic computation.  The stopwatch is now moved inside the
diagnostic branch and clean runs return zero direct-defect time.

Performance remains the first objective.  The widened-domain variant is retained because
it is now close to the 24.0561 s baseline and materially tighter.  The next optimisation
target, after obtaining a clean profile with corrected timing, is the 2.68 s validated
widened-defect construction, especially avoiding redundant materialisation of the
derivative operand if the already materialised centre polynomial can provide equivalent
validated data without reintroducing separate restriction.


### 9.113 Clean profile establishes a new best production baseline (2026-09-25)

With the direct-defect diagnostic and its stopwatch completely absent from clean runs, the
widened-domain subtract-before-restrict variant completes the Van der Pol benchmark at:

```
elapsed_seconds       22.7361
achieved_final_error  4.7221144502652554e-8
final_radius          0.0403
reach_sets            2000
dense_calls           985603
```

This is the best measured production point so far. Relative to the previous
behaviour-preserving bit-equivalent baseline (24.0561 s,
8.5378288508794491e-8), runtime improves by about 5.49% while final Taylor-model error
improves by about 44.7%. Relative to the first roughly 55.5 s final-sweep implementation,
the runtime reduction is about 59.0%.

The corrected clean profile confirms `direct_defect_seconds=0`. At 2000 recurrence
calls the main centre-polynomial costs are:

```
graded_flow_iterate              5.70563 s
validated widened defect         2.68550 s
final Procedure g(P)             2.47454 s
centre polynomial flow_function  1.75418 s
centre total                    12.7763 s
```

Carried-state cumulative composition costs are:

```
flowpipe_compose_seconds  4.59934
endpoint_compose_seconds  1.93101
state_compose_seconds     1.90520
```

The widened-domain construction is therefore retained as the new baseline. Performance
remains the first objective.

Next optimisation target: decompose the 2.6855 s validated widened-defect construction
into the following phase costs before changing its algorithm:

1. materialising the widened derivative model;
2. materialising the widened recurrence-field model;
3. subtracting the two model vectors;
4. restricting the combined residual once.

Instrument only phase boundaries. The 22.7361 s clean run remains the performance
reference; any instrumented wall time is diagnostic only.


### 9.114 Phase profile for the validated widened defect (2026-09-25)

The 22.7361 s clean baseline leaves 2.6855 s in the validated widened-defect construction.
Before changing that algorithm, the production path is instrumented only at phase
boundaries.  The new cumulative marker is:

```
[WidenedDefectCostProfile]
```

It reports four mutually sequential phases:

```
derivative_materialise_seconds
field_materialise_seconds
subtract_seconds
restrict_seconds
```

corresponding to:

```
make_taylor_function_model(dP/dt, widened_domain)
make_taylor_function_model(g(P), widened_domain)
wide_derivative - wide_field
restriction(wide_defect, forward_domain)
```

No arithmetic, enclosure, cutoff, or step-selection semantics are changed.  The
instrumented wall time is diagnostic only; the clean performance reference remains
22.7361 s with final error 4.7221144502652554e-8.


### 9.115 Widened-defect phase profile: materialisation dominates (2026-09-25)

The four-phase profile completes all 2000 Van der Pol reach sets with unchanged final
error and radius.  The instrumented wall time is 21.7731 s, but the clean reference
remains 22.7361 s because timing noise and instrumentation make this run unsuitable as a
new performance baseline.

Cumulative widened-defect costs at 2000 calls are:

```
derivative materialisation  0.987009 s
field materialisation       1.39446 s
model subtraction           0.0242431 s
combined restriction        0.0853471 s
total measured phases       2.49106 s
outer defect timer          2.56990 s
```

Thus about 95.6% of the measured phase cost is in the two
`make_taylor_function_model` calls.  Field materialisation is the largest individual
phase (about 56.0% of the four-phase total), derivative materialisation is second (about
39.6%), while subtraction and the single final restriction together are only about
4.4%.

This rejects further optimisation of the final restriction as the immediate target.
The next experiment should target derivative materialisation first because an already
materialised centre polynomial exists, making that operand the one with the clearest
potential redundancy.  However, replacing the widened derivative by differentiating the
already restricted centre polynomial would reintroduce the restriction-order problem and
is not equivalent.

Before changing production semantics, add a diagnostic equivalence/cost experiment that
constructs the derivative Taylor model directly from the already materialised widened
centre model, before its forward restriction.  If coefficient and Error semantics match
the current `make_taylor_function_model(derivative_dphi, wide_domain)` operand, retain
the derived model and eliminate the redundant ~0.99 s materialisation.  If they do not
match, keep the current validated path and move to field-materialisation optimisation.
