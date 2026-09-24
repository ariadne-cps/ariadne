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
