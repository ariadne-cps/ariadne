# Preconditioned Graded Taylor Series Integrator: Investigation Log

**Branch:** `solvers-integrator#357`  
**Last updated:** 2026-09-23  
**Current HEAD when this log was created:** `cb1496eb436a1d4ed226554a4f18eaa4da39f29a`  
**Latest analysed investigation HEAD:** `639782197e01889cad88e0b7ee9e66260f8ceb5d`

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

