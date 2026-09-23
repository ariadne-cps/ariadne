# Preconditioned Graded Taylor Series Integrator: Investigation Log

**Branch:** `solvers-integrator#357`  
**Last updated:** 2026-09-23  
**Current HEAD when this log was created:** `cb1496eb436a1d4ed226554a4f18eaa4da39f29a`  
**Latest analysed investigation HEAD:** `a0f4842485b8dddd58fde132203f10603071a458`

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

