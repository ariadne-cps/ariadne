# Preconditioned Graded Taylor Series Integrator: Investigation Log

**Branch:** `solvers-integrator#357`  
**Last updated:** 2026-09-23  
**Current HEAD when this log was created:** `cb1496eb436a1d4ed226554a4f18eaa4da39f29a`

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

