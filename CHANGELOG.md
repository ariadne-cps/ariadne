## Changelog

Each release lists the issues it addresses. The issue number is followed by the kind of feature (in parentheses) and a short description. 

Legenda for the issue kind:

- N: new feature
- A: addition to a feature
- C: change to the behavior of a feature
- F: fixed feature
- R: removed feature

### 2.2

*Date: 25/04/2021*

- [#441](https://github.com/ariadne-cps/ariadne/issues/441) (N) Add support for Gnuplot output, including animated gif plot of sets and tridimensional plots for PDEs
- [#507](https://github.com/ariadne-cps/ariadne/issues/507) (N) Add a simulator for vector field dynamics
- [#514](https://github.com/ariadne-cps/ariadne/issues/514) (N) Add Python examples in python/examples
- [#509](https://github.com/ariadne-cps/ariadne/issues/509) (A) Additionally support a set as input in simulators, using the midpoint as the effective point
- [#513](https://github.com/ariadne-cps/ariadne/issues/513) (A) Add missing Python bindings for verify_safety in (Hybrid)ReachabilityAnalyser
- [#516](https://github.com/ariadne-cps/ariadne/issues/516) (A) Add missing Python bindings for Real predicates to be used in automata specification
- [#518](https://github.com/ariadne-cps/ariadne/issues/518) (A) Add missing Python bindings for evolver configuration and initial set assignment
- [#520](https://github.com/ariadne-cps/ariadne/issues/520) (A) Add missing Python bindings for plotting using HybridFigure
- [#543](https://github.com/ariadne-cps/ariadne/issues/543) (A) Add missing Python bindings for iterating across ListSet of Enclosure classes
- [#527](https://github.com/ariadne-cps/ariadne/issues/527) (A) Allow to draw a Labelled/Hybrid orbit directly to a Labelled/Hybrid figure
- [#492](https://github.com/ariadne-cps/ariadne/issues/492) (C) Modify SFINAE code to use C++20 concepts, currently preventing AppleClang compilation under macOS until the compiler supports Concepts
- [#529](https://github.com/ariadne-cps/ariadne/issues/529) (C) Disallow construction of VectorField and IteratedMap from a Function, since it was broken
- [#533](https://github.com/ariadne-cps/ariadne/issues/533) (C) Map<K,V> now checks for existing key using ARIADNE_ASSERT, yielding errors also for Release builds
- [#447](https://github.com/ariadne-cps/ariadne/issues/447) (F) Check that a VectorField is defined with dynamics for all involved variables, fixes a segfault within evolution
- [#532](https://github.com/ariadne-cps/ariadne/issues/532) (F) RealExpressionBoundedConstraintSet could be constructed in an incoherent way, due to missing checks
- [#211](https://github.com/ariadne-cps/ariadne/issues/211) (R) Remove various deprecated functions

### 2.1 

*Date: 09/03/2021*

- [#412](https://github.com/ariadne-cps/ariadne/issues/412) (N) Supply Aptitude and Homebrew packages for installation
- [#267](https://github.com/ariadne-cps/ariadne/issues/267) (A) Complete Python bindings for dynamics and hybrid modules
- [#477](https://github.com/ariadne-cps/ariadne/issues/477) (C) Make the tutorials the same for C++ and Python

### 2.0

*Date: 17/04/2020*

First release.