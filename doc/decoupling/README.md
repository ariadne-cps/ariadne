# Decoupling Ariadne components

This directory collects the architectural analysis and the evolving plan for
separating the components currently stored under `source/` so that they can be
built, tested, versioned, and eventually developed in separate repositories.

The documents describe dependencies in the direction **provider → dependent**.
An arrow therefore points to the component that consumes the other component.

## Documents

- [Coupling analysis](coupling-analysis.md): baseline at commit
  `949d7e044ae65837fc02e6387701b10e7c15ddc6`, scoring method, complete
  dependency matrix, architectural interpretation, and supporting evidence.
- [Working plan](work-plan.md): living document for milestones, current work,
  decisions, risks, and validation results.
- [Dependency graph](dependency-graph.svg): full directed graph. A V arrowhead
  denotes low coupling, an open triangle medium coupling, and a filled triangle
  high coupling. Two arrowheads denote a mutual dependency; each end retains
  its own coupling level.
- [Dependency summary](dependency-summary.csv): one row for every direct
  inter-directory dependency.
- [Include evidence](include-evidence.csv): source location for every include
  directive used by the analysis.

## Current conclusion

The baseline at `949d7e...` contained one strongly connected component across
all ten source directories, with 55 direct dependency relations: 9 high,
35 medium, and 11 low. Those measurements remain the historical reference.

The current `decoupling` branch has progressed beyond that baseline. Foundation
has been made dependent only on Utility, extracted to
`ariadne-cps/foundation`, validated by standalone Unix/Windows/Coverage CI,
and re-integrated into Ariadne as a pinned submodule. The two small backward
includes from geometry/dynamics into hybrid have also been removed, and Tensor
graphics support has been moved out of algebra.

The remaining work follows the same pattern: remove backward implementation
dependencies, make target boundaries explicit, validate components independently,
and extract repositories only after those boundaries are enforced.

