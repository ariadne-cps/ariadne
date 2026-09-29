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

All ten directories form one strongly connected component. Moving each
directory to a repository now would preserve the cycles while adding version
coordination overhead. The first stage is therefore to make dependencies
explicit in CMake, add component-level tests, and remove or invert the small
backward dependencies that close the cycles. Repository extraction follows
only after those boundaries are verified in the monorepo.

The baseline contains 55 direct dependency relations: 9 high, 35 medium, and
11 low. The most extensive relation is `algebra → function`: `function`
contains 142 includes of `algebra`, spread over 41 files and 31 headers.

