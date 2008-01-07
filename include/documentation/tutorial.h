/*!
 
\page tutorial Tutorial (Incomplete)


\section creating_functions Creating functions

Currently, user-defined functions cannot be created from Python functions. 
Instead, user-defined functions 

\section creating_sets Creating sets

The most general way of creating subsets of Euclidean space in Ariadne is as a constraint. 

\code

circles = InterpolatedFunction("1-x[0]*x[0]-x[1]*x[1]")
constraint = Constraint(circles)
set = ConstraintSet(constraint)

\endcode

Abstract sets can be intersected and joined
\code

set = intersection(join(set1,set2),set3)

\endcode


Built-in abstract sets RectangularSet and PolyhedralSet are provided for convenience.





*/
