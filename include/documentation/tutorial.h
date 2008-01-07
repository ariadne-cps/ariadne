/*!
 
\page tutorial Tutorial 

In this page we present a quick-start guide to %Ariadne, with an overview of the core types and structures, a guide on how to build and analyse dynamic systems.

\section tutorial_numeric Numerical Types

%Ariadne supports the numerical types \c Integer, \c Rational, \c Float and \c Interval. 

The \c Integer and \c Rational types support exact arithmetic operations. Unlike C++ and the Python 2, integer division is performed exactly and returns a rational.
The operations \c quot(Integer,Integer) and \c rem(Integer,Integer) can be used to perform integer division. 

The \c Float class represents floating-point numbers. Since most arithmetic operations on floating-point numbers can only be performed approximately, %Ariadne uses <em>interval arithmetic</em> to represent the results of floating-point computations. The result of any floating-point computation is represented as an interval \f$[l,u]\f$ enclosing the exact value of the result. In this way, round-off errors can be propagated automatically. 

%Ariadne numerical types can be constructed by conversion from built-in Python types or from string literals. Note that a string literal representing a \c Float must be exacly representable on the machine. Hence <c>Float(3.3)</c> and <c>Float("3.25")</c> are both valid (the former has a value of \f$3.2999999999999998224\ldots\f$) but <c>Float("3.3")</c> is an error. Note that <c>Interval(3.3)</c> yields the singleton interval \f$[3.2999999999999998224,3.2999999999999998224]\f$ (the constant is first interpreted by Python to give a Python \c float, whereas <c>Interval("3.3")</c> yields the interval \f$[3.2999999999999998224,3.3000000000000002665]\f$ enclosing \f$3.3\f$.

Comparison tests on \c Interval use the idea that an interval represents a single number with an unknown value. Hence the result is of type \c tribool, which can take values { \c True, \c False, \c Indeterminate }.  Hence a test \f$[l_1,u_1]\leq [l_2,u_2]\f$ returns \c True if \f$u_1\leq u_2\f$, since in this case \f$x_1\leq x_2\f$ whenever \f$x_1\in[l_1,u_2]\f$ and \f$x_2\in[l_2,u_2]\f$, \c False if \f$l_1>u_2\f$, since in this case we know \f$x_1>x_2\f$, and \c Indeterminate otherwise, since in this case we can find \f$x_1,x_2\f$ making the result either true or false. In the case of equality, the comparison \f$[l_1,u_1]\f$==\f$[l_2,u_2]\f$ only returns \c True if both intervals are singletons, since otherwise we can find values making the result either true of false.

To obtain the lower and upper bounds of an interval, use \c ivl.lower() and \c ivl.upper(). 
To obtain the midpoint and radius, use \c ivl.midpoint() and \c ivl.radius().
Alternatives \c midpoint(ivl) and \c radius(ivl) are also provided. 
Note that \c midpoint and \c radius return approximations to the true midpoint and radius of the interval. If \f$m\f$ and \f$r\f$ are the returned midpoint and radius of the interval \f$[l,u]\f$, the using exact arithmetic, we guarentee \f$m-r\leq l\f$ and \f$m+r\geq u\f$

To test if an interval contains a point or another interval, use \c encloses(Interval,Float) or \c encloses(Interval,Interval). 
The test \c refines(Interval,Interval) can also be used. 

In the C++ kernel, the \c Float class is a generic (templated) class parameterised by the base number type \c T, which can either be double precision (\c Float64) or multiple precision (\c FloatMP). The Interval class is also a generic class parameterised by the base number type \c R, which can be an \c Integer, \c Rational or \c Float<T>. 
In the Python interface, the C++ type \c FloatPy used for floating-point numbers must be chosen at compile-time, and the interval type is then \c Interval<FloatPy>.

\sa  Numeric sub-module, Ariadne::Numeric::Integer, Ariadne::Numeric::Rational, Ariadne::Numeric::Float, Ariadne::Numeric::Interval

<br><br>



\section tutorial_linear_algebra LinearAlgebra

%Ariadne supports three basic types for linear algebra, a \c Vector (column vector) class, a \c Covector (row vector) class and a \c Matrix class, with the standard operations of addition, subtraction, scalar multiplication and matrix multiplication. 
The standard classes \c Vector class has values of type \c Float, but \c RationalVector and \c IntervalVector versions are also supported; similiarly for \c Covector and \c Matrix. 
As in the Numeric module, operations on floating-point types return interval types.

The inverse of a floating-point matrix can be computed using interval arithmetic using the function <c>inverse(Matrix)</c>. 
Additional implementations of this functionality may also be provided. 
Systems of linear equations can be solved using <c>solve(Matrix)</c>. 
Approximate matrix factorisations can be computed using <c>(P,L,U)=plu_approx(A)</c> and <c>(Q,R)=qr_approx(A)</c>.

Vectors, covectors and matrices can be constructed by using Python list objects. Note that a string literal representing a Float must be possible to represent exactly. For example
\code
v=Vector( [1,2.2,Float(3.3),"4.5"] )
v=IntervalVector( [1,2.2,Float(3.3),"4.4",Interval(5.5,5.6) ] )
A=Matrix( [ [1,2.2,Float(3.3)], ["4.5",Float("5.75"),6.875] ] )
\endcode



\subsection tutorial_linear_algebra_slicing Slicing (Future addition)

In the future, it will be possible to create \c VectorSlice and \c MatrixSlice objects using Python's slicing operatior. As in Matlab and other C++ libraries, and unlike Python sequences, these operations return a "view" into the structure, and not a copy, so can be used for in-place modification.

The semantics of this operator have not been determined, since Python has different conventions from C++ libraries such as the STL \c valarray class, glas and uBlas \c vector classes and Matlab classes. For example, if \c v is the vector <c>[0,1,2,3,4,5,6,7,8,9]</c>, then different possibilities are

 - Python: <c> v[start:finish:stride] </c>  e.g. <c>v[1:6:2]</c> = <c>v[1:7:2]</c> = <c>[1,3,5]</c>; 

 - C++ STL \c valarray : <c> v[slice(start,size,stride)]</c> e.g.  <c> v[slice(1,3,2)]</c>\c  = <c>[3,5,7]</c>; 

 - C++ uBlas: \c <c> project(v,slice(start,stride,size)) </c>  e.g.  <c>project(v,slice(1,2,3))</c> = <c>[1,3,5]</c>; 

 - C++ glas: \c <c> slice(start,finish,stride)) </c>  e.g.  <c>v[slice(1,6,2)])</c> = <c>v[slice(1,7,2)])</c> = <c>[1,3,5]</c>; 

 - Matlab: <c> v(first:stride:last) </c> with one-based indexing; e.g. <c>v(2:2:6)</c> = <c>v(2:2;7)</c> = <c>[1,3,5]</c>;  

The main design decision is whether to use a single "Ariadne" slicing semantics which is used for any and all interfaces, or whether to use the natural semantics for the given environment. For the C++ API, we could provide STL \c valarray semantics <c> v[slice(start,size,stride=1)] </c> and Python semantics <c> v[range(start,finish,stride=1)] </c>.

\sa  \ref LinearAlgebra sub-module, Ariadne::LinearAlgebra::Vector, Ariadne::LinearAlgebra::Covector, Ariadne::LinearAlgebra::Matrix

<br><br>


  

\section tutorial_function Functions

Currently %Ariadne provides a \c FunctionInterface class which specifies functions \f$f:\R^n\rightarrow\R^m\f$, i.e. functions of a single vector-valued variable. 

There are currently two ways of creating user-defined functions, either via the Python interface or in C++.

\subsection tutorial_function_python Python functions

An %Ariadne function can be constructed directly from a Python function. 
\code
def henon_function(x):
    return [ 1.5-x[0]*x[0]+0.375*x[1], x[0] ]
\endcode
Note that there must be a single argument of sequence type (i.e. supporting subscripting), and that the result must be a Python \c list object.

To prepare the Python function for use in %Ariadne, some meta-information is required.
\code
henon_function.result_size = 2
henon_function.argument_size = 2
\endcode

Finally, the function can be exported to the C++ wrapper.
\code
henon = AriadneFunction(henon_function)
\endcode
The newly-created function can be used exactly as a native C++ function, and supports automatic differentiation.
\code
x=Vector([1,0])
print henon.evaluate(x)
print henon.jacobian(x)
print henon.derivative(x,3)
\endcode

An %Ariadne function can actually be constructed from any Python object with methods \c "result_size", \c "argument_size" and the special method \c "__call__". For example, to create a parameterised version of the Henon map, use:
\code
class HenonFunction:
    def __init__(self,p):
        self.p=p
        self.result_size=2
        self.argument_size=2
    def __call__(self,x):
        return [ p[0]-x[0]*x[0]+p[1]*x[1], x[0] ]

henon = AriadneFunction( HenonFunction( [1.5,0.375] ) )
\endcode   

Unfortunately, since all calls to the function require the Python interpreter, function evaluations will be much less efficient than a native C++

\subsection tutorial_function_python_cpp C++ functions

To construct a function in C++, a template definition of the function must be given. To facilitate the instantiation of the function with multiple result types, the result of the evaluation is given by a non-constant reference paramter.

\code
template<class R, class A>
void henon_function(R& r, const A& x) {
  r[0] = 1.5-x[0]*x[0]+0.375*x[1];
  r[1] = x[0];
}
\endcode
A function object with a method 
\code 
template<class R, class A> void operator() (R&, const A&)
\endcode
may be used instead.

The convenience macro \c ARIADNE_BUILD_FUNCTION is defined in the header file \c "macros/build_function.h", and can be used to create a function given some meta-information. For example
\code
int henon_result_size = 2;
int henon_argument_size = 2;
int henon_smoothness = MAXIMUM_SMOOTHNESS; 

// Maybe "smoothness" should not be part of the interface.
ARIADNE_BUILD_FUNCTION("Henon",henon_function,henon_result_size,henon_argument_size,henon_smoothness);
ARIADNE_BUILD_FUNCTION("Henon",henon_function,henon_result_size,henon_argument_size);
\endcode

\subsection tutorial_function_interpreted Interpreted functions (Deprecated)

%Ariadne has its own \c InterpretedFunction class, which can be created from a string literal in a Modelica-like syntax. However, the use of this class is deprecated now the \c AriadneFunction class is available from the Python interface.

\subsection tutorial_function_linking Automatic Linking for C++ functions (Possible future addition)

It would be nice to have a script in which a function could be specified in C++ (or some other function description language) and directly linked into the Python library. Ideally it would be possible to make a user-specific library of C++ functions and alongside this a separate user-specific library of Python bindings. However, for reasons relating to how Python performs symbol table lookup, it may be impossible to import user-defined C++ functions into Python without re-linking the Python interface. Help on this issue would be appreciated.

\subsection tutorial_function_multivariable Multivariable Functions (Future addition)

In the future, it may be possible to create functions with multiple vector-valued arguments, \f$f:\R^{n_1}\times\cdots\times\R^{n_k}\rightarrow\R^m\f$. This is useful in specifying functions with different argument types such as parameters, space and time, or input and output signals.
For example, a parameterised function could be created (in C++) as follows:
\code
template<class R, class A, class P>
void parameterised_henon_function(R& r, const A& x, const P& p) {
  r[0] = p[0]-x[0]*x[0]+p[1]*x[1];
  r[1] = x[0];
}

ARIADNE_BUILD_BINARY_FUNCTION("Henon", parameterised_henon_function, henon_result_size, henon_argument_size, henon_number_of_parameters)

\endcode

\sa \ref Function, Ariadne::Function::FunctionInterface

<br><br>





\section tutorial_geometry Sets and geometric operations

\subsection tutorial_geometry_box Boxes

The most fundamental set type in %Ariadne is \c Box, which is a coordinate-aligned orthogonal subset of Euclidean space.
The easiest way to construct a box is as a list of intervals:
\code
box = Box( [ [0,2], [2,3] ] )
box = Box( Interval(0,1)]*d )

a=Float(2); b=Interval(2.5,3.3)
box = Box( [ [0,2], Interval(a,3), b ] )
\endcode

\subsection tutorial_geometry_abstract Sets defined by functions

To specify systems, it is possible to define sets either as the image or preimage of a box using the \c ConstraintSet and \c ImageSet classes:
\code
constraint_set = ConstraintSet(function, box)
image_set = ImageSet(box, function) 
\endcode
For example, a disc can be described as the preimage of a halfspace.
\code
def radius_function(x):
     return [ x[0]*x[0]+x[1]*x[1] ]

radius_function.result_size=2
radius_function.argument_size=2
radius_function = AriadneFunction(radius_function)

codomain = Box( [Interval(-inf(),1.0)] ) 

# create the disc -infinity < x^2+y^2 < 1
disc = ConstraintSet(radius_function,codomain) # create 
\endcode 

In general, the \c ConstraintSet class is useful for specifying constraints on the evolution, such as guard sets, and the \c ImageSet class may be more useful for creating initial state sets. In general, a \c ConstraintSet need not be bounded, and even if it is, it may be impossible to find a bounding box. An \c ImageSet may have an empty interior, and is generally not suitable for specifying domains for a flow.

\endcode

Built-in constraint sets \c RectangularSet and \c PolyhedralSet are provided for convenience. For example, the call
\code
r=RectangularSet( [ [-1,1], [0,1] ] )
p=PolyhedralSet(A,b)
\endcode
will create the box \f$[-1,1]\times[0,1]\f$ and the polyhedron \f$ Ax \leq b \f$.

\subsection tutorial_geometry_denotable Sets defined as unions of boxes

When performing concrete computations, it is useful to represent sets as unions of boxes. %Ariadne provides a number of these sets, currently \c BoxListSet, \c GridCellListSet, \c GridMaskSet and \c GridTreeSet.

\note GridTreeSet is currently called \c PartitionTreeSet.

\subsection tutorial_geometry_approximation Approximation of sets

%Ariadne provides a number of approximation concepts relating abstract sets (described by functions) with concrete sets (described by unions of boxes). 
We have the following concepts relating an approximation \f$ A =\bigcup_{j=0}^{k} I_j\f$ with a concrete set \c S.
 - Inner approximation: \f$\overline{A}\subset S^\circ \f$
 - Outer approximation: \f$A^\circ \supset \overline{S} \f$.
 - \f$\epsilon\f$-lower approximation: \f$ N_\epsilon(A) \supset S \f$

An approximation can also be specified by  a \e cover \f$\mathcal{A}=\{ I_1,\ldots,I_k\}\f$, in which not only the union of the cells is important, but also the cells themselves.
 - Lower cover: \f$ \forall I \in \mathcal{A},\ I\cap S\neq\emptyset \f$.

\sa Geometry 

<br><br>




\section tutorial_system Describing systems

\subsection tutorial_system_dynamic Dynamical Systems

To define a simple dynamical system in discrete or continuous time, we need only give the function defining the evolution.
\code 
henon_map = Map(henon_function)
lorenz_system = VectorField(lorenz_function)
\endcode
Special constructors are given for important subclasses.
\code 
# Define the affine map x -> Ax+b.
affine_map = AffineMap(A,b)
\endcode 

\note
In some future version, a dynamic system may be described by a binary function \f$f(x,p)\f$, with the second argument representing parameters. This will allow aspects such as bifurcations and parameter sensitivity to be studied. Time-dependent systems may also be supported explicitly using \f$f(x,t,p)\f$.


\subsection tutorial_system_hybrid Hybrid Automata

See the documentation for the Ariadne::System::SetBasedHybridAutomaton.

\note The current distinction between "SetBased" and "ConstraintBased" hybrid systems may at some point be eliminated.


\sa \ref System, Ariadne::System::Map, Ariadne::System::VectorField, Ariadne::System::SetBasedHybridAutomaton

\note 
The rationale for having separate "Function" and "System" classes is that the "Function" classes are supposed to be elementary building blocks which can be used for a variety of purposes, whereas a "Map" or "VectorField" object is tagged as specifying a dynamic system Further, extra information such as variable and parameter names may be attached to system classes.

<br><br>



\section tutorial_evaluation Evolution and analysis of systems

%Ariadne can perform a number of analyses of system behaviour. The following core functionality is provided:
 - <c>evolve(system, initial_set, time)</c>: Compute the set of states which the system can reach at the specified time from the given initial set.
 - <c>reach(system, initial_set, time)</c>: Compute the set of states which the system can reach at times up to and including the specified time from the given initial set.

 - <c>reach(system, initial_set)</c>: Compute the set of all states which the system can reach from the given initial set.

 - <c>viable(system, bounding_set)</c>: Compute the set of states for which the evolution remains in the bounding set for all times.

 - <c>verify(system, initial_set, bounding_set)</c>: Attempts to decide whether any evolution starting in the initial set remains in the bounding set for all times.

The \c evolve, \c reach and \c viable operators each come in two flavours \c lower_ and \c upper_ , representing different "semantics of evolution". Hence \c lower_evolve computes a lower-approximation to the evolution, and \c upper_evolve computes an outer approximation. An outer approximation to the reachable set is also called \c chain_reach or \c chainreach, since the approximations need not converge to the reachable set itself. 

Currently, only lower and upper approximations to the evolution can be computed. In a future version, we hope to combine these computations to give an approximation to the evolved set with a known accuracy.

To compute the chain reachable set of the Henon map, the following code can be used:
\code
# Construct the henon map (builtin)
henon = HenonMap(1.4,0.3)

# Specify the initial set as the origin
initial_set = RectangularSet( [0,0] )
initial_set = ConstraintSet( IdentityFunction(2), Box([0,0]) )

# Use default evolution parameters
paramrters=EvolutionParameters() 
print parameters 

#Use Kuhn's algorithm for computing the image of sets
applicator=KuhnApplicator(3)

# Construct an object which can compute the evolution
evolver = MapEvolver(parameters,applicator)

# Compute the chain-reachable set
chain_reachable_set = evolver.chain_reach(henon,initial_set)

# Write the set to a postscript output 
eps = EpsPlot()
eps.open( "henon_chainreach.eps", chain_reachable_set.bounding_box() )
eps.set_fill_colour(cyan)
eps.write(chain_reachable_set)
eps.close()
\endcode


\note 
 In a future version of %Ariadne, it may be possible to use a default MapEvolver object, and set the accuracy with a \c set_accuracy method.

\note 
 An alternative way of defining the semantics of evolution could be <c>evolve(system,initial_set,time,semantics)</c>, where \c semantics can take the constant  \c Lower_Semantics or \c Upper_Semantics, if this would be preferable to users.

\note 
The finite-time \c evolve and \c reach routines could be combined into a single routine, with signature <c>evolve(system, initial_set, time_interval)</c> or <c>evolve(system, initial_set, lower_time, upper_time)</c> if this would be preferable to users.

\sa Ariadne::Evaluation::MapEvolver, Ariadne::Evaluation::VectorFieldEvolver, Ariadne::Evaluation::SetBasedHybridEvolver




\subsection tutorial_evaluation_simulation Simulation

Currently there are no facilities for simulation, but an approximation of the trajectory can be computed using the method
 - <c>lower_reach(system, point, time)</c>

\subsection tutorial_evaluation_parameters Evolution parameters and plugins

The \c Evolver objects are function objects which support a wide variety of operations. The main rationale for having these operations as methods of some class, rather than as stand-along functions, is so that the parameters used to govern the accuracy can be stored in the data of the class and thus hidden from the caller. The novice user should be able to use the classes with "default parameters", and merely need to select a single "accuracy level" but the experienced user should have the flexibility to change the parameters describing the evolution method, or even the evolution method itself.

All evolution methods rely on a number of core components:
 - A fundamental data type, called a <em>basic set</em> which is updated at single step of the evolution.
 - An approximation scheme for storing information about evolved sets.
 - Methods for accurately computing the image of a basic set under a map or flow.
 - Methods for checking satisfiability of constraints and for detecting crossings with guard sets (for hybrid systems).
 - Methods for performing geometric operations on basic sets, such as subdividing into smaller pieces, or approximating basic sets which have a complex or numerically ill-conditioned description by simpler sets.
 - Methods for converting between basic sets and finite approximations based on grids, polyhedral partitions or covers 
Additionally, for hybrid systems, we may wish to use different evolution methods, approximations or even different representations of sets in different modes of the system. 

In the current version of %Ariadne, different "plugin" classes can be used to provide different services. Each plugin is packaged to provide related services, which are combined in the evolution. The exact services provided by the plugins are not currently stable.

Core evolution paramters which are common to (most) evolution routines are given in the class \c EvolutionParameters.

\sa Ariadne::Evaluation::EvolutionParameters


*/
