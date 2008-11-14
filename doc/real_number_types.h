/***************************************************************************
 *            real_number_types.h
 *
 *  Copyright  2004-7  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */


/*! 

\file real_number_types.h
\brief Documentation on real number types



\page realnumbertypes Real Number Types (Outdated)

\section Introduction

%Ariadne currently supports three different real number types, \a Float64, \a FloatMP and \a Rational.
The %Float64 type is a \em finite-precision floating-point type, the %FloatMP type is a \em multiple-precision floating-point type
and %Rational is an \em exact real number type. Future realeases will also support an \em arbitrary-precision
real number type \a Real. Additionally, an \a Integer type is supported.

\section floatingpoint Floating-Point Numbers

Real numbers are traditionally described by floating-point types \a float and \a double.
The set of elements which can be represented by these typesis <em>finite</em>.
Arithmetic operations can only be performed approximately to maximal precision determined by the data type.
However, these types have the advantage of requiring a known amount of memory, which means they can be statically
allocated, and having hardware-supported arithmetical approximations. This makes them especially suitable in
situations where execution speed and memory usage are more important than a knowledge of the computational accuracy.

\section finiteprecision Finite-Precision Real Numbers

For a package such as %Ariadne, in which we are concerned with keeping track of numerical errors,
we cannot use these types directly. Instead, we provide four different arithmetic operations for
the same type. We can specify that the result can be rounded up or down, or that the result
should be rounded to the nearest representable value. However, the default option is to
<em>return an interval</em> containing the exact value of the operation. 
To summarize, we have the following 
operations.
 - Exact arithmetic where possible:
   <dl><dt><dd><tt>Float neg_exact(Float);</tt>
       <dt><dd><tt>Float neg(Float);</tt>
       <dt><dd><tt>Float operator-(Float);</tt></dl>

 - Approximate and rounded arithmetic:
   <dl><dt><dd><tt>Float add_approx(Float);</tt>
       <dt><dd><tt>Float add_down(Float);</tt>
       <dt><dd><tt>Float add_up(Float);</tt></dl>

 - Interval arithmetic:
   <dl><dt><dd><tt>Interval<Float> operator+(Interval<Float>,Interval<Float>);</tt></dl>

 - Automatic conversion to interval arithmetic:
   <dl><dt><dd><tt>Interval<Float> operator+(Float,Float);</tt>
       <dt><dd><tt>Interval<Float> operator+(Interval<Float>,Float);</tt>
       <dt><dd><tt>Interval<Float> operator+(Float,Interval<Float>);</tt></dl>
 .
For types based on floating point numbers, the arithmetical operators should follow IEEE standards.

Algebraic and transcendental functions are provided with both rounded and interval versions.
 - Approximate and rounded algebraic transcendental functions:
   <dl><dt><dd><tt>Float exp_approx(Float);</tt>
       <dt><dd><tt>Float exp_down(Float);</tt>
       <dt><dd><tt>Float exp_up(Float);</tt></dl>

 - Interval algebraic and transcendental functions:
   <dl><dt><dd><tt>Interval<Float> exp(Interval<Float>);</tt>
       <dt><dd><tt>Interval<Float> exp(Float);</tt></dl>

The precision of the algebraic and transcendental functions is not guarenteed,
except that the functions must give values in the correct range.

The upper and lower rounded versions satisfy the mathematical postcondition 
\f$ \underline{f}(x) \leq f(x) \leq \overline{f}(x) \f$
and the implementation postcondition 
<tt>f_down(x) <= f_approx(x) <= f_up(x)</tt> .

Note that no other conditions are required of the implementation.
Hence a valid (but useless) implementation of \f$ \sin(x) \f$ is
<dl><dt><dd><tt>Float sin_approx(Float x) { return 0.0; }</tt>
    <dt><dd><tt>Float sin_down(Float x) { return -1.0; }</tt>
    <dt><dd><tt>Float sin_up(Float x) { return 1.0; }</tt>
</dl>


\section multipleprecision Multiple-precision types.

If we have a problem for which a fixed-precision type is not sufficient to obtain an accurate answer,
we can switch to a \em multiple-precision type. The semantics of arithmetic on multiple-precision types 
is mostly the same as that of a fixed-precision type; arithmetic is approximate, and the default is to
return an interval. However, a multiple-precision type has a 
<dl><dt><dd><tt>MPReal::set_precision(unsigned int)</tt></dl>
method, which sets the precision to which the type can store its result.
Using a higher precision yields a more accurate answer. Further, the result of an
arithmetic operation is guarenteed to converge as the precision is increased.

For efficiency, all elements of an array of a multiple precision type 

\section arbitraryprecision Arbitrary-precision types.

An arbitrary-precision type stores a number in a form so that it can be
recovered to any desired precision. This is typically acheived by expressing
the number as a formula in terms of other arbitrary-precision or exact number
types. However, the high computational overhead of such numbers makes them
impractical for describing higher-order types, such as matrices or sets.

Arbitrary-precision numbers may be used to store constants occurring in a 
system definition.

\section exactarithmetic Exact arithmetic types.

Elements of a countable set of numbers, such as integer, dyadic, rational or algebraic numbers,
may be stored using a finite amount of data, though the size of the data depends on the 
element used. Hence in an array of an exact arithmetic type, each element needs to be
dynamically allocated, which is expensive in terms of spacial overhead. 
Since these sets are typically closed under arithmetical operations, 
arithmetic can be performed exactly for these types. Hence these types are
appropriate where time and space overhead are not at a premium.

Since these types require arbitrarily
large amounts of memory which is typically dynamically allocated, and hardware support for arithmetic does not
exist, they are typically less efficient than the fixed-precision and multiple-precision types.

In %Ariadne, we support arithmetic on the Rational number type, but no algebraic or transcendental functions.
This type is primarily useful for testing.

Exact arithmetic operators:
<dl><dt><dd><tt>Rational operator-(Rational);</tt>
    <dt><dd><tt>Rational operator+(Rational,Rational);</tt>
</dl>


Unlike a finite-precision type, arbitrary-precision types are intended for use when precise error specifications are required.
Hence, whenever a function returns an approximation, the suffix \c _approx is \em always added to the function name.
This ensures that the user is always aware of the use of approximations.

The Dyadic type does not in general support exact division.
The exception is that division by a power of \f$2\f$ yields an exact answer.
Generic code intended for use with all exact arithmetic types should \em not use division, except by a power of two.



\section interval Interval Arithmetic

Interval arithmetic is commonly used in finite-precision computations to obtain guarenteed bounds for the result.
However, it can also be used in arbtitrary-precision computations to compute the image of basic sets.
The semantics is as follows:

\subsection finite_precision_interval Finite-precision interval arithmetic

Finite precision functions are of the form
<dl><dt><dd><tt>Interval<Float> f(Interval<Float>);</tt></dl>
and satisfy the mathematical postcondition \f$ \forall x\in I,\ f(x)\in\f$\c f(I),
and the implementation postcondition <tt>I.contains(x)</tt> implies <tt>f(I).contains(f_approx(x))</tt>.

\subsection multiple_precision_interval Multiple-precision interval arithmetic

Multiple-precision interval functions follow the same conditions as fixed-precision
interval functions, together with the convergence criterion
that the length of <tt>f(I)</tt> approaches \f$0\f$ as the length of <tt>I</tt> approaches \f$0\f$.

*/
