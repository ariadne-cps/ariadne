/***************************************************************************
 *            function_model_archetype.h
 *
 *  Copyright  2008  Pieter Collins
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
01234567890123456789012345678901234567890123456789012345678901234567890123456789
/*! \file function_model_archetype.h
 *  \brief Archetype for function models.
 *
 * An <em>archetype</em> is a class which declares precisely the methods for 
 * a particular <em>concept</em>. The methods do not need to be implemented.
 * An archetype can be use for documenting the concept, or for ensuring 
 * correctness of a concept checking class or algorithm at compile-time. 
 */
 
#ifndef ARIADNE_FUNCTION_MODEL_ARCHETYPE_H
#define ARIADNE_FUNCTION_MODEL_ARCHETYPE_H

#include "boost/concept_check.hpp"
//#include "boost/concept_archetype.hpp"
//#include "boost/concept/assert.hpp"
//#include "boost/concept/usage.hpp"

namespace Ariadne {

    template<class R> class Interval;
    template<class X> class Vector;
    template<class X> class Matrix;
    template<class R> class FunctionInterface;

    /*!\ingroup Function
     * \brief Archetype class for function models.
     * 
     * A function model is an approximation to a function on a 
     * restricted domain that can be manipulated numerically. 
     * Function model classes are the main workhorses of the 
     * Ariadne numerical engine. They contain methods for
     * evaluation, differentiation, solving implicit function
     * systems and integration.
     */
    template<class Real>
    class FunctionModelArchetype {
     public:
      /*! \brief The real number type. */
      typedef typename R::real_type real_type;
     private:
      typedef Interval<R> I;
     public:
      /// Basic constructor for a function model.
      ///
      /// \parameter function  The function to be approximated
      /// \parameter domain  The domain on which the model approximates the function.
      /// \parameter  accuracy  The accuracy of the approximation, given by the 
      ///    the order of an approximating series.
      /// \parameter smoothness  The maximum order to which the model gives bounds 
      ///    on the derivatives (optional; may be limited by the implementation)
      FunctionModelArchetype(Box<Real> domain, 
                             FunctionInterface<Real> function,
                             Integer accuracy,
                             Integer smoothness);
                             
      /// Basic constructor for a function model. (Optional)
      ///
      /// \parameter  accuracy  The accuracy of the approximation, given by 
      ///    the maximum error in the supremum norm.
      FunctionModelArchetype(Box<Real> domain, 
                             FunctionInterface<Real> function,
                             Real accuracy,
                             Integer smoothness);
                             
      /// The domain on which the model yields a valid approximation.
      Box<Real> domain() const;

      /// An over-approximation to the range of the model. 
      Box<Real> codomain() const;

      /// The maximal derivative which may be computed by the model.
      Integer smoothness() const;

      /// Evaluate the model on the interval vector \a X.
      /// 
      /// The result is guarenteed to contain \f$f(x)\f$ for all 
      /// \f$x\f$ in the box \f$X\f$, where \f$f\f$ is the function
      /// approximated by the model.
      Vector<Interval> evaluate(Vector<Interval> X) const;

      /// Compute the jacobian on the interval vector \a X.
      /// 
      /// \precondition  The smoothness is at least 1.
      ///
      /// The result is guarenteed to \f$Df(x)\f$ for all 
      /// \f$x\f$ in the box \f$X\f$.
      Matrix<Interval> jacobian(Vector<Interval> X) const;

      /// Compute the function value and all its derivatives \f$D^{\alpha}f\f$
      /// for multi-indices \f$\alpha\f$ with \f$|\alpha|\leq s\f$
      /// on the interval vector \a X.
      /// 
      /// \precondition  \a s is at most the smoothness of the model.
      Differential<Interval> expansion(Vector<Interval> X, Integer s) const;

      
      /// Add a constant vector.
      FunctionModel<Real>& operator+=(Vector<Interval> c);
      /// Multiply by a scalar.
      FunctionModel<Real>& operator*=(Interval c);

      /// Multiply by a matrix.
      friend FunctionModel<Real> operator*(Matrix<Interval> A,
                                           FunctionModel<Real> f);

      /// A function model representing the same function as \a f, but 
      /// only valid on the restricted domain \a d.
      ///
      /// \prerequisite \a D is a subset of the domain of \a f.
      friend FunctionModel<Real> restrict(FunctionModel<Real> f,
                                          Box<Real> D);
                                         
      /// Embed the function model \a f in a domain \a D of larger dimension
      /// by means of a projection map \a p. The result is a model for the
      /// composition \f$f\circ p\f$.
      ///
      /// \remark The syntax of this method is not entirely fixed. The
      /// method may be superceded by use of a compose.
      friend FunctionModel<Real> embed(FunctionModel<Real> f,
                                       Projection p,
                                       Box<Real> D);
                                         
      /// Partially evaluate the function model by setting the value of
      /// \f$x_j\f$ equal to \f$w_i\f$ if \f$j=p(i)\f$. 
      ///
      /// If the domain of \f$f\f$ has dimension \f$n\f$ and \f$w\f$ has 
      /// dimension \f$m\f$, then the domain of the result has dimension
      /// \f$n-m\f$.
      ///
      /// \remark The syntax of this method is not entirely fixed. The 
      /// map \a p may instead be an embedding (the adjoint of a projection). 
      /// The method may be superceded by use of a compose.
      friend FunctionModel<Real> evaluate(FunctionModel<Real> f,
                                          Vector<Interval> w,
                                          Projection p);
                                         
      /// Construct the function model \f$(f,g)\f$.
      ///
      /// \prerequisite \a f and \a g have the same domain.
      friend Vector<Interval> join(FunctionModel<Real> f,
                                   FunctionModel<Real> g);

      /// Computes a model for the composition \f$f\circ g\f$ over the domain
      /// of \f$g\f$.
      ///
      /// \prerequisite  The codomain of \a g is a subset of the domain of \a f.
      ///
      /// The composition must contain any function \a h which could 
      /// be the result of any function represented by \a f with any function 
      /// represented by \a g.
      ///
      friend FunctionModel<Real> compose(FunctionModel<Real> f,
                                         FunctionModel<Real> g);
                                         
      /// Computes a model for the composition \f$f\circ g\f$ over the domain
      /// of \f$g\f$. (Optional)
      ///
      /// The composition must contain any function \a h which could 
      /// be the result of composing \a f with any function 
      /// represented by \a g.
      friend FunctionModel<Real> compose(FunctionInterface<Real> f,
                                         FunctionModel<Real> g);
                                         
      /// Computes a model for the composition \f$f\circ g\f$ over the domain
      /// of \f$D\f$. (Optional)
      ///
      /// \prerequisite The image of \a D under \a g must be a subset of the 
      /// domain of \a f.      
      ///
      /// \remark The status and syntax of this method is still unclear. 
      /// It may be preferable to (partially) replace this method with
      /// specialist methods to perform embeddings and partial evaluations.
      friend FunctionModel<Real> compose(FunctionModel<Real> f,
                                         FunctionInterface<Real> g,
                                         Box<R> D);
                                         
      /// Solve the equation \f$f(x)=0\f$.
      friend Vector<Interval> solve(FunctionModel<Real> f);
       
      /// Computes a model for the function \f$h\f$ defined by \f$f(x,h(x))=0\f$.
      ///
      /// The equation should have a unique solution for all \f$x\f$ in the
      /// projection of the domain of \f$f\f$ on the first coordinates.
      ///
      /// If a solution exists but is not unique, then the result is 
      /// implementation-dependent. If a solution only exists on a subset of
      /// the domain, then the implementation may return a solution on a 
      /// reduced domain.
      ///
      /// The function model \a f must represent a 
      /// function \f$f:\R^m\rightarrow \R^n\f$ with \f$m\geq n\f$.
      /// The result is a model for \f$h:\R^{m-n}\rightarrow \R^n\f$.
      friend FunctionModel<Real> implicit(FunctionModel<Real> f);
                                         
      /// Construct the function model representing
      /// \f$\partial f/\partial x_i\f$.
      ///
      /// \prerequisite The smoothness of \a f is at least 1.
      friend FunctionModel<Real> derivative(FunctionModel<Real> f,
                                            Integer i);

      /// Construct the function model representing a function \a g such that
      /// \f$f = \partial g/\partial x_i\f$.
      friend FunctionModel<Real> derivative(FunctionModel<Real> f,
                                            Integer i);

      /// Computes a model for the flow of the vector field 
      /// defined by \a f for initial conditions in the box \a X
      /// for times between 0 and \a t, assuming the solution remains within
      /// the domain \a D.
      ///
      /// (Optionally the time may be an interval.)
      /// (Optionally the domain may be taken to be the domain of \a f.
      /// 
      /// The function model \a f must represent a 
      /// function \f$f:\R^n\rightarrow \R^n\f$. The result is a function
      /// \f$\phi:\R^{n+1}\rightarrow\R^{n}\f$ with arguments \f$\phi(x,t)\f$
      /// with \f$x\in\R^n\f$ and \f$t\in\R\f$.
      friend FunctionModel<Real> flow(FunctionModel<Real> f,
                                      Box<Real> X,
                                      Real t,
                                      Box<R> D);
                                         
      

    };
  
} // namespace Ariadne

#endif /* ARIADNE_FUNCTION_MODEL_ARCHETYPE_H */
