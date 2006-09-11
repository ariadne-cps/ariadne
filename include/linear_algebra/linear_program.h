/***************************************************************************
 *            linear_program.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it Pieter.Collins@cwi.nl
 ****************************************************************************/
/*
 * Based on the linear programming algorithms in PPL-0.8
 *   Copyright (C) 2001-2006 Roberto Bagnara <bagnara@cs.unipr.it>
 */

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
 
/*! \file linear_program.h
 *  \brief Linear programming problems.
 */

#ifndef _ARIADNE_LINEAR_PROGRAM_H
#define _ARIADNE_LINEAR_PROGRAM_H

#include <iosfwd>
#include <map>


#include "../linear_algebra/matrix.h"

namespace Ariadne {
  namespace LinearAlgebra {

    /*! \ingroup LinearAlgebra
     *  \brief Linear programming problems
     *
     *  FIXME: The dyadic version doesn't work properly due to divisions. 
     */
    template<typename R>
    class LinearProgram {
     public:
      /*! \brief The type of denotable real number. */
      typedef R real_type;
      /*! \brief The type of matrix used for the tableau. */
      typedef Matrix<R> matrix_type;
      /*! \brief The type of vector used to represent the constraint values and costs. */
      typedef Vector<R> vector_type;
     
      /*! \brief Destructor. */
      ~LinearProgram();
    
      /*! \brief Default constructor: builds a trivial LP problem. 
       *
       * The trivial LP problem requires to maximize the objective function
       * \f$0\f$ on the zero-dimensional vector space under no constraints
       * at all: the origin of the vector space is the optimal solution.
       */
      LinearProgram();

      /*! \brief
       * Builds an LP problem from the constraint system \f$x\geq0\f$ and \f$Ax\leq b\f$ and the objective function
       * \f$\max c^T x\f$.
       *
       * \param A
       * The matrix on the left-hand side of the system of linear inequalities.
       *
       * \param b
       * The vector describing the right-hand side of the system of linear inequalities.
       * \param c
       * The objective function for the LP problem (optional argument with
       * default value \f$0\f$).
       *
       * \exception std::invalid_argument
       */
      explicit LinearProgram(const Matrix<R>& A, const Vector<R>& b, const Vector<R>& c);

      /*! \brief Builds an LP problem from the _tableau T. 
       *
       * The tableau T is a bordered matrix of the form \f$\begin{array}{c|c} A&b\\\hline -c^T&v \end{array}\f$
       * with \f$ b\geq0 \f$. 
       * The optimization problem is 
       * \f[ \max c^T x +v \quad \textrm{s.t.} \quad Ax+s=b,\ x,s\geq 0 \f]
       * The variables \f$s\f$ are the <em>slack variables</em>, and a feasible 
       * solution is given \f[ \max c^T x +v \quad \textrm{s.t.} \quad Ax+s=b,\ x,s\geq 0 \f]
       */
      explicit LinearProgram(const Matrix<R>& T);
    
      /*! \brief Ordinary copy-constructor. */
      LinearProgram(const LinearProgram& y);
    
      /*! \brief Assignment operator. */
      LinearProgram& operator=(const LinearProgram& y);
    
      /*! \brief Returns the number of constraints. */
      size_type number_of_constraints() const;
    
      /*! \brief Returns the number of non-basic variables. */
      size_type number_of_free_variables() const;
    
      /*! \brief Returns the number of variables and constraints. */
      size_type number_of_variables() const;
    
      /*! \brief The current working tableau. */
      const Matrix<R>& tableau() const;

      /*! \brief The current working tableau. */
      Matrix<R>& tableau();
   
      /*! \brief Checks satisfiability.
       *
       * \return
       * <CODE>true</CODE> if and only if the LP problem is satisfiable.
       */
      bool is_satisfiable() const;
    
      /*! \brief Optimizes the current LP problem using the primal simplex algorithm.
       */
      void solve() const;
    
      /*! \brief Evaluates the objective function at \p p. */
      real_type objective_function(const Vector<R>& p) const;
    
      /*! \brief Returns a feasible point, if it exists.
       *
       * \exception std::domain_error
       * Thrown if the LP problem is not satisfiable.
       */
      Vector<R> feasible_point() const;
    
      /*! \brief Returns an optimal point, if it exists.
       *
       * \exception std::domain_error
       * Thrown if \p *this doesn't not have an optimizing point, i.e.,
       * if the LP problem is unbounded or not satisfiable.
       */
      Vector<R> optimizing_point() const;
    
      /*! \brief Computes and returns the optimal value, if it exists.
       *
       * \exception std::domain_error
       * Thrown if \p *this doesn't not have an optimizing point, i.e.,
       * if the LP problem is unbounded or not satisfiable.
       */
      real_type optimal_value() const;
      
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;

     private:
      // An enumerated type describing the internal status of the LP problem.
      enum Status {
        //! The LP problem has not been solved yet.
        UNSOLVED,
        //! The LP problem is unsatisfiable.
        UNSATISFIABLE,
        //! The LP problem is satisfiable; a feasible solution has been computed.
        SATISFIABLE,
        //! The LP problem is unbounded; a feasible solution has been computed.
        UNBOUNDED,
        //! The LP problem is optimized; an optimal solution has been computed.
        OPTIMIZED,
        /*! \brief
          The feasible region of the LP problem has been changed by adding
          new constraints; a feasible solution for the old constraints has
          been computed.
        */
        PARTIALLY_SATISFIABLE
      };

     private:
      //  The matrix encoding the current feasible region in _tableau form.
      mutable Matrix<R> _tableau;
      // The current basic variables.
      mutable std::vector<size_type> _variable_indices;
       
      // Current status of the problem
      mutable Status _status;
    
     private:
      /*! \brief
        Optimizes the current LP problem using the second phase of the
        primal simplex algorithm.
      */
      void solve_tableau() const;
       
      /* 
       * Checks for optimality and, if it does not hold, computes the column
       * index of the variable entering the base of the LP problem.
       * Implemented with anti-cycling rule.
       *
       * Returns the column index of the variable that enters the base. If no such
       * variable exists, optimality was achieved and n is retuned.
       */
      size_type compute_entering_variable_index() const;
    
      /*
       * Computes the row index of the variable exiting the base
       * of the LP problem. Implemented with anti-cycling rules.
       *
       * Returns the row index of the variable exiting the base.
       */
      size_type compute_exiting_base_index(size_type entering_variable_index) const;
    
   
      /*
       * Swaps two variables in base during the simplex algorithm,
       * performing the needed linear combinations.
       *
       * entering_var_index is the index of the variable entering the base.
       * exiting_base_index is the index of the row exiting the base.
      */
      void swap_base(const size_type entering_var_index,
                     const size_type exiting_base_index) const;
   
      /* \brief
       * Checks for optimality and, if it does not hold, computes the column
       * index of the variable entering the base of the LP problem.
       *
       * \return
       * The column index of the variable that enters the base. If no such
       * variable exists, optimality was achieved and <CODE>0</CODE> is retuned.
       *
       * To compute the entering_index, the steepest edge algorithm chooses
       * the index `j' such that \f$\frac{d_{j}}{\|\Delta x^{j} \|}\f$ is the
       * largest in absolute value, where
       * \f[
       *   \|\Delta x^{j} \|
       *     = \left(
       *         1+\sum_{i=1}^{m} \alpha_{ij}^2
       *       \right)^{\frac{1}{2}}.
       * \f]
       * Recall that, due to the Integer implementation of the algorithm, our
       * _tableau doesn't contain the ``real'' \f$\alpha\f$ values, but these
       * can be computed dividing the value of the cofficient by the value of
       * the variable in base. Obviously the result may not be an Integer, so
       * we will proceed in another way: the following code will compute the
       * lcm of all the variables in base to get the good ``weight'' of each
       * Coefficient of the _tableau.
       */
      size_type steepest_edge() const;

    };

    template<typename R>
    inline
    std::ostream& 
    operator<<(std::ostream& os, const LinearProgram<R>& lp) {
      return lp.write(os);
    }
    
    
    
    
    template<typename R>
    inline 
    const Matrix<R>&
    LinearProgram<R>::tableau() const 
    {
      return _tableau;
    }
    
    template<typename R>
    inline 
    Matrix<R>&
    LinearProgram<R>::tableau() 
    {
      return _tableau;
    }
    
    template<typename R>
    inline 
    size_type
    LinearProgram<R>::number_of_constraints() const 
    {
      return _tableau.size1()-1;
    }
    
    template<typename R>
    inline 
    size_type
    LinearProgram<R>::number_of_free_variables() const 
    {
      return _tableau.size2()-1;
    }
    
    template<typename R>
    inline 
    size_type
    LinearProgram<R>::number_of_variables() const 
    {
      return _tableau.size1()+_tableau.size2()-2;
    }
    

    
    template<typename R>
    inline 
    void
    LinearProgram<R>::solve() const 
    {
      this->solve_tableau();
    }
    
    
    template<typename R>
    inline
    bool 
    LinearProgram<R>::is_satisfiable() const 
    {
      assert(false);
      return true;
    }
    
/*
    template<typename R>
    inline 
    Vector<R>
    LinearProgram<R>::feasible_point() const 
    {
      if (is_satisfiable()) {
        return last_point;
      }
      throw std::domain_error("PPL::LinearProgram::feasible_point():\n"
                              "*this is not satisfiable.");
    }
*/    
    
  }
}


#endif /* _ARIADNE_LINEAR_PROGRAM_H */
