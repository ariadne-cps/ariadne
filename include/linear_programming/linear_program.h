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

#ifndef ARIADNE_LINEAR_PROGRAM_H
#define ARIADNE_LINEAR_PROGRAM_H

#include <iosfwd>
#include <map>


#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

namespace Ariadne {
  namespace LinearProgramming {

    //! \ingroup LinearProgramming
    //! Solve the linear programming problem \f$\max c^Tx\f$ such that \f$Ax=b;\ x\geq0\f$.
    template<class R> R solve(const LinearAlgebra::Matrix<R>& A, const LinearAlgebra::Vector<R>& b, const LinearAlgebra::Vector<R>& c);
    //! Solve the linear programming problem \f$\max c^Tx\f$ such that \f$Ax=b;\ l \leq x\leq u\f$.
    template<class R> R solve(const LinearAlgebra::Matrix<R>& A, 
                              const LinearAlgebra::Vector<R>& b, const LinearAlgebra::Vector<R>& c,
                              const LinearAlgebra::Vector<R>& l, const LinearAlgebra::Vector<R>& u);

    //! Tests feasibility of the primal linear programming problem \f$Ax=b;\ x\geq0\f$. 
    template<class R> tribool feasible(const LinearAlgebra::Matrix<R>& A, const LinearAlgebra::Vector<R>& b);
    //! Tests feasibility of the constrained primal linear programming problem \f$Ax\leq b;\ l\leq x\leq u\f$. 
    template<class R> tribool feasible(const LinearAlgebra::Matrix<R>& A, const LinearAlgebra::Vector<R>& b,
                                       const LinearAlgebra::Vector<R>& l, const LinearAlgebra::Vector<R>& u);
    //! Tests feasibility of the dual feasibility problem \f$A^Ty\leq c\f$. 
    template<class R> tribool dual_feasible(const LinearAlgebra::Matrix<R>& A, const LinearAlgebra::Vector<R>& c);


    enum Comparison { equal, greater, less, greater_or_equal, less_or_equal };
    
    class Variable;
    template<class R> class LinearTerm;
    template<class R> class LinearCombination;
    template<class R1, class R2> class LinearConstraint;
    template<class R> class LinearProgram;
   
    class Variable {
     public:
      Variable(size_type begin, size_type size, int increment=1) 
        : _size(size), _begin(begin), _increment(increment) { }
      size_type operator[] (size_type i) { 
        return this->_begin+this->_increment*i; }
      size_type size() const { return this->_size; }
     private:
      size_type _size;
      size_type _begin;
      int _increment;
    };
   
    
    template<class R>
    class LinearTerm {
      friend class LinearCombination<R>;
     public:
      LinearTerm(const LinearAlgebra::Matrix<R>& m, const Variable& s)
        : _coefficient(1), _matrix(&m), _variable(&s) { }
      LinearTerm(const R& c, const LinearAlgebra::Matrix<R>& m, const Variable& s)
        : _coefficient(c), _matrix(&m), _variable(&s) { }
     private:
      const R _coefficient;
      const LinearAlgebra::Matrix<R>* _matrix;
      const Variable* _variable;
    };
    
    template<class R>
    class LinearCombination {
      template<class R1, class R2> friend class LinearConstraint;
     public:
      LinearCombination(const LinearTerm<R>& lt) 
        : _coefficients(1,lt._coefficient), _matrices(1,lt._matrix), _variables(1,lt._variable) { }
      LinearCombination(const LinearCombination<R>& lc) 
        : _coefficients(lc._coefficients), _matrices(lc._matrices), _variables(lc._variables) { }
      void add(const LinearTerm<R>& lt) {
        this->_coefficients.push_back(lt._coefficient); 
        this->_matrices.push_back(lt._matrix); 
        this->_variables.push_back(lt._variables); }
      void subtract(const LinearTerm<R>& lt) {
        this->_coefficients.push_back(-lt._coefficient); 
        this->_matrices.push_back(lt._matrix); 
        this->_variables.push_back(lt._variables); }
     private:
      std::vector<R> _coefficients;
      std::vector<const LinearAlgebra::Matrix<R>*> _matrices;
      std::vector<const Variable*> _variables;
    };
    
    template<class R1, class R2>
    class LinearConstraint {
      friend class LinearProgram<typename Numeric::traits<R1,R2>::arithmetic_type>;
     public:
      LinearConstraint(const LinearCombination<R1>& lc, const Comparison& cmp, const LinearAlgebra::Vector<R2>& v) 
        : _coefficients(lc._coefficients), _matrices(lc._matrices), _variables(lc._variables),
          _comparison(cmp), _vector(&v) { }
     private:
      std::vector<R1> _coefficients;
      std::vector<const LinearAlgebra::Matrix<R1>*> _matrices;
      std::vector<const Variable*> _variables;
      Comparison _comparison;
      const LinearAlgebra::Vector<R2>* _vector;
    };
     
    
    template<class R> inline
    LinearTerm<R> operator*(const LinearAlgebra::Matrix<R>& m, const Variable& s) {
      return LinearTerm<R>(m,s);
    }
    
    template<class R> inline
    LinearTerm<R> operator*(const R& c, const LinearTerm<R>& lt) {
      return LinearTerm<R>(c*lt._coefficient, lt._matrix, lt._variable);
    }
    
    template<class R> inline
    LinearCombination<R> operator+(const LinearTerm<R>& lt1, const LinearTerm<R>& lt2) {
      return LinearCombination<R>(lt1).add(lt2); }
        
    template<class R> inline
    LinearCombination<R> operator+(const LinearCombination<R>& lc1, const LinearTerm<R>& lt2) {
      return LinearCombination<R>(lc1).add(lt2); }
    
    template<class R> inline
    LinearCombination<R> operator-(const LinearCombination<R>& lc1, const LinearTerm<R>& lt2) {
      return LinearCombination<R>(lc1).subtract(lt2); }
    
    template<class R1, class R2> inline
    LinearConstraint<R1,R2> operator==(const LinearCombination<R1>& lc, const LinearAlgebra::Vector<R2>& v) {
      return LinearConstraint<R1,R2>(lc,equal,v); }
      
    template<class R1, class R2> inline
    LinearConstraint<R1,R2> operator==(const LinearTerm<R1>& lt, const LinearAlgebra::Vector<R2>& v) {
      return LinearConstraint<R1,R2>(LinearCombination<R1>(lt),equal,v); }
      
    template<class R1, class R2> inline
    LinearConstraint<R1,R2> operator<=(const LinearCombination<R1>& lc, const LinearAlgebra::Vector<R2>& v) {
      return LinearConstraint<R1,R2>(lc,less_or_equal,v); }
      
    template<class R1, class R2> inline
    LinearConstraint<R1,R2> operator>=(const LinearCombination<R1>& lc, const LinearAlgebra::Vector<R2>& v) {
      return LinearConstraint<R1,R2>(lc,greater_or_equal,v); }
      
    
      
    /*!\ingroup LinearProgramming
     * \brief Linear programming problems
     *
     * Solve a linear program in the standard form \f$\text{minimize } cx \text{ subject to } Ax=b;\ x\geq0\f$.
     *
     * The problem is stored as a tableau 
     * \f$\left(\begin{array}{c|c} A&b\\\hline c^T&d \end{array}\right)\f$
     * where the constraints are \f$Ax+y = b\f$ and the cost of putting variable
     * \f$x_i\f$ into the basis is given by \f$c_i\f$, and the current value is \f$d\f$.
     *
     * For a two-stage problem, \f$c^T\f$ stores both the actual value function
     * and the value function associated with the constraints.
     *
     * The standard MATLAB constructor is x = linprog(c ;A; b; Aeq; beq; l; u) 
     * for the linear program \f$\text{minimize } c^Tx \text{ subject to } Ax\leq b;\ A_{\mathrm{eq}}x=b_{\mathrm{eq}};\ l\leq b\leq u\f$.
     */
    template<class R>
    class LinearProgram {
     public:
      /*! \brief The type of denotable real number. */
      typedef R real_type;
      /*! \brief The type of matrix used for the tableau. */
      typedef LinearAlgebra::Matrix<R> matrix_type;
      /*! \brief The type of vector used to represent the constraint values and costs. */
      typedef LinearAlgebra::Vector<R> vector_type;
     
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
      explicit LinearProgram(const LinearAlgebra::Matrix<R>& A, const LinearAlgebra::Vector<R>& b, const LinearAlgebra::Vector<R>& c);

      /*! \brief Builds an LP problem from the tableau \a T. 
       *
       * The tableau T is a bordered matrix of the form \f$\begin{array}{c|c} A&b\\\hline -c^T&v \end{array}\f$
       * with \f$ b\geq0 \f$. 
       * The optimization problem is 
       * \f[ \max c^T x +v \quad \textrm{s.t.} \quad Ax+s=b,\ x,s\geq 0 \f]
       * The variables \f$s\f$ are the <em>slack variables</em>, and a feasible 
       * solution is given \f[ \max c^T x +v \quad \textrm{s.t.} \quad Ax+s=b,\ x,s\geq 0 \f]
       */
      explicit LinearProgram(const LinearAlgebra::Matrix<R>& T);
    
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
      const LinearAlgebra::Matrix<R>& tableau() const;

      /*! \brief The current working tableau. */
      LinearAlgebra::Matrix<R>& tableau();
   
      /*! \brief The variable indices. */
      array<int>& variable_indices();

      /*! \brief Checks satisfiability.
       *
       * \return
       * <CODE>true</CODE> if and only if the LP problem is satisfiable.
       */
      tribool is_feasible() const;
    
      /*! \brief Optimizes the current LP problem using the primal simplex algorithm.
       */
      void solve() const;
    
      /*! \brief Evaluates the objective function at \p p. */
      real_type objective_function(const LinearAlgebra::Vector<R>& p) const;
    
      /*! \brief Returns a feasible point, if it exists.
       *
       * \exception std::domain_error
       * Thrown if the LP problem is not satisfiable.
       */
      LinearAlgebra::Vector<R> feasible_point() const;
    
      /*! \brief Returns an optimal point, if it exists.
       *
       * \exception std::domain_error
       * Thrown if \p *this doesn't not have an optimizing point, i.e.,
       * if the LP problem is unbounded or not satisfiable.
       */
      LinearAlgebra::Vector<R> optimizing_point() const;
    
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
      };

     private:
      //  The matrix encoding the current feasible region in _tableau form.
      mutable LinearAlgebra::Matrix<R> _tableau;
      // The current basic variables.
      mutable std::vector<size_type> _variable_indices;
       
      // Current status of the problem
      mutable Status _status;
    
     private:
      /*! \brief
        Optimizes the current LP problem using the second phase of the
        primal simplex algorithm.
      */
      void compute_feasible_point() const;
       
      /*! \brief
        Optimizes the current LP problem using the second phase of the
        primal simplex algorithm.
      */
      void compute_optimizing_point() const;
       
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

    template<class R>
    inline
    std::ostream& 
    operator<<(std::ostream& os, const LinearProgram<R>& lp) {
      return lp.write(os);
    }
    
    
    
    
    template<class R>
    inline 
    const LinearAlgebra::Matrix<R>&
    LinearProgram<R>::tableau() const 
    {
      return _tableau;
    }
    
    template<class R>
    inline 
    LinearAlgebra::Matrix<R>&
    LinearProgram<R>::tableau() 
    {
      return _tableau;
    }
    
    template<class R>
    inline 
    size_type
    LinearProgram<R>::number_of_constraints() const 
    {
      return _tableau.number_of_rows()-2;
    }
    
    template<class R>
    inline 
    size_type
    LinearProgram<R>::number_of_free_variables() const 
    {
      return _tableau.number_of_columns()-1;
    }
    
    template<class R>
    inline 
    size_type
    LinearProgram<R>::number_of_variables() const 
    {
      return _tableau.number_of_rows()+_tableau.number_of_columns()-2;
    }
    

    
    
/*
    template<class R>
    inline 
    LinearAlgebra::Vector<R>
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


#endif /* ARIADNE_LINEAR_PROGRAM_H */
