/***************************************************************************
 *            float_concept.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
/*! \file float_concept.h
 *  \brief Specification of the floating-point concept
 */

namespace Ariadne {
  namespace Numeric {

#ifdef DOXYGEN
    /*!\ingroup Numeric
     * \brief The specification of the floating-point concept.
     *
     * \internal This specification is not finished; some functions may be added or removed. 
     */
    class Float {
     public:
      //@{ 
      //! \name Constructors and assignment operators.
      /*! \brief Default constructor. Yields zero. */
      Float();
      /*! \brief Construct from an integer. */
      Float(int n);
      /*! \brief Construct from a double-precision number. */
      Float(double x);
      /*! \brief Construct from a string literal. Throws an exception if the numerical value of the string is not exactly representable. */
      Float(const std::string& str);
      /*! \brief Copy constructor. */
      Float(const Float& x);
      
      /*! \brief Assign from an integer. */
      Float& operator=(const int& n);
      /*! \brief Assign from a double-precision number. */
      Float& operator=(const double& x);
      /*! \brief Assign from another floating-point number. Throws an exception if \a x has higher precision and the assignment cannot be performed exactly. */
      Float& operator=(const Float& x);
      //@}

      //@{ 
      //! \name Conversion operators.
      /*! \brief Convert to a rational. (Recommended); */
      operator Rational (); const;
      //@}

      //@{ 
      //! \name Precision
      /*! \brief . */
      friend unsigned int precision(const Float& x);
      /*! \brief (Optional) */
      friend void set_precision(const Float& x, unsigned int p);
      //@}


      //@{ 
      //! \name Comparison operators
      /*! \brief . */
      friend bool operator==(const Float& x1, const Float& x2);
      /*! \brief . */
      friend bool operator!=(const Float& x1, const Float& x2);
      /*! \brief . */
      friend bool operator<=(const Float& x1, const Float& x2);
      /*! \brief . */
      friend bool operator>=(const Float& x1, const Float& x2);
      /*! \brief . */
      friend bool operator< (const Float& x1, const Float& x2);
      /*! \brief . */
      friend bool operator> (const Float& x1, const Float& x2);
      //@}

      //@{ 
      //! \name Mixed comparison operators
      /*! \brief . */
      friend bool operator==(const Float& x1, const int& x2);
      /*! \brief . */
      friend bool operator!=(const Float& x1, const int& x2);
      /*! \brief . */
      friend bool operator<=(const Float& x1, const int& x2);
      /*! \brief . */
      friend bool operator>=(const Float& x1, const int& x2);
      /*! \brief . */
      friend bool operator< (const Float& x1, const int& x2);
      /*! \brief . */
      friend bool operator> (const Float& x1, const int& x2);
   

      /*! \brief . */
      friend bool operator==(const Float& x1, const double& x2);
      /*! \brief . */
      friend bool operator!=(const Float& x1, const double& x2);
      /*! \brief . */
      friend bool operator<=(const Float& x1, const double& x2);
      /*! \brief . */
      friend bool operator>=(const Float& x1, const double& x2);
      /*! \brief . */
      friend bool operator< (const Float& x1, const double& x2);
      /*! \brief . */
      friend bool operator> (const Float& x1, const double& x2);
   

      /*! \brief . */
      friend bool operator==(const Float& x1, const Rational& x2);
      /*! \brief . */
      friend bool operator!=(const Float& x1, const Rational& x2);
      /*! \brief . */
      friend bool operator<=(const Float& x1, const Rational& x2);
      /*! \brief . */
      friend bool operator>=(const Float& x1, const Rational& x2);
      /*! \brief . */
      friend bool operator< (const Float& x1, const Rational& x2);
      /*! \brief . */
      friend bool operator> (const Float& x1, const Rational& x2);
       //@}
  
      
      //@{ 
      //! \name Exact arithmetic
      /*! \brief . */
      friend Float min(const Float& x1, const Float& x2);
      /*! \brief . */
      friend Float max(const Float& x1, const Float& x2);
      /*! \brief . */
      friend Float neg(const Float& x);
      /*! \brief . */
      friend Float abs(const Float& x);
      /*! \brief . */
      friend Float floor(const Float& x); 
      /*! \brief . */
      friend Float ceil(const Float& x); 
      
      //@}


      //@{ 
      //! \name Conversion operators
      /*! \brief . */
      friend double conv_approx(const Float& x);

      /*! \brief . */
      friend Float conv_exact(const Float& x);
      /*! \brief . */
      friend Float conv_approx(const Float& x);
      /*! \brief . */
      friend Float conv_down(const Float& x);
      /*! \brief . */
      friend Float conv_up(const Float& x);

      /*! \brief . */
      friend Float conv_exact(const int& n);
      /*! \brief . */
      friend Float conv_approx(const int& n);
      /*! \brief . */
      friend Float conv_down(const int& n);
      /*! \brief . */
      friend Float conv_up(const int& n);

      /*! \brief . */
      friend Float conv_exact(const double& x);
      /*! \brief . */
      friend Float conv_approx(const double& x);
      /*! \brief . */
      friend Float conv_down(const double& x);
      /*! \brief . */
      friend Float conv_up(const double& x);

      /*! \brief . */
      friend Float conv_approx(const Rational& x);
      /*! \brief . */
      friend Float conv_down(const Rational& x);
      /*! \brief . */
      friend Float conv_up(const Rational& x);
      
      /*! \brief . */
      friend int int_down(const Float& x);
      /*! \brief . */
      friend int int_up(const Float& x);
      //@}


 
   
      //@{ 
      //! \name Approximate arithmetic
      /*! \brief . */
      friend Float add_down(const Float& x1,const Float& x2); 
      /*! \brief . */
      friend Float add_up(const Float& x1,const Float& x2); 
      /*! \brief . */
      friend Float add_approx(const Float& x1,const Float& x2); 
    
      /*! \brief . */
      friend Float sub_down(const Float& x1,const Float& x2); 
      /*! \brief . */
      friend Float sub_up(const Float& x1,const Float& x2); 
      /*! \brief . */
      friend Float sub_approx(const Float& x1,const Float& x2); 
    
      /*! \brief . */
      friend Float mul_down(const Float& x1,const Float& x2); 
      /*! \brief . */
      friend Float mul_up(const Float& x1,const Float& x2); 
      /*! \brief . */
      friend Float mul_approx(const Float& x1,const Float& x2); 
    
      /*! \brief . */
      friend Float div_down(const Float& x1,const Float& x2); 
      /*! \brief . */
      friend Float div_up(const Float& x1,const Float& x2); 
      /*! \brief . */
      friend Float div_approx(const Float& x1,const Float& x2); 
    
      /*! \brief . */
      friend Float med_approx(const Float& x1,const Float& x2); 
    
      /*! \brief . */
      friend Float pow_down(const Float& x1,const int& n2);
      /*! \brief . */
      friend Float pow_up(const Float& x1,const int& n2);
      /*! \brief . */
      friend Float pow_approx(const Float& x1, const int& n2); 
      //@}

      //@{
      //! \name Approximate mixed arithmetic
      /*! \brief . */
      friend Float mul_approx(const int& n, const Float& x); 
      /*! \brief . */
      friend Float mul_approx(const Float& x, const int& n); 
      /*! \brief . */
      friend Float div_approx(const Float& x, const int& n); 
      
      /*! \brief . */
      friend Float mul_approx(const uint& n, const Float& x); 
      /*! \brief . */
      friend Float mul_approx(const Float& x, const uint& n); 
      /*! \brief . */
      friend Float div_approx(const Float& x, const uint& n); 
      
      /*! \brief . */
      friend Float mul_approx(const double& d, const Float& x); 
      /*! \brief . */
      friend Float mul_approx(const Float& x, const double& d); 
      //@}

      
      //@{
      //! \name Approximate algebraic functions
      /*! \brief . */
      friend Float sqrt_down(const Float& x); 
      /*! \brief . */
      friend Float sqrt_up(const Float& x); 
      /*! \brief . */
      friend Float sqrt_approx(const Float& x); 

      /*! \brief . */
      friend Float hypot_down(const Float& x1,const Float& x2);
      /*! \brief . */
      friend Float hypot_up(const Float& x1,const Float& x2);
      /*! \brief . */
      friend Float hypot_approx(const Float& x1,const Float& x2); 
      //@}

      //@{
      //! \name Approximate constants
      /*! \brief . */
      friend Float pi_approx(); 
      /*! \brief . */
      friend Float pi_down(); 
      /*! \brief . */
      friend Float pi_up();
      //@}

      //@{
      //! \name Approximate transcendental functions
      /*! \brief . */
      friend Float exp_approx(const Float& x); 
      /*! \brief . */
      friend Float exp_down(const Float& x); 
      /*! \brief . */
      friend Float exp_up(const Float& x); 
    
      /*! \brief . */
      friend Float log_approx(const Float& x); 
      /*! \brief . */
      friend Float log_down(const Float& x); 
      /*! \brief . */
      friend Float log_up(const Float& x); 
      //@}


      //@{
      //! \name Interval arithmetic
      /*! \brief . */
      friend Interval<Float> add(const Float& x1, const Float& x2); 
      /*! \brief . */
      friend Interval<Float> add(const Float& x1, const Interval<Float>& x2); 
      /*! \brief . */
      friend Interval<Float> add(const Interval<Float>& x1, const Float& x2); 
      /*! \brief . */
      friend Interval<Float> add(const Interval<Float>& x1, const Interval<Float>& x2); 
      /*! \brief . */
      friend Interval<Float> sub(const Float& x1, const Float& x2); 
      /*! \brief . */
      friend Interval<Float> sub(const Float& x1, const Interval<Float>& x2); 
      /*! \brief . */
      friend Interval<Float> sub(const Interval<Float>& x1, const Float& x2); 
      /*! \brief . */
      friend Interval<Float> sub(const Interval<Float>& x1, const Interval<Float>& x2); 
      /*! \brief . */
      friend Interval<Float> mul(const Float& x1, const Float& x2); 
      /*! \brief . */
      friend Interval<Float> mul(const Float& x1, const Interval<Float>& x2); 
      /*! \brief . */
      friend Interval<Float> mul(const Interval<Float>& x1, const Float& x2); 
      /*! \brief . */
      friend Interval<Float> mul(const Interval<Float>& x1, const Interval<Float>& x2); 
      /*! \brief . */
      friend Interval<Float> div(const Float& x1, const Float& x2); 
      /*! \brief . */
      friend Interval<Float> div(const Float& x1, const Interval<Float>& x2); 
      /*! \brief . */
      friend Interval<Float> div(const Interval<Float>& x1, const Float& x2); 
      /*! \brief . */
      friend Interval<Float> div(const Interval<Float>& x1, const Interval<Float>& x2); 
      /*! \brief . */
      friend Interval<Float> pow(const Float& x1, const int& n2); 
      /*! \brief . */
      friend Interval<Float> pow(const Interval<Float>& x1, const int& n2); 
      //@}

      //@{
      //! \name Interval algebraic functions
      /*! \brief . */
      friend Interval<Float> sqrt(const Float& x); 
      /*! \brief . */
      friend Interval<Float> sqrt(const Interval<Float>& x); 
      /*! \brief . */
      friend Interval<Float> hypot(const Float& x1, const Float& x2); 
      /*! \brief . */
      friend Interval<Float> hypot(const Interval<Float>& x1, const Interval<Float>& x2); 
      //@}

      //@{
      //! \name Interval constants
      /*! \brief . */
      friend Interval<Float> pi();
      //@}

      //@{
      //! \name Interval transcendental functions
      /*! \brief . */
      friend Interval<Float> exp(const Interval<Float>& x); 

      /*! \brief . */
      friend Interval<Float> log(const Interval<Float>& x); 

      /*! \brief . */
      friend Interval<Float> sin(const Interval<Float>& x); 

      /*! \brief . */
      friend Interval<Float> cos(const Interval<Float>& x); 

      /*! \brief . */
      friend Interval<Float> tan(const Interval<Float>& x); 

      /*! \brief . */
      friend Interval<Float> sinh(const Interval<Float>& x); 

      /*! \brief . */
      friend Interval<Float> cosh(const Interval<Float>& x); 

      /*! \brief . */
      friend Interval<Float> tanh(const Interval<Float>& x); 

      /*! \brief . */
      friend Interval<Float> asin(const Interval<Float>& x); 

      /*! \brief . */
      friend Interval<Float> acos(const Interval<Float>& x); 

      /*! \brief . */
      friend Interval<Float> atan(const Interval<Float>& x); 

      /*! \brief . */
      friend Interval<Float> asinh(const Interval<Float>& x); 

      /*! \brief . */
      friend Interval<Float> acosh(const Interval<Float>& x); 

      /*! \brief . */
      friend Interval<Float> atanh(const Interval<Float>& x); 
      //@}
      
    
      //@{ 
      //! \name Exact arithmetic
      /*! \brief . */
      friend Float operator+(const Float& x);
      /*! \brief . */
      friend Float operator-(const Float& x);
      //@}
      
    
      //@{ 
      //! \name Interval arithmetic
      /*! \brief . */
      friend Interval<Float> operator+(const Float& x1, const Float& x2);
      /*! \brief . */
      friend Interval<Float> operator-(const Float& x1, const Float& x2);
      /*! \brief . */
      friend Interval<Float> operator*(const Float& x1, const Float& x2);
      /*! \brief . */
      friend Interval<Float> operator/(const Float& x1, const Float& x2);
      //@}
      

      //@{ 
      //! \name Stream input / output
      /*! \brief . */
      friend std::ostream& operator<<(std::ostream& os, const Float& x);
      /*! \brief . */
      friend std::istream& operator>>(std::istream& is, Float& x);
      //@}

    /*! \brief An interval of a floating-point type. */
    template< > class Interval<Float>
    {
   };
    };
      

  }
}
#endif
