/***************************************************************************
 *            constraint.h
 *
 *  Thu Feb  3 09:31:28 2005
 *  Copyright  2005  Alberto Casagrande
 *  casagrande@dimi.uniud.it
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
 
/*! \file constraint.h
 *  \brief Linear inequalities and constraints.
 */
 
#ifndef _ARIADNE_CONSTRAINT_H
#define _ARIADNE_CONSTRAINT_H

#include <algorithm>
#include <ppl.hh>

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

namespace Ariadne {
namespace LinearAlgebra {

/*! \brief A collection of linear equalities, inequalities and strict inequalities. */
  
template <typename R> 
class ConstraintSystem {
 public:
  
  enum RelType{
    EQUALITY,
    INEQUALITY,
    STRICT_INEQUALITY
  };
  
  typedef ::Ariadne::LinearAlgebra::matrix<R> matrix_type;
  typedef ::Ariadne::LinearAlgebra::vector<R> vector_type;

  typedef std::vector<RelType> Rel_vector_type;
  
  ConstraintSystem() { }
  
  ConstraintSystem(const Parma_Polyhedra_Library::Constraint_System& ppl_cs);
  
  ConstraintSystem(matrix_type new_C, vector_type new_d)
    : _C(new_C), _d(new_d), _rel_types(new_C.size1()) 
  { 
    for(size_type i=0; i!=_rel_types.size(); ++i) {
      _rel_types[i]=INEQUALITY;
    }
  }

  ConstraintSystem(matrix_type new_C, vector_type new_d , Rel_vector_type new_r_type)
    : _C(new_C), _d(new_d), _rel_types(new_r_type) 
  { }

  Parma_Polyhedra_Library::Constraint_System ppl_constraint_system() const;
  Parma_Polyhedra_Library::Constraint_System open_ppl_constraint_system() const;
  operator Parma_Polyhedra_Library::Constraint_System() const { return ppl_constraint_system(); }

  inline size_t space_dimension() const {
    return this->_C.size2();
  }
  
  inline size_t number_of_constraints() const {
    return this->_C.size1();
  }
  
  inline size_t number_of_equalities() const {
    return std::count(_rel_types.begin(), _rel_types.end(),EQUALITY);
  }
  
  inline bool empty() const{
    return (this->_C.size1()==0);
  }
  
  
  inline bool already_at_precision(const R& delta) const;
  inline void reduce_precision_to_expanding(const R& delta);
  inline void reduce_precision_to_shrinking(const R& delta);
  inline void divide_all_by(const R& delta);
  inline void expand_by(const R& delta);
  inline void reset_dimensions(const size_t& dim, const size_t& nb_constraints);
  inline void expand_equality_by(const R &delta);

  inline bool is_equality(const size_t& i) const { 
    return(this->_rel_types[i]==EQUALITY); 
  }
  
  inline void set_as_equality(const size_t &i) {
    this->_rel_types[i]=EQUALITY;
  }
  
    
  inline ConstraintSystem<R>& operator=(const ConstraintSystem<R>& cs) {    
    if(this != &cs) {
      this->_C = cs._C;
      this->_d = cs._d;
      this->_rel_types = cs._rel_types;
    }
    return *this;
  }
 private:
  matrix_type _C;
  vector_type _d;
  Rel_vector_type _rel_types;
 private:  
  inline bool any_equality() const{
    for (size_t i=0; i< this->number_of_constraints(); i++) {
      if (this->_rel_types[i]==EQUALITY) {
        return true;
      }
    }    
    return false;
  }
};





template<typename R>
bool 
ConstraintSystem<R>::already_at_precision(const R& delta) const 
{
  R min;
  R too_big=delta*delta*delta;
  
  for (size_t i=0; i<this->number_of_constraints(); ++i) {
    min=0.0;
    if ( abs(this->_d(i)) != 0 ) {
      min=abs(this->_d(i));
    }
    for (size_t j=0; j<this->space_dimension(); ++j) {
      if ( (abs(this->_C(i,j))!=0) && ((min==0)||(abs(this->_C(i,j))<min)) ) {
        min=abs(this->_C(i,j));
      }
    }
    if (min>too_big) {
      return false;
    }
  }
  return true;
}


template<typename R>
void 
ConstraintSystem<R>::reduce_precision_to_expanding(const R &delta)
{
  R too_big=delta*delta*delta;
  R min;

  for (size_t i=0; i<this->number_of_constraints(); ++i) {
    min=0.0;
    if (abs(this->_d(i))!=0) {
      min=abs(this->_d(i));
    }
    for (size_t j=0; j< this->space_dimension(); ++j) {
      if  ((abs(this->_C(i,j))!=0)&&((min==0)||(abs(this->_C(i,j))<min))) {
        min=abs(this->_C(i,j));
      }
    }
    if (min>too_big) {
      while (min>too_big) {
        min=numerator(min)/numerator(delta);
        for (size_t j=0; j<this->space_dimension(); ++j) {
          this->_C(i,j)=numerator(this->_C(i,j))/numerator(delta);
        }
        this->_d(i)=numerator(this->_d(i))/numerator(delta);
      }
      this->_d(i)-=1;
    }
  }
}

template<typename R>
void 
ConstraintSystem<R>::reduce_precision_to_shrinking(const R &delta)
{  
  size_t i,j;
  R min;
  
  for (j=0; j<this->number_of_constraints(); j++) {
    min=0.0;
    if (abs(this->_d(j))!=0) {
      min=abs(this->_d(j));
    }
    for (i=0; i< this->space_dimension(); i++) {
      if  ((abs(this->_C(j,i))!=0)&&((min==0)||(abs(this->_C(j,i))<min))) {
        min=abs(this->_C(j,i));
      }
    }
    if (min>delta) {
      while (min>delta) {
        min=numerator(min)/numerator(delta);
        for (i=0; i< this->space_dimension(); i++) {
          this->_C(j,i)=numerator(this->_C(j,i))/numerator(delta);
        }
        this->_d(j)=numerator(this->_d(j))/numerator(delta);
      }
      this->_d(j)--;
    }
  }
}


template<typename R>
void 
ConstraintSystem<R>::divide_all_by(const R& delta)
{
  for (size_t i=0; i<this->number_of_constraints(); ++i) {
    for (size_t j=0; j<this->space_dimension(); ++j) {
      this->_C(i,j)=numerator(this->_C(i,j))/numerator(delta);
    }
    this->_d(i)=numerator(this->_d(i))/numerator(delta);
  }
}


template<typename R>
void
ConstraintSystem<R>::expand_by(const R& delta) 
{
  if (delta==0) {
    return;
  }
  
  const size_t dim=this->space_dimension();

  if (!this->any_equality()) {
    for (size_t j=0; j< this->number_of_constraints(); j++) {
      this->_d(j)+=delta;
    }
    return;
  }

  if (delta<0) {
    matrix_type new_C(1, dim);
    vector_type new_d(1);
    this->_rel_types.resize(1);
    for (size_t i=0; i< dim; i++) {
      new_C(1,i)=0;
    }
    new_d(1)=-1;
    this->_C=new_C;
    this->_d=new_d;
  
    return;
  } 
  else {
    this->expand_equality_by(delta);
    return;
  }
}


template<typename R>
inline 
void 
ConstraintSystem<R>::reset_dimensions(const size_t& dim, 
                                      const size_t& nb_constraints) 
{
  matrix_type new_C(nb_constraints,dim);
  vector_type new_d(nb_constraints);
  
  Rel_vector_type new_r_type(nb_constraints);
  
  _C=new_C;
  _d=new_d;
  _rel_types = new_r_type;
}


template<typename R>
void 
ConstraintSystem<R>::expand_equality_by(const R& delta) 
{
  const size_t dim=this->space_dimension();
  size_t new_constr=this->number_of_equalities()+this->number_of_constraints();
  size_t new_i=0;
  
  R den=denominator(delta), num=numerator(delta);
  
  matrix_type new_C(new_constr, dim);
  vector_type new_d(new_constr);
  
  /* TODO: REIMPLEMENT */
  for (size_t i=0; i<this->number_of_constraints(); ++i) {
    if (!this->is_equality(i)) {
      for (size_t j=0; j<dim; ++j) {
        new_C(new_i,j)=this->_C(i,j);
      }
      new_d(new_i)=this->_d(i);
      ++new_i;
    } 
    else {
      for (size_t j=0; j<dim; ++j) {
        new_C(new_i,j)=-this->_C(i,j)*den;
        new_C(new_i+1,j)=this->_C(i,j)*den;
      }
      new_d(new_i)=-this->_d(i)*den-num;
      new_d(new_i+1)=this->_d(i)*den-num;
      new_i+=2;
    }
  }
  
  this->_C=new_C;
  this->_d=new_d;
  
  this->_rel_types.resize(new_constr);
  
  for (size_t i=0; i<dim; ++i) {
    this->_rel_types[i]=INEQUALITY;
  }

}
  

template <typename R>
ConstraintSystem<R>::ConstraintSystem(const Parma_Polyhedra_Library::Constraint_System& ppl_cs)
  : _C(), _d(), _rel_types()
{
  Parma_Polyhedra_Library::Constraint_System::const_iterator j_cs, begin, end;
  
  begin=ppl_cs.begin();
  end=ppl_cs.end();
            
  size_t space_dim=ppl_cs.space_dimension();
  size_t nb_constr=0,i,j;
            
  for (j_cs=begin; j_cs!=end; j_cs++) {
    nb_constr++;
  }
  
  this->reset_dimensions(space_dim, nb_constr);
  
  j=0;
  
  for (j_cs=begin; j_cs!=end; j_cs++) {
    for (i=0; i< space_dim; ++i) {
        this->_C(j,i)=j_cs->coefficient(Parma_Polyhedra_Library::Variable(i));
    }
    this->_d(j)=-j_cs->inhomogeneous_term();    
    if (j_cs->is_equality()) {
      this->_rel_types[j]=EQUALITY;
    } 
    else {
      if (j_cs->is_strict_inequality()) {
        this->_rel_types[j]=STRICT_INEQUALITY;
      } 
      else {
        this->_rel_types[j]=INEQUALITY;
      }
    }
    ++j;
  }
}


template <typename R>
Parma_Polyhedra_Library::Constraint_System
ConstraintSystem<R>::ppl_constraint_system() const
{
  Parma_Polyhedra_Library::Constraint_System ppl_cs;
  Parma_Polyhedra_Library::Linear_Expression ppl_lin_expr;
  
  for (size_t i=0; i<this->number_of_constraints(); ++i) {
    ppl_lin_expr=Rational(this->_C(i,0))*Parma_Polyhedra_Library::Variable(0);
    for (size_t j=1; j<this->space_dimension(); ++j) {
      ppl_lin_expr+=Rational(this->_C(i,j))*Parma_Polyhedra_Library::Variable(j);
    }
    if (this->is_equality(i)) {    
      ppl_cs.insert(ppl_lin_expr == Rational(this->_d(i)));
    } 
    else {
      ppl_cs.insert(ppl_lin_expr >= Rational(this->_d(i)));
    }
  }
  return ppl_cs;
}


template <typename R>
inline 
Parma_Polyhedra_Library::Constraint_System
ConstraintSystem<R>::open_ppl_constraint_system() const
{
  const ConstraintSystem<R>& cs(*this);
  
  /* if there is an equality, the open set is empty */
  if (cs.any_equality()) {
    return Parma_Polyhedra_Library::Constraint_System();
  }
   
  Parma_Polyhedra_Library::Constraint_System ppl_cs(cs);
  Parma_Polyhedra_Library::Linear_Expression ppl_lin_expr;
   
 
  for (size_t i=0; i<cs.number_of_constraints(); ++i) {
    ppl_lin_expr=Rational(cs._C(i,0))*Parma_Polyhedra_Library::Variable(0);
    for (size_t j=1; j<cs.space_dimension(); ++j) {
      ppl_lin_expr += Rational(cs._C(i,j))*Parma_Polyhedra_Library::Variable(j);
    }
    ppl_cs.insert(ppl_lin_expr > Rational(cs._d(i)));
  }
  
  return ppl_cs;
}


}
}

#endif /* _ARIADNE_CONSTRAINT_H */
