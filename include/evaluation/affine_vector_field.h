/***************************************************************************
 *            affine_vector_field.h
 *
 *  Fri Feb  4 08:57:39 2005
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
 
#ifndef _AFFINE_VECTOR_FIELD_H
#define _AFFINE_VECTOR_FIELD_H


#include "../vector_field.h"
#include "../affine_map.h"
#include "../approx_type.h"

namespace Ariadne {
namespace VectorField{

enum AffineKind {
  HOMOGENEOUS,
  NON_HOMOGENEOUS,
  TRASLATION
};


template<typename R>
integrate(const AffineVectorField<R>& vf, const Geometry::Polyhedron<R> p, R t, R h) {
  
}


template <typename R>
class AffineVectorField : public VectorField<R,Geometry::Point> 
{
  typedef typename Geometry::Polyhedron<R> Polyhedron;

 public:
  typedef R Real;
  typedef typename Geometry::Point<Real> State;
  
  typedef typename boost::numeric::ublas::matrix<Real> Matrix;
  typedef typename boost::numeric::ublas::vector<Real> Vector;

  AffineVectorField(const AffineVectorField<R>& T) : _map(T.A(),T.b()) { }
  AffineVectorField(const Matrix &A, const Vector &b) : _map(A,b) { }

  State operator() (const State& s);
  Polyhedron operator() (const Polyhedron& p) { return apply(_map,p); }

  inline const Matrix& A() const { return this->_map.A(); }
  inline const Vector& b() const { return this->_map.b(); }
  
  inline size_t dimension() const {
    return (this->_map).dimension();
  }
  
  /*! Deprecated. */
  inline size_t dim() const {
    return this->dimension();
  }

 private:
  Map::AffineMap<Real> _map;
};
 


template<class R, template<typename> class BS> class AffineMapIntegrator;

template<class R> 
class AffineMapIntegrator<R,Geometry::Polyhedron> 
{
  typedef R Real;
  typedef Geometry::Point<Real> State;
  typedef Geometry::Polyhedron<Real> BasicSet;
  typedef Geometry::Polyhedron<Real> Polyhedron;
  typedef Geometry::ListSet<Real,Polyhedron> DenotableSet;
  typedef AffineVectorField<Real> VectorField;
  typedef typename AffineVectorField<Real>::Matrix Matrix;
  typedef typename AffineVectorField<Real>::Vector Vector;
  
  AffineIntegrator() { }
  AffineIntegrator(const int& n) : _terms(n) { };

  AffineMap solution_map(const VectorField& vf, const Real& t) {
    Matrix expA = Ariadne::LinearAlgebra::exp_Ah(vf.A(),t,n);
    Matrix expb = Ariadne::LinearAlgebra::exp_b(vf.A(),vf.b(),h,n);
    return AffineMap(expA,expB);
  }
    
  inline BasicSet integrate(const VectorField& vf,
                            const BasicSet& A, 
                            const Real& t) 
  {
    return solution_map(vf,h)(A); 
  }
      
  inline DenotableSet integrate(const VectorField& vf,
                                const DenotableSet& A, 
                                const Real& t)
  { 
    DenotableSet result;
    AffineMap solution=solution_map(vf,t);
    for(uint i=0; i< A.size(); i++) {
      result.inplace_union(solution(A[i]));
    }
    return result;
  }

  inline DenotableSet integrate(const VectorField& vf,
                                const DenotableSet& A, 
                                const Real& t1,
                                const Real& t2)
  {
    assert(0<=t1 && t1<=t2);
    
    DenotableSet result;
    DenotableSet B;

    if(t1==0) {
      B=A;
    } 
    else {
      B=integrate(vf,A,t);
    }
    
    if(t2==t1) {
      return B;
    }

    AffineMap affmap=solution_map(vf,t2-t1);
    
    for(uint i=0; i!=B.size(); ++i) {
      BasicSet bs=B[i];
      result.inplace_union(Ariadne::Geometry::convex_hull(bs,affmap(bs)));
    }
    return result;
  }
 private:
  Real _time_step;
  unsigned int _terms;
};
  



template<typename R>
ListSet<R,Geometry::Polyhedron>
AffineIntegrator<R,Geometry::Polyhedron>::
integrate(const VectorField& vf,
          const Polyhedron& p,
          const Real& h)
{
  ListSet<R,Geometry::Polyhedron> result;

  const Real& h=_time_step;
  const unsigned int& n=_terms;
  const Matrix& A=vf.A();
  const Matrix& b=vf.b();
  
  Matrix expA = Ariadne::LinearAlgebra::exp_Ah(A,h,n);
  Matrix expb = Ariadne::LinearAlgebra::exp_b(A,b,h,n);
  AffineMap solution(expA,expB);
  
  Real t=0;
  Polyhedron q=p;
  
  while(t<t_1) {
    q=solution(q);
  }
  

}
  
  inline BasicSet get_flow_tube_from_to(const BasicSet& A, 
                                        const BasicSet& B, const SolutionMap& sol_map,
                                        const Ariadne::Geometry::ApproxKind& atype) {
    
    return (this->_bs_i).get_flow_tube_to(A,B,
                                          sol_map,atype);
  }
  
  inline DenotableSet get_flow_tube_from_to(const DenotableSet& A, 
                                            const DenotableSet& B, const SolutionMap& sol_map, 
                                            const Ariadne::Geometry::ApproxKind& atype) {
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
    DenotableSet flow_tube(A.dimension());
    size_t i;
    
    for (i=0; i< A.size(); i++) {
      flow_tube.inplace_union(
                              (this->_bs_i).get_flow_tube_from_to(A[i], B[i],
                                                                  sol_map,atype));
    }
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
    return flow_tube;
  }
  
  inline AffineIntegrator& operator()(const  VectorFieldMap& vfield) {
    
    this->_vf_map=vfield;
    
    return *this;
  }
  
  inline VectorFieldMap& vector_field() const {
    return _vf_map;
  }
  
 private:
  VectorFieldMap _vf_map;
  
  BasicSetIntegrator _bs_i;
 
}}
#endif /* _AFFINE_VECTOR_FIELD_H */
