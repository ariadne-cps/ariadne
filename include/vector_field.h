/***************************************************************************
 *            vector_field.h
 *
 *  Thu Feb  3 21:06:54 2005
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
 
#ifndef _VECTOR_FIELD_H
#define _VECTOR_FIELD_H

namespace Ariadne {	
namespace VectorField{
	
enum VectorFieldKind {
	LINEAR,
	AFFINE,
	MULTIVALUE,
	GENERAL
};

}
}

#include <affine_vector_field.h>
#include <approx_type.h>

namespace Ariadne {	
  namespace VectorField{

    template <typename I>
    class Integrator {	
     public:
      typedef typename I::VectorFieldMap VectorFieldMap;
      typedef typename I::SolutionMap SolutionMap;
      typedef typename SolutionMap::DenotableSet DenotableSet;
      typedef typename DenotableSet::BasicSet BasicSet;
      typedef typename SolutionMap::Real Real;
      
     private:
      inline void _change_delta_step(const Real &delta) {
        
        this->_sol_map=(this->_integrator).get_solution(delta);
        
        this->_delta_step=delta;
      }
      
     public:
      Integrator() {}
      
      
      Integrator(const Integrator<I>& orig): 
        _integrator(orig._integrator),
        _sol_map(orig._sol_map),_delta_step(orig._delta_step) {}
      
      Integrator(const VectorFieldMap &vfield): 
        _integrator(vfield),_delta_step(0){}
      
      Integrator(const VectorFieldMap &vfield, 
                 const Real &delta): _integrator(vfield), _delta_step(delta) {
        
        this->_sol_map=(this->_integrator).get_solution(delta);
      }
      
      template < typename VF >
      Integrator(const VF &vfield): _delta_step(0){
        
        this->_integrator(vfield.vector_field());	
      }
      
      template < typename VF >
      Integrator(const VF &vfield, const Real &delta ): 
        _delta_step(delta){
        
        this->_integrator(vfield.vector_field());
      }
      
      inline BasicSet integrate_from_for(const BasicSet &A, 
                                         const Real &delta) { 		
        
        if (delta!=this->_delta_step) {
          this->_change_delta_step(delta);
        }
        
        return _sol_map(A); 
        
      }
      
      inline DenotableSet integrate_from_for(const DenotableSet &A, 
                                             const Real &delta) { 
        
        if (delta!=this->_delta_step) {
          this->_change_delta_step(delta);
        }
        
        return _sol_map(A);
      }
      
      inline BasicSet get_flow_tube_from_for_to(const BasicSet &A, 
                                                const Real &delta, BasicSet &B, 
                                                const Ariadne::Geometry::ApproxKind &atype) {	
        
        B=this->integrate_from_for(A,delta);
        
        return (this->_integrator).get_flow_tube_from_to(A,B,
                                                         this->_sol_map,atype); 
      }
      
      inline DenotableSet get_flow_tube_from_for_to(const DenotableSet &A, 
                                                    const Real &delta, DenotableSet &B, 
                                                    const Ariadne::Geometry::ApproxKind &atype) {
        
#ifdef DEBUG
        std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
	
        B=this->integrate_from_for(A,delta);
        
        
#ifdef DEBUG
        
        DenotableSet out=(this->_integrator).get_flow_tube_from_to(A,B,
                                                                   this->_sol_map,atype);
        
        std::cout << __FILE__ << ":" << __LINE__ << std::endl;
        
        return out;
#endif
	
        
        return (this->_integrator).get_flow_tube_from_to(A,B,
                                                         this->_sol_map,atype);
      }
      
      inline Integrator& operator()(const  VectorFieldMap &vfield) {
        
        this->_integrator(vfield);
        
        this->_delta_step=0;
        
        return *this;			
      }
      
      inline VectorFieldMap& vector_field() const {
        return (this->_integrator).vector_field();
      }
      
     private:
      I _integrator;
      
      SolutionMap _sol_map;
      
      Real _delta_step;
    };
    
    
    template <typename VF>
    class VectorField {
      
     public:
      typedef typename VF::VectorFieldMap VectorFieldMap;
      typedef typename VF::DenotableSet DenotableSet;
      typedef typename DenotableSet::BasicSet BasicSet;
      typedef typename BasicSet::State State;
      typedef typename State::Real Real;
      
      VectorField(const VF &T):_vector_field(T) {}
      
      inline VectorFieldKind vector_field_kind() {
        return _vector_field.vector_field_kind();}
      
      inline VectorFieldMap &vector_field() const{
        return T._vector_field();
      }
      
      inline size_t dimension() const {
        return (this->_vector_field).dimension();
      }
      
     private:
      
      VF _vector_field;
    };
    
  }
}

#endif /* _VECTOR_FIELD_H */
