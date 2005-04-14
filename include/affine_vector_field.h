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


#include <affine_map.h>
#include <approx_type.h>

namespace Ariadne {	
namespace VectorField{
namespace Affine{

enum AffineKind {
	HOMOGENEOUS,
	NON_HOMOGENEOUS,
	TRASLATION
};


template <typename BS_I>
class AriadneAffineIntegrator {	
	
	public:
		typedef BS_I BasicSetIntegrator;
		typedef typename BasicSetIntegrator::Map Map;
		typedef typename Ariadne::Map::Affine::AriadneAffineMap< Map > SolutionMap;
		typedef typename Ariadne::Map::Affine::AriadneAffineMap< Map > VectorFieldMap;
		typedef typename Map::DenotableSet DenotableSet;
		typedef typename DenotableSet::BasicSet BasicSet;
		typedef typename SolutionMap::Real Real;
	
		typedef typename SolutionMap::Matrix Matrix;
		typedef typename SolutionMap::Vector Vector;
	
		AriadneAffineIntegrator() {}
		
		AriadneAffineIntegrator(const AriadneAffineIntegrator<BS_I>& orig): 
				_vf_map(orig._vf_map) {}
	
		AriadneAffineIntegrator(const VectorFieldMap &vfield): 
				_vf_map(vfield) {}
		
		AriadneAffineIntegrator(const Matrix &A): 
				_vf_map(A) {}
					
		AriadneAffineIntegrator(const Vector &b): 
				_vf_map(b) {}
					
		AriadneAffineIntegrator(const Matrix &A, const Vector &b): 
				_vf_map(A,b) {}
					
		inline SolutionMap get_solution(const Real& delta) const { 	
			
			Matrix A_sol=
				Ariadne::LinearAlgebra::exp_Ah((this->_vf_map).A() , delta,5);	
			
			Vector b_sol=
				Ariadne::LinearAlgebra::exp_b((this->_vf_map).A() , 
								(this->_vf_map).b(), delta,5);
			
			SolutionMap sol(A_sol,b_sol);
			
			return sol;
		}
		
		inline BasicSet get_flow_tube_from_to(const BasicSet &A, 
				const BasicSet &B, const SolutionMap &sol_map,
				const Ariadne::Geometry::ApproxKind &atype) {
					
				return (this->_bs_i).get_flow_tube_to(A,B,
								sol_map,atype);				
		}
		
		inline DenotableSet get_flow_tube_from_to(const DenotableSet &A, 
				const DenotableSet &B, const SolutionMap &sol_map, 
				const Ariadne::Geometry::ApproxKind &atype) {
				
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
					
			DenotableSet flow_tube(A.dim());
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
		
		inline AriadneAffineIntegrator& operator()(const  VectorFieldMap &vfield) {
			
			this->_vf_map=vfield;
		
			return *this;			
		}
		
		inline VectorFieldMap& vector_field() const {
			return 	_vf_map;
		}
		
	private:
		VectorFieldMap _vf_map;
	
		BasicSetIntegrator _bs_i;
};

template <typename BS_MAP>
class AriadneAffineVectorField {
	
	public:
		typedef BS_MAP VectorFieldMap;
		typedef typename VectorFieldMap::DenotableSet DenotableSet;
		typedef typename DenotableSet::BasicSet BasicSet;
		typedef typename BasicSet::State State;
		typedef typename State::Real Real;
	
		typedef typename boost::numeric::ublas::matrix<Real> Matrix;
		typedef typename boost::numeric::ublas::vector<Real> Vector;
	
		AriadneAffineVectorField(const AriadneAffineVectorField<BS_MAP> &T):
				 _vf_map(T._vf_map){}
	
		AriadneAffineVectorField(const Matrix &A): _vf_map(A) {}
	
		AriadneAffineVectorField(const Vector &b): _vf_map(b) {}
		
		AriadneAffineVectorField(const Matrix &A, const Vector &b):
			_vf_map(A,b) {}
		
		inline AffineKind affine_kind() {
			return (this->vector_field()).affine_kind();
		}
		
		inline const Matrix& A() const { return _vf_map.A(); }
		
		inline const Vector& b() const { return _vf_map.b(); }
		
		inline VectorFieldKind vector_field_kind() { return AFFINE; }
		
		inline const VectorFieldMap &vector_field() const {
			return (this->_vf_map);
		}
		
		inline size_t dim() const {
			return (this->_vf_map).dim();
		}
		
	private:
	
		VectorFieldMap _vf_map;
};

}}}
#endif /* _AFFINE_VECTOR_FIELD_H */
