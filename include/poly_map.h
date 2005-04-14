/***************************************************************************
 *            poly_map.h
 *
 *  Tue Feb  1 12:18:22 2005
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
 
#ifndef _POLY_MAP_H
#define _POLY_MAP_H

#include <linear_algebra.h>
#include <map.h>
#include <constraint.h>
#include <generator.h>
#include <polyhedron.h>

namespace Ariadne {	
namespace Map{
namespace Affine {
	
/*! \brief Affine trasformation for polyhedron.
 *
 * This class represents AriadnePolyhedron's affine trasformations. 
 * Hence all the objects of this class are of the type 
 * \f$T(x)=A\vec{x}+\vec{b}\f$. We can distringuish two cases: \f$A\f$ is 
 * invertible and \f$A\f$ is not invertible.
 *
 * If \f$A\f$ invertible, to apply a transformation \f$T\f$ to a polyhedron 
 * \f$P\f$, we should consider the polyhedron constraint
 * representation: \f$\vec{x} \in P\f$ if and only if 
 * \f$C\vec{x}\geq \vec{d}\f$, where \f$C\f$ and \f$\vec{d}\f$ are 
 * an integer matrix and an integer vector
 * respectevely. The set \f$T(P)\f$ is the set vector \f$x'\f$  such that 
 * \f$C A^{-1} (x' - \vec{b}) \geq \vec{d}\f$. Thus \f$T(P)\f$ is the
 * set of costraints \f$C A^{-1} x' \geq C A^{-1}\vec{b} + \vec{d}\f$.
 * Let \f$E\f$ and \f$n_3\f$ be respectevely an integer matrix and a 
 * natural value 
 * such that \f$n_3 A^{-1}=\bar{A}\f$, then \f$T(P)\f$ is 
 * the set of costraints \f$C \bar{A} x' \geq C \bar{A} \vec{b} + 
 * n_3 \vec{d}\f$.
 * Moreover, if we assume \f$E=n_2 \bar{A}\f$, 
 * \f$n_2 \bar{A} \vec{b}=\vec{f}\f$ and \f$n_1= n_2 n_3\f$, where 
 * \f$\vec{f}\f$ and \f$n_2\f$ are respectevely an integer vector and 
 * a natural value, \f$T(P)\f$ is defined by
 * the costraint system \f$ C E x' \geq C \vec{f} + n_1 \vec{d}\f$.
 *
 * Otherwise, if \f$A\f$ is not invertible, to apply a transformation, 
 * we should consider the polyhedron generator representation i.e.~ the 
 * set of points \f$\mathcal{p}\f$, closure point\f$\mathcal{cp}\f$, 
 * ray \f$\mathcal{r}\f$ and line \f$\mathcal{l}\f$ defining the polyhedron.
 * If \f$P = (\mathcal{p},\mathcal{cp},\mathcal{r},\mathcal{l})\f$ then
 * \f$T(P) = (A\mathcal{p} + b,A\mathcal{cp} + b, A\mathcal{r}, 
 * A\mathcal{l})\f$.
 */
template <typename S>
class AriadnePolyAffineMap {
	
	public:	
		typedef typename S::Real Real;
		typedef S State;
		typedef typename Ariadne::Geometry::AriadnePolyhedron< S > BasicSet;
		typedef typename Ariadne::Geometry::AriadnePolyhedron< S > Polyhedron;
		typedef typename Ariadne::Geometry::AriadneDenotableSet< Polyhedron > DenotableSet;
	
		typedef typename boost::numeric::ublas::identity_matrix<Real> IdentityMatrix;
		typedef typename boost::numeric::ublas::matrix<Real> Matrix;
		typedef typename boost::numeric::ublas::vector<Real> Vector;
	
		AriadnePolyAffineMap() {}
	
		AriadnePolyAffineMap(const AriadnePolyAffineMap<S> &map): 
				_A(map._A), _A_invertible(map._A_invertible), _b(map._b), 
				_map_type(map._map_type), _E(map._E), _n1(map._n1), _f(map._f){}
	
		AriadnePolyAffineMap(const Matrix &A):_A(A), _A_invertible(true),
				_b(A.size1()), _map_type(HOMOGENEOUS) {
			if ((A.size1()==0)||(A.size2()==0))
				throw std::invalid_argument("The matrix of affine trasformation should have at least 1 row and 1 column.");			
			
			this->_get_solution(A);
		}
	
		AriadnePolyAffineMap(const Vector &b):_A(b.size(),b.size()), 
				_A_invertible(true),_b(b), _map_type(TRASLATION) {
			if (b.size()==0)
				throw std::invalid_argument("The vector of affine trasformation should have at least 1 row.");
					
			this->_get_solution(b);
		}
		
		AriadnePolyAffineMap(const Matrix &A, const Vector &b):
			_A(A), _A_invertible(true), _b(b), _map_type(NON_HOMOGENEOUS) {
			
			if ((A.size1()==0)||(A.size2()==0))
				throw std::invalid_argument("The matrix of affine trasformation should have at least 1 row and 1 column.");
			
			if (b.size()==0)
				throw std::invalid_argument("The vector of affine trasformation should have at least 1 row.");

			this->_get_solution(A,b);		
		}
		
		inline size_t column_nb() const {
				return this->_A.size2();			
		}
		
		inline size_t row_nb() const {
			switch(this->_map_type) {
				case TRASLATION:
					return this->_b.size();
					break;
				default:
					return this->_A.size1();
			}				
		}
		
		inline size_t dim() const {
				return this->column_nb();			
		}
		
		inline AriadnePolyAffineMap<S>  &operator=(
					const AriadnePolyAffineMap<S> &orig) {
			
			this->_A=orig._A;
			this->_A_invertible=orig._A_invertible;
			this->_b=orig._b;
			this->_map_type=orig._map_type;
			this->_E=orig._E;
			this->_n1=orig._n1;
			this->_f=orig._f;
						
			return *this;

		}
		
		inline Polyhedron operator()(const Polyhedron &A) const {
			
			if (this->_A_invertible) {
				return _A_not_invertible_apply(A);
			}
			else {
				return _A_not_invertible_apply(A);
			}
		}
		
		inline const Matrix& A() const { return _A; }
		
		inline const Vector& b() const { return _b; }
		
		inline AffineKind affine_kind() { return _map_type; }
		
		inline bool invertible() const {return this->_A_invertible;}
		
	private:
		/*! \brief The matrix \f$A\f$ */
		Matrix _A;
	
		/*! \brief The invertibility flag of matrix \f$A\f$ */
		bool _A_invertible;
		
		/*! \brief The vector \f$b\f$ */
		Vector _b;
	
		/*! \brief Type of affine trasformation.*/
		AffineKind _map_type;
	
		/*! \brief The matrix \f$E\f$.*/	
		Matrix _E;
	
		/*! \brief The value \f$n_1\f$. */
		Real _n1;
		
		/*! \brief The vector \f$\vec{f}\f$. */	
		Vector _f;
		
		inline Polyhedron _A_not_invertible_apply(const Polyhedron &A) const {
	
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif	
			
			Ariadne::LinearAlgebra::AriadneGeneratorSystem<Real> gen;	
			
			Ariadne::Geometry::_extract_generators_from_polyhedron(A, gen);

			gen.point = prod (this->_A, gen.point);
			gen.sum_vector_to_all_points(this->_b);
			
			gen.ray = prod (this->_A, gen.ray);
			
			Polyhedron new_poly=Ariadne::Geometry::_create_polyhedron_from_generators(
											A,gen);
											
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif		
			
			return new_poly;
		}
	
		
		inline Polyhedron _A_invertible_apply(const Polyhedron &A) const {
			
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif	
	
			Ariadne::LinearAlgebra::AriadneConstrainSystem<Real> cs;	
			
			Ariadne::Geometry::_extract_constraints_from_polyhedron(A, cs);
			
			switch(this->_map_type) {
				case HOMOGENEOUS:
				
					this->_apply_homogeneous(cs);
				
					break;
				
				case NON_HOMOGENEOUS:
					
					this->_apply_non_homogeneous(cs);				
				
					break;
				
				case TRASLATION:	
					
					this->_apply_traslation(cs);
				
					break;
				
				default:
					throw std::invalid_argument("There is a problem in the application of AriadnePolyAffineMap.");
				
			}	

			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif	
			
			return Ariadne::Geometry::_create_polyhedron_from_constraints(
											A,cs);
		}
	
		inline void _apply_homogeneous(
			Ariadne::LinearAlgebra::AriadneConstrainSystem<Real> &cs) const{
			
			cs.d=this->_n1 * cs.d;
				
			cs.C=boost::numeric::ublas::prod(cs.C, this->_E);
			
		}
		
		inline void _apply_non_homogeneous(
			Ariadne::LinearAlgebra::AriadneConstrainSystem<Real> &cs) const{
			
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
			
			cs.d=this->_n1 * cs.d;

			/* an optimizated version of d= C f + d */
			//boost::numeric::ublas::axpy_prod(cs.C, this->_f,cs.d, false);
				
			cs.d+= boost::numeric::ublas::prod(cs.C, this->_f);
				
			cs.C=boost::numeric::ublas::prod(cs.C, this->_E);
				
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
		}
		
		inline void _apply_traslation(
			Ariadne::LinearAlgebra::AriadneConstrainSystem<Real> &cs) const {
		
			
			cs.d=this->_n1 * cs.d;
			
			/* an optimizated version of d= C f + d */
			//boost::numeric::ublas::axpy_prod(cs.C, this->_f,cs.d, false);
				
			cs.d+= boost::numeric::ublas::prod(cs.C, this->_f);
							
		}
			
		inline bool _solution_calculated() {
			return ((this->_E.size1()!=0)||(!this->_A_invertible));
		}
	
		inline void _empty_solution() {
			
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
			
			Matrix empty_matrix(0,0);
			Vector empty_vector(0);
				
			this->_E=empty_matrix;
			this->_f=empty_vector;
			
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
			
		}
	
		inline void _get_solution(const Matrix &A) {
			
			try {
				
				Matrix inv_A=Ariadne::LinearAlgebra::invert_matrix(A);
				
				Real n3=Ariadne::LinearAlgebra::find_matrix_denumerator(inv_A);
			
				this->_E=n3*inv_A;

				this->_n1=n3;

				/* useless, but safer */
				this->_f=Ariadne::LinearAlgebra::null_vector<Real>(A.size1());;
						
			} catch (std::exception &e) {
				
				this->_A_invertible=false;

			}
			
			
		}
		
		inline void _get_solution(const Vector &b) {
			
			Real n2=Ariadne::LinearAlgebra::find_vector_denumerator(b);
			
			this->_f=n2 * b;
			
			this->_n1=n2;
			
			/* E=n2*I    useless, but safer */
			this->_E=n2*Ariadne::LinearAlgebra::identity_matrix<Real>(b.size());	
	
		}
		
		inline void _get_solution(const Matrix &A, const Vector &b) {
			
			try {
				
				Matrix inv_A=Ariadne::LinearAlgebra::invert_matrix(A);
				
				Real n3=Ariadne::LinearAlgebra::find_matrix_denumerator(inv_A);
				Real n2=Ariadne::LinearAlgebra::find_vector_denumerator(b);
			
				this->_n1=n2 * n3;
				this->_E=this->_n1*inv_A;
				this->_f= prod(this->_E, b);
				
			} catch (std::exception &e) {
				this->_A_invertible=false;
			}
			
		}	

};

}
}
}
#endif /* _POLY_MAP_H */
