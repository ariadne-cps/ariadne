/***************************************************************************
 *            polyhedron.h
 *
 *  Thu Jan 27 10:26:36 2005
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
 
#ifndef _POLYHEDRON_H
#define _POLYHEDRON_H

#include<ppl.hh>
#include<iostream>
#include<vector>

#include <rectangle.h>


namespace Ariadne {	
namespace Geometry {
	
template <typename S> class Polyhedron;

}
}

#include <poly_map.h>
#include <constraint.h>
#include <generator.h>
#include <set_const.h>
#include <polyhedron_io.h>


namespace Ariadne {	
namespace Geometry {

template <typename S>
inline void _extract_constraints_from_polyhedron(
			const Polyhedron<S> &poly, 
			Ariadne::LinearAlgebra::ConstrainSystem<typename S::Real> &cs){
				
	Parma_Polyhedra_Library::Constraint_System::const_iterator j_cs, begin, end;
			
	begin=((poly._poly).constraints()).begin();
	end=((poly._poly).constraints()).end();
						
	size_t space_dim=poly.dim();
	size_t nb_constr=0,i,j;
						
	for (j_cs=begin; j_cs!=end; j_cs++) {
		nb_constr++;
	}
	
	cs.reset_dimensions(space_dim, nb_constr);
	
	j=0;
	
	for (j_cs=begin; j_cs!=end; j_cs++) {
					
		for (i=0; i< space_dim; i++) {
				cs.C(j,i)=j_cs->coefficient(Parma_Polyhedra_Library::Variable(i));
		}
					
		cs.d(j)=-j_cs->inhomogeneous_term();
		
		if (j_cs->is_equality()) {
			cs.r_type[j]=Ariadne::LinearAlgebra::ConstrainSystem<typename S::Real>::EQUALITY;
		} else {
			if (j_cs->is_strict_inequality()) {
				cs.r_type[j]=Ariadne::LinearAlgebra::ConstrainSystem<typename S::Real>::STRICT_INEQUALITY;
			} else {
				cs.r_type[j]=Ariadne::LinearAlgebra::ConstrainSystem<typename S::Real>::INEQUALITY;
			}
		}
					
		j++;
	}

}

template <typename S>
inline void _extract_generators_from_polyhedron(
			const Polyhedron<S> &poly, 
			Ariadne::LinearAlgebra::GeneratorSystem<typename S::Real> &gen){

	#ifdef DEBUG
		std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	#endif		

	typedef typename S::Real Real;
				
	Parma_Polyhedra_Library::Generator_System::const_iterator j_gen, begin, end;
				
	begin=((poly._poly).generators()).begin();
	end=((poly._poly).generators()).end();
						
	size_t space_dim=poly.dim();
	size_t point_nb=0,ray_nb=0, i,j_ray=0, j_point=0;
						
	for (j_gen=begin; j_gen!=end; j_gen++) {
		
		if ((j_gen->is_line())||(j_gen->is_ray())) {
			ray_nb++;
		} else {
			point_nb++;	
		}
	}
	
	gen.reset_dimensions(point_nb,ray_nb,space_dim);

	for (j_gen=begin; j_gen!=end; j_gen++) {
		
		
		if (j_gen->is_line()) {
		
			for (i=0; i< space_dim; i++) {
				gen.ray(i,j_ray)=j_gen->coefficient(Parma_Polyhedra_Library::Variable(i));
			}
					
			gen.set_line(j_ray);
			
			j_ray++;
			
		} else {
			if (j_gen->is_ray()) {
		
				for (i=0; i< space_dim; i++) {
					gen.ray(i,j_ray)=j_gen->coefficient(Parma_Polyhedra_Library::Variable(i));
				}
					
				gen.set_ray(j_ray);
			
				j_ray++;
				
			} else {
				if (j_gen->is_point()) {
		
					Real den=j_gen->divisor();
					
					for (i=0; i< space_dim; i++) {
						gen.point(i,j_point)=j_gen->coefficient(Parma_Polyhedra_Library::Variable(i))/den;
					}
					
					gen.set_point(j_point);
			
					j_point++;
				
				} else { /*j_gen is a closure point */
					
					Real den=j_gen->divisor();
					
					for (i=0; i< space_dim; i++) {
						gen.point(i,j_point)=j_gen->coefficient(Parma_Polyhedra_Library::Variable(i))/den;
					}
					
					gen.set_closure_point(j_point);
			
					j_point++;
					
				}
			}
		}
				
	}
	
	#ifdef DEBUG
		std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	#endif	
}


template <typename S>
inline Polyhedron<S> _create_polyhedron_from_generators(
			const Polyhedron<S> &poly,
			const Ariadne::LinearAlgebra::GeneratorSystem<typename S::Real> &gen){
	
	#ifdef DEBUG
		std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	#endif			

	typedef typename S::Real Real;
				
	Real den;
				
	Parma_Polyhedra_Library::Generator_System poly_gen;
	Parma_Polyhedra_Library::Linear_Expression exp;
	size_t j,i;
					
	if (poly.dim()!=gen.space_dim()) 
		throw std::invalid_argument("The generator system's matrix and the variable vector should have the same dimensions.");
	
	for (j=0; j< gen.point_nb(); j++) {
		
		den=denumerator(gen.point(0,j));
		
		for (i=1; i< poly.dim(); i++) {
			den=lcm(numerator(den),denumerator(gen.point(i,j)));
		}
		
		exp=(den*gen.point(0,j))*(Parma_Polyhedra_Library::Variable(0));

		for (i=1; i< poly.dim(); i++) {
			exp+=(den*gen.point(i,j))*(Parma_Polyhedra_Library::Variable(i));
		}
		
		if (gen.is_point(j)) {		
			poly_gen.insert(Parma_Polyhedra_Library::Generator::point(exp,den));
		} else {
			poly_gen.insert(Parma_Polyhedra_Library::Generator::closure_point(exp,den));
		}
	}
	
	for (j=0; j< gen.ray_nb(); j++) {
		
		exp=gen.ray(0,j)*(Parma_Polyhedra_Library::Variable(0));

		for (i=1; i< poly.dim(); i++) {
			exp+=gen.ray(i,j)*(Parma_Polyhedra_Library::Variable(i));
		}
		
		if (gen.is_ray(j)) {		
			poly_gen.insert(Parma_Polyhedra_Library::Generator::ray(exp));
		} else {
			poly_gen.insert(Parma_Polyhedra_Library::Generator::line(exp));
		}
	}
	
	Polyhedron<S> new_poly(poly_gen);
	
	#ifdef DEBUG
		std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	#endif

	return new_poly;
}


template <typename Real >
inline Parma_Polyhedra_Library::NNC_Polyhedron 
		_create_open_PPL_poly_from_constraints(
			const Ariadne::LinearAlgebra::ConstrainSystem<Real> &cs){
	
	Parma_Polyhedra_Library::Constraint_System poly_cs;
	Parma_Polyhedra_Library::Linear_Expression exp;
	
	size_t i,j;
				
	/* if there is an equality, the open set is empty */
	if (cs.any_equality()) {
		Parma_Polyhedra_Library::NNC_Polyhedron 
			open_poly(cs.space_dim(),Parma_Polyhedra_Library::Polyhedron::EMPTY);
		
		return open_poly;
	}
	
	for (j=0; j< cs.nb_constraints(); j++) {
				
		exp=cs.C(j,0)*(Parma_Polyhedra_Library::Variable(0));

		for (i=1; i< cs.space_dim(); i++) {
			exp+=cs.C(j,i)*(Parma_Polyhedra_Library::Variable(i));
		}
		
		poly_cs.insert(exp > cs.d(j));
		
	}
	
	Parma_Polyhedra_Library::NNC_Polyhedron poly(poly_cs);
	
	return poly;
}

template <typename Real >
inline Parma_Polyhedra_Library::NNC_Polyhedron 
		_create_closed_PPL_poly_from_constraints(
			const Ariadne::LinearAlgebra::ConstrainSystem<Real> &cs){
	
	Parma_Polyhedra_Library::Constraint_System poly_cs;
	Parma_Polyhedra_Library::Linear_Expression exp;
	
	size_t i,j;
	
	for (j=0; j< cs.nb_constraints(); j++) {
				
		exp=cs.C(j,0)*(Parma_Polyhedra_Library::Variable(0));

		for (i=1; i< cs.space_dim(); i++) {
			exp+=cs.C(j,i)*(Parma_Polyhedra_Library::Variable(i));
		}
		
		if (cs.is_equality(j)) {
			poly_cs.insert(exp == cs.d(j));
		} else {
			poly_cs.insert(exp >= cs.d(j));
		}
	}
	
	Parma_Polyhedra_Library::NNC_Polyhedron poly(poly_cs);
	
	return poly;
}

template <typename S>
inline Polyhedron<S> _create_polyhedron_from_constraints(
			const Polyhedron<S> &poly,
			const Ariadne::LinearAlgebra::ConstrainSystem<typename S::Real> &cs){
	
	#ifdef DEBUG
		std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	#endif			
				
	Parma_Polyhedra_Library::Constraint_System poly_cs;
	Parma_Polyhedra_Library::Linear_Expression exp;
	size_t j,i;
					
	if (poly.dim()!=cs.space_dim()) 
		throw std::invalid_argument("The constraint system's matrix and the variable vector should have the same dimensions.");
	
	for (j=0; j< cs.nb_constraints(); j++) {
				
		exp=cs.C(j,0)*(Parma_Polyhedra_Library::Variable(0));

		for (i=1; i< poly.dim(); i++) {
			exp+=cs.C(j,i)*(Parma_Polyhedra_Library::Variable(i));
		}
		
		if (cs.is_equality(j)) {		
			poly_cs.insert(exp == cs.d(j));
		} else {
			poly_cs.insert(exp >= cs.d(j));
		}
	}
			
	Polyhedron<S> new_poly(poly_cs);

	#ifdef DEBUG
		std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	#endif

	return new_poly;
}


template <typename S>
inline Polyhedron<S> convex_hull(const Polyhedron<S> &A, 
						const Polyhedron<S> &B){
			
	#ifdef DEBUG
		std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	#endif
							
	Polyhedron<S> chull(A);
							
	(chull._poly).poly_hull_assign_and_minimize(B._poly);
	chull._evaluate_interior();
							
	#ifdef DEBUG
		std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	#endif
							
	return chull;
}
	
/*! \brief Transforms an State into a Parma_Polyhedra_Library::C_Polyhedron */
template <typename S>
inline Parma_Polyhedra_Library::NNC_Polyhedron _from_State_to_PPL_Polihedron(const S &s){
	
	Parma_Polyhedra_Library::Constraint_System cs;
	Rational num;
	
	for (size_t i=0; i<s.dim(); i++) {
			num=transform_into_rational(s[i]);
			
			cs.insert(Parma_Polyhedra_Library::Variable(i) * 
							denumerator(num) == numerator(num));
	}
	
	Parma_Polyhedra_Library::NNC_Polyhedron p(cs);
	
	return p;
	
}


/*! \brief Transforms an Rectangle into a Parma_Polyhedra_Library::C_Polyhedron */
template <typename R>
inline Parma_Polyhedra_Library::NNC_Polyhedron _from_Rectangle_to_closed_PPL_Polihedron(const R &r){
	
	typedef typename R::State State;
	
	State u_corner(r.upper_corner()), l_corner(r.lower_corner());
	
	Parma_Polyhedra_Library::Constraint_System cs;
  	std::vector<Parma_Polyhedra_Library::Variable *> x;
	Rational num;
	
	for (size_t i=0; i<r.dim(); i++) {
			x.push_back(new Parma_Polyhedra_Library::Variable(i));
	}
	
	for (size_t i=0; i<r.dim(); i++) {
			num=transform_into_rational(u_corner[i]);
			cs.insert((Parma_Polyhedra_Library::Variable(i))*denumerator(num) <= 
						numerator(num));
		
			num=transform_into_rational(l_corner[i]);
			cs.insert((Parma_Polyhedra_Library::Variable(i))*denumerator(num) >=
						numerator(num));
	}
	
	Parma_Polyhedra_Library::NNC_Polyhedron p(cs);
	
	return p;
	
}


template <typename S>
inline bool disjoint(const Polyhedron<S> &A, 
				const Polyhedron<S> &B){
							
	return (A._poly).is_disjoint_from(B._poly);
									
}

template <typename S>
inline bool intersects_interior(const Polyhedron<S> &A, 
				const Polyhedron<S> &B){
	
	#ifdef DEBUG
		std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	#endif
					
	return !((A._poly).is_disjoint_from(B._interior_poly));				
}

template <typename S >
inline bool interiors_intersect(const Rectangle< S > &rect, 
				const Polyhedron< S > &A) {
	
	return interiors_intersect(A,rect);
}

template <typename S >
inline bool interiors_intersect(const Polyhedron< S > &A,
				const Rectangle< S > &rect) {
	
	Parma_Polyhedra_Library::NNC_Polyhedron P_rect=_from_Rectangle_to_open_PPL_Polihedron(rect);

	return !((A._poly).is_disjoint_from(P_rect));
									
}

template <typename S>
inline bool subset_of_interior(const Polyhedron<S> &A, 
				const Polyhedron<S> &B){				
		
	return (B._interior_poly).contains(A._poly);
}

template <typename S >
inline bool interior_subset_of_interior(const Rectangle< S > &rect, 
				const Polyhedron< S > &A) {
	
	Parma_Polyhedra_Library::NNC_Polyhedron P_rect=_from_Rectangle_to_open_PPL_Polihedron(rect);					
										
	return (A._interior_poly).contains(P_rect);
}

template <typename S >
inline bool interior_subset_of_interior(const Polyhedron< S > &A,
					const Rectangle< S > &rect) {

	Parma_Polyhedra_Library::NNC_Polyhedron P_rect=_from_Rectangle_to_open_PPL_Polihedron(rect);					
										
	return (P_rect).contains(A._interior_poly);
}


template <typename S>
inline bool subset_of_open_cover(const Polyhedron<S> &A, 
				const std::vector< Polyhedron<S> > &vector){
					
		
	typename std::vector< Polyhedron <S> >::const_iterator i;

	/* TO REIMPLEMENT */
	for (i = vector.begin(); i != vector.end(); i++) {

		if (subset_of_interior(A,*i)) return true;
		
	}

	return false;
}


template <typename S>
inline Polyhedron<S> closure_of_intersection_of_interior (
				const Polyhedron<S> &A, const Polyhedron<S> &B) {
					
	#ifdef DEBUG
		std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	#endif
					
	Polyhedron<S> new_poly(A);				
		
	(new_poly._poly).intersection_assign(B._poly);
	(new_poly._interior_poly).intersection_assign(B._interior_poly);

	#ifdef DEBUG
		std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	#endif
	
	return (new_poly);
}

/* TO REIMPLEMENT */
template <typename S>
inline Polyhedron<S> minkowski_sum(const Polyhedron<S> &A,
			const Polyhedron<S> &B) {
		
	typedef typename S::Real Real;
	
	#ifdef DEBUG
		std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	#endif
	
	if ((!(A._poly).is_bounded()) || (!(B._poly).is_bounded())) {
		throw std::domain_error("Minkosky sum is implemented for bounded polyhedra only");
	}	
				
	Ariadne::LinearAlgebra::GeneratorSystem<Real> gen_A, gen_B;
	_extract_generators_from_polyhedron(A, gen_A);
				
	Polyhedron<S> sum=B, new_poly;
				
	for (size_t k=0; k<gen_A.point_nb(); k++) {
					
		_extract_generators_from_polyhedron(B, gen_B);
			
		for (size_t i=0; i<gen_B.space_dim(); i++) {
				
			for (size_t j=0; j<gen_B.point_nb() ; j++) {
						
				gen_B.point(i,j)+=gen_A.point(i,k);
			}
		}
			
		new_poly=_create_polyhedron_from_generators(B,gen_B);
				
		sum=convex_hull(sum,new_poly);
	}
	
	#ifdef DEBUG
		std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	#endif
	
	return sum;
}

/*! \brief The polyhedron class.
 *
 * This is a wrapper from Parma Polyhedra Library Polyhedron class
 * to Ariadne representation.
 */ 
template <typename S>
class Polyhedron {
	
	public:
		typedef S State;
		typedef typename State::Real Real;
	
		/*! \brief Tests disjointness. */
		friend bool disjoint <> (const Polyhedron< S > &A, 
				const Polyhedron< S > &B);

		/*! \brief Tests intersection of interiors. */
		friend bool intersects_interior <> (const Polyhedron< S > &A, 
				const Polyhedron< S > &B);
	
		/*! \brief Tests intersection of interiors. */	
		//friend bool intersects_interior <> (const Rectangle< S > &rect, 
		//					const Polyhedron<S> &A);
		
		/*! \brief Tests intersection of interiors. */	
		friend bool interiors_intersect <> (const Polyhedron<S> &A,
							const Rectangle< S > &rect);
		
		/*! \brief Tests inclusion of interiors. */
		friend bool subset_of_interior <> (const Polyhedron< S > &A, 
				const Polyhedron< S > &B);
	
		/*! \brief Tests inclusion of interiors. */
		friend bool interior_subset_of_interior <> (const Rectangle< S > &rect, 
				const Polyhedron< S > &A);

		/*! \brief Tests inclusion of interiors. */
		friend bool interior_subset_of_interior <> (const Polyhedron< S > &A,
				const Rectangle< S > &rect);
				
		/*! \brief Tests inclusion. */
		friend bool subset_of_open_cover <> (const Polyhedron< S > &A, 
				const std::vector< Polyhedron< S > > &vector);

		/*! \brief Makes intersection of interiors. */
		friend Polyhedron< S > closure_of_intersection_of_interior <> (
				const Polyhedron< S > &A, const Polyhedron< S > &B);

		/*! \brief Prints polyhedron. */
		template <typename SType >
		friend std::ostream& IO_Operators::operator<<(std::ostream &os, 
										const Polyhedron< SType > &r);

		template <typename STATE>
		friend class IO_Operators::PolyhedronMatlabExporter;
		
		template <typename STATE>
		friend class IO_Operators::PolyhedronOneDimMatlabExporter;
	
		friend void _extract_constraints_from_polyhedron <> (
						const Polyhedron< S > &poly, 
						Ariadne::LinearAlgebra::ConstrainSystem< typename S::Real > &cs);
		
		friend void _extract_generators_from_polyhedron <> (
						const Polyhedron< S > &poly, 
						Ariadne::LinearAlgebra::GeneratorSystem< typename S::Real > &gen);
		
		friend Polyhedron<S> _create_polyhedron_from_constraints <> (
			const Polyhedron<S> &poly,
			const Ariadne::LinearAlgebra::ConstrainSystem<typename S::Real> &cs);
		
		friend Polyhedron<S> _create_polyhedron_from_generators <> (
			const Polyhedron<S> &poly,
			const Ariadne::LinearAlgebra::GeneratorSystem<typename S::Real> &gen);
		
		friend Polyhedron< S > convex_hull <> (
						const Polyhedron< S > &A, 
						const Polyhedron< S > &B);
		
		friend Polyhedron< S > minkowski_sum <> (
						const Polyhedron< S > &A, 
						const Polyhedron< S > &B);
						
		friend class Ariadne::Map::Affine::PolyAffineMap<S>;
		
		inline Polyhedron<S> _interior() const {
			
			if (this->_interior_poly.space_dimension()!=0) {
				Polyhedron<S> interior(*this);
				
				interior._poly=this->_interior_poly;
				
				return interior;
			}
			
			Ariadne::LinearAlgebra::ConstrainSystem<Real> cs;
			
			_extract_constraints_from_polyhedron(*this, cs);

			// FIXME! This method doesn't exist!			
			return _create_open_polyhedron_from_constraints(*this, cs);
			
		}
			
	private:
		
		inline void _evaluate_interior(){			
			
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
			
			Ariadne::LinearAlgebra::ConstrainSystem<Real> cs;
			
			_extract_constraints_from_polyhedron(*this, cs);
			
			this->_interior_poly=_create_open_PPL_poly_from_constraints(cs);
			
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
		}			
		
		inline void _set_empty() {
	
			Parma_Polyhedra_Library::NNC_Polyhedron 
							poly(this->dim(),
							Parma_Polyhedra_Library::Polyhedron::EMPTY);
	
			this->_poly=poly;
			this->_interior_poly=poly;
		}
		
		/*! \brief The polyhedron data structure. */
		Parma_Polyhedra_Library::NNC_Polyhedron _poly; //(thanks Parma, I LOVE you! :-) )
	
		Parma_Polyhedra_Library::NNC_Polyhedron _interior_poly;
	public:
		/*! \brief A polyhedron constructor. 
		 *
		 * Builds a polyhedron.
		 * \param kind is the kind of the new polyhedron.
		 */
		Polyhedron(DegenerateSetKind kind = EMPTY): 
			_poly(Parma_Polyhedra_Library::Polyhedron::EMPTY),
			_interior_poly(Parma_Polyhedra_Library::Polyhedron::EMPTY) {
			
			if (kind==UNIVERSE) {
					Parma_Polyhedra_Library::NNC_Polyhedron 
							poly(Parma_Polyhedra_Library::Polyhedron::UNIVERSE);
				
					this->_poly=poly;
					this->_interior_poly=poly;
			}
			

		}			
		
		/*! \brief A polyhedron constructor. 
		 *
		 * Builds a polyhedron with dimension \a dim.
		 * \param dim is the new polyhedron's number of dimensions.
		 * \param kind is the kind of the new polyhedron.
		 */
		Polyhedron(size_t dim, DegenerateSetKind kind = EMPTY): 
			_poly(dim, Parma_Polyhedra_Library::Polyhedron::EMPTY), 
			_interior_poly(dim, Parma_Polyhedra_Library::Polyhedron::EMPTY) {
				
			if (kind==UNIVERSE) {
					Parma_Polyhedra_Library::NNC_Polyhedron 
							poly(Parma_Polyhedra_Library::Polyhedron::UNIVERSE);
				
					this->_poly=poly;
					this->_interior_poly=poly;
			}	
		}
			
		/*! \brief A polyhedron constructor. 
		 *
		 * Builds a hypercube of dimension \f$dim\f$ centered in the 
		 * origin.
		 * \param dim is the new polyhedron's space dimension.
		 * \param size is the dimension per space dimension.
		 */
		Polyhedron(const size_t &dim, const Real &size) {
			
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
			
			Parma_Polyhedra_Library::Constraint_System cs;
			Real num=numerator(abs(size/2)),den=denumerator(abs(size/2));
			
			for (size_t i=0; i< dim; i++) {
				cs.insert( den* Parma_Polyhedra_Library::Variable(i) >=  -num );
				cs.insert( den* Parma_Polyhedra_Library::Variable(i) <=  num );
			}
			
			Parma_Polyhedra_Library::NNC_Polyhedron new_poly(cs);
			
			this->_poly=new_poly;	
			this->_evaluate_interior();

			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif

		}
		
		/*! \brief A polyhedron constructor. 
		 *
		 * Builds a polyhedron from an other polyhedron.
		 * \param original is the original polyhedron.
		 */
		Polyhedron(const Polyhedron<State> &original): 
			_poly(original._poly), _interior_poly(original._interior_poly){}
				
		/*! \brief A polyhedron constructor. 
		 *
		 * Builds a polyhedron from constraint system.
		 * \param cs is the constraint system.
		 */
		Polyhedron(Parma_Polyhedra_Library::Constraint_System &cs): 
			_poly(cs){
				
			this->_evaluate_interior();
		}
		
		/*! \brief A polyhedron constructor. 
		 *
		 * Builds a polyhedron from generator system.
		 * \param gen is the generator system.
		 */
		Polyhedron(Parma_Polyhedra_Library::Generator_System &gen): 
			_poly(gen){
				
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
				
			this->_evaluate_interior();
				
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif	
		}		
		
		/*! \brief A polyhedron constructor. 
		 *
		 * Builds a polyhedron from a rectangle.
		 * \param rect is a rectangle.
		 */
		template<typename R>
		Polyhedron(const R& rect){
			_poly=_from_Rectangle_to_closed_PPL_Polihedron(rect);
			
			this->_evaluate_interior();
		}

		/*! \brief Returns polyhedron's space dimensions.
		 *
		 * \return The polyhedron's space dimensions.
		 */
		inline size_t dim() const {
			return ((size_t)this->_poly.space_dimension());
		}
		
		/*! \brief Checks the emptyness.
		 *
		 * \return \a true if the polyhedron is empty,
		 * \a false otherwise.		
		 */
		inline bool empty() const {
			return (this->_poly.is_empty());
		}
		
		/*! \brief Tests if a state is included into a polyhedron.
		 *
		 * \param state is a state in the polyhedron's space.
		 * \return \a true if the state is contained into 
		 * the current polyhedron, \a false otherwise.
		 */
		inline bool contains(const State& state) const {
			
			if (state.dim()!=this->dim()) 
				throw std::domain_error("This object and parameter have different space dimentions");
			
			Parma_Polyhedra_Library::NNC_Polyhedron p(this->_poly);
			
			p.topological_closure_assign();
			
			return p.contains(_from_State_to_PPL_Polihedron(state));
			
		}
		
		/*! \brief Tests if a state is included into the interior of a polyhedron.
		 *
		 * \param state is a state in the polyhedron's space.
		 * \return \a true if the state is contained into the interior of 
		 * the current polyhedron, \a false otherwise.
		 */
		inline bool interior_contains(const State& state) const {
			
			if (state.dim()!=this->dim()) 
				throw std::domain_error("This object and parameter have different space dimentions");
			
			return this->_poly.contains(_from_State_to_PPL_Polihedron(state));
			
		}
		
		
		
		inline Polyhedron<State> operator+(const Polyhedron<State>& A) const{
			
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif		
			
			return minkowski_sum(*this,A);
	
		}
		
		inline Polyhedron<State> &expand_by(const Real &delta) {
			
			Ariadne::LinearAlgebra::ConstrainSystem<typename S::Real> cs;
			
			_extract_constraints_from_polyhedron(*this, cs);
			
			if (cs.any_equality()) {
				this->_set_empty();
				return *this;				
			}
			
			cs.expand_by(delta);
			
			this->_poly=_create_closed_PPL_poly_from_constraints(cs);
			this->_interior_poly=_create_open_PPL_poly_from_constraints(cs);
			
			return *this;
		}
		
		/*! \brief Copies a polyhedron. 
		 *
		 * \param original is the original polyhedron.
		 * \return A reference to the current object.
		 */
		inline Polyhedron<State> &operator=(
					const Polyhedron<State> &original) {				
						
			this->_poly=original._poly;
			this->_interior_poly=original._interior_poly;
						
			return *this;
		}
		
		inline Polyhedron<State> &set_precision_to_upperapproximating(const Real &delta) {
		
			Real denum=denumerator(delta);

			Ariadne::LinearAlgebra::ConstrainSystem<typename S::Real> cs;
			
			_extract_constraints_from_polyhedron(*this, cs);
			
			if (cs.already_at_precision(denum)) return *this;
		
			cs.reduce_precision_to_expanding(denum);
			
			Parma_Polyhedra_Library::NNC_Polyhedron new_poly=_create_open_PPL_poly_from_constraints(cs);
			
			if (!new_poly.is_empty()) {
					
				this->_poly=_create_closed_PPL_poly_from_constraints(cs);
				this->_interior_poly=new_poly;
			}
			
			return *this;
						
		}
	
		inline Polyhedron<State> &set_precision_to_upperapproximating_for_output(const Real &delta) {
		
			Real denum=denumerator(delta);

			Ariadne::LinearAlgebra::ConstrainSystem<typename S::Real> cs;
			
			_extract_constraints_from_polyhedron(*this, cs);
			
			if (cs.already_at_precision(denum)) return *this;
		
			cs.reduce_precision_to_expanding(denum);
			
			this->_poly=_create_closed_PPL_poly_from_constraints(cs);
			this->_interior_poly=_create_open_PPL_poly_from_constraints(cs);
			
			return *this;
						
		}

		
		inline Polyhedron<State> project_on_dimentions(const std::vector<uint> &dims) const {
			
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif		
			
			if (dims.size()==0) {
				throw "I can not project on zero dimensions.";	
			}
			if (dims.size()>this->dim()) {
				throw "I can not project on more dimensions than the polyhedron ones.";
			}
			
			
			boost::numeric::ublas::matrix<Real> projection_map=
									Ariadne::LinearAlgebra::zero_matrix<Real>(this->dim());
			
			for (size_t i=0; i< dims.size(); i++) {
				projection_map(dims[i],dims[i])=1.0;
			}
			
			Ariadne::Map::Affine::PolyAffineMap<State> projection(projection_map);
			
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif	
			
			return projection(*this);
			
		}
};

}
}
#endif /* _POLYHEDRON_H */
