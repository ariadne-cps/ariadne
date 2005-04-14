/***************************************************************************
 *            denotable_set.h
 *
 *  Thu Sep 24 17:01:05 2004
 *  Copyright  2004  Alberto Casagrande, Pieter Collins
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
 
#ifndef _DENOTABLE_SET_H
#define _DENOTABLE_SET_H

#include <vector>

#include <set_const.h>
#include <denotableset_io.h>

namespace Ariadne {
namespace Geometry {

template <typename BS, uint BS_PER_BLOCK >
bool disjoint (const AriadneDenotableSet< BS , BS_PER_BLOCK > &A, 
				const AriadneDenotableSet< BS , BS_PER_BLOCK > &B){

	if (A.dim()!=B.dim()) 
		throw std::invalid_argument("The two denotable set have different space dimensions.");
	
					
	size_t i,j;
					
	for (i=0; i<A.size() ; i++) {
	
		for (j=0; j<B.size() ; j++) {
	
				if (!disjoint(A[i],B[j])) return false;
		
		}	
		
	}
	
	return true;
						
}

template <typename BS, uint BS_PER_BLOCK >
bool intersects_interior(const AriadneDenotableSet< BS , BS_PER_BLOCK > &A, 
				const AriadneDenotableSet< BS , BS_PER_BLOCK > &B){
	
	if (A.dim()!=B.dim()) 
		throw std::invalid_argument("The two denotable set have different space dimensions.");
	
					
	size_t i,j;
					
	for (i=0; i<A.size() ; i++) {
	
		for (j=0; j<B.size() ; j++) {
	
			if (intersects_interior(A[i],B[j])) return true;
		
		}	
		
	}				
	
	return false;
}
	
template <typename BS, uint BS_PER_BLOCK >
bool intersects_interior(const AriadneRectangle< typename BS::State > &rect, 
				const AriadneDenotableSet< BS , BS_PER_BLOCK > &A){
		
	if (A.dim()!=rect.dim()) 
		throw std::invalid_argument("The denotable set and the rectangle have different space dimensions.");
	
					
	for (size_t i=0; i<A.size() ; i++) {				

		if (intersects_interior(rect,A[i])) return true;
	}
	
	return false;
}
		
template <typename BS, uint BS_PER_BLOCK >
bool intersects_interior(const AriadneDenotableSet< BS , BS_PER_BLOCK > &A,
				const AriadneRectangle< typename BS::State > &rect){
		
	if (A.dim()!=rect.dim()) 
		throw std::invalid_argument("The denotable set and the rectangle have different space dimensions.");
	
					
	for (size_t i=0; i<A.size() ; i++) {				

		if (intersects_interior(A[i],rect)) return true;
	}
	
	return false;
}

template <typename BS, uint BS_PER_BLOCK >
bool subset_of_interior(const AriadneDenotableSet< BS , BS_PER_BLOCK > &A, 
				const AriadneDenotableSet< BS , BS_PER_BLOCK > &B){
	
	if (A.dim()!=B.dim()) 
		throw std::invalid_argument("The two denotable set have different space dimensions.");
					
					
	for (size_t i=0; i<A.size() ; i++) {
	
		if (!subset_of_open_cover(A[i], B._vector)) return false;
	}				
	
	return true;
						
}
	
template <typename BS, uint BS_PER_BLOCK >
bool subset_of_interior(const AriadneRectangle<typename BS::State > &rect, 
				const AriadneDenotableSet< BS , BS_PER_BLOCK > &A){

	#ifdef DEBUG
		std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	#endif
					
	if (A.dim()!=rect.dim()) 
		throw std::invalid_argument("The denotable set and the rectangle have different space dimensions.");
	

	/* TO REIMPLEMENT */
	for (size_t i=0; i<A.size() ; i++) {

		if (subset_of_interior(rect,A[i])) {
		
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif			
			
			return true;
		}
		
	}
	
	#ifdef DEBUG
		std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	#endif
	
	return false;					
}

template <typename BS, uint BS_PER_BLOCK >
bool subset_of_interior(const AriadneDenotableSet< BS , BS_PER_BLOCK > &A,
				const AriadneRectangle<typename BS::State > &rect){
	
	#ifdef DEBUG
		std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	#endif
					
	if (A.dim()!=rect.dim()) 
		throw std::invalid_argument("The denotable set and the rectangle have different space dimensions.");
	
					
	for (size_t i=0; i<A.size() ; i++) {

		if (!subset_of_interior(A[i], rect)) {
		
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
			
			return false;
		}
		
	}

	#ifdef DEBUG
		std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	#endif
	
	return true;					
}


template <typename BS, uint BS_PER_BLOCK >
AriadneDenotableSet< BS , BS_PER_BLOCK > join(
				const AriadneDenotableSet< BS , BS_PER_BLOCK > &A, 
				const AriadneDenotableSet< BS , BS_PER_BLOCK > &B) {

	#ifdef DEBUG
		std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	#endif
					
	if (A.dim()!=B.dim()) 
		throw std::invalid_argument("The two denotable set have different space dimensions.");
	
					
	AriadneDenotableSet< BS , BS_PER_BLOCK > ds_union(A);
					
	ds_union.inplace_union(B);

	#ifdef DEBUG
		std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	#endif
	
	return ds_union;
					
}

template <typename BS, uint BS_PER_BLOCK >
AriadneDenotableSet< BS , BS_PER_BLOCK > closure_of_intersection_of_interior(
				const AriadneDenotableSet< BS , BS_PER_BLOCK > &A, 
				const AriadneDenotableSet< BS , BS_PER_BLOCK > &B) {

	AriadneDenotableSet< BS , BS_PER_BLOCK > ds_inter;
	size_t i,j,&k=ds_inter._basic_sets;
	size_t new_size=((A.size()*B.size())/BS_PER_BLOCK)+1;
	std::vector<BS> &vector=ds_inter._vector;
	
	#ifdef DEBUG
		std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	#endif
					
	if (A.dim()!=B.dim()) 
		throw std::invalid_argument("The two denotable set have different space dimensions.");
			
					
	(vector).resize(new_size);
					
	for (i=0; i<A.size(); i++) {
	
		for (j=0; j<B.size(); j++) {

			if (intersects_interior(A[i],B[j])) {
				
				vector[k]=closure_of_intersection_of_interior(A[i],B[j]);
				
				k++;
			}
		}
	}
	
	ds_inter._dim=A._dim;

	#ifdef DEBUG
		std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	#endif
	
	return ds_inter;		
					
}
	
/*! \class DenotableSet
 * \brief It is a list of basic set. 
 */
template <class BS , size_t BS_PER_BLOCK= 2000>
class AriadneDenotableSet{
		
	private:
		/*! \brief The basic set vector */
		std::vector<BS> _vector;
	
		size_t _basic_sets;
	
		size_t _dim;
		
	public:
		typedef BS BasicSet;
		typedef typename BS::State State;
		typedef typename BS::State::Real Real;

		typedef typename std::vector<BS>::const_iterator const_iterator;
		typedef typename std::vector<BS>::iterator iterator;
	
		/*! \brief A denotable set constructor. */
		AriadneDenotableSet(): _basic_sets(0), _dim(0) {
			this->_vector.resize(BS_PER_BLOCK);	
		}
			
		/*! \brief A denotable set constructor. */
		AriadneDenotableSet(const BasicSet &A): _basic_sets(0), _dim(0) {
			
			this->_vector.resize(BS_PER_BLOCK);
			this->_dim=A.dim();
			
			if (A.empty()) return;
			
			this->_vector[0]=A;
			this->_basic_sets++;
			
		}
			
		/*! \brief A denotable set constructor. */
		AriadneDenotableSet(const AriadneDenotableSet< BasicSet , BS_PER_BLOCK> &A): 
					_basic_sets(0), _dim(A.dim()) {
						
			if (A.empty()) {
				(this->_vector).resize(BS_PER_BLOCK);
				
				return;
			}
			
			(this->_vector).resize((A._vector).size());
			
			for (size_t i=0; i< A.size() ; i++) {
				this->_vector[i]=A._vector[i];
			}
			
			this->_dim=A._dim;
			this->_basic_sets=A._basic_sets;
						
		}
		
		/*! \brief A denotable set constructor. */
		AriadneDenotableSet(size_t dim, 
				DegenerateSetKind kind = EMPTY):_basic_sets(0),_dim(dim) {
					
			if (kind!=EMPTY) {
				BasicSet bs(dim,kind);
			
				this->_vector.resize(BS_PER_BLOCK);

				this->_vector[0]=bs;
				this->_basic_sets++;
			}				
					
		}

		/*! \brief A denotable set constructor. */
		AriadneDenotableSet(const size_t &dim, const Real &size):_basic_sets(0),_dim(dim) {
					
			BasicSet bs(dim,size);
			
			this->_vector.resize(BS_PER_BLOCK);
			
			this->_vector[0]=bs;
			this->_basic_sets++;	
					
		}
		
		/*! \brief A denotable set destructor. */
		~AriadneDenotableSet() {
			this->_vector.clear();
		}
	
		/*! \brief Return the number of basic sets forming this object. 
		 *
		 * \return The number of basic sets forming this object.
		 */
		inline const size_t &size() const {
			return (this->_basic_sets);
		}
		
		/*! \brief Return the denotable set's space dimention. 
		 *
		 * \return The space dimention of the DenotableSet.
		 */
		inline const size_t &dim() const{
			return this->_dim;
		}
		
		/*! \brief Accesses the i-th BasicSet.
		 *
		 * \param index is the index of the returned basic set.
		 * \return The i-th basic set maitained by the DenotableSet.
		 */
		inline const BasicSet &operator[](size_t index) const{
			
			if (this->size()<=index) 
				throw std::invalid_argument("Index overlaps vector bounds.");
			
			return this->_vector[index];
		}
		
		inline AriadneDenotableSet<BS> operator+(const AriadneDenotableSet<BS>& A) const{
		
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
			
			AriadneDenotableSet<BS> sum(A.dim());
			
			for (size_t i=0; i< this->size(); i++) {
				for (size_t j=0; j< A.size(); j++) {
					sum.inplace_union((this->_vector[i])+A[j]);
				}
			}
			
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
			
			return sum;
		}
		
		inline AriadneDenotableSet<BS> &expand_by(const Real &delta) {
			
			for (size_t i=0; i< this->size(); i++) {
				(this->_vector[i]).expand_by(delta);
			}
	
			return *this;
		}
		
		inline AriadneDenotableSet<BS> &set_precision_to_upperapproximating(const Real &delta) {
			
			for (size_t i=0; i< this->size(); i++) {
				(this->_vector[i]).set_precision_to_upperapproximating(delta);
			}	
			
			if (this->empty()) {
				(this->_vector).resize(BS_PER_BLOCK);
				this->_basic_sets=0;
			}
			
			return *this;	
		}
		
		/*! \brief Copy the denotable set. */
		inline const AriadneDenotableSet<BS> &operator=(
					const AriadneDenotableSet<BS> &A) {
		
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
			
			(this->_vector).resize((A._vector).size());
			
			for (size_t i=0; i< A.size() ; i++) {
				this->_vector[i]=A._vector[i];
			}
			
			this->_dim=A._dim;
			this->_basic_sets=A._basic_sets;
			
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
			
			return *this;
		}
		
		/*! \brief Checks if a denotable set includes a state.
		 *
		 * This method checks whenever the current denotable set
		 * includes the state \a s.
		 * \param s is the state of which inclusion 
		 * into the current denotable set should be tested.
		 * \return  \a true, if \a s is contained into the 
		 * current set, \a false otherwise.
		 */
		inline bool contains(const State &s) const {
			
			for (size_t i=0; i<this->size(); i++) {
				if ((this->_vector[i]).contains(s)) 
					return true;
			}

			return false;			
		}
		
		/*! \brief Checks if the interior of the denotable set includes a state.
		 *
		 * This method checks whenever the interior of the current denotable set
		 * includes the state \a s.
		 * \param s is the state of which inclusion 
		 * into the current denotable set should be tested.
		 * \return  \a true, if \a s is contained into the 
		 * current set, \a false otherwise.
		 */
		inline bool interior_contains(const State & state) const {
			
			for (size_t i=0; i<this->size(); i++) {
				if ((this->_vector[i]).interior_contains(s)) 
					return true;
			}

			return false;
			
		}
		
		/*! \brief Checks the emptyness.
		 *
		 * \return \a true if the denotable set is empty,
		 * \a false otherwise.		
		 */
		inline bool empty() const {
						
			for (size_t i=0; i<this->size(); i++) {
				if (!(this->_vector[i]).empty()) 
					return false;
			}
		
			return true;
		}
		
		/*! \brief Returns the begin of the basic set vector.
	     *
		 * \return The begin of the maintained basic set vector.
		 */
		inline const_iterator begin() const {
			return _vector().begin();
		}
		
		/*! \brief Returns the end of the maintained basic set vector.
	     *
		 * \return The end of the maintained basic set vector.
		 */
		inline const_iterator end() const {
			return (_vector().begin()+this->_basic_sets);
		}
		
		/*! \brief Makes the union a denotable set and a basic set.
		 *
		 * Makes the union of the current denotable set with a basci set.
		 * The result is stored into the current object.
		 * \param A is a basic set.
		 */
		inline void inplace_union(const AriadneDenotableSet<BasicSet>& A) {	
			
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
			
			if (this->dim()==0) this->_dim=A.dim();
				
			if (A.dim()!=this->dim()) {
				throw std::invalid_argument("The two denotable set have different space dimensions.");
			}
			
			if (A.empty()) {
								
				#ifdef DEBUG
					std::cout << __FILE__ << ":" << __LINE__ << std::endl;
				#endif
				
				return;
			}
			
			size_t new_size= (A.size() + this->size());
			
			if (new_size > (this->_vector).size() ) {
				
				this->_vector.resize(new_size);
			}
	
			for (size_t i = 0; i < A.size(); i++) {
				this->_vector[this->_basic_sets]=A._vector[i];
				this->_basic_sets++;
			}
			
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
			
		}
		
		/*! \brief Makes the union of two denotable set.
		 *
		 * Makes the union of the current denotable set with an other one.
		 * The result is stored into the current object.
		 * \param A is a basic set.
		 */
		inline void inplace_union(const BasicSet &A) {

			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
			
			if (this->dim()==0) this->_dim=A.dim();
				
			if (A.dim()!=this->dim()) {
				throw std::invalid_argument("The denotable set the basic set have different space dimensions.");
			}
			
			if (A.empty()) return;
			
			if (this->dim()+1 > (this->_vector).size()) {
				
				this->_vector.resize((this->_vector).size()+BS_PER_BLOCK);
			}
			
			this->_vector[this->_basic_sets]=A;
			
			this->_basic_sets++;	
			
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
			
		}
		
		/*! \brief Prints the ptree on a ostream */
		template <typename BSet , uint BSet_PER_BLOCK>
		friend std::ostream& IO_Operators::operator<<(std::ostream &os, 
				const AriadneDenotableSet< BSet , BSet_PER_BLOCK > &r);
		
		/*! \brief Tests disjointness. 
	 	 *
	     	 * Tests disjointness of two denotable sets.
		 * \param A is a denotable set.
		 * \param B is a denotable set.
		 * \return \a true if A and B are disjoint, \a false otherwise.
	     	 */
		friend bool disjoint <> (
				const AriadneDenotableSet< BasicSet , BS_PER_BLOCK > &A, 
				const AriadneDenotableSet< BasicSet , BS_PER_BLOCK > &B);

		/*! \brief Tests intersection of interiors. 
		 *
	     	 * Tests intersection of a denotable set with the interior of an other
		 * denotable set.
		 * \param A is a denotable set.
		 * \param B is a denotable set.
		 * \return \a true if A intersects the interior of B, 
		 * \a false otherwise.
	     	 */
		friend bool intersects_interior <> (
				const AriadneDenotableSet< BasicSet , BS_PER_BLOCK > &A, 
				const AriadneDenotableSet< BasicSet , BS_PER_BLOCK > &B);
	
		/*! \brief Tests intersection of interiors. 
		 *
	     	 * Tests intersection of a denotable set with the interior of a rectangle.
		 * \param A is a rectangle.
		 * \param B is a denotable set.
		 * \return \a true if A intersects the interior of B, 
		 * \a false otherwise.
	     	 */	
		friend bool intersects_interior <> (
				const AriadneRectangle< State > &rect, 
				const AriadneDenotableSet< BasicSet , BS_PER_BLOCK > &A);
		
		/*! \brief Tests inclusion of interiors. 
		 *
	     	 * Tests inclusion of a denotable set into the interior of a denotable set.
		 * \param A is a denotable set.
		 * \param B is a denotable set.
		 * \return \a true if A is a subset of the interior of B, 
		 * \a false otherwise.
	     	 */
		friend bool subset_of_interior <> (
				const AriadneDenotableSet< BasicSet , BS_PER_BLOCK > &A, 
				const AriadneDenotableSet< BasicSet , BS_PER_BLOCK > &B);
	
		/*! \brief Tests inclusion of interiors. 
		 *
	    	 * Tests inclusion of a rectangle into the interior of a denotable set.
		 * \param A is a rectangle.
		 * \param B is a denotable set.
		 * \return \a true if A is a subset of the interior of B, 
		 * \a false otherwise.
	    	 */
		friend bool subset_of_interior <> (
				const AriadneRectangle< State > &rect, 
				const AriadneDenotableSet< BasicSet , BS_PER_BLOCK > &A);

		/*! \brief Tests inclusion of interiors. 
		 *
		 * Tests inclusion of a denotable set into the interior of a rectangle.
		 * \param A is a denotable set.
		 * \param B is a rectangle.
		 * \return \a true if A is a subset of the interior of B, 
		 * \a false otherwise.
		 */
		friend bool subset_of_interior <> (
				const AriadneDenotableSet< BasicSet , BS_PER_BLOCK > &A,
				const AriadneRectangle< State > &rect);

		/*! \brief Makes union of two interiors. 
		 *
	     	 * Evalutates the union of two denotable sets.
		 * \param A is a denotable set.
		 * \param B	is a denotable set.
		 * \return The union of A and B.
	     	 */
		friend AriadneDenotableSet< BasicSet , BS_PER_BLOCK > join <> (
				const AriadneDenotableSet< BasicSet , BS_PER_BLOCK > &A, 
				const AriadneDenotableSet< BasicSet , BS_PER_BLOCK > &B);

		/*! \brief Makes intersection of two interiors. 
		 *
	     	 * Evalutates the closure of the intersection of a denotable set 
		 * with the interiors of an other denotable sets.
		 * \param A is a denotable set.
		 * \param B is a denotable set.
		 * \return The closure of the intersection of A with the interiors of B.
	    	 */
		friend AriadneDenotableSet< BasicSet , BS_PER_BLOCK > 
			closure_of_intersection_of_interior <> (
				const AriadneDenotableSet< BasicSet , BS_PER_BLOCK > &A, 
				const AriadneDenotableSet< BasicSet , BS_PER_BLOCK > &B);
};
	
}
}

#endif /* _DENOTABLE_SET_H */
