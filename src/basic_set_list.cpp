/***************************************************************************
 *            basic_set_list.cpp
 *
 *  Thu Aug 19 10:57:59 2004
 *  Copyright  2004  Alberto Casagrande
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
 
#include <list>

#include "basic_set_list.h"
#include "basic_set.h"
#include "set.h"

using namespace Ariadne;
	
void BasicSetList::clear() {
	std::list<BasicSet *>::iterator i;
	
	/* delete all the list's sets. */
	for ( i=(this->l).begin(); i!=(this->l).end(); i++ ) {
		delete (*i);
		(this->l).pop_front();
	}
	
	/* clear the list */
	(this->l).clear();
}

BasicSetList& BasicSetList::operator=(const BasicSetList& tl){
	
	/* clear the object */
	this->clear();

	/* join a copy of tl on the object */
	this->attach_copy(tl);

	return(*this);
}

BasicSetList::BasicSetList() {
	this->clear();
}

BasicSetList::BasicSetList(const BasicSetList& tl) {
	this->clear();

	/* attach a copy of tl on the object */
	this->attach_copy(tl);
}


BasicSetList::~BasicSetList() {
	this->clear();
}

BasicSetList BasicSetList::join(const BasicSetList& l1) const {
	
	BasicSetList bsl;
	
	/* copy into \a bsl all the current object's basic sets and
	 * attach a copy of that of {\a l1}. */
	bsl=*this;
	bsl.attach_copy(l1);

	return(bsl);
}

BasicSetList BasicSetList::intersect(const BasicSetList& l1,
		const ApproxType atype) const {
	
	std::list<BasicSet *>::const_iterator i,j;
	BasicSet *ij;
	BasicSetList bsl;

	/*for each basic set \a i in the current object ...*/
	for ( i=(this->l).begin(); i!=(this->l).end(); i++ ) {

		/* for each basic set \a j in {\a l1}... */
		for ( j=(l1.l).begin(); j!= (l1.l).end(); j++ ) {

			if ((*i)->does_intersect(*(*j))) {
			
				/* intersect \a i and \a j and create the 
				 * object \a *ij */
				ij=(*i)->intersect(*(*j),atype);
		
				/* include the intersection of i and j into
				 * the list of this object */
				(bsl.l).push_back(ij);
			}
		}
	}

	return(bsl);
}

BasicSetList BasicSetList::intersect(const BasicSet& bset,
		const ApproxType atype) const {
	
	std::list<BasicSet *>::const_iterator i,j;
	BasicSet *ij;
	BasicSetList bsl;
	
	/*for each basic set i in the current object ...*/
	for ( i=(this->l).begin(); i!= (this->l).end(); i++ ) {

		if ((*i)->does_intersect(bset)) {
			/* intersect i and j and create the object *ij */
			ij=(*i)->intersect(bset,atype);
		
			/* delete the basic set \a i */
			delete(*i);
		
			/* include the intersection of i and j into
			 * the list of this object */
			(bsl.l).push_back(ij);
		}
	}

	return(bsl);
}

#ifdef FASTER_BUT_DIRTIER
BasicSetList& BasicSetList::attach(BasicSetList& l1) {
	
	/* attach l1 to the object */
	while (!(l1.empty())) {
			
		(this->l).push_back((l1.l).front());
		(l1.l).pop_front();
	}

	return(*this);
}
#endif


BasicSetList& BasicSetList::attach_copy(const BasicSetList& l1) {
	std::list<BasicSet *>::const_iterator i;
	
	/* join a copy of l1 with the object */
	for ( i=(l1.l).begin(); i!= (l1.l).end(); i++ ) {
		(this->l).push_back((*i)->copy());	
	}

	return(*this);
}

BasicSetList& BasicSetList::remove_doubly_included_set(){
	std::list<BasicSet *>::iterator i,j;

	/* for each basic set \a i in the object ... */
	for ( i=(this->l).begin(); i!= (this->l).end(); i++ ) {

		/* for each basic set \a j following \a i in 
		 * this object */
		for ( j=i; j!= (this->l).end(); ++j ) { // TODO: check this code

			if ((*j)->is_subset_of(*(*i))) {
				/* delete the basic set \a j and copy 
				 * the basic set's pointer in front of 
				 * the list on the position \a j */
				delete(*j);
				(*j)=(this->l).front();

				/* remove the element in front of the list */
				(this->l).pop_front();
			}
		}
	}

	return(*this);
}



#ifdef FAST_BUT_DIRTIER
BasicSetList& BasicSetList::join_on_obj(const BasicSetList& l1) {
	
	this->attach_copy(l1);

	return(*this);
}

BasicSetList& BasicSetList::intersect_on_obj(const BasicSet& bset,
		const ApproxType atype) {
	
	/*for each basic set i in the current object ...*/
	for ( i=(this->l).begin(); i!= (this->l).end(); i++ ) {

		if (bset.does_intersect(*i)) {
			/* intersect i and j and create the object *ij */
			ij=new BasicSet((*i)->intersect(bset,atype));
		
			/* delete the basic set \a i */
			delete(*i);
		
			/* include the intersection of i and j into
			 * the list of this object */
			(*i)=ij;
		}
	}

	return(*this);
}
#endif

bool BasicSetList::is_subset_of(const BasicSetList& l1) const {
	
	/* WARNING: verify the code. I should prove the corectness. */
	BasicSetList *l2, *l3, *sub;
	std::list<BasicSet *>::const_iterator i;
	BasicSet *bs;
	
	l2=new BasicSetList(*this);
	
	/* for each basic set \a i in \a l1 ... */
	for ( i=(l1.l).begin(); i!= (l1.l).end(); i++ ) {
				
		l3=l2;
		l2=new BasicSetList();
	
		while (!((l3->l).empty())) {
	
			bs=(l3->l).front();
			(l3->l).pop_front();
			
			if ((*i)->does_intersect(*bs)) {
				sub=bs->subtract(*(*i),OVER);
				
#ifdef FASTER_BUT_DIRTIER
				l2->attach(*sub);
#else
				l2->attach_copy(*sub);
#endif
				delete sub;
			}
		}

		delete l3;
	}

	
	/* if \a l2 is empty, \a l1 contains completly 
	 * the current object then return {\a true}, otherwise
	 * return {\a false}. */
	if (l2->empty()) {
		delete l2;
		return true;
	} else {
		delete l2;
		return false;
	}
}

bool BasicSetList::does_intersect(const BasicSet& bset) const {
	std::list<BasicSet *>::const_iterator i;	
	
	/* for each basic set \a i in the object ... */
	for ( i=(this->l).begin(); i!= (this->l).end(); i++ ) {
		/* check if \a i intersects \a bset */
		if ((*i)->does_intersect(bset)) {
			return true;
		}
	}

	return false;
}

bool BasicSetList::does_intersect(const BasicSetList& l1) const {
	std::list<BasicSet *>::const_iterator i;	
	
	/* for each basic set \a i in the object ... */
	for ( i=(this->l).begin(); i!= (this->l).end(); i++ ) {
		
		/* check if \a l1 intersects \a i and, if it does,
		 * return {\a true}. */
		if (l1.does_intersect(*(*i))) {
			return true;	
		}
	}

	/* NOTE: if the method's control passes through this point
	 * none of the object's basic sets does not intersect 
	 * \a l1 and then the method should return {\a false}. */
	return false;

}

BasicSetList& BasicSetList::include(const BasicSet &bset){
	
	/* if \a bset is not empty ...*/
	if (!bset.empty()) {
		
		/* create a new object, copy bset on it 
		 * and include it on the list. */
		
		BasicSet *aus=bset.copy();
		this->include(aus);
	}

	return(*this);
}

BasicSetList& BasicSetList::include(BasicSet *bset){
	
	/* if \a bset is not a null pointer and 
	 * it points to a non-empty basic set ...*/
	if ((bset!=NULL)&&(!(bset->empty()))) {

		/* attach it to the list */
		(this->l).push_back(bset);
	}
		
	return(*this);
}


bool BasicSetList::contains(const State &s) const {

	std::list<BasicSet *>::const_iterator i;

	/* for each basic set \a i in the object ... */
	for ( i=(this->l).begin(); i!= (this->l).end(); i++ ) {

		/* if \a i contains {\a s}, return \a true */
		if ((*i)->contains(s)) {
			return true;
		}
	}

	/* none of the basic set into the object contians \a s and then
	 * \a s is not into the current object. */
	return false;
	
}

bool BasicSetList::contains(const BasicSet &bset) const {
	
	BasicSetList l1;
	bool result;
	
	l1.include(bset);

	result=l1.is_subset_of(*this);

	l1.clear();

	return result;
}

BasicSetList& BasicSetList::compact_list(const unsigned int max_elem, 
				ApproxType atype) {

	std::list<BasicSet *>::iterator i,j, idx_i, idx_j;
	double dist_value, max_dv; 
	
	this->remove_doubly_included_set();

	while ((this->l).size() > max_elem) {
		
		max_dv=+INFINITY; 
		idx_i=(this->l).begin();
		idx_j=idx_i;
		idx_j++;
		
		/* for each basic set \a i in the object ... */
		for ( i=(this->l).begin(); i!= (this->l).end(); i++ ) {

			/* for each basic set \a j following i into this
			 * object ...*/
			for ( j=i; j!= (this->l).end(); ++j ) {
				
				/* evaluate the \a distance between 
				 * \a i and \a j */
				dist_value=(*i)->distance(*(*j));

				if (dist_value<max_dv) {
					max_dv=dist_value;
					idx_i=i;
					idx_j=j;
				}
			}
		}
#ifdef FASTER_BUT_DIRTIER
		/* make the union of \a i and \a j */
		(*idx_i)->join_on_obj(*(*idx_j),atype);
#else 
		/* make the union of \a i and \a j */
		BasicSet *ij=(*idx_i)->join(*(*idx_j),atype);

		/* delete the basic set \a i and copy the pointer of
		 * the union on the position \a idx_i */
		delete (*idx_i);
		(*idx_i)=ij;
#endif
		/* delete the basic set \a j and copy the basic set's pointer 
		 * in front of the list on the position \a idx_j */
		delete(*idx_j);
		(*idx_j)=(this->l).front();

		/* remove the element in front of the list */
		(this->l).pop_front();
	}

	return(*this);
}
		
BasicSetList& BasicSetList::compact_list(const unsigned int max_elem){

	return this->compact_list(max_elem,OVER);

}

bool BasicSetList::empty() {
	std::list<BasicSet *>::iterator i;

	/* for each basic set \a i in the object ... */
	for ( i=(this->l).begin(); i!= (this->l).end(); i++ ) {

		/* if \a i is not empty return \a false */
		if (!((*i)->empty())) {
			return false;
		}
		else {
			/* if \a *i is empty, the method can delete
			 * it from the list*/
			delete (*i);

			/* copy the pointer in the list's front in
			 * position \a i and pop the element in 
			 * the front of the list. */
			(*i)=(this->l).front();
			
			(this->l).pop_front();
		}
	}

	/* if, for all \a i in the object, \a i is empty return
	 * {\a true}. */
	return true;
}

BasicSetList BasicSetList::apply(const Map &M, const ApproxType atype) const {
	
	BasicSet *bs;
	BasicSetList out;
	
	std::list<BasicSet *>::const_iterator i;

	/* for each basic set \a i in the object ... */
	for ( i=(this->l).begin(); i!= (this->l).end(); i++ ) {

		/* insert into the output list M(i). */
		bs= (*i)->apply(M,atype);

		(out.l).push_back(bs);
	}

	return out;	
}

BasicSetList* BasicSetList::apply_p(const Map &M, 
		const ApproxType atype) const {
	
	BasicSet *bs;
	BasicSetList *out= new BasicSetList();
	
	std::list<BasicSet *>::const_iterator i;

	/* for each basic set \a i in the object ... */
	for ( i=(this->l).begin(); i!= (this->l).end(); i++ ) {

		/* insert into the output list M(i). */
		bs= (*i)->apply(M,atype);

		(out->l).push_back(bs);
	}

	return out;	
}

#ifdef FASTER_BUT_DIRTIER
BasicSetList& BasicSetList::apply_on_obj(const Map &M, const ApproxType atype) {
	
	std::list<BasicSet *>::iterator i;

	/* for each basic set \a i in the object ... */
	for ( i=(this->l).begin(); i!= (this->l).end(); i++ ) {

		/* apply M to *i */
		(*i)->apply_on_obj(M,atype);
	}

	return(*this);	
}
#endif

BasicSetList* BasicSetList::flow_slice_p(const VectorField &field,
				const double delta,
				const ApproxType atype) const {

		
	BasicSet *bs;
	BasicSetList *out= new BasicSetList();
	
	std::list<BasicSet *>::const_iterator i;

	/* for each basic set \a i in the object ... */
	for ( i=(this->l).begin(); i!= (this->l).end(); i++ ) {

		/* insert into the output list the slice of \a i . */
		bs= (*i)->flow_slice(field,delta,atype);

		(out->l).push_back(bs);
	}

	return out;	
	
}

#ifdef FASTER_BUT_DIRTIER
BasicSetList& BasicSetList::flow_slice_on_obj(const VectorField &field,
				const double delta,
				const ApproxType atype) {
		
	BasicSet *bs;
	
	std::list<BasicSet *>::iterator i;

	/* for each basic set \a i in the object ... */
	for ( i=(this->l).begin(); i!= (this->l).end(); i++ ) {

		/* insert into the output list the slice of \a i . */
		bs= (*i)->flow_slice(field,delta,atype);

		delete (*i);

		(*i)=bs;
	}

	return(*this);	
}
#endif 

BasicSetList BasicSetList::expand_of(const double delta) const {

		
	BasicSet *bs;
	BasicSetList out;
	
	std::list<BasicSet *>::const_iterator i;

	/* for each basic set \a i in the object ... */
	for ( i=(this->l).begin(); i!= (this->l).end(); i++ ) {

		/* insert into the output list the expanded set. */
		bs= (*i)->expand_of(delta);

		(out.l).push_back(bs);
	}

	return out;	
	
}

BasicSetList* BasicSetList::expand_of_p(const double delta) const {
		
	BasicSet *bs;
	BasicSetList *out= new BasicSetList();
	
	std::list<BasicSet *>::const_iterator i;

	/* for each basic set \a i in the object ... */
	for ( i=(this->l).begin(); i!= (this->l).end(); i++ ) {

		/* insert into the output list the expanded set. */
		bs= (*i)->expand_of(delta);

		(out->l).push_back(bs);
	}

	return out;	
	
}

#ifdef FASTER_BUT_DIRTIER
BasicSetList& BasicSetList::expand_of_on_obj(const double delta){
		
	BasicSet *bs;
	
	std::list<BasicSet *>::iterator i;

	/* for each basic set \a i in the object ... */
	for ( i=(this->l).begin(); i!= (this->l).end(); i++ ) {

		/* insert into the output list the expanded set. */
		bs= (*i)->expand_of(delta);

		delete (*i);

		(*i)=bs;
	}

	return(*this);	
}
#endif

BasicSetList BasicSetList::shrink_of(const double delta) const {

		
	BasicSet *bs;
	BasicSetList out;
	
	std::list<BasicSet *>::const_iterator i;

	/* for each basic set \a i in the object ... */
	for ( i=(this->l).begin(); i!= (this->l).end(); i++ ) {

		/* insert into the output list the shrinked set. */
		bs= (*i)->shrink_of(delta);

		(out.l).push_back(bs);
	}

	return out;	
	
}


BasicSetList* BasicSetList::shrink_of_p(const double delta) const {

		
	BasicSet *bs;
	BasicSetList *out= new BasicSetList();
	
	std::list<BasicSet *>::const_iterator i;

	/* for each basic set \a i in the object ... */
	for ( i=(this->l).begin(); i!= (this->l).end(); i++ ) {

		/* insert into the output list the shrinked set. */
		bs= (*i)->shrink_of(delta);

		(out->l).push_back(bs);
	}

	return out;	
	
}

#ifdef FASTER_BUT_DIRTIER
BasicSetList& BasicSetList::shrink_of_on_obj(const double delta){
		
	BasicSet *bs;
	
	std::list<BasicSet *>::iterator i;

	/* for each basic set \a i in the object ... */
	for ( i=(this->l).begin(); i!= (this->l).end(); i++ ) {

		/* insert into the output list the shrinked set. */
		bs= (*i)->shrink_of(delta);

		delete (*i);

		(*i)=bs;
	}

	return(*this);	
}
#endif


In_BSL& In_BSL::compact_list(Ex_BSL &ex, const unsigned int max_elem,
		const ApproxType atype) {
	std::list<BasicSet *>::iterator i;

	/* compact the object's list */
	((BasicSetList *)this)->compact_list(max_elem,atype);
	
	/* for each basic set \a i in the object ... */
	for ( i=(ex.l).begin(); i!= (ex.l).end(); i++ ) {

		/* if the object's list does not intersect \a i ... */
		if (!(this->does_intersect(*(*i)))) {
			
			/* delete the basic set \a i and copy the basic set's 
			 * pointer in front of the list on the 
			 * position \a i */
			delete(*i);
			(*i)=(ex.l).front();

			/* remove the element in front of the list */
			(ex.l).pop_front();
		}

	}

	/* if a over-approximation is needed, the methods 
	 * should under-approximates the \a ExBSL list
	 * [ (A/B) => (A/C) <==> B <= C ]. Otherwise if a 
	 * under-approximation is needed, the methods 
	 * should over-approximates the \a ExBSL list
	 * [ (A/B) <= (A/C) <==> B => C ]. */
	((BasicSetList &)ex).compact_list(max_elem,invert_approx(atype));
	
	
	/* if \a this is a subset of ex, \f$this\ex=\emptyset\f$,
	 * thus the method can clear both \a this and {\a ex}. */
	if (this->is_subset_of(ex)) {
		this->clear();
		ex.clear();
	}
	
	return(*this);
}

In_BSL& In_BSL::operator=(const BasicSetList& tl){
	
	((BasicSetList)(*this))=tl;

	return(*this);
}

Ex_BSL& Ex_BSL::compact_list(In_BSL &in, const unsigned int max_elem,
		const ApproxType atype) {

	in.compact_list(*this, max_elem, atype);

	return(*this);
}

Ex_BSL& Ex_BSL::operator=(const BasicSetList& tl){
	
	((BasicSetList)(*this))=tl;

	return(*this);
}

