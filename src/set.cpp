/***************************************************************************
 *            set.cpp
 *
 *  Fri Aug  6 10:57:59 2004
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

#include "arc.h"
#include "set.h"
#include "basic_set.h"
#include "vectorfield.h"

using namespace Ariadne;
	
	
ASet::ASet() {
	this->max_basic_sets=DEFAULT_MAX_OBJECTS;
	this->atype=NONE;
}
		
ASet::ASet(const ASet &t) {

	/* copy the {\a t}'s lists */
	this->in_list=t.in_list;
	this->ex_list=t.ex_list;

	/* copy the others {\a t}'s members */
	this->max_basic_sets=t.max_basic_sets;
	this->atype=t.atype;
}
		
ASet::ASet(unsigned int max_basic_sets) {
	this->max_basic_sets=max_basic_sets;
	this->atype=NONE;
}
		
ASet::ASet(ApproxType type) {
	this->max_basic_sets=DEFAULT_MAX_OBJECTS;
	this->atype=type;
}
		
ASet::ASet(unsigned int max_basic_sets, ApproxType type) {
	this->max_basic_sets=max_basic_sets;
	this->atype=type;
}

ASet::~ASet() {

	/* useless, but... "paranoia is good" */
	
	/* clear the \a in_list list's */
	(this->in_list).clear();

	/* clear the \a ex_list list's */
	(this->ex_list).clear();
}
	
ASet& ASet::set_max_basic_sets(unsigned int max_basic_sets) {
	this->max_basic_sets=max_basic_sets;

	return(*this);
}
	
unsigned int ASet::get_max_basic_sets() {
	return this->max_basic_sets;
}
		
ASet& ASet::include(const BasicSet &bset) {

	/* if \a bset is not empty... */
	if (!bset.empty()) {

		/* create a new basic set using as template \a bset and
		 * include the new basic set from the current object */
		this->include(bset.copy());
	}

	return(*this);
}

ASet& ASet::include(BasicSet *bset) {

	/* if \a bset is not null... */
	if (bset!=NULL) {

		/* include \a bset in the \a in_list list */
		(this->in_list).include(bset);

		/* compact lists */
		this->compact_lists();

	}
	
	return(*this);
}

ASet& ASet::exclude(const BasicSet &bset) {

	/* if \a bset is not empty... */
	if (!bset.empty()) {

		/* create a new basic set using as template \a bset 
		 * and exclude the new basic set from the current object */
		this->exclude(bset.copy());
		
	}

	return(*this);
}

ASet& ASet::exclude(BasicSet *bset) {

	/* if \a bset is not null... */
	if (bset!=NULL) {
		
		/* include \a bset in the \a ex_list list */
		(this->ex_list).include(bset);
		
		/* compact lists */
		this->compact_lists();
	}
	
	return(*this);
}

ASet& ASet::operator=(const ASet& ts) {
	
	this->max_basic_sets=ts.max_basic_sets;
	this->atype=ts.atype;

	this->in_list=ts.in_list;
	this->ex_list=ts.ex_list;
	
	return(*this);
}

ASet ASet::operator&&(const BasicSet &bset) const {
	ASet out;
		
	out.max_basic_sets=this->max_basic_sets;
	out.atype=this->atype;
	
	out.in_list=(this->in_list).intersect(bset, this->atype);
	
	out.ex_list=this->ex_list;
	
	out.compact_lists();

	return(out);	
}
		
ASet* ASet::intersect(const BasicSet &bset) const{

#ifdef FASTER_BUT_DIRTIER
	ASet *out=new ASet(*this);
		
	(out->in_list).intersect_on_obj(bset, this->atype);
#else
	ASet *out=new ASet(this->max_basic_sets, this->atype);
	
	out->in_list=(this->in_list).intersect(bset, this->atype);
	out->ex_list=this->ex_list;
	
#endif	
	out->compact_lists();

	return(out);	
}


ASet& ASet::compact_lists() {
	
	/* compact in_list and ex_list */
	(this->in_list).compact_list(this->ex_list,
				     this->max_basic_sets,
				     this->atype);

	return(*this);
}
		
bool ASet::does_intersect(const ASet &R) const {
	
	BasicSetList in_inter, ex_union;

	/* evaluate intersection of the in_list of this object and 
	 * of \a R */
	in_inter=(this->in_list).intersect(R.in_list,OVER);
	
	/* if \a in_inter is empty this object and \a R do not 
	 * intersect */
	if (in_inter.empty()) {
		return false;
	}
		
	/* evaluate union of the ex_list of this object and 
	 * of \a R  */
	ex_union=(this->ex_list).join(R.ex_list);

	/* if \a in_inter is a subset of \a ex_union this object 
	 * and \a R do not intersect ( see "prove paper 1" for
	 * motivations) */
	if (in_inter.is_subset_of(ex_union)) {
		return false;
	} 
	
	/* if \a in_inter is not a subset of \a ex_union this object 
	 * and \a R intersects ( see "prove paper 1" for motivations) */
	return true;
}

bool ASet::empty() {

	/* if \a in_list list are empty or in_list is a subset of
	 * \a ex_list return true, otherwise return false. */
	
	return (((this->in_list).empty())||
			((this->in_list).is_subset_of(ex_list)));
}

bool ASet::is_subset_of(const ASet& R) const {
	BasicSetList ie_inter, ei_union;

	/* evaluate intersection of the current object's \a in_list 
	 * and the \a ex_list of \a R */
	ie_inter=(this->in_list).intersect(R.ex_list,OVER);
	
	/* evaluate union of the current object's \a ex_list 
	 * and the \a in_list of \a R */
	ei_union=(this->ex_list).join(R.ex_list);

	
	/* If \a ie_inter is not a subset of \a this->ex_list
	 * or \a this->in_list is not a subset of \a ei_union
	 * this is not a subset of R. ( see "prove paper 1" 
	 * for motivations ) */
	if ((!ie_inter.is_subset_of(this->ex_list))||
			(!(this->in_list).is_subset_of(ei_union))) {
		return false;
	}

	return true;
}
		
bool ASet::is_disjoint_from(const ASet& R) const {
	
	/* this object is disjoint from \a R if
	 * is does not intersect {\a R}. */
	return (!(this->does_intersect(R)));
}
		
bool ASet::contains(const State &s) const {

	/* this object contains a state \a s if
	 * some basic set in \a in_list
	 * contain \a s and \a s is not contained 
	 * by any basic set of 
	 * \a ex_list */

	if ((this->ex_list).contains(s)) {
		return false;
	}
	
	/* NOTE: if the method's control passes through this 
	 * point \a s is not contained into the \a ex_list
	 * list. */

	if ((this->in_list).contains(s)) {
		return true;
	}

	/* NOTE: if the method's control passes through this 
	 * point \a s is not contained into the \a in_list
	 * list and then the method should return false.*/
	
	return false;

}

ASet ASet::operator&&(const ASet &R) const {

	/* TODO: reimplement this method: can be
	 * implemented with faster code. */
	
	ASet *aus=this->intersect(R);

	ASet out=*aus;

	delete aus;

	return out;
}

ASet* ASet::intersect(const ASet &R) const {
	
	/* create the output ASet */
	ApproxType type=this->atype;
	
	if (type==NONE) {
		type=R.atype;
	} else {
		if (type!=R.atype) {
			type=OVER;
		}
	}

	return (this->intersect(R,type));
}

ASet* ASet::intersect(const ASet &R, const ApproxType type) const {

	/* create the output ASet */
	unsigned int max_basic=this->max_basic_sets;
	
	if (R.max_basic_sets>max_basic) max_basic=R.max_basic_sets;

	ASet *out=new ASet(max_basic, type);

	/* evaluate the intersection of the \a in_list of this object 
	 * and of \a R and assign it to the \a in_list of the output
	 * object. */
	out->in_list=(this->in_list).intersect(R.in_list,type);
	
	/* evaluate the union of the \a ex_list of this object and 
	 * of \a R and assign it to the \a ex_list of the output
	 * object. */
#ifdef FASTER_BUT_DIRTIER
	out->ex_list=this->ex_list;
	(out->ex_list).join_on_obj(R.ex_list);
#else
	out->ex_list=(this->ex_list).join(R.ex_list);
#endif

	/* compact lists \a in_list and \a ex_list */ 	
	out->compact_lists();
	

	/* out=(this->in_list && R.in_list ) \ (this->ex_list !! R.ex_list )
	 * (see "prove paper 1" for motivations). */
	return out;
}
	
ASet ASet::operator||(const ASet &R) const {
	
	/* TODO: reimplement this method: can be
	 * implemented with faster code. */
	
	ASet *aus=this->join(R);

	ASet out=*aus;

	delete aus;

	return out;
}

ASet* ASet::join(const ASet &R) const {

	/* find the needed approximation type */
	ApproxType type=this->atype;
	
	if (type==NONE) {
		if (R.atype==NONE) {

			/* if ( Ex_A && In_B ) = ( Ex_B && In_A ) = emptyset 
			 * then the type approximation of union could be
			 * NONE (see "prove paper 1" for motivations). */
			if (!(((this->ex_list).does_intersect(R.in_list))||
				((this->in_list).does_intersect(R.ex_list)))){

				type=NONE;
			} else {
				type=OVER;
			}
		} else {
			type=R.atype;
		}
	} else {
		if (type!=R.atype) {
			type=OVER;
		}
	}
	

	return (this->join(R,type));
}

ASet* ASet::join(const ASet &R, const ApproxType type) const {
	
	/* Create the output ASet */
	unsigned int max_basic=this->max_basic_sets;
	
	if (R.max_basic_sets>max_basic) max_basic=R.max_basic_sets;
	
	ASet *out=new ASet(max_basic, type);
	
	/* the union of two ASets is approximated 
	 * by the formulae presented in "prove paper 1".*/ 

	/* include from the new object all the basic sets
	 * included by this object and by \a R */
#ifdef FASTER_BUT_DIRTIER
	out->in_list=this->in_list;
	(out->in_list).join_on_obj(R.in_list);
#else
	out->in_list=(this->in_list).join(R.in_list);
#endif
	
	switch(type) {
		case OVER:
			/* if type is OVER, A->join(B) = 
			 * ( In_A || In_B , Ex_A && Ex_B ) 
			 * (see "prove paper 1" for motivations). */
			out->ex_list=(this->ex_list).intersect(R.ex_list,
					UNDER);
	
			break;
		case UNDER:
		case NONE:
			/* if \a comput is UNDER, A->join(B) = 
			 * ( In_A || In_B , Ex_A || Ex_B ) 
			 * (see "prove paper 1" for motivations). */
			
			/* exclude from the new object all the basic sets
			 * excluded by this object and by \a R */
#ifdef FASTER_BUT_DIRTIER
			out->ex_list=this->ex_list;
			(out->ex_list).join_on_obj(R.ex_list);
#else
			out->ex_list=(this->ex_list).join(R.ex_list);
#endif

			break;
			
		default:
			/* NOTE: if the method's control arrives here
			 * there are some problems. */
			break;
	}

	out->compact_lists();
	
	return out;
}

Set ASet::operator&&(const Set &R) const{
	return (R&&(*this));
}
		
Set* ASet::intersect(const Set &R) const{
	return (R.intersect(*this));
}

Set ASet::operator||(const Set &R) const {
	return (R||(*this));
}

Set* ASet::join(const Set &R) const {
	return (R.join(*this));
}

HSet ASet::operator&&(const HSet &R) const{
	return (R&&(*this));
}
		
HSet* ASet::intersect(const HSet &R) const{
	return (R.intersect(*this));
}

HSet ASet::operator||(const HSet &R) const {
	return (R||(*this));
}

HSet* ASet::join(const HSet &R) const {
	return (R.join(*this));
}

ASet* ASet::apply(const Map &M) const {

#ifdef FASTER_BUT_DIRTIER

	/* create the output object */
	ASet *out=new ASet(*this);

	/* apply M to the output object's lists */
	(out->in_list).apply_on_obj(M,this->atype);
	(out->ex_list).apply_on_obj(M,invert_approx(this->atype));

#else
	
	/* create the output object */
	ASet *out=new ASet(this->max_basic_sets,this->atype);

	/* apply M to the lists */
	(out->in_list)=(this->in_list).apply(M,this->atype);
	(out->ex_list)=(this->ex_list).apply(M,invert_approx(this->atype));

#endif
	return out;
}

ASet ASet::operator<<(const Map &M) const{

#ifdef FASTER_BUT_DIRTIER

	/* create the output object */
	ASet out(*this);

	/* apply M to the output object's lists */
	(out.in_list).apply_on_obj(M,this->atype);
	(out.ex_list).apply_on_obj(M,invert_approx(this->atype));

#else
	
	/* create the output object */
	ASet out(this->max_basic_sets,this->atype);

	/* apply M to the lists */
	out.in_list=(this->in_list).apply(M,this->atype);
	out.ex_list=(this->ex_list).apply(M,invert_approx(this->atype));

#endif
	return out;
}

#ifdef FASTER_BUT_DIRTIER
ASet &ASet::apply_on_obj(const Map &M) {

	/* apply M to the object's lists */
	(this->in_list).apply_on_obj(M,this->atype);
	(this->ex_list).apply_on_obj(M,invert_approx(this->atype));

	return (*this);
}
#endif

ASet* ASet::flow_slice_p(const VectorField &field,const double delta) const {
#ifdef FASTER_BUT_DIRTIER

	/* create the output object */
	ASet *out=new ASet(*this);

	/* evaluate the slice lists of the object's lists */
	(out->in_list).flow_slice_on_obj(field,delta,this->atype);
	(out->ex_list).flow_slice_on_obj(field,delta,
					 invert_approx(this->atype));

#else
	
	/* create the output object */
	ASet *out=new ASet(this->max_basic_sets,this->atype);

	/* evaluate the slice lists of the object's lists */
	BasicSetList *in=(this->in_list).flow_slice_p(field,delta,this->atype);
	
	BasicSetList *ex=(this->ex_list).flow_slice_p(field,delta,
					 invert_approx(this->atype));
	
	(out->in_list)=*in;
	(out->ex_list)=*ex;

	delete in;
	delete ex;

	in=NULL;
	ex=NULL;
#endif
	return out;	
}

ASet ASet::flow_slice(const VectorField &field,const double delta) const {
#ifdef FASTER_BUT_DIRTIER

	/* create the output object */
	ASet out(*this);

	/* evaluate the slice lists of the object's lists */
	(out.in_list).flow_slice_on_obj(field,delta,this->atype);
	(out.ex_list).flow_slice_on_obj(field,delta,
					 invert_approx(this->atype));

#else
	
	/* create the output object */
	ASet out(this->max_basic_sets,this->atype);

	/* evaluate the slice lists of the object's lists */
	BasicSetList *in=(this->in_list).flow_slice_p(field,delta,this->atype);
	
	BasicSetList *ex=(this->ex_list).flow_slice_p(field,delta,
					 invert_approx(this->atype));
	
	(out.in_list)=*in;
	(out.ex_list)=*ex;

	delete in;
	delete ex;

	in=NULL;
	ex=NULL;
#endif
	return out;	
}

#ifdef FASTER_BUT_DIRTIER
ASet& ASet::flow_slice_on_obj(const VectorField &field,
		const double delta) {

	/* evaluate the slice lists of the object's lists */
	(this->in_list).flow_slice_on_obj(field,delta,this->atype);
	(this->ex_list).flow_slice_on_obj(field,delta,
					  invert_approx(this->atype));

	return(*this);
}
#endif


ASet* ASet::flow_p(const VectorField &field,const double delta) const {

	/* TODO: try an adaptive numerical solution */
	/* TODO: error should be a map rather than a scalar
	 * value. As a matter of facts, it depends on the
	 * current object.*/

	Map *f_sol;
	double error=0.0;
	
	/* evaluate the flow map */
	if (field.exactly_solvable()) {
		f_sol=field.get_solution(delta);
	} else {
		f_sol=field.get_numerical_solution(*this, delta, &error);
	}

	/* apply the map to this object */
	ASet *out=this->apply(*f_sol);

	out->take_in_account_error(error);
	
	return out;
}

ASet ASet::flow(const VectorField &field,const double delta) const {

	/* TODO: try an adaptive numerical solution */
	/* TODO: error should be a map rather than a scalar
	 * value. As a matter of facts, it depends on the
	 * current object.*/
	
	const Map *f_sol;
	double error=0.0;
	
	/* evaluate the flow map */
	if (field.exactly_solvable()) {
		f_sol=field.get_solution(delta);
	} else {
		f_sol=field.get_numerical_solution(*this, delta, &error);
	}

#ifdef FASTER_BUT_DIRTIER
	/* apply the map to this object */
	ASet out(*this);
	
	out.apply_on_obj(*f_sol);
#else 
	/* apply the map to this object */
	ASet out=(*this)<<(*f_sol);
#endif 
	out.take_in_account_error(error);
	
	return out;
}

#ifdef FASTER_BUT_DIRTIER
ASet& ASet::flow_on_obj(const VectorField &field,
		const double delta) {

	/* TODO: try an adaptive numerical solution */
	/* TODO: error should be a map rather than a scalar
	 * value. As a matter of facts, it depends on the
	 * current object.*/
	
	const Map *f_sol;
	double error=0.0;
	
	/* evaluate the flow map */
	if (field.exactly_solvable()) {
		f_sol=field.get_solution(delta);
	} else {
		f_sol=field.get_numerical_solution(*this, delta, &error);
	}

	/* apply the map to this object */
	this->apply_on_obj(*f_sol);

	this->take_in_account_error(error);
	
	return (*this);
}
#endif


ASet& ASet::take_in_account_error(const double e) {

	switch (this->atype) {
		case OVER:
		case NONE:
#ifdef FASTER_BUT_DIRTIER
			(this->in_list).expand_of_on_obj(e);
			(this->ex_list).shrink_of_on_obj(e);
#else
			(this->in_list)=(this->in_list).expand_of(e);
			(this->ex_list)=(this->ex_list).shrink_of(e);
#endif
			break;
		case UNDER:
#ifdef FASTER_BUT_DIRTIER
			(this->ex_list).expand_of_on_obj(e);
			(this->in_list).shrink_of_on_obj(e);
#else
			(this->ex_list)=(this->ex_list).expand_of(e);
			(this->in_list)=(this->in_list).shrink_of(e);
#endif
			break;
		default:
			/* TODO: code to manage errors */
			break;
	}

	return(*this);
}

OverSet::OverSet(const ASet &t): ASet(t) {
	this->atype=OVER;
}

OverSet::OverSet(unsigned int max_basic_sets) {
	this->max_basic_sets=max_basic_sets;
	this->atype=OVER;
}

OverSet::OverSet() {
	this->max_basic_sets=DEFAULT_MAX_OBJECTS;
	this->atype=OVER;
}

UnderSet::UnderSet(const ASet &t): ASet(t) {
	this->atype=UNDER;
}

UnderSet::UnderSet(unsigned int max_basic_sets) {
	this->max_basic_sets=max_basic_sets;
	this->atype=UNDER;
}

UnderSet::UnderSet() {
	this->max_basic_sets=DEFAULT_MAX_OBJECTS;
	this->atype=UNDER;
}


Set::Set(){
	this->over_appr=new OverSet();
	this->under_appr=new UnderSet();
}
		
Set::Set(const Set &t){
	this->over_appr=new OverSet(*t.over_appr);
	this->under_appr=new UnderSet(*t.under_appr);
}
		
Set::Set(const ASet &over,const ASet &under){
	this->over_appr=new OverSet(over);
	this->under_appr=new UnderSet(under);
}

Set::Set(ASet *over, ASet *under){
	this->over_appr=(OverSet *)over;
	this->under_appr=(UnderSet *)under;
}

Set::Set(unsigned int max_basic_sets) {
			
	this->over_appr=new OverSet(max_basic_sets);
	this->under_appr=new UnderSet(max_basic_sets);
}
		
Set::~Set() {
	delete (this->under_appr);
	delete (this->over_appr);

	/* useless, but I prefer to do it */
	this->under_appr=NULL;
	this->over_appr=NULL;
}

Set& Set::operator=(const Set& ts) {

	(*this->under_appr)=(*ts.under_appr);
	(*this->over_appr)=(*ts.over_appr);
	
	return (*this);
}

void Set::set_max_basic_sets(unsigned int max_basic_sets) {

	(this->under_appr)->set_max_basic_sets(max_basic_sets);
	(this->over_appr)->set_max_basic_sets(max_basic_sets);
}
	
unsigned int Set::get_max_basic_sets() const {

	/* Under and over approximations use the same number of
	 * basic sets. */
	return ((this->under_appr)->get_max_basic_sets());
}
		
Set& Set::include(BasicSet *set) {
	
	(this->under_appr)->include(set);
	(this->over_appr)->include(set);

	return (*this);
}	
	  	
Set& Set::exclude(BasicSet *set) {
	(this->under_appr)->exclude(set);
	(this->over_appr)->exclude(set);

	return (*this);
}

HSet* Set::intersect(const HSet &set) const {
	return set.intersect(*this);
}

HSet Set::operator&&(const HSet &set) const {
	return (set&&(*this));
}

HSet* Set::join(const HSet &set) const {
	return set.join(*this);
}

HSet Set::operator||(const HSet &set) const {
	return (set||(*this));
}

Set* Set::intersect(const BasicSet &set) {
	
	ASet *under=(this->under_appr)->intersect(set);
	ASet *over=(this->over_appr)->intersect(set);

	Set *out= new Set(over,under);
	
	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;
}

Set Set::operator&&(const BasicSet &set) {

	ASet *under=(this->under_appr)->intersect(set);
	ASet *over=(this->over_appr)->intersect(set);

	Set out(over,under);
	
	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;
}


/*
 * ASet::join(const BasicSet &set) is NOT
 * implemented (I should use the 
 * set subtraction!!!)

Set* Set::join(const BasicSet &set) {
	
	ASet *under=(this->under_appr)->join(set);
	ASet *over=(this->over_appr)->join(set);

	Set *out= new Set(over,under);
	
	// useless, but... "paranoia is good"
	over=NULL;
	under=NULL;
	
	return out;
}

Set& Set::operator||(const BasicSet &set) {

	ASet *under=(this->under_appr)->join(set);
	ASet *over=(this->over_appr)->join(set);

	Set out(over,under);
	
	// useless, but... "paranoia is good" 
	over=NULL;
	under=NULL;
	
	return out;
}

*/

bool Set::does_intersect(const ASet &set, ApproxType atype) const {
	
	if (atype==UNDER) {
		return ((this->under_appr)->does_intersect(set));

	} else {
		return ((this->over_appr)->does_intersect(set));
	}	
	
}

bool Set::does_intersect(const Set &set, ApproxType this_type, 
		ApproxType set_type) const {

	if (set_type==UNDER) {
		return (this->does_intersect(*set.under_appr,this_type));

	} else {
		return (this->does_intersect(*set.over_appr,this_type));
	}	
}
	
bool Set::does_intersect(const Set &set) const {
	return (this->does_intersect(set,OVER,OVER));
}

bool Set::does_intersect(const ASet &set) const {
	return (this->does_intersect(set,OVER));
}

bool Set::is_subset_of(const ASet &set, ApproxType atype) const {
	
	if (atype==UNDER) {
		return ((this->under_appr)->is_subset_of(set));

	} else {
		return ((this->over_appr)->is_subset_of(set));
	}	
	
}

bool Set::is_subset_of(const Set &set, ApproxType this_type, 
		ApproxType set_type) const {

	if (set_type==UNDER) {
		return (this->is_subset_of(*set.under_appr,this_type));

	} else {
		return (this->is_subset_of(*set.over_appr,this_type));
	}	
}

bool Set::is_subset_of(const Set& set) const {
	return (this->is_subset_of(set, OVER, OVER));
}

bool Set::is_subset_of(const ASet& set) const {
	return (this->is_subset_of(set, OVER));
}


bool Set::is_disjoint_from(const Set& set, ApproxType this_type,
				ApproxType set_type) const {
	
	return (!(this->does_intersect(set, this_type,set_type)));
}


bool Set::is_disjoint_from(const Set& set) const {
	return (!(this->does_intersect(set, OVER, OVER)));
}
	
bool Set::contains(const State& s, ApproxType type) const {

	if (type==UNDER) {
		return ((this->under_appr)->contains(s));
	} else {
		return ((this->over_appr)->contains(s));
	}
}

bool Set::contains(const State &s) const {
	return ((this->over_appr)->contains(s));
}
	
Set* Set::intersect(const Set &set) const {

	/* evalutate the intersections between
	 * over-approximations and between
	 * under-approxiamtions (the type 
	 * of approximation done during the
	 * computation depends on the type
	 * of the ASet).
	 */
	ASet *over=(this->over_appr)->
		intersect(*(set.over_appr));
	
	ASet *under=(this->under_appr)->
		intersect(*(set.under_appr));

	Set *out=new Set(over,under);

	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;
}

Set Set::operator&&(const Set &set) const {

	/* evalutate the intersections between
	 * over-approximations and between
	 * under-approxiamtions (the type 
	 * of approximation done during the
	 * computation depends on the type
	 * of the ASet).
	 */
	ASet *over=(this->over_appr)->
		intersect(*(set.over_appr));
	
	ASet *under=(this->under_appr)->
		intersect(*(set.under_appr));

	Set out(over,under);

	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;

}

Set* Set::join(const Set &set) const {

	/* evalutate the unions between
	 * over-approximations and between
	 * under-approxiamtions (the type 
	 * of approximation done during the
	 * computation depends on the type
	 * of the ASet).
	 */
	ASet *over=(this->over_appr)->
		join(*(set.over_appr));
	
	ASet *under=(this->under_appr)->
		join(*(set.under_appr));

	Set *out=new Set(over,under);

	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;
}

Set Set::operator||(const Set &set) const {

	/* evalutate the unions between
	 * over-approximations and between
	 * under-approxiamtions (the type 
	 * of approximation done during the
	 * computation depends on the type
	 * of the ASet).
	 */
	ASet *over=(this->over_appr)->
		join(*(set.over_appr));
	
	ASet *under=(this->under_appr)->
		join(*(set.under_appr));

	Set out(over,under);

	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;

}

Set* Set::intersect(const ASet &set) const {

	/* evalutate the intersections between
	 * over-approximation and between
	 * under-approxiamtion (the type 
	 * of approximation done during the
	 * computation depends on the type
	 * of the ASet).
	 */
	ASet *over=(this->over_appr)->intersect(set);
	
	ASet *under=(this->under_appr)->intersect(set);

	Set *out=new Set(over,under);

	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;
}

Set Set::operator&&(const ASet &set) const {

	/* evalutate the intersections between
	 * over-approximation and between
	 * under-approxiamtion (the type 
	 * of approximation done during the
	 * computation depends on the type
	 * of the ASet).
	 */
	ASet *over=(this->over_appr)->
		intersect(set);
	
	ASet *under=(this->under_appr)->
		intersect(set);

	Set out(over,under);

	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;

}

Set* Set::join(const ASet &set) const {

	/* evalutate the unions between
	 * over-approximation and between
	 * under-approxiamtion (the type 
	 * of approximation done during the
	 * computation depends on the type
	 * of the ASet).
	 */
	ASet *over=(this->over_appr)->join(set);
	
	ASet *under=(this->under_appr)->join(set);

	Set *out=new Set(over,under);

	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;
}

Set Set::operator||(const ASet &set) const {

	/* evalutate the unions between
	 * over-approximation and between
	 * under-approxiamtion (the type 
	 * of approximation done during the
	 * computation depends on the type
	 * of the ASet).
	 */
	ASet *over=(this->over_appr)->join(set);
	
	ASet *under=(this->under_appr)->join(set);

	Set out(over,under);

	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;

}

Set* Set::apply(const Map &M) const {
	/* evalutate the application of the 
	 * map \a M (the approximation type 
	 * done during the computation 
	 * depends on the type of the 
	 * ASet).
	 */
	ASet *over=(this->over_appr)->apply(M);
	
	ASet *under=(this->under_appr)->apply(M);

	Set *out=new Set(over,under);

	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;
}

Set Set::operator<<(const Map &M) const {
	/* evalutate the application of the 
	 * map \a M (the approximation type 
	 * done during the computation 
	 * depends on the type of the 
	 * ASet).
	 */
	ASet *over=(this->over_appr)->apply(M);
	
	ASet *under=(this->under_appr)->apply(M);

	Set out(over,under);

	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;
}
		
#ifdef FASTER_BUT_DIRTIER
Set& Set::apply_on_obj(const Map &M) {
	/* evalutate the application of the 
	 * map \a M (the approximation type 
	 * done during the computation 
	 * depends on the type of the 
	 * ASet).
	 */
	(this->over_appr)->apply_on_obj(M);
	
	(this->under_appr)->apply_on_obj(M);

	return (*this);
}
#endif

Set* Set::flow_slice_p(const VectorField &field, const double delta) const {
	
	/* evalutate the slices of both over
	 * and under approximations (the 
	 * approximation type done during 
	 * the computation depends on the 
	 * type of the ASet).
	 */
	ASet *over=(this->over_appr)->flow_slice_p(field, delta);
	
	ASet *under=(this->under_appr)->flow_slice_p(field, delta);

	Set *out=new Set(over,under);

	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;
}

Set Set::flow_slice(const VectorField &field, const double delta) const {
	
	/* evalutate the slices of both over
	 * and under approximations (the 
	 * approximation type done during 
	 * the computation depends on the 
	 * type of the ASet).
	 */
	ASet *over=(this->over_appr)->flow_slice_p(field, delta);
	
	ASet *under=(this->under_appr)->flow_slice_p(field, delta);

	Set out(over,under);

	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;
}

#ifdef FASTER_BUT_DIRTIER
Set& Set::flow_slice_on_obj(const VectorField &field,
		const double delta) {

	(this->over_appr)->flow_slice_on_obj(field, delta);
	(this->under_appr)->flow_slice_on_obj(field, delta);

	return (*this);
}
#endif
		 
Set* Set::flow_p(const VectorField &field, const double delta) const {
	
	/* evalutate the flow from both 
	 * over and under approximations
	 * (the approximation type done 
	 * during the computation depends 
	 * on the type of the ASet).
	 */
	ASet *over=(this->over_appr)->flow_p(field, delta);
	
	ASet *under=(this->under_appr)->flow_p(field, delta);

	Set *out=new Set(over,under);

	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;
}

Set Set::flow(const VectorField &field, const double delta) const {
	
	/* evalutate the flow from both 
	 * over and under approximations
	 * (the approximation type done 
	 * during the computation depends 
	 * on the type of the ASet).
	 */
	ASet *over=(this->over_appr)->flow_p(field, delta);
	
	ASet *under=(this->under_appr)->flow_p(field, delta);

	Set out(over,under);

	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;
}

#ifdef FASTER_BUT_DIRTIER
Set& Set::flow_on_obj(const VectorField &field,
		const double delta) {

	(this->over_appr)->flow_on_obj(field, delta);
	(this->under_appr)->flow_on_obj(field, delta);

	return (*this);
}
#endif

		
ASet* Set::over_approximation() const {
	return this->over_appr;
}
		 
ASet* Set::under_approximation() const {
	return this->under_appr;
}
		 
HSet::HSet(const HSet &ts): Set((Set)ts){
	this->location=ts.location;
}

HSet::HSet(Location *l): Set(){
	this->location=l;
}
		
HSet::HSet(const Set &set, Location *l) : Set(set) {

	this->location=l;
}

HSet::HSet(ASet *over, ASet *under, Location *l): Set(over,under) {

	this->location=l;
}

	
HSet::HSet(int max_basic_sets, Location *l) : Set() {

	this->location=l;
	this->set_max_basic_sets(max_basic_sets);
}
		
HSet::~HSet(){
	
	delete (this->under_appr);
	delete (this->over_appr);

	/* useless, but I prefer to do it */
	this->under_appr=NULL;
	this->over_appr=NULL;
	this->location=NULL;
}
	
HSet& HSet::operator=(const HSet& ts) {

	(*this->under_appr)=(*ts.under_appr);
	(*this->over_appr)=(*ts.over_appr);

	this->location=ts.location;
	
	return (*this);
}

bool HSet::same_location(const HSet &hset) const{
	return (this->location==hset.location);
}
	
HSet* HSet::apply(const Map &M) const {
	/* evalutate the application of the 
	 * map \a M (the approximation type 
	 * done during the computation 
	 * depends on the type of the 
	 * ASet).
	 */
	ASet *over=(this->over_appr)->apply(M);
	
	ASet *under=(this->under_appr)->apply(M);

	HSet *out=new HSet(over,under,this->location);

	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;
}

HSet HSet::operator<<(const Map &M) const {
	/* evalutate the application of the 
	 * map \a M (the approximation type 
	 * done during the computation 
	 * depends on the type of the 
	 * ASet).
	 */
	ASet *over=(this->over_appr)->apply(M);
	
	ASet *under=(this->under_appr)->apply(M);

	HSet out(over,under,this->location);

	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;
}
		
#ifdef FASTER_BUT_DIRTIER
HSet& HSet::apply_on_obj(const Map &M) {
	/* evalutate the application of the 
	 * map \a M (the approximation type 
	 * done during the computation 
	 * depends on the type of the 
	 * ASet).
	 */
	(this->over_appr)->apply_on_obj(M);
	
	(this->under_appr)->apply_on_obj(M);

	return (*this);
}
#endif

HSet* HSet::intersect(const BasicSet &set) {
	
	ASet *under=(this->under_appr)->intersect(set);
	ASet *over=(this->over_appr)->intersect(set);

	HSet *out= new HSet(over,under,this->location);
	
	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;
}

HSet HSet::operator&&(const BasicSet &set) {

	ASet *under=(this->under_appr)->intersect(set);
	ASet *over=(this->over_appr)->intersect(set);

	HSet out(over,under,this->location);
	
	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;
}

HSet* HSet::intersect(const ASet &set) const {

	/* evalutate the intersections between
	 * over-approximation and between
	 * under-approxiamtion (the type 
	 * of approximation done during the
	 * computation depends on the type
	 * of the ASet).
	 */
	ASet *over=(this->over_appr)->intersect(set);
	
	ASet *under=(this->under_appr)->intersect(set);

	HSet *out=new HSet(over,under,this->location);

	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;
}

HSet HSet::operator&&(const ASet &set) const {

	/* evalutate the intersections between
	 * over-approximation and between
	 * under-approxiamtion (the type 
	 * of approximation done during the
	 * computation depends on the type
	 * of the ASet).
	 */
	ASet *over=(this->over_appr)->
		intersect(set);
	
	ASet *under=(this->under_appr)->
		intersect(set);

	HSet out(over,under,this->location);

	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;

}

HSet* HSet::join(const ASet &set) const {

	/* evalutate the unions between
	 * over-approximation and between
	 * under-approxiamtion (the type 
	 * of approximation done during the
	 * computation depends on the type
	 * of the ASet).
	 */
	ASet *over=(this->over_appr)->join(set);
	
	ASet *under=(this->under_appr)->join(set);

	HSet *out=new HSet(over,under,this->location);

	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;
}

HSet HSet::operator||(const ASet &set) const {

	/* evalutate the unions between
	 * over-approximation and between
	 * under-approxiamtion (the type 
	 * of approximation done during the
	 * computation depends on the type
	 * of the ASet).
	 */
	ASet *over=(this->over_appr)->join(set);
	
	ASet *under=(this->under_appr)->join(set);

	HSet out(over,under,this->location);

	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;

}

HSet* HSet::intersect(const Set &set) const {
	
	/* evalutate the intersections between
	 * over-approximations and between
	 * under-approxiamtions (the type 
	 * of approximation done during the
	 * computation depends on the type
	 * of the ASet).
	 */
	ASet *over=(this->over_appr)->
		intersect(*(set.over_appr));
	
	ASet *under=(this->under_appr)->
		intersect(*(set.under_appr));

	HSet *out=new HSet(over,under,this->location);

	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;
}



HSet HSet::operator&&(const Set &set) const {

	/* evalutate the intersections between
	 * over-approximations and between
	 * under-approxiamtions (the type 
	 * of approximation done during the
	 * computation depends on the type
	 * of the ASet).
	 */
	ASet *over=(this->over_appr)->
		intersect(*(set.over_appr));
	
	ASet *under=(this->under_appr)->
		intersect(*(set.under_appr));

	HSet out(over,under,this->location);

	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;

}

HSet* HSet::join(const Set &set) const {

	/* evalutate the unions between
	 * over-approximations and between
	 * under-approxiamtions (the type 
	 * of approximation done during the
	 * computation depends on the type
	 * of the ASet).
	 */
	ASet *over=(this->over_appr)->
		join(*(set.over_appr));
	
	ASet *under=(this->under_appr)->
		join(*(set.under_appr));

	HSet *out=new HSet(over,under,this->location);

	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;
}

HSet HSet::operator||(const Set &set) const {

	/* evalutate the unions between
	 * over-approximations and between
	 * under-approxiamtions (the type 
	 * of approximation done during the
	 * computation depends on the type
	 * of the ASet).
	 */
	ASet *over=(this->over_appr)->
		join(*(set.over_appr));
	
	ASet *under=(this->under_appr)->
		join(*(set.under_appr));

	HSet out(over,under,this->location);

	/* useless, but... "paranoia is good" */
	over=NULL;
	under=NULL;
	
	return out;

}


HSet* HSet::intersect(const HSet &set) const {
	if (this->same_location(set)) {
		return (this->intersect((Set)set));
	} else {
		return (new HSet(this->location));
	}
}

HSet HSet::operator&&(const HSet &set) const {
	if (this->same_location(set)) {
		return ((*this)&&((Set)set));
	} else {
		HSet out(this->location);
		
		return (out);
	}
}
	
HSet* HSet::join(const HSet &set) const {
	if (this->same_location(set)) {
		return (this->join((Set)set));
	} else {
		/* TODO: manage error code */

		return NULL;
	}
}
		
HSet HSet::operator||(const HSet &set) const {
	if (this->same_location(set)) {
		return ((*this)||((Set)set));
	} else {
		/* TODO: manage error code */

		HSet out(this->location);
		
		return (out);
	}
}


HSet* HSet::flow_slice_p(const double delta) const {
	
	const VectorField *field=(this->location)->get_vectorfield();

	/* evalutate the slices of both over
	 * and under approximations (the 
	 * approximation type done during 
	 * the computation depends on the 
	 * type of the ASet).
	 */
	ASet *over=(this->over_appr)->flow_slice_p(*field, delta);
	
	ASet *under=(this->under_appr)->flow_slice_p(*field, delta);

	HSet *out=new HSet(over,under,this->location);

	/* useless, but... "paranoia is good" */
	field=NULL;
	over=NULL;
	under=NULL;
	
	return out;
}

HSet HSet::flow_slice(const double delta) const {
	
	const VectorField *field=(this->location)->get_vectorfield();

	/* evalutate the slices of both over
	 * and under approximations (the 
	 * approximation type done during 
	 * the computation depends on the 
	 * type of the ASet).
	 */
	ASet *over=(this->over_appr)->flow_slice_p(*field, delta);
	
	ASet *under=(this->under_appr)->flow_slice_p(*field, delta);

	HSet out(over,under,this->location);

	/* useless, but... "paranoia is good" */
	field=NULL;
	over=NULL;
	under=NULL;
	
	return out;
}

#ifdef FASTER_BUT_DIRTIER
HSet& HSet::flow_slice_on_obj(const double delta) {

	const VectorField *field=(this->location)->get_vectorfield();
	
	(this->over_appr)->flow_slice_on_obj(*field, delta);
	(this->under_appr)->flow_slice_on_obj(*field, delta);
	
	/* useless, but... "paranoia is good" */
	field=NULL;

	return (*this);
}
#endif
		 
HSet* HSet::flow_p(const double delta) const {
	
	const VectorField *field=(this->location)->get_vectorfield();
	
	/* evalutate the flow from both 
	 * over and under approximations
	 * (the approximation type done 
	 * during the computation depends 
	 * on the type of the ASet).
	 */
	ASet *over=(this->over_appr)->flow_p(*field, delta);
	
	ASet *under=(this->under_appr)->flow_p(*field, delta);

	HSet *out=new HSet(over,under,this->location);

	/* useless, but... "paranoia is good" */
	field=NULL;
	over=NULL;
	under=NULL;
	
	return out;
}

HSet HSet::flow(const double delta) const {
	
	const VectorField *field=(this->location)->get_vectorfield();
	
	/* evalutate the flow from both 
	 * over and under approximations
	 * (the approximation type done 
	 * during the computation depends 
	 * on the type of the ASet).
	 */
	ASet *over=(this->over_appr)->flow_p(*field, delta);
	
	ASet *under=(this->under_appr)->flow_p(*field, delta);

	HSet out(over,under,this->location);

	/* useless, but... "paranoia is good" */
	field=NULL;
	over=NULL;
	under=NULL;
	
	return out;
}

#ifdef FASTER_BUT_DIRTIER
HSet& HSet::flow_on_obj(const double delta) {

	const VectorField *field=(this->location)->get_vectorfield();
	
	(this->over_appr)->flow_on_obj(*field, delta);
	(this->under_appr)->flow_on_obj(*field, delta);

	/* useless, but I prefer to do it */
	field=NULL;

	return (*this);
}
#endif

const Location *HSet::get_location() const {
	return this->location;
}
		 
Set* HSet::get_Set() const {
	Set *set=new Set(this->over_appr,this->under_appr); 

	return set;
}
		 
HSet& HSet::set_location(Location *l) {
	this->location=l;

	return(*this);
}

bool HSet::does_activate(const LeavingArc &e) const {
	
	const ASet *activ=e.get_activation();

	/* check if the current set intersects the
	 * activation regione of \a e . */
	
	/* NOTE: The method must NOT destroy the
	 * \a activ object because it is a member
	 * of \a e . */
	
	return (this->does_intersect(*activ));
	
}

bool HSet::does_activate(const Arc &e) const {

	if ((*(e.get_source()))!=(*(this->location))) {
		return false;
	}
	
	const ASet *activ=e.get_activation();

	/* check if the current set intersects the
	 * activation regione of \a e . */

	/* NOTE: The method must NOT destroy the
	 * \a activ object because it is a member
	 * of \a e . */

	return (this->does_intersect(*activ));
	
}

HSet* HSet::jump_p(const LeavingArc &e) const {
	
	if (!(this->does_activate(e))) {
		return (new HSet(*this));
	}
	
	/* get the activation region and the 
	 * reset of \a e */
	const ASet *activ=e.get_activation();
	const Map *reset=e.get_reset();

	/* evaluate the intersection of
	 * the leaving arc's activation
	 * region and the current set */
	HSet *a_set=this->intersect(*activ);


#ifdef FASTER_BUT_DIRTIER
	/* evaluate the reset */
	a_set->apply_on_obj(*reset);

	HSet *out=a_set;
#else
	/* evaluate the reset */
	HSet *out=a_set->apply(*reset);

	delete a_set;
#endif
	/* useless, but I prefer to do it */
	a_set=NULL;
	activ=NULL;
	reset=NULL;

	/* set as new location the leaving arc's
	 * destination */
	out->set_location(e.get_destination());

	return out;
}

HSet* HSet::jump_p(const Arc &e) const {
	
	if (!(this->does_activate(e))) {
		return (new HSet(*this));
	}

	return (this->jump_p((LeavingArc)e));
}

HSet HSet::jump(const LeavingArc &e) const {
	
	if (!(this->does_activate(e))) {
		return (HSet(*this));
	}
	
	/* get the activation region and the 
	 * reset of \a e */
	const ASet *activ=e.get_activation();
	const Map *reset=e.get_reset();

	/* evaluate the intersection of
	 * the leaving arc's activation
	 * region and the current set */
	HSet *a_set=this->intersect(*activ);


#ifdef FASTER_BUT_DIRTIER
	/* evaluate the reset */
	a_set->apply_on_obj(*reset);

	HSet out(*a_set);
#else
	/* evaluate the reset */
	HSet out((*a_set)<<(*reset));
#endif
	delete a_set;
	
	/* useless, but I prefer to do it */
	a_set=NULL;
	activ=NULL;
	reset=NULL;

	/* set as new location the leaving arc's
	 * destination */
	out.set_location(e.get_destination());

	return out;
}

HSet HSet::jump(const Arc &e) const {
	
	if (!(this->does_activate(e))) {
		return (HSet(*this));
	}

	return (this->jump((LeavingArc)e));
}

