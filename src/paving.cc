/***************************************************************************
 *            paving.cc
 *
 *  Copyright 2008  Pieter Collins
 * 
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
 
#include <cassert>
#include "paving.h" 

namespace Ariadne {

typedef Vector<Interval> Box;

//--------------------------------------------------------------------------
// Misc Utility functions

bool intersection(Interval& ivlr, 
                  const Interval& ivl1 , 
                  const Interval& ivl2)
{
    double upb  = min(ivl1.upper(),ivl2.upper());
        double lowb = max(ivl1.lower(),ivl2.lower());
            if (lowb>upb) { return false; }
                ivlr = Interval(lowb,upb);
                    return true;
                        }

bool intersection(IntervalVector& r, 
                  const IntervalVector& a,
                  const IntervalVector& b)
{
    assert(a.size()==b.size());
        r = a;
            for (uint i=0; i!=a.size(); ++i) {
                if(!intersection(r[i],a[i],b[i])) { return false; } 
                    }
                return true;
                    }

IntervalVector split(const IntervalVector& iv, int j, bool lu) 
{
    IntervalVector r=iv;
        if(lu==0) {
            r[j] = Interval(iv[j].lower(),iv[j].midpoint());
                } else {
            r[j] = Interval(iv[j].midpoint(),iv[j].upper());
                }    
            return r;
                }



double diam (const IntervalVector& x, int& c)
{
    double res = diam(x[0]); 
        c = 0;             
            
            for (uint i=1; i!=x.size(); ++i) {
                if (diam(x[i])>res) { 
                    res = diam(x[i]); c = i; } }    
                
                return res;
                    }


double diam (const IntervalVector& x)
{
    double res = diam(x[0]); 
        int c = 0;             
            
            for (uint i=1; i!=x.size(); ++i) {
                if (diam(x[i])>res) { 
                    res = diam(x[i]); c = i; } }    
                
                return res;
                    }



//--------------------------------------------------------------------------
// Node class


Node::Node(const IntervalVector& v, tribool tb) {
    this->theBox=v;
        this->theValue=tb;
            this->leftChild=NULL; 
                this->rightChild=NULL;
                    }

Node::Node(const Node& n)
{ 
    this->theBox=n.theBox;
        this->theValue=n.theValue;
            //recursion on the children
            if (n.leftChild) { this->leftChild=new Node(*n.leftChild); }
            else this->leftChild=NULL;
                if (n.rightChild) { this->rightChild=new Node(*n.rightChild); }
                else this->rightChild=NULL; 
                    }

Node::~Node() {
    delete leftChild; 
        delete rightChild;
            }

std::ostream& operator<< (std::ostream& os, const Subpaving a)
{
    static uint indent=0;
        if (a==NULL) { return os << "NULL" << std::flush; } 
            for(uint i=0; i!=indent; ++i) { os << " "; }
                os << (a->theValue)<<": "<<(a->theBox) << std::endl; 
                    //os << int(a->theValue)<<": "<<(a->theBox) << endl; 
                    if (is_leaf(a)) { }
                    else { ++indent; os << (a->leftChild); os << (a->rightChild); --indent; }
                        return os;
                            }



const IntervalVector& box(const Subpaving a) { 
    return a->theBox; 
        }    


bool is_empty(const Subpaving a) {
    assert(a!=NULL);
        return a->theValue == false;
            }


bool is_full(const Subpaving a) {
    assert(a!=NULL);
        return a->theValue == true;
            }


bool is_leaf(const Subpaving a) {
    assert(a!=NULL);
        assert((a->leftChild==NULL) == (a->rightChild==NULL));
            return(!a->leftChild && !a->rightChild);
                }

bool is_null(const Subpaving a) {
    return a==NULL;
        }


uint size (const Subpaving a)
{
    assert(a!=NULL);
        if (is_leaf(a)) { 
            assert(a->theValue!=indeterminate);
                return (a->theValue==true) ? 1u : 0u; 
                    } else {
            return size(a->leftChild) + size(a->rightChild);
                }
            }


double volume (const Subpaving a)
{
    assert(a!=NULL);
        
        double vol=0;
            if (is_empty(a)) { return 0; }
                
                if (is_full(a)) {
                    return volume(box(a));
                        }
                    
                    assert(!is_leaf(a));
                        vol += volume(a->leftChild);
                            vol += volume(a->rightChild);
                                
                                return (vol);
                                    }


tribool inside(const IntervalVector& z, const Subpaving X)
{
    // v is assumed not to be empty
    assert(X!=NULL);
        
        if (is_empty(X)) { return false; }
            
            IntervalVector r;
                
                if (is_leaf(X)) {
                    if (inside(z,box(X))) { return true; }
                        if (!intersection(r,z,box(X))) { return false; }
                            }
                    
                    tribool Ltest=true; 
                        tribool Rtest=true;
                            
                            IntervalVector Lz,Rz;
                                
                                if (!is_empty(X->leftChild)&&!is_empty(X->rightChild)) {
                                    if (intersection(Lz,z,box(X->leftChild))) {
                                        Ltest = inside(Lz,X->leftChild);
                                            if (intersection(Rz,z,box(X->rightChild))) {
                                                Rtest = inside(Rz,X->rightChild);
                                                    if (Ltest==Rtest) {
                                                        return Ltest; 
                                                            } else {
                                                        return indeterminate; 
                                                            }
                                                        } else {
                                                return Ltest; 
                                                    }
                                                } else {
                                        if (intersection(Rz,z,box(X->rightChild))) {
                                            return inside(Rz,X->rightChild);
                                                } else { 
                                            return true;
                                                }
                                            }
                                        }
                                else if (!is_empty(X->leftChild)) {
                                    if (!intersection(Lz,z,box(X->leftChild))) {
                                        return false; }
                                        Ltest = inside(Lz,X->leftChild);
                                            if (!(Lz==z)) { 
                                                Rtest = false; }
                                            else { 
                                                return Ltest; }
                                                }
                                else {
                                    if (!intersection(Rz,z,box(X->rightChild))) {
                                        return false; }
                                        Rtest = inside(Rz,X->rightChild);
                                            if (!(Rz==z)) {
                                                Ltest = false; }
                                            else {
                                                return Rtest; }
                                                }
                                    
                                    if (Ltest==Rtest) {
                                        return Ltest; 
                                            } else {
                                        return indeterminate; 
                                            }
                                        }


void expand(Subpaving A)
{
    assert(A && is_leaf(A));
        
        int comp; diam(box(A),comp);
                      A->leftChild = new Node(split(box(A),comp,0),A->theValue);
                          A->rightChild = new Node(split(box(A),comp,1),A->theValue);
                              }


Subpaving reunite(Subpaving lChild, Subpaving rChild, const IntervalVector& x)
{
    assert((lChild!=NULL) && (rChild!=NULL));
        
        {
            int maxdiamcomp; diam(x,maxdiamcomp);
                                 assert(box(lChild)==split(x,maxdiamcomp,0)); 
                                     assert(box(rChild)==split(x,maxdiamcomp,1));
                                         }
            
            if(is_leaf(lChild) && is_leaf(rChild) && lChild->theValue == rChild->theValue) {
                Subpaving new_node = new Node(x,lChild->theValue);
                    delete lChild; 
                        delete rChild; 
                            return new_node; 
                                }
                
                Subpaving result = new Node(x,indeterminate);
                    result->leftChild = lChild; 
                        result->rightChild = rChild;
                            return result;
                                }


//--------------------------------------------------------------------------
// sivia components

Subpaving sivia (IntervalPredicate pred, const IntervalVector& box, double eps, tribool aprx)
{
    tribool test = pred(box);
        //cout << "test=" << test << " box="<<box << flush;
        
        if ( test==true ) { return new Node(box,false); }
            if ( test==false ) { return new Node(box,true); }
                
                int subdv; double diam = Ariadne::diam(box,subdv);
                               if ( diam<eps ) { 
                                   return (aprx==true) ? new Node(box,true)
                                       : new Node(box,false); 
                                       }
                                   
                                   return reunite(sivia(pred,split(box,subdv,0),eps,aprx),
                                                  sivia(pred,split(box,subdv,1),eps,aprx),box);
                                       }




std::ostream& write(std::ostream& os, const Node* n) {
    return os << n->theValue << " " << n->theBox << " "
              << static_cast<const void*>(n->leftChild) << " "
              << static_cast<const void*>(n->rightChild) << std::endl;
        }



//--------------------------------------------------------------------------
// image components

void mince(Subpaving A, double eps)
{
    assert(A!=NULL);
        if (is_empty(A)) { return; }
            if (is_leaf(A)) {
                if(diam(box(A))<eps) { return; }
                    expand(A); 
                        }
                mince(A->leftChild,eps);
                    mince(A->rightChild,eps);
                        }


void flatten(const Subpaving a, ImageList& list) 
{
    if(!is_leaf(a)) { 
        flatten(a->leftChild,list); 
            flatten(a->rightChild,list); 
                } else if(is_full(a)) { 
        list.push_back(box(a)); 
            } 
        return;
            }


void evaluate(IntervalFunction f, 
              Subpaving A, 
              ImageList& list,
              IntervalVector& hull)
{
    if (is_leaf(A)) {
        IntervalVector image(f(box(A)));
            
            if (list.empty()) { hull = image; }
            else { hull = Ariadne::hull(hull,image); }
                
                list.push_back(Box(image));
                    return;
                        }
        
        evaluate(f, A->leftChild, list, hull);
            evaluate(f, A->rightChild, list, hull);
                }

void contract (Subpaving a)
{
    if(is_leaf(a)) { return; }
        contract(a->leftChild);
            contract(a->rightChild);
                if(is_leaf(a->leftChild) && is_leaf(a->rightChild)
                   && a->leftChild->theValue == a->rightChild->theValue) 
                    {
                        a->theValue = a->leftChild->theValue;
                            delete a->leftChild;
                                delete a->rightChild;
                                    a->leftChild=NULL;
                                        a->rightChild=NULL;
                                            }
                    }

Subpaving regularize (const ImageList& list, const IntervalVector& hull, double eps)
{
    if (list.empty()) { return new Node(hull,false); }
        
        if (inside(hull,list.front())) { return new Node(hull,true); }
            
            int maxdiamcomp;
                if (diam(hull,maxdiamcomp)<eps) { return new Node(hull,true); }
                    
                    IntervalVector lefthull = split(hull,maxdiamcomp,0);
                        IntervalVector righthull = split(hull,maxdiamcomp,1);
                            
                            ImageList leftlist,rightlist;
                                
                                IntervalVector Inter;
                                    
                                    ImageList::const_iterator curr=list.begin();
                                        while (curr!=list.end()) {
                                            if (intersection(Inter,*curr,lefthull)) {
                                                leftlist.push_back( Box(Inter) );
                                                    }
                                                if (intersection(Inter,*curr,righthull)) {
                                                    rightlist.push_back( Box(Inter) );
                                                        }
                                                    ++curr;
                                                        }
                                            
                                            return reunite(regularize(leftlist,lefthull,eps),
                                                           regularize(rightlist,righthull,eps),hull);
                                                }


Subpaving image(IntervalFunction f, Subpaving A, double eps)
{
    
    ImageList images;
        IntervalVector hull;
            
            mince(A,eps);
                evaluate(f,A,images,hull);
                    Subpaving result = regularize(images, hull, eps);
                        return result;
                            }





Paving mince(const Paving& pv, double eps) {
    Paving result(pv);
        mince(result.root(),eps);
            return result;
                }

Paving contract(const Paving& pv) {
    Paving result(pv);
        contract(result.root());
            return result;
                }

ImageList flatten(const Paving& pv) {
    ImageList list;
        flatten(pv.root(),list);
            return list;
                }



Paving inner_set (IntervalPredicate pred, const IntervalVector& box, double eps)
{
    return Paving(sivia(pred,box,eps,false));
        }

Paving outer_set (IntervalPredicate pred, const IntervalVector& box, double eps)
{
    return Paving(sivia(pred,box,eps,true));
        }

Paving outer_set(const ImageList& set, const IntervalVector& box, double eps)
{
    return Paving(regularize(set,box,eps)); 
        }

Paving image(IntervalFunction f, const Paving& pav, double eps) {
    return Paving(image(f,new Node(*pav.root()),eps));
        }

} // namespace Ariadne
