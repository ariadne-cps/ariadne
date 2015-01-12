/***************************************************************************
 *            polynomial.tcc
 *
 *  Copyright 2008-15  Pieter Collins
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

namespace Ariadne {

template<class X> Polynomial<X>::Polynomial(SizeType as) : _expansion(as) { }

template<class X> template<class XX> Polynomial<X>::Polynomial(const Polynomial<XX>& p) : _expansion(p._expansion) { }

template<class X> template<class XX> Polynomial<X>::Polynomial(const Expansion<XX>& e) : _expansion(e) { }




template<class X> Polynomial<X> Polynomial<X>::create_zero() const {
    return Polynomial<X>(this->argument_size()); }


template<class X> Polynomial<X> Polynomial<X>::constant(SizeType as, const X& c) {
    Polynomial<X> r(as); r[MultiIndex::zero(as)]=c; return r; }

template<class X> Polynomial<X> Polynomial<X>::variable(SizeType as, SizeType j) {
    ARIADNE_ASSERT(j<as); Polynomial<X> r(as); r[MultiIndex::unit(as,j)]=1; return r; }
template<class X> Polynomial<X> Polynomial<X>::coordinate(SizeType as, SizeType j) {
    ARIADNE_ASSERT(j<as); Polynomial<X> r(as); r[MultiIndex::unit(as,j)]=1; return r; }


template<class X> Vector<Polynomial<X>> Polynomial<X>::variables(SizeType as) {
    Vector<Polynomial<X>> r(as); for(SizeType i=0; i!=as; ++i) { r[i]=variable(as,i); } return r; }
template<class X> Polynomial<X>& Polynomial<X>::operator=(const X& x) { this->_expansion.clear(); this->_expansion.append(MultiIndex(this->argument_size()),x); return *this; }





template<class X> template<class XX> Bool Polynomial<X>::operator==(const Polynomial<XX>& p) const {
    const_cast<Polynomial<X>*>(this)->cleanup(); const_cast<Polynomial<XX>&>(p).cleanup();
    return this->_expansion==p._expansion; }

template<class X> template<class XX> Bool Polynomial<X>::operator!=(const Polynomial<XX>& p) const {
    return !(*this==p); }






template<class X> SizeType Polynomial<X>::argument_size() const { return this->_expansion.argument_size(); }

template<class X> SizeType Polynomial<X>::number_of_nonzeros() const { return this->_expansion.number_of_nonzeros(); }

template<class X> SizeType Polynomial<X>::degree() const { return this->_expansion.degree(); }

template<class X> const X& Polynomial<X>::value() const { return this->_expansion[MultiIndex::zero(this->argument_size())]; }

template<class X> X& Polynomial<X>::operator[](const MultiIndex& a) { return this->_expansion.at(a); }

template<class X> const X& Polynomial<X>::operator[](const MultiIndex& a) const { return this->_expansion.get(a); }

template<class X> const Expansion<X>& Polynomial<X>::expansion() const { return this->_expansion; }

template<class X> Expansion<X>& Polynomial<X>::expansion() { return this->_expansion; }






template<class X> typename Polynomial<X>::Iterator Polynomial<X>::begin() { return this->_expansion.begin(); }

template<class X> typename Polynomial<X>::Iterator Polynomial<X>::end() { return this->_expansion.end(); }

template<class X> typename Polynomial<X>::Iterator Polynomial<X>::find(const MultiIndex& a) { return this->_expansion.find(a); }

template<class X> typename Polynomial<X>::ConstIterator Polynomial<X>::begin() const { return this->_expansion.begin(); }

template<class X> typename Polynomial<X>::ConstIterator Polynomial<X>::end() const { return this->_expansion.end(); }

template<class X> typename Polynomial<X>::ConstIterator Polynomial<X>::find(const MultiIndex& a) const { return this->_expansion.find(a); }







template<class X> Void Polynomial<X>::append(const MultiIndex& a, const X& c) { this->_expansion.append(a,c); }

template<class X> Void Polynomial<X>::insert(const MultiIndex& a, const X& c) { this->_expansion.insert(a,c); }

template<class X> Void Polynomial<X>::reserve(SizeType n) { this->_expansion.reserve(n); }

template<class X> Void Polynomial<X>::erase(Iterator iter) { this->_expansion.erase(iter); }

template<class X> Void Polynomial<X>::clear() { this->_expansion.clear(); }







template<class X>
Polynomial<X>::Polynomial(InitializerList< PairType<InitializerList<Int>,X> > lst)
    : _expansion(lst)
{
    this->cleanup();
}

template<class X>
Polynomial<X>::Polynomial(SizeType as, DegreeType deg, InitializerList<X> lst)
    : _expansion(as,deg,lst)
{
    this->cleanup();
}

template<class X>
Polynomial<X>&
Polynomial<X>::differentiate(SizeType j) {
    for(typename Polynomial<X>::Iterator iter=this->begin(); iter!=this->end(); ++iter) {
        MultiIndex& a=iter->key();
        X& c=iter->data();
        c*=a[j];
        if(a[j]!=0u) { ++a[j]; }
    }
    return *this;
}

template<class X>
Polynomial<X>&
Polynomial<X>::antidifferentiate(SizeType j) {
    for(typename Polynomial<X>::Iterator iter=this->begin(); iter!=this->end(); ++iter) {
        MultiIndex& a=iter->key();
        X& c=iter->data();
        ++a[j];
        c/=a[j];
    }
    return *this;
}

template<class X>
Polynomial<X>& Polynomial<X>::truncate(DegreeType d) {
    Polynomial<X> r(this->argument_size());
    for(typename Polynomial<X>::ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        if(iter->key().degree()<=d && iter->data()!=0) {
            r.append(iter->key(),iter->data());
        }
    }
    this->swap(r);
    return *this;
}

template<class X>
Void Polynomial<X>::cleanup()
{
    Polynomial<X>* self=const_cast<Polynomial<X>*>(this);
    self->_expansion.reverse_lexicographic_sort();
    Iterator new_end=unique_key(self->_expansion.begin(), self->_expansion.end(), std::plus<X>());
    self->_expansion.resize(new_end-self->_expansion.begin());
}

template<class X> inline
Void Polynomial<X>::check() const
{
    this->_expansion.check();
}


template<class X> inline Polynomial<X> operator+(const Polynomial<X>& p) {
    return p; }

template<class X> inline Polynomial<X> operator-(const Polynomial<X>& p) {
    Polynomial<X> r(p.argument_size());  typedef typename Polynomial<X>::ConstIterator Iter;
    for(Iter iter=p.begin(); iter!=p.end(); ++iter) { r[iter->key()]-=iter->data(); } return r; }


template<class X> inline Polynomial<X> operator+(const Polynomial<X>& p, const X& c) {
    Polynomial<X> r(p); r[MultiIndex(p.argument_size())]+=c; return r; }

template<class X> inline Polynomial<X> operator-(const Polynomial<X>& p, const X& c) {
    Polynomial<X> r(p); r[MultiIndex(p.argument_size())]-=c; return r; }

template<class X> inline Polynomial<X> operator*(const Polynomial<X>& p, const X& c) {
    if(c==0) { return Polynomial<X>(p.argument_size()); }
    Polynomial<X> r(p); typedef typename Polynomial<X>::Iterator Iter;
    for(Iter iter=r.begin(); iter!=r.end(); ++iter) { iter->data()*=c; } return r; }

template<class X> inline Polynomial<X> operator/(const Polynomial<X>& p, const X& c) {
    Polynomial<X> r(p); typedef typename Polynomial<X>::Iterator Iter;
    for(Iter iter=r.begin(); iter!=r.end(); ++iter) { iter->data()/=c; } return r; }


template<class X> inline Polynomial<X> operator+(const X& c, const Polynomial<X>& p) {
    return p+c; }

template<class X> inline Polynomial<X> operator-(const X& c, const Polynomial<X>& p) {
    return (-p)+c; }

template<class X> inline Polynomial<X> operator*(const X& c, const Polynomial<X>& p) {
    return p*c; }


template<class X> inline Polynomial<X> operator+(const Polynomial<X>& p1, const Polynomial<X>& p2) {
    ARIADNE_ASSERT(p1.argument_size()==p2.argument_size());
    Polynomial<X> r(p1);  typedef typename Polynomial<X>::ConstIterator Iter;
    for(Iter iter=p2.begin(); iter!=p2.end(); ++iter) { r[iter->key()]+=iter->data(); } return r; }

template<class X> inline Polynomial<X> operator-(const Polynomial<X>& p1, const Polynomial<X>& p2) {
    ARIADNE_ASSERT(p1.argument_size()==p2.argument_size());
    Polynomial<X> r(p1);  typedef typename Polynomial<X>::ConstIterator Iter;
    for(Iter iter=p2.begin(); iter!=p2.end(); ++iter) { r[iter->key()]-=iter->data(); } return r; }

template<class X> inline Polynomial<X> operator*(const Polynomial<X>& p1, const Polynomial<X>& p2) {
    ARIADNE_ASSERT(p1.argument_size()==p2.argument_size());
    Polynomial<X> r(p1.argument_size());
    typedef typename Polynomial<X>::ConstIterator Iter;
    for(Iter iter1=p1.begin(); iter1!=p1.end(); ++iter1) {
        for(Iter iter2=p2.begin(); iter2!=p2.end(); ++iter2) {
            MultiIndex a=iter1->key()+iter2->key();
            r[a]+=iter1->data()*iter2->data(); } } return r; }


template<class X> inline Polynomial<X> sqr(const Polynomial<X>& p) {
    return p*p; }


template<class X> inline Polynomial<X> pow(const Polynomial<X>& p, Nat m) {
    Polynomial<X> r=Polynomial<X>::constant(p.argument_size(),1.0); Polynomial<X> q(p);
    while(m) { if(m%2) { r=r*q; } q=q*q; m/=2; } return r; }


template<class X, class XX> inline Polynomial<X>& operator+=(Polynomial<X>& p, const Polynomial<XX>& q) {
    ARIADNE_ASSERT(p.argument_size()==q.argument_size());
    typedef typename Polynomial<X>::ConstIterator Iter;
    for(Iter iter=q.begin(); iter!=q.end(); ++iter) { p[iter->key()]+=iter->data(); } return p; }

template<class X, class XX> inline Polynomial<X>& operator-=(Polynomial<X>& p, const Polynomial<XX>& q) {
    ARIADNE_ASSERT(p.argument_size()==q.argument_size());
    typedef typename Polynomial<X>::ConstIterator Iter;
    for(Iter iter=q.begin(); iter!=q.end(); ++iter) { p[iter->key()]-=iter->data(); } return p; }


template<class X> inline Polynomial<X>& operator+=(Polynomial<X>& p, const X& c) {
    p[MultiIndex(p.argument_size())]+=c; return p; }

template<class X> inline Polynomial<X>& operator-=(Polynomial<X>& p, const X& c) {
    p[MultiIndex(p.argument_size())]-=c; return p; }

template<class X> inline Polynomial<X>& operator*=(Polynomial<X>& p, const X& c) {
    typedef typename Polynomial<X>::Iterator Iter;
    if(c==0) { p.clear(); }
    for(Iter iter=p.begin(); iter!=p.end(); ++iter) { iter->data()*=c; } return p; }

template<class X> inline Polynomial<X>& operator/=(Polynomial<X>& p, const X& c) {
    typedef typename Polynomial<X>::Iterator Iter;
    for(Iter iter=p.begin(); iter!=p.end(); ++iter) { iter->data()/=c; } return p; }


template<class X> inline Polynomial<X>& operator*=(Polynomial<X>& p, const Monomial<X>& m) {
    typedef typename Polynomial<X>::Iterator Iter;
    if(m.data()==0) { p.clear(); }
    for(Iter iter=p.begin(); iter!=p.end(); ++iter) { iter->key()+=m.key(); iter->data()*=m.data(); } return p; }

template<class X> inline Polynomial<X> operator+(const Polynomial<X>& p, double c) { return p+X(c); }
template<class X> inline Polynomial<X> operator-(const Polynomial<X>& p, double c) { return p-X(c); }
template<class X> inline Polynomial<X> operator*(const Polynomial<X>& p, double c) { return p*X(c); }
template<class X> inline Polynomial<X> operator/(const Polynomial<X>& p, double c) { return p/X(c); }
template<class X> inline Polynomial<X> operator+(double c, const Polynomial<X>& p) { return X(c)+p; }
template<class X> inline Polynomial<X> operator-(double c, const Polynomial<X>& p) { return X(c)-p; }
template<class X> inline Polynomial<X> operator*(double c, const Polynomial<X>& p) { return X(c)*p; }
template<class X> inline Polynomial<X>& operator+=(Polynomial<X>& p, double c) { return p+=X(c); }
template<class X> inline Polynomial<X>& operator-=(Polynomial<X>& p, double c) { return p-=X(c); }
template<class X> inline Polynomial<X>& operator*=(Polynomial<X>& p, double c) { return p*=X(c); }
template<class X> inline Polynomial<X>& operator/=(Polynomial<X>& p, double c) { return p/=X(c); }

template<class X> inline Polynomial<MidpointType<X>> midpoint(const Polynomial<X>& p) {
    Polynomial<MidpointType<X>> r(p.argument_size());
    for(typename Polynomial<X>::ConstIterator iter=p.begin(); iter!=p.end(); ++iter) {
        r.append(iter->key(),static_cast<MidpointType<X>>(midpoint(iter->data()))); }
    return r;
}

template<class X> inline Vector< Polynomial<MidpointType<X>> > midpoint(const Vector<Polynomial<X>>& p) {
    Vector< Polynomial<MidpointType<X>> > r(p.size());
    for(Nat i=0; i!=p.size(); ++i) {
        r[i]=midpoint(p[i]); }
    return r;
}


template<class X, class Y> inline
Y evaluate(const Polynomial<X>& p, const Vector<Y>& x)
{
    return evaluate(p.expansion(),x);
}


template<class X, class Y> inline
Vector<Y> evaluate(const Vector<Polynomial<X>>& p, const Vector<Y>& x)
{
    ARIADNE_ASSERT(p.size()>0 && p.zero_element().argument_size()==x.size());
    Y zero = x.zero_element(); zero*=0;
    Vector<Y> r(p.size(),zero);
    for(Nat i=0; i!=p.size(); ++i) {
        r[i]=evaluate(p[i],x);
    }
    return r;
}


template<class X>
Polynomial<X>
partial_evaluate(const Polynomial<X>& x, Nat k, const X& c)
{
    Polynomial<X> r(x.argument_size()-1);
    MultiIndex ra(r.argument_size());
    if(c==0) {
        for(typename Polynomial<X>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
            const MultiIndex& xa=xiter->key();
            MultiIndex::IndexType xak=xa[k];
            if(xak==0) {
                const X& xv=xiter->data();
                for(Nat i=0; i!=k; ++i) { ra[i]=xa[i]; }
                for(Nat i=k; i!=ra.size(); ++i) { ra[i]=xa[i+1]; }
                r.expansion().append(ra,xv);
            }
        }
    } else if(c==1) {
        Polynomial<X> s(x.argument_size()-1);
        Array< Polynomial<X> > p(x.degree()+1,Polynomial<X>(x.argument_size()-1));

        for(typename Polynomial<X>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
            const MultiIndex& xa=xiter->key();
            const X& xv=xiter->data();
            MultiIndex::IndexType xak=xa[k];
            for(Nat i=0; i!=k; ++i) { ra[i]=xa[i]; }
            for(Nat i=k; i!=ra.size(); ++i) { ra[i]=xa[i+1]; }
            assert(ra.degree()+xak==xa.degree());
            p[xak].expansion().append(ra,xv);
        }

        r=p[0];
        for(Nat i=1; i!=p.size(); ++i) {
            r+=p[i];
        }
    } else {
        Polynomial<X> s(x.argument_size()-1);
        Array< Polynomial<X> > p(x.degree()+1,Polynomial<X>(x.argument_size()-1));

        Array<X> cpowers(x.degree()+1);
        cpowers[0]=static_cast<X>(1); cpowers[1]=c;
        if(x.degree()>=2) { cpowers[2]=sqr(c); }
        for(Nat j=3; j<=x.degree(); ++j) {
            cpowers[j]=cpowers[j-2]*cpowers[2];
        }

        for(typename Polynomial<X>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
            const MultiIndex& xa=xiter->key();
            const X& xv=xiter->data();
            MultiIndex::IndexType xak=xa[k];
            for(Nat i=0; i!=k; ++i) { ra[i]=xa[i]; }
            for(Nat i=k; i!=ra.size(); ++i) { ra[i]=xa[i+1]; }
            assert(ra.degree()+xak==xa.degree());
            p[xak].expansion().append(ra,xv);
        }
        for(Nat i=1; i!=p.size(); ++i) {
            p[i]*=cpowers[i];
        }

        r=p[0];
        for(Nat i=1; i!=p.size(); ++i) {
            r+=p[i];
        }
    }

    return r;
}


template<class X> inline
Polynomial<X> compose(const Polynomial<X>& p, const Vector<Polynomial<X>>& q)
{
    return evaluate(p,q);
}

template<class X> inline
Vector<Polynomial<X>> compose(const Vector<Polynomial<X>>& p, const Vector<Polynomial<X>>& q)
{
    return evaluate(p,q);
}


template<class X> inline
Polynomial<X> embed(SizeType before_size, const Polynomial<X>& x, SizeType after_size)
{
    return Polynomial<X>(embed(before_size,x.expansion(),after_size));
}

template<class X>
Polynomial<X>
derivative(const Polynomial<X>& p, Nat j) {
    Polynomial<X> r(p.argument_size());
    MultiIndex ar(p.argument_size());
    for(typename Polynomial<X>::ConstIterator iter=p.begin(); iter!=p.end(); ++iter) {
        const MultiIndex& ap=iter->key();
        if(ap[j]>0) {
            const X& val=iter->data();
            ar=ap;
            ar[j]-=1;
            r.append(ar,val*static_cast<X>(ap[j]));
        }
    }
    return r;
}

template<class X>
Vector<Polynomial<X>>
derivative(const Polynomial<X>& p) {
    Vector<Polynomial<X>> r(p.argument_size(), Polynomial<X>(p.argument_size()) );
    for(Nat j=0; j!=r.size(); ++j) {
        r[j]=derivative(p,j);
    }
    return r;
}

template<class X>
Vector<Polynomial<X>>
derivative(const Vector<Polynomial<X>>& p, Nat j) {
    Vector<Polynomial<X>> r(p.size());
    for(Nat i=0; i!=p.size(); ++i) {
        r[i]=derivative(p[i],j);
    }
    return r;
}

template<class X>
Polynomial<X>
antiderivative(const Polynomial<X>& p, Nat j) {
    Polynomial<X> r(p.argument_size());
    MultiIndex ar(p.argument_size());
    for(typename Polynomial<X>::ConstIterator iter=p.begin(); iter!=p.end(); ++iter) {
        const MultiIndex& ap=iter->key();
        const X& val=iter->data();
        ar=ap;
        ++ar[j];
        r.append(ar,val/ar[j]);
    }
    return r;
}

template<class X>
Vector<Polynomial<X>>
antiderivative(const Vector<Polynomial<X>>& p, Nat j) {
    Vector<Polynomial<X>> r(p.size());
    for(Nat i=0; i!=p.size(); ++i) {
        r[i]=antiderivative(p[i],j);
    }
    return r;
}

template<class X>
Polynomial<X>
truncate(const Polynomial<X>& p, unsigned short d) {
    Polynomial<X> r(p.argument_size());
    for(typename Polynomial<X>::ConstIterator iter=p.begin(); iter!=p.end(); ++iter) {
        if(iter->key().degree()<=d && iter->data()!=0) {
            r.append(iter->key(),iter->data());
        }
    }
    return r;
}

template<class X>
Vector<Polynomial<X>>
truncate(const Vector<Polynomial<X>>& p, Nat d) {
    Vector<Polynomial<X>> r(p.size());
    for(Nat i=0; i!=p.size(); ++i) {
        r[i]=truncate(p[i],d);
    }
    return r;
}



template<class X>
Vector<Polynomial<X>> flow(const Vector<Polynomial<X>>& p, Nat order)
{
    Nat n=p.size();
    Vector<Polynomial<X>> p0(n,Polynomial<X>(n+1));
    for(Nat i=0; i!=n; ++i) { p0[i][MultiIndex::unit(n+1,i)]=1.0; }

    Vector<Polynomial<X>> r=p0;
    for(Nat k=0; k!=order; ++k) {
        r=compose(p,r);
        for(Nat i=0; i!=n; ++i) {
            r[i].cleanup();
            r[i]=truncate(r[i],k+1);
            r[i]=antiderivative(r[i],n);
            r[i]+=p0[i];
        }
    }

    return r;
}


template<class X, class Y> Vector<Polynomial<X>> operator*(const Polynomial<X>& p, const Vector<Y> e) {
    Vector<Polynomial<X>> r(e.size(),Polynomial<X>(p.argument_size()));
    for(Nat i=0; i!=r.size(); ++i) { r[i]=p; r[i]*=X(e[i]); }
    return r;
}



/*
template<class X>
OutputStream& operator<<(OutputStream& os, const Polynomial<X>& p) {
    if(p.begin()==p.end()) {
        return os << "{"<<MultiIndex::zero(p.argument_size())<<":0.0}"; }

    os << "{";
    for(typename Polynomial<X>::ConstIterator iter=p.begin(); iter!=p.end(); ++iter) {
        os << (iter==p.begin() ? "" : ",");
        for(SizeType i=0; i!=iter->key().size(); ++i) {
            os << (i==0?" ":",") << Int(iter->key()[i]); }
        os << ":" << iter->data(); }
    return os << " }";
}
*/

template<class F> struct NamedArgumentRepresentation {
    const F& function; const std::vector<String>& argument_names;
};
template<class F> NamedArgumentRepresentation<F> named_argument_repr(const F& function, const std::vector<String>& argument_names) {
    NamedArgumentRepresentation<F> r={function,argument_names}; return r; }

template<class X>
OutputStream& operator<<(OutputStream& os, const Polynomial<X>& q) {
    Bool first_term=true;
    Bool identically_zero=true;

    //os <<"[P"<<q.argument_size()<<"]";
    Polynomial<X> p=q;
    p.expansion().graded_sort();
    for(typename Polynomial<X>::ConstIterator iter=p.begin(); iter!=p.end(); ++iter) {
        MultiIndex a=iter->key();
        X v=iter->data();
        if(decide(v!=0)) {
            identically_zero=false;
            Bool first_factor=true;
            if(decide(v>0) && !first_term) { os << "+"; }
            first_term=false;
            if(decide(v==1)) { }
            else if (decide(v==-1)) { os << '-'; }
            else { os << v; first_factor=false; }
            for(Nat j=0; j!=a.size(); ++j) {
                if(a[j]!=0) {
                    if(first_factor) { first_factor=false; } else { os <<"*"; }
                    os<<"x"<<j; if(a[j]!=1) { os<<"^"<<Int(a[j]); } }
            }
            if(first_factor) { os << '1'; }
        }
    }
    if(identically_zero) { os << "0"; }
    return os;
}

template<class X>
OutputStream& operator<<(OutputStream& os, const NamedArgumentRepresentation< Polynomial<X> >& repr) {
    const Polynomial<X>& p=repr.function;
    const std::vector<String>& n=repr.argument_names;
    Bool first_term=true;
    Bool identically_zero=true;
    for(typename Polynomial<X>::ConstIterator iter=p.begin(); iter!=p.end(); ++iter) {
        MultiIndex a=iter->key();
        X v=iter->data();
        if(v!=0) {
            identically_zero=false;
            Bool first_factor=true;
            if(v>0 && !first_term) { os << "+"; }
            first_term=false;
            if(v==1) { } else if (v==-1) { os << '-'; }
            else { os << 'v'; first_factor=false; }
            for(Nat j=0; j!=a.size(); ++j) {
                if(a[j]!=0) {
                    if(first_factor) { first_factor=false; } else { os <<"*"; }
                    os<<n[j]; if(a[j]!=1) { os<<"^"<<Int(a[j]); } }
            }
            if(first_factor) { os << v; }
        }
    }
    if(identically_zero) { os << "0"; }
    return os;
}




} // namespace Ariadne
