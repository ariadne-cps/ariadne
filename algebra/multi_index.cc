/***************************************************************************
 *            multi_index.cc
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

#include "multi_index.h"

namespace Ariadne {

OutputStream& operator<<(OutputStream& os, const MultiIndex& a) {
    //os << "("<<Int(a.degree());
    //for(MultiIndex::SizeType i=0; i!=a.size(); ++i) { os << (i==0?';':',') << Int(a[i]); }
    if(a.size()==0) { return os << "()"; }
    os << '(' << a.degree();
    for(MultiIndex::SizeType i=0; i!=a.size(); ++i) { os << (i==0?';':',') << Int(a[i]); }
    return os << ')';
}

Bool graded_less(const MultiIndex& a1, const MultiIndex& a2) {
    if(a1.size()!=a2.size()) {
        throw std::runtime_error("graded_less(MultiIndex,MultiIndex): number of variables must match");
    }

    if(a1.degree()!=a2.degree()) {
        return a1.degree()<a2.degree();
    } else {
        for(SizeType i=0; i!=a1.size(); ++i) {
            if(a1[i]!=a2[i]) {
                return a1[i]>a2[i];
            }
        }
        return false;
        //for(size_type j=0; j!=a1.word_size(); ++j) {
        for(SizeType j=a1.size(); j!=0u;) {
            --j;
            if(a1[j]!=a2[j]) {
                return a1[j]<a2[j];
            }
        }
        return false;
    }
}

Bool lexicographic_less(const MultiIndex& a1, const MultiIndex& a2) {
    if(a1.size()!=a2.size()) {
        throw std::runtime_error("lexicographic_less(MultiIndex,MultiIndex): number of variables must match");
    }

    for(SizeType i=0; i!=a1.size(); ++i) {
        if(a1[i]!=a2[i]) {
            return a1[i]<a2[i];
        }
    }
    return false;
}

Bool reverse_lexicographic_less(const MultiIndex& a1, const MultiIndex& a2) {
    if(a1.size()!=a2.size()) {
        throw std::runtime_error("reverse_lexicographic_less(MultiIndex,MultiIndex): number of variables must match");
    }

    SizeType i=a1.size();
    while(i!=0) {
        --i;
        if(a1[i]!=a2[i]) {
            return a1[i]>a2[i];
        }
    }
    return false;
}



typedef SizeType BigSizeType;

BigSizeType MultiIndex::number() const
{
    BigSizeType result=fac(this->degree());
    for(SizeType k=0; k!=this->size(); ++k) {
        result/=fac((*this)[k]);
    }
    return result;
}

BigSizeType MultiIndex::factorial() const
{
    BigSizeType result=1;
    for(SizeType k=0; k!=this->size(); ++k) {
        result*=fac((*this)[k]);
    }
    return result;
}

BigSizeType MultiIndex::position() const
{
    DegreeType deg=this->degree()-1;
    DegreeType nvar=this->size();
    BigSizeType result=bin(deg+nvar,nvar);
    for(SizeType k=0; k!=this->size()-1; ++k) {
        --nvar;
        deg-=(*this)[k];
        result+=bin(deg+nvar,nvar);
    }
    return result;
}


MultiIndex& MultiIndex::operator++()
{
    //std::cerr<<"MultiIndex::operator++() with *this="<<*this<<" "<<std::flush;
    assert(_n>0);

    SizeType const n=this->_n;
    IndexType* const p=reinterpret_cast<IndexType*>(this->_p);

    if(n==1) {
        ++p[0];
        ++p[n];
        return *this;
    }
    if(p[n-2]!=0) {
        --p[n-2];
        ++p[n-1];
        return *this;
    } else {
        IndexType li=p[n-1];
        p[n-1]=0;
        for(SizeType k=n-1; k!=0; --k) {
            if(p[k-1]!=0) {
                --p[k-1];
                p[k]=li+1;
                return *this;
            }
        }
        p[0]=li+1;
        ++p[n];
    }
    return *this;
}



Array<SizeType> complement(SizeType nmax, Array<SizeType> vars) {
    Array<SizeType> cmpl(nmax-vars.size());
    SizeType kr=0; SizeType kv=0;
    for(SizeType j=0; j!=nmax; ++j) {
        if(kv==vars.size() || j!=vars[kv]) {
            cmpl[kr]=j; ++kr;
        } else {
            ++kv;
        }
    }
    return cmpl;
}


} // namespace Ariadne
