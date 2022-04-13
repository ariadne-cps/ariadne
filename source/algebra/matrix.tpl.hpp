/***************************************************************************
 *            matrix.tpl.hpp
 *
 *  Copyright  2005-20  Alberto Casagrande, Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "algebra/vector.hpp"
#include "algebra/covector.hpp"

namespace Ariadne {

template<class X> inline X create_zero();

template<class X> Matrix<X>::Matrix(SizeType m, SizeType n, const X* p)
    : _zero(nul(p[0])), _rs(m), _cs(n), _ary(p,p+m*n) {
}


template<class X> Matrix<X>::Matrix(InitializerList<InitializerList<X>> lst)
    : _zero(nul(*lst.begin()->begin())), _rs(lst.size()), _cs(lst.begin()->size()), _ary(_rs*_cs,_zero)
{
    typename InitializerList<InitializerList<X>>::const_iterator row_iter=lst.begin();
    for(SizeType i=0; i!=this->row_size(); ++i, ++row_iter) {
        ARIADNE_PRECONDITION(row_iter->size()==this->column_size());
        typename InitializerList<X>::const_iterator col_iter=row_iter->begin();
        for(SizeType j=0; j!=this->column_size(); ++j, ++col_iter) {
            this->at(i,j)=*col_iter;
        }
    }
}

template<class X> Matrix<X> join(Matrix<X> const& A1, Matrix<X> const& A2) {
    ARIADNE_PRECONDITION(A1.column_size()==A2.column_size());
    Matrix<X> R(A1.row_size()+A2.row_size(),A1.column_size(),A1.zero_element());
    const SizeType m1=A1.row_size();
    for(SizeType i=0; i!=A1.row_size(); ++i) {
        for(SizeType j=0; j!=A1.column_size(); ++j) {
            R[i][j]=A1[i][j];
        }
    }
    for(SizeType i=0; i!=A2.row_size(); ++i) {
        for(SizeType j=0; j!=A2.column_size(); ++j) {
            R[m1+i][j]=A2[i][j];
        }
    }
    return R;
}

template<class X> Matrix<X> join(Matrix<X> const& A1, Covector<X> const& u2) {
    ARIADNE_PRECONDITION(A1.column_size()==u2.size());
    Matrix<X> R(A1.row_size()+1u,A1.column_size(),A1.zero_element());
    const SizeType m1=A1.row_size();
    for(SizeType i=0; i!=A1.row_size(); ++i) {
        for(SizeType j=0; j!=A1.column_size(); ++j) {
            R[i][j]=A1[i][j];
        }
    }
    for(SizeType j=0; j!=u2.size(); ++j) {
        R[m1][j]=u2[j];
    }
    return R;
}

template<class X> Matrix<X> join(Covector<X> const& u1, Matrix<X> const& A2) {
    ARIADNE_PRECONDITION(u1.size()==A2.column_size());
    Matrix<X> R(1u+A2.row_size(),A2.column_size(),A2.zero_element());
    for(SizeType j=0; j!=u1.size(); ++j) {
        R[0u][j]=u1[j];
    }
    for(SizeType i=0; i!=A2.row_size(); ++i) {
        for(SizeType j=0; j!=A2.column_size(); ++j) {
            R[1u+i][j]=A2[i][j];
        }
    }
    return R;
}

template<class X> Matrix<X> cojoin(Matrix<X> const& A1, Matrix<X> const& A2) {
    ARIADNE_PRECONDITION(A1.row_size()==A2.row_size());
    Matrix<X> R(A1.row_size(),A1.column_size()+A2.column_size(),A1.zero_element());
    const SizeType n1=A1.column_size();
    for(SizeType i=0; i!=A1.row_size(); ++i) {
        for(SizeType j=0; j!=A1.column_size(); ++j) {
            R[i][j]=A1[i][j];
        }
        for(SizeType j=0; j!=A2.column_size(); ++j) {
            R[i][n1+j]=A2[i][j];
        }
    }
    return R;
}

template<class X> Matrix<X> cojoin(Matrix<X> const& A1, Vector<X> const& v2) {
    ARIADNE_PRECONDITION(A1.row_size()==v2.size());
    Matrix<X> R(A1.row_size(),A1.column_size()+1u,A1.zero_element());
    const SizeType n1=A1.column_size();
    for(SizeType i=0; i!=A1.row_size(); ++i) {
        for(SizeType j=0; j!=A1.column_size(); ++j) {
            R[i][j]=A1[i][j];
        }
        R[i][n1]=v2[i];
    }
    return R;
}

template<class X> Matrix<X> cojoin(Vector<X> const& v1, Matrix<X> const& A2) {
    ARIADNE_PRECONDITION(v1.size()==A2.row_size());
    Matrix<X> R(A2.row_size(),1u+A2.column_size(),A2.zero_element());
    for(SizeType i=0; i!=A2.row_size(); ++i) {
        R[i][0]=v1[i];
        for(SizeType j=0; j!=A2.column_size(); ++j) {
            R[i][1u+j]=A2[i][j];
        }
    }
    return R;
}


#ifdef ARIADNE_OMIT
template<class X> InputStream& Matrix<X>::read(InputStream& is) {
    Matrix<X>& A=*this;
    char c;
    is >> c;
    is.putback(c);
    if(c=='[') {
        is >> c;
        /* Representation as a literal [a11,a12,...,a1n; a21,a22,...a2n; ... ; am1,am2,...,amn] */
        std::vector< std::vector<X> > v;
        X x;
        c=';';
        while(is && c==';') {
            v.push_back(std::vector<X>());
            c=',';
            while(is && c==',') {
                is >> x;
                v.back().push_back(x);
                is >> c;
            }
        }
        if(is) {
            A=Matrix<X>(v.size(),v.front().size());
            for(SizeType i=0; i!=A.row_size(); ++i) {
                if(v[i].size()!=A.column_size()) {
                    // TOOD: Add exception
                    assert(false);
                    //ARIADNE_THROW(InvalidInput,"Matrix::read(istream&)","row[0].size()="<<v[0].size()<<", row["<<i<<"].size()="<<v[i].size());
                }
                for(SizeType j=0; j!=A.column_size(); ++j) {
                    A[i][j]=v[i][j];
                }
            }
        }
    }
    else {
        // TOOD: Add exception
        assert(false);
        //ARIADNE_THROW(InvalidInput,"Matrix::read(istream&)"," separator c="<<c);
    }
    return is;
}
#endif

template<class X> decltype(mag(declval<X>())) sup_norm(const Matrix<X>& A)
{
    typedef decltype(mag(declval<X>())) R;
    R zero=mag(A.zero_element());
    R result=zero;
    for(SizeType i=0; i!=A.row_size(); ++i) {
        R row_sum=zero;
        for(SizeType j=0; j!=A.column_size(); ++j) {
            row_sum+=mag(A[i][j]);
        }
        // NOTE: The arguments must be this way round to propagate a nan row_sum
        result=max(row_sum,result);
    }
    return result;
}

template<class X> decltype(declval<X>()+mag(declval<X>())) log_norm(Matrix<X> const& A)
{
    ARIADNE_PRECONDITION(A.row_size()==A.column_size());
    typedef decltype(declval<X>()+mag(declval<X>())) R;
    R r=A.zero_element(); r=-inf;
    for(SizeType i=0; i!=A.row_size(); ++i) {
        R t=A[i][i];
        for(SizeType j=0; j!=A.column_size(); ++j) {
            if(j!=i) {
                t+=mag(A[i][j]);
            }
        }
        r=max(r,t);
    }
    return r;
}


} // namespace Ariadne
