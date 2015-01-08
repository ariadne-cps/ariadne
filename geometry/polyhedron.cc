/***************************************************************************
 *            polyhedron.cc
 *
 *  Copyright  2006-8  Alberto Casagrande, Pieter Collins
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
 *  along with this program; if not, write to bouthe Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "function/functional.h"
#include "config.h"

#include <iostream>
#include <string>
#include <sstream>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "geometry/polyhedron.h"

#include "utility/macros.h"
#include "utility/exceptions.h"

#include "numeric/numeric.h"

#include "algebra/vector.h"
#include "algebra/matrix.h"

#include "function/function.h"

#include "geometry/box.h"
#include "geometry/polytope.h"



namespace Ariadne {


extern Int global_verbosity;
Int verbosity=global_verbosity;


Void ddconv(std::vector< Vector<Float> >&, const std::vector< Vector<Float> >&);

Polyhedron polyhedron(const ExactBox& bx);
Polyhedron polyhedron(const Polytope& p);
Polytope polytope(const Polyhedron& p);



Polyhedron::Polyhedron()
    : _A(), _b()
{
}


Polyhedron::Polyhedron(Nat d)
    : _A(0,d), _b(0)
{
}


Polyhedron::Polyhedron(const Matrix<Float>& A, const Vector<Float>& b)
    : _A(A), _b(b)
{
    ARIADNE_ASSERT_MSG(A.column_size()==b.size(),"Invalid sizes of A="<<A<<" and b="<<b<<" for polyhedral constraints");
}




Polyhedron::Polyhedron(const ExactBox& bx)
    : _A(bx.dimension()*2,bx.dimension()), _b(bx.dimension()*2)
{
    const Nat n=bx.dimension();
    for(Nat i=0; i!=n; ++i) {
        _A[i][i]=-1;
        _b[i]=-bx[i].lower();
        _A[i+n][i]=+1;
        _b[i+n]=+bx[i].upper();
    }
}

Polyhedron::Polyhedron(const Polytope& p)
{
    ARIADNE_NOT_IMPLEMENTED;
}

Polyhedron* Polyhedron::clone() const
{
    return new Polyhedron(*this);
}



Matrix<Float>
Polyhedron::A() const
{
    return this->_A;
}

Vector<Float>
Polyhedron::b() const
{
    return this->_b;
}



Nat
Polyhedron::dimension() const
{
    return this->_A.column_size();
}

Tribool
Polyhedron::empty() const
{
    ARIADNE_NOT_IMPLEMENTED;
}


Tribool
Polyhedron::bounded() const
{
    ARIADNE_NOT_IMPLEMENTED;
}


Polyhedron
Polyhedron::halfspace(SizeType i) const
{
    ARIADNE_NOT_IMPLEMENTED;
}



Tribool Polyhedron::overlaps(const ExactBox& bx) const {
    ARIADNE_NOT_IMPLEMENTED;
}

Tribool Polyhedron::covers(const ExactBox& bx) const {
    ARIADNE_NOT_IMPLEMENTED;
}

Tribool Polyhedron::separated(const ExactBox& bx) const {
    ARIADNE_NOT_IMPLEMENTED;
}



Polyhedron
intersection(const Polyhedron& plhd1, const Polyhedron& plhd2)
{
    ARIADNE_ASSERT(plhd1.dimension()==plhd2.dimension());
    Nat d=plhd1.dimension();
    SizeType nc1=plhd1.number_of_constraints();
    SizeType nc2=plhd2.number_of_constraints();
    Matrix<Float> A(nc1+nc2,d);
    project(A,range(0,nc1),range(0,d)) = plhd1.A();
    project(A,range(nc1,nc1+nc2),range(0,d)) = plhd1.A();
    Vector<Float> b=join(plhd1.b(),plhd2.b());
    return Polyhedron(A,b);
}





Polyhedron
polyhedron(const ExactBox& bx)
{
    return Polyhedron(bx);
}

Polyhedron
polyhedron(const Polytope& pltp)
{
    ARIADNE_NOT_IMPLEMENTED;
}



Polytope
polytope(const Polyhedron& pltp)
{
    ARIADNE_NOT_IMPLEMENTED;
}







OutputStream&
Polyhedron::write(OutputStream& os) const
{
    //return os << "Polyhedron( A=" << this->A() << ", b=" << this->b() << " )";
    const Matrix<Float> A=this->A();
    const Vector<Float> b=this->b();
    os << "Polyhedron( constraints=";
    Nat d=this->dimension();
    SizeType nc=this->number_of_constraints();
    for(SizeType i=0; i!=nc; ++i) {
        os << ( i==0 ? "[" : "," );
        for(SizeType j=0; j!=d; ++j) {
            os << ( j==0 ? "(" : ",");
            os << A[i][j];
        }
        os << ":" << b[i] << ")";
    }
    os << "] )";
    return os;
}

InputStream&
operator>>(InputStream& is, Polyhedron& p)
{
    std::vector< std::vector<Float> > Alst;
    std::vector< Float > Blst;

    std::vector<Float> a;
    Float b;

    char c;
    is >> c;
    assert(c=='[');

    c=is.peek();
    while(c=='[') {
        // Read constraint ax<=b in form [a_1,a_2,...,a_n;b];
        read_sequence(is,a,'[',';',',');
        is >> b;
        is >> c;
        assert(c==']');
        Alst.push_back(a);
        Blst.push_back(b);
    }

    SizeType m=Alst.size();
    SizeType n=Alst[0].size();
    Matrix<Float> A(m,n);
    Vector<Float> B(m);
    for(Nat i=0; i!=m; ++i) {
        for(SizeType j=0; j!=n; ++j) {
            A[i][j]=Alst[i][j];
        }
        B[i]=Blst[i];
    }

    p=Polyhedron(A,B);

    return is;
}




} // namespace Ariadne

