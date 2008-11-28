/***************************************************************************
 *            utilities.cc
 *
 *  Copyright  2005-8  Alberto Casagrande, Pieter Collins
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

#include "config.h"

#include "utilities.h"

#include "numeric.h"

namespace Ariadne {

void read_interval(Interval& ivl, const boost::python::list& pair);
void read_interval(Interval& ivl, const boost::python::tuple& pair);

template<class X> void read_array(array<X>& ary, const boost::python::object& obj);


void
read_scalar(bool& n, const boost::python::object& obj)
{
    n=boost::python::extract<bool>(obj);
}

void
read_scalar(int& n, const boost::python::object& obj)
{
    n=boost::python::extract<int>(obj);
}

void
read_scalar(long int& n, const boost::python::object& obj)
{
    n=boost::python::extract<long int>(obj);
}

// Read a scalar variable of type X from a Python object
void
read_scalar(unsigned int& n, const boost::python::object& obj)
{
    n=boost::python::extract<unsigned int>(obj);
}

void
read_scalar(unsigned long int& n, const boost::python::object& obj)
{
    n=boost::python::extract<unsigned long int>(obj);
}


#ifdef HAVE_GMPXX_H

typedef mpz_class Integer;

/*


void
read_scalar(Integer& z, const boost::python::object& obj)
{
    boost::python::extract<std::string> sz(obj);
    boost::python::extract<int> nz(obj);
    boost::python::extract<Integer> zz(obj);
    if(sz.check()) {
        z=static_cast<Integer>(sz());
    } else if(nz.check()) {
        z=static_cast<Integer>(nz());
    } else {
        z=zz();
    }
}
 
*/

void
read_scalar(Rational& q, const boost::python::object& obj)
{
    boost::python::extract<std::string> sq(obj);
    boost::python::extract<int> nq(obj);
    boost::python::extract<double> dq(obj);
    boost::python::extract<Integer> zq(obj);
    boost::python::extract<Rational> qq(obj);
    if(sq.check()) {
        q=static_cast<Rational>(sq());
    } else if(nq.check()) {
        q=static_cast<Rational>(nq());
    } else if(dq.check()) {
        q=static_cast<Rational>(dq());
    } else if(zq.check()) {
        q=static_cast<Rational>(zq());
    } else {
        q=static_cast<Rational>(qq());
    }
}

#endif


void
read_scalar(Float& x, const boost::python::object& obj)
{
    boost::python::extract<int> nx(obj);
    boost::python::extract<double> dx(obj);
    boost::python::extract<Float> rx(obj);
    if(nx.check()) {
        x=static_cast<Float>(nx());
    } else if(dx.check()) {
        x=static_cast<Float>(dx());
    } else {
        x=rx();
    }
}


void
read_scalar(Interval& x, const boost::python::object& obj)
{
    // The calls lx.check() and tx.check() produce compiler warnings 
    //   "dereferencing type-punned pointer will break strict-aliasing rules"
    // but this code works, and I don't how to change the extract type
    // to disable the warning without making the check fail.
    boost::python::extract<boost::python::list> lx(obj);
    boost::python::extract<boost::python::tuple> tx(obj);
    boost::python::extract<int> nx(obj);
    boost::python::extract<double> dx(obj);
    boost::python::extract<Float> rx(obj);
    boost::python::extract<Interval> ix(obj);

    if(lx.check()) {
        read_interval(x,lx());
    } else if(tx.check()) {
        read_interval(x,tx());
    } else if(nx.check()) {
        x=static_cast<Interval>(nx());
    } else if(dx.check()) {
        x=static_cast<Interval>(dx());
    } else if(rx.check()) {
        x=static_cast<Interval>(rx());
    } else if(ix.check()) {
        x=ix();
    }
  
}


inline
void
read_interval(Interval& ivl, const boost::python::list& pair)
{
    if(boost::python::len(pair)!=2) {
        throw std::runtime_error("Interval must be list of pairs representing intervals");
    }
    Float& l=ivl.l; Float& u=ivl.u;
    read_scalar(l,pair[0]); read_scalar(u,pair[1]);
}

inline
void
read_interval(Interval& ivl, const boost::python::tuple& pair)
{
    if(boost::python::len(pair)!=2) {
        throw std::runtime_error("Interval must be list of pairs representing intervals");
    }
    Float& l=ivl.l; Float& u=ivl.u;
    read_scalar(l,pair[0]); read_scalar(u,pair[1]);
}


} // namespace Ariadne
