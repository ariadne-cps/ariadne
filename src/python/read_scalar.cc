/***************************************************************************
 *            python/read_scalar.cc
 *
 *  Copyright  2007   Pieter Collins
 *  Pieter.Collins@cwi.nl
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

#include "base/exceptions.h"

#include "numeric/integer.h"
#include "numeric/rational.h"
#include "numeric/float.h"
#include "numeric/interval.h"

#include "python/float.h"
#include "python/read_scalar.h"

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>

namespace Ariadne {
namespace Python {

template<class R> void read_interval(Interval<R>& ivl, const boost::python::list& pair);
template<class R> void read_interval(Interval<R>& ivl, const boost::python::tuple& pair);
template<class R> void read_interval(Interval<R>& ivl, const boost::python::dict& pair);

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
 
template<class T>
void
read_scalar(Float<T>& x, const boost::python::object& obj)
{
  typedef Float<T> R;
  boost::python::extract<std::string> sx(obj);
  boost::python::extract<int> nx(obj);
  boost::python::extract<double> dx(obj);
  boost::python::extract<R> rx(obj);
  if(sx.check()) {
    x=static_cast<R>(sx());
  } else if(nx.check()) {
    x=static_cast<R>(nx());
  } else if(dx.check()) {
    x=static_cast<R>(dx());
  } else {
    x=rx();
  }
}

template<class T>
void
read_scalar(ApproximateFloat<T>& a, const boost::python::object& obj)
{
  Float<T> x; 
  read_scalar(x,obj);
  a=ApproximateFloat<T>(x);
}


template<class R> 
void
read_scalar(Interval<R>& x, const boost::python::object& obj)
{
  // The calls lx.check() and tx.check() produce compiler warnings 
  //   "dereferencing type-punned pointer will break strict-aliasing rules"
  // but this code works, and I don't how to change the extract type
  // to disable the warning without making the check fail.
  typedef Interval<R> I;
  boost::python::extract<std::string> sx(obj);
  boost::python::extract<boost::python::list> lx(obj);
  boost::python::extract<boost::python::tuple> tx(obj);
  boost::python::extract<boost::python::dict> mx(obj);
  boost::python::extract<int> nx(obj);
  boost::python::extract<double> dx(obj);
  boost::python::extract<R> rx(obj);
  boost::python::extract<I> ix(obj);

  if(sx.check()) {
    x=static_cast<I>(sx());
  } else if(lx.check()) {
    read_interval(x,lx());
  } else if(tx.check()) {
    read_interval(x,tx());
  } else if(mx.check()) {
    read_interval(x,mx());
  } else if(nx.check()) {
    x=static_cast<I>(nx());
  } else if(dx.check()) {
    x=static_cast<I>(dx());
  } else if(rx.check()) {
    x=static_cast<I>(rx());
  } else if(ix.check()) {
    x=ix();
  }
  
}


template void read_scalar(ApproximateFloatPy&, const boost::python::object&);
template void read_scalar(FloatPy&, const boost::python::object&);
template void read_scalar(IntervalPy&, const boost::python::object&);
template void read_scalar(Interval<Rational>&, const boost::python::object&);



template<class R> inline
void
read_interval(Interval<R>& ivl, const boost::python::list& pair)
{
  if(boost::python::len(pair)!=2) {
    throw std::runtime_error("Interval must be list of pairs representing intervals");
  }
  R& l=ivl._lower; R& u=ivl._upper;
  read_scalar(l,pair[0]); read_scalar(u,pair[1]);
}

template<class R> inline
void
read_interval(Interval<R>& ivl, const boost::python::tuple& pair)
{
  if(boost::python::len(pair)!=2) {
    throw std::runtime_error("Interval must be list of pairs representing intervals");
  }
  R& l=ivl._lower; R& u=ivl._upper;
  read_scalar(l,pair[0]); read_scalar(u,pair[1]);
}

template<class R> inline
void
read_interval(Interval<R>& ivl, const boost::python::dict& pair)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
  //if(boost::python::len(pair)!=1) {
  //  throw std::runtime_error("Interval must be a pair [a,b] or {a:b}.");
  //}
  //R& l=ivl._lower; R& u=ivl._upper;
  //read_scalar(l,pair[0].key()); read_scalar(u,pair[1].value());
}



}
}
