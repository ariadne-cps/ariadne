/****************************************************************************
 *            texstream.template.h
 *
 *  Copyright  2007  Pieter Collins
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

#include "numeric/integer.h"
#include "numeric/rational.h"
#include "numeric/float64.h"
#include "numeric/floatmp.h"

namespace Ariadne {

inline Output::texstream::texstream() : _os_ptr(&std::cout) { }
inline Output::texstream::~texstream() { }
inline Output::texstream::texstream(std::ostream& os) : _os_ptr(&os) { }
inline void Output::texstream::redirect(std::ostream& os) { this->_os_ptr=&os; }
  
inline Output::texfstream::texfstream() : texstream(), _ofs_ptr(new std::ofstream())  { 
  texstream::redirect(*this->_ofs_ptr); }
inline void Output::texfstream::open(const char* fn, const char* preamble) { 
  this->_ofs_ptr->open(fn); *this->_ofs_ptr << "\\documentclass{article}\n" << preamble << "\n\\begin{document}\n"; }
inline void Output::texfstream::close() { 
  *this->_ofs_ptr << "\n\\end{document}\n"; this->_ofs_ptr->close(); }
inline Output::texfstream::~texfstream() { 
  this->_ofs_ptr->close(); delete this->_ofs_ptr; }

inline Output::texstream& Output::operator<<(texstream& txs, const char& c) { *txs._os_ptr << c; return txs; }
inline Output::texstream& Output::operator<<(texstream& txs, const char* s) { *txs._os_ptr << s; return txs; }
inline Output::texstream& Output::operator<<(texstream& txs, const int& n) { *txs._os_ptr << n; return txs; }
inline Output::texstream& Output::operator<<(texstream& txs, const double& x) { *txs._os_ptr << x; return txs; }

inline Output::texstream& Output::operator<<(texstream& txs, const Numeric::Integer& z) { *txs._os_ptr << z; return txs; }
inline Output::texstream& Output::operator<<(texstream& txs, const Numeric::Rational& q) { *txs._os_ptr << q; return txs; }
inline Output::texstream& Output::operator<<(texstream& txs, const Numeric::Float64& x) { *txs._os_ptr << x; return txs; }
inline Output::texstream& Output::operator<<(texstream& txs, const Numeric::FloatMP& x) { *txs._os_ptr << x; return txs; }

template<class R> 
Output::texstream& 
Output::operator<<(texstream& txs, const Numeric::Interval<R>& ivl) 
{
  return txs << "[" << ivl.lower() << ":" << ivl.upper() << "]";
}

template<class R> 
Output::texstream& 
Output::operator<<(texstream& txs, const LinearAlgebra::Vector<R>& v) 
{
  txs << "\\left(\\begin{array}{c}";
  for(size_type i=0; i!=v.size(); ++i) {
    if(i!=0) { txs << "\\\\"; }
    txs << v[i];
  } 
  txs << "\\end{array}\\right)";
  return txs;
}

template<class R> 
Output::texstream& 
Output::operator<<(texstream& txs, const LinearAlgebra::Matrix<R>& A) 
{
  txs << "\\left(\\begin{array}{ccccccccccc}";
  for(size_type i=0; i!=A.number_of_rows(); ++i) {
    if(i!=0) { txs << "\\\\"; }
    for(size_type j=0; j!=A.number_of_columns(); ++j) {
      if(j!=0) { txs << "&"; }
      txs << A(i,j);
    }
  } 
  txs << "\\end{array}\\right)";
  return txs;
}

template<class R> 
Output::texstream& 
Output::operator<<(texstream& txs, const Geometry::Rectangle<R>& r) 
{
  for(size_type i=0; i!=r.dimension(); ++i) {
    if(i!=0) { txs << "\\times"; }
    Numeric::Interval<R> ivl=r[i];
    txs << ivl;
    //    txs << Numeric::Interval<R>(r[i]);
  } 
  return txs;
}

template<class R> 
Output::texstream& 
Output::operator<<(texstream& txs, const Function::Polynomial<R>& p) 
{
  /*
  assert(p.result_size()==1);
  for(size_type i=0; i!=p.number_of_terms(); ++i) {
    txs << p[i];
  } 
  */
  return txs;
}

  
}

