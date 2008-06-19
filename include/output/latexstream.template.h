/****************************************************************************
 *            latexstream.template.h
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

inline latexstream::latexstream() : _os_ptr(&std::cout) { }
inline latexstream::~latexstream() { }
inline latexstream::latexstream(std::ostream& os) : _os_ptr(&os) { }
inline void latexstream::redirect(std::ostream& os) { this->_os_ptr=&os; }
  
inline latexfstream::latexfstream() : latexstream(), _ofs_ptr(new std::ofstream())  { 
  latexstream::redirect(*this->_ofs_ptr); }
inline void latexfstream::open(const char* fn, const char* preamble) { 
  this->_ofs_ptr->open(fn); *this->_ofs_ptr << "\\documentclass{article}\n" << preamble << "\n\\begin{document}\n"; }
inline void latexfstream::close() { 
  *this->_ofs_ptr << "\n\\end{document}\n"; this->_ofs_ptr->close(); }
inline latexfstream::~latexfstream() { 
  this->close(); delete this->_ofs_ptr; }

inline latexstream& operator<<(latexstream& txs, const char& c) { *txs._os_ptr << c; return txs; }
inline latexstream& operator<<(latexstream& txs, const char* s) { *txs._os_ptr << s; return txs; }
inline latexstream& operator<<(latexstream& txs, const int& n) { *txs._os_ptr << n; return txs; }
inline latexstream& operator<<(latexstream& txs, const long int& n) { *txs._os_ptr << n; return txs; }
inline latexstream& operator<<(latexstream& txs, const unsigned int& n) { *txs._os_ptr << n; return txs; }
inline latexstream& operator<<(latexstream& txs, const unsigned long int& n) { *txs._os_ptr << n; return txs; }
inline latexstream& operator<<(latexstream& txs, const double& x) { *txs._os_ptr << x; return txs; }

inline latexstream& operator<<(latexstream& txs, const Integer& z) { *txs._os_ptr << z; return txs; }
inline latexstream& operator<<(latexstream& txs, const Rational& q) { *txs._os_ptr << q; return txs; }
inline latexstream& operator<<(latexstream& txs, const Float64& x) { *txs._os_ptr << x; return txs; }
inline latexstream& operator<<(latexstream& txs, const FloatMP& x) { *txs._os_ptr << x; return txs; }

template<class R> 
latexstream& 
operator<<(latexstream& txs, const Interval<R>& ivl) 
{
  return txs << "[" << ivl.lower() << ":" << ivl.upper() << "]";
}

template<class R> 
latexstream& 
operator<<(latexstream& txs, const Vector<R>& v) 
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
latexstream& 
operator<<(latexstream& txs, const Matrix<R>& A) 
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
latexstream& 
operator<<(latexstream& txs, const Box<R>& r) 
{
  for(size_type i=0; i!=r.dimension(); ++i) {
    if(i!=0) { txs << "\\times"; }
    Interval<R> ivl=r[i];
    txs << ivl;
    //    txs << Interval<R>(r[i]);
  } 
  return txs;
}

  
} // namespace Ariadne


